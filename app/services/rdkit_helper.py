from __future__ import annotations

import importlib
import logging
import math
from typing import Any

from app.config import get_settings

logger = logging.getLogger(__name__)

BOLTZMANN_R_KCAL = 0.0019872041  # kcal / (mol * K)
BODY_TEMPERATURE_K = 298.15

# RDKit's monolithic ``AllChem`` (and ``Draw.rdMolDraw2D``) pull in ``rdMolAlign``,
# which Windows Application Control blocks on this machine. Importing the narrow
# submodules instead lets all the other chemistry keep working.
_SUBMODULES = {
    "rdMolDescriptors": "rdkit.Chem.rdMolDescriptors",
    "rdDepictor": "rdkit.Chem.rdDepictor",
    "rdmolfiles": "rdkit.Chem.rdmolfiles",
    "rdDistGeom": "rdkit.Chem.rdDistGeom",
    "rdForceFieldHelpers": "rdkit.Chem.rdForceFieldHelpers",
}

_CAPABILITY_CACHE: dict[str, Any] | None = None


def rdkit_compat(reload: bool = False) -> dict[str, Any]:
    """Import only the RDKit pieces this OS policy allows, degrading per-module.

    On the first successful import the result is cached for the rest of the process
    (the import is expensive and stable). A failure is deliberately NOT cached: a
    transient problem — a module not yet on ``sys.path``, or a slow first import —
    can succeed on the next request, so every call retries the import while RDKit is
    unavailable. The exception type and message are logged and carried in the
    returned ``note`` so callers can report *why* it is unavailable.
    """
    global _CAPABILITY_CACHE
    if _CAPABILITY_CACHE is not None and not reload:
        return _CAPABILITY_CACHE
    mods: dict[str, Any] = {"Chem": None}
    try:
        from rdkit import Chem

        mods["Chem"] = Chem
    except Exception as exc:
        detail = f"{type(exc).__name__}: {exc}"
        logger.warning("RDKit import failed (%s); not caching this failure", detail)
        return {**mods, "note": detail}
    for name, path in _SUBMODULES.items():
        try:
            mods[name] = importlib.import_module(path)
        except Exception as exc:
            logger.warning("RDKit submodule %s import failed (%s: %s)", name, type(exc).__name__, exc)
            mods[name] = None
    _CAPABILITY_CACHE = mods
    return mods


def rdkit_unavailable_reason(mods: dict[str, Any] | None = None) -> str:
    """One sentence explaining why RDKit is unavailable, including the import error.

    The exception type and message captured by :func:`rdkit_compat` are folded in,
    so an "unavailable" reason never hides the actual cause.
    """
    caps = mods if mods is not None else rdkit_compat()
    detail = caps.get("note")
    if detail:
        return f"RDKit is not importable in this environment ({detail})"
    return "RDKit is not importable in this environment"


def rdkit_capabilities() -> dict[str, bool]:
    mods = rdkit_compat()
    return {
        "rdkit": mods.get("Chem") is not None,
        "descriptors": mods.get("rdMolDescriptors") is not None,
        "conformer_3d": mods.get("rdDistGeom") is not None,
        "force_field": mods.get("rdForceFieldHelpers") is not None,
        "depiction_2d": mods.get("rdDepictor") is not None,
    }


def normalize_smiles_safe(smiles: str) -> str | None:
    s = smiles.strip()
    if not s:
        return None
    try:
        from rdkit import Chem
    except Exception:
        return s

    mol = Chem.MolFromSmiles(s)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol, isomericSmiles=True)


def _calc_crippen_logp(mol: Any, rd_mol_descriptors: Any) -> float | None:
    """Crippen logP, tolerating RDKit builds/OS policies where it is unavailable."""
    fn = getattr(rd_mol_descriptors, "CalcCrippenMolLogP", None)
    if fn is not None:
        try:
            return float(fn(mol))
        except Exception:
            pass
    try:
        from rdkit.Chem import Crippen

        return float(Crippen.MolLogP(mol))
    except Exception:
        return None


def rdkit_properties_from_smiles(smiles: str) -> dict[str, Any]:
    mods = rdkit_compat()
    Chem = mods.get("Chem")
    rd_mol_descriptors = mods.get("rdMolDescriptors")
    if Chem is None or rd_mol_descriptors is None:
        raise RuntimeError(rdkit_unavailable_reason(mods))

    # rdFreeSASA (and rdPartialCharges, which backs Crippen logP) can be blocked
    # by OS application-control policy. Neither is required for the core
    # descriptors, so treat them as optional rather than failing the whole call.
    try:
        from rdkit.Chem import rdFreeSASA
    except Exception:
        rdFreeSASA = None

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("SMILES could not be parsed by RDKit")

    mw = rd_mol_descriptors.CalcExactMolWt(mol)
    logp = _calc_crippen_logp(mol, rd_mol_descriptors)
    # These two fields feed Lipinski's rule of five (the violation checks below),
    # and that rule counts every N and O atom, so use the Lipinski variants.
    # CalcNumHBA/CalcNumHBD implement a stricter definition that under-counts:
    # on RDKit 2026.03.6 CalcNumHBA(aspirin) is 3, while the molecule really has
    # four acceptor oxygens (two carbonyl, one ester, one carboxyl). Note that
    # Chem.Lipinski.NumHAcceptors is an alias for the strict count and is also 3
    # on this build; only CalcNumLipinskiHBA returns 4.
    h_donors = rd_mol_descriptors.CalcNumLipinskiHBD(mol)
    h_acceptors = rd_mol_descriptors.CalcNumLipinskiHBA(mol)
    # CalcNumRotatableBonds defaults to RDKit's *strict* definition, which also
    # hides the ester/carboxyl C-O bond next to a carbonyl: it reports 2 rotatable
    # bonds for aspirin and 1 for ethyl acetate. This field feeds Veber's rule,
    # and Veber uses the conventional count (a single, non-ring bond between two
    # heavy atoms that each have another heavy neighbour), so ask for NonStrict
    # explicitly. Verified on RDKit 2026.03.6 against ethanol, butane, toluene,
    # propanoic acid, ethyl acetate and aspirin.
    rot_bonds = rd_mol_descriptors.CalcNumRotatableBonds(
        mol, rd_mol_descriptors.NumRotatableBondsOptions.NonStrict
    )
    formula = rd_mol_descriptors.CalcMolFormula(mol)
    tpsa = rd_mol_descriptors.CalcTPSA(mol)
    exact_mass = float(mw)
    heavy_atoms = mol.GetNumHeavyAtoms()
    rings = rd_mol_descriptors.CalcNumRings(mol)
    aromatic_rings = rd_mol_descriptors.CalcNumAromaticRings(mol)
    fraction_csp3 = rd_mol_descriptors.CalcFractionCSP3(mol)

    violations = 0
    if mw > 500:
        violations += 1
    if logp is not None and logp > 5:
        violations += 1
    if h_donors > 5:
        violations += 1
    if h_acceptors > 10:
        violations += 1

    # Veber's oral bioavailability rules (rotatable bonds <= 10, TPSA <= 140).
    veber_pass = (rot_bonds <= 10) and (tpsa <= 140.0)

    sasa = None
    if rdFreeSASA is not None:
        try:
            sasa = rdFreeSASA.CalcSASA(mol)
        except Exception:
            sasa = None

    return {
        "molecular_weight": round(float(mw), 2),
        "logp": round(float(logp), 2) if logp is not None else None,
        "h_bond_donors": int(h_donors),
        "h_bond_acceptors": int(h_acceptors),
        "rotatable_bonds": int(rot_bonds),
        "formula": formula,
        "exact_mass": exact_mass,
        "tpsa": round(float(tpsa), 2),
        "sasa": round(float(sasa), 2) if sasa is not None else None,
        "heavy_atom_count": int(heavy_atoms),
        "ring_count": int(rings),
        "aromatic_ring_count": int(aromatic_rings),
        "fraction_csp3": round(float(fraction_csp3), 3),
        "veber": {
            "pass": bool(veber_pass),
            "rotatable_bonds": int(rot_bonds),
            "tpsa": round(float(tpsa), 2),
        },
        "lipinski_rule_of_five": {
            "pass": (violations == 0) if logp is not None else None,
            "violations": int(violations),
        },
    }


def relative_energies(energies: list[float]) -> list[float]:
    """Shift a set of conformer energies so the lowest minimum is 0.0 kcal/mol.

    Absolute MMFF94 energies are not comparable between runs: they depend on the
    conformer pool (``numConfs`` and the pruning threshold), the embedding seed,
    the force-field variant and the RDKit build. Only the differences inside one
    ensemble carry meaning, so those are what get reported.
    """
    if not energies:
        return []
    lowest = min(energies)
    return [round(e - lowest, 3) for e in energies]


def boltzmann_populations(
    energies: list[float],
    temperature_k: float = BODY_TEMPERATURE_K,
) -> list[float]:
    """Boltzmann populations for a set of conformer energies in kcal/mol.

    These are estimates, not measured populations: they treat the MMFF94
    electronic energies as if they were free energies, which ignores entropy and
    any solvent effect. Use them to say which shape dominates, not to quote an
    exact fraction.
    """
    if not energies:
        return []
    rt = BOLTZMANN_R_KCAL * float(temperature_k)
    lowest = min(energies)
    weights = [math.exp(-(e - lowest) / rt) for e in energies]
    total = sum(weights)
    if total <= 0:
        return []
    return [round(w / total, 4) for w in weights]


POPULATION_ESTIMATE_NOTE = (
    "Populations are Boltzmann estimates computed from MMFF94 energies at 298.15 K. "
    "They treat electronic energy as free energy, so they ignore entropy and solvent "
    "and should be read as which conformer dominates, not as an exact fraction."
)

ABSOLUTE_ENERGY_NOTE = (
    "Absolute MMFF94 energies are not comparable between runs. The same molecule "
    "reports a different absolute minimum when the conformer pool changes (numConfs, "
    "RMS pruning threshold), when the embedding seed changes, or on a different RDKit "
    "build. Only relative_energy_kcal_per_mol, which is measured within a single "
    "ensemble, is meaningful."
)


def _prepare_mol(Chem: Any, smiles: str) -> Any:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("SMILES could not be parsed by RDKit")
    return mol


def _molblock_from_conformer(Chem: Any, rdmolfiles: Any, mol: Any, conf_id: int = 0) -> str | None:
    if rdmolfiles is None:
        return None
    target = mol
    if mol.GetNumConformers() == 0:
        return None
    try:
        if conf_id != 0:
            target = Chem.Mol(mol)
            conf = mol.GetConformer(conf_id)
            new_conf = Chem.Conformer(conf)
            target.RemoveAllConformers()
            target.AddConformer(new_conf, assignId=True)
        return rdmolfiles.MolToMolBlock(target, confId=0)
    except Exception:
        try:
            return rdmolfiles.MolToMolBlock(mol)
        except Exception:
            return None


def _is_conformer_3d(Chem: Any, mol: Any) -> bool:
    try:
        conf = mol.GetConformer(0)
    except Exception:
        return False
    is_3d = getattr(conf, "Is3D", None)
    if callable(is_3d):
        try:
            return bool(is_3d())
        except Exception:
            pass
    zs = [conf.GetAtomPosition(i).z for i in range(mol.GetNumAtoms())]
    return any(abs(z) > 1e-3 for z in zs)


def _optimise_conformers(ff: Any, mol: Any, max_iters: int) -> tuple[list[Any], str] | tuple[None, None]:
    """Minimise every conformer; prefer MMFF94, fall back to UFF."""
    if ff is None:
        return None, None
    try:
        results = list(ff.MMFFOptimizeMoleculeConfs(mol, maxIters=max_iters, mmffVariant="MMFF94"))
        if results:
            return results, "MMFF94"
    except Exception:
        pass
    try:
        results = list(ff.UFFOptimizeMoleculeConfs(mol, maxIters=max_iters))
        if results:
            return results, "UFF"
    except Exception:
        pass
    return None, None


def smiles_to_conformer_ensemble(
    smiles: str,
    num_conformers: int = 8,
    max_iters: int = 400,
    temperature_k: float = BODY_TEMPERATURE_K,
) -> dict[str, Any]:
    """Generate and minimise a conformer ensemble, then Boltzmann-weight it.

    This is a real molecular-mechanics simulation: ETKDG distance-geometry
    embedding followed by force-field energy minimisation. The populations are
    the equilibrium Boltzmann weights of those minima at ``temperature_k``.
    """
    mods = rdkit_compat()
    Chem = mods.get("Chem")
    rd_dist_geom = mods.get("rdDistGeom")
    ff = mods.get("rdForceFieldHelpers")
    rdmolfiles = mods.get("rdmolfiles")

    if Chem is None:
        return {"available": False, "note": "RDKit is not importable in this environment"}
    if rd_dist_geom is None:
        return {
            "available": False,
            "note": "3D conformer generation (rdDistGeom) is not available in this environment",
        }

    mol = _prepare_mol(Chem, smiles)
    mol = Chem.AddHs(mol)

    try:
        params = rd_dist_geom.ETKDGv3()
        params.randomSeed = 0xC0FFEE
        params.pruneRmsThresh = 0.35
        cids = list(
            rd_dist_geom.EmbedMultipleConfs(mol, numConfs=max(1, int(num_conformers)), params=params)
        )
    except Exception as exc:
        return {
            "available": False,
            "note": f"Conformer embedding failed ({type(exc).__name__})",
        }

    if not cids:
        return {"available": False, "note": "Conformer embedding produced no conformers"}

    energy_results, force_field = _optimise_conformers(ff, mol, max_iters)

    conformers: list[dict[str, Any]] = []
    if energy_results and len(energy_results) == len(cids):
        for cid, (not_converged, energy) in zip(cids, energy_results):
            conformers.append(
                {
                    "index": int(cid),
                    "energy_kcal_per_mol": round(float(energy), 3),
                    "minimised": not bool(not_converged),
                    "molblock": _molblock_from_conformer(Chem, rdmolfiles, mol, int(cid)),
                }
            )
    else:
        force_field = None
        for cid in cids:
            conformers.append(
                {
                    "index": int(cid),
                    "energy_kcal_per_mol": None,
                    "minimised": False,
                    "molblock": _molblock_from_conformer(Chem, rdmolfiles, mol, int(cid)),
                }
            )

    energies = [c["energy_kcal_per_mol"] for c in conformers if c["energy_kcal_per_mol"] is not None]
    populations: list[float] | None = None
    if energies:
        minimised = [c for c in conformers if c["energy_kcal_per_mol"] is not None]
        populations = boltzmann_populations(energies, temperature_k) or None
        relatives = relative_energies(energies)
        for conf, pop, rel in zip(minimised, populations or [], relatives):
            conf["population"] = pop
            conf["relative_energy_kcal_per_mol"] = rel
        lowest = min(
            minimised,
            key=lambda c: c["relative_energy_kcal_per_mol"],
        )
    else:
        lowest = conformers[0]

    if force_field:
        method = (
            f"ETKDGv3 distance-geometry embedding of {len(conformers)} conformers + "
            f"{force_field} energy minimisation. Energies are reported relative to the "
            f"lowest minimum in this run; populations are Boltzmann estimates at "
            f"{temperature_k:.2f} K."
        )
    else:
        method = (
            f"ETKDGv3 distance-geometry embedding of {len(conformers)} conformers. "
            "No force field was available, so energies and populations are not reported."
        )

    return {
        "available": True,
        "method": method,
        "force_field": force_field,
        "temperature_k": float(temperature_k),
        "conformer_count": len(conformers),
        # Headline number: the spread inside this single ensemble. Comparable.
        "energy_span_kcal_per_mol": round(max(energies) - min(energies), 3) if energies else None,
        # Kept for transparency but explicitly not comparable across runs.
        "lowest_absolute_energy_kcal_per_mol": round(min(energies), 3) if energies else None,
        "energy_reference": "relative_to_lowest_conformer_in_this_run",
        "absolute_energies_comparable": False,
        "absolute_energy_note": ABSOLUTE_ENERGY_NOTE,
        "populations_are_estimates": True,
        "population_note": POPULATION_ESTIMATE_NOTE,
        "conformers": conformers,
        "lowest_energy_index": lowest.get("index"),
        "lowest_energy_relative_kcal_per_mol": lowest.get("relative_energy_kcal_per_mol"),
        "lowest_energy_molblock": lowest.get("molblock"),
        "note": None,
    }


def smiles_to_3d_structure(
    smiles: str,
    num_conformers: int = 8,
    max_iters: int = 400,
) -> dict[str, Any]:
    """Single best 3D conformer for viewing, with an honest 2D fallback."""
    mods = rdkit_compat()
    Chem = mods.get("Chem")
    if Chem is None:
        return {
            "available": False,
            "is_3d": False,
            "molblock": None,
            "note": rdkit_unavailable_reason(mods),
        }

    ensemble = smiles_to_conformer_ensemble(smiles, num_conformers=num_conformers, max_iters=max_iters)
    if ensemble.get("available"):
        molblock = ensemble.get("lowest_energy_molblock")
        mol = Chem.MolFromMolBlock(molblock, removeHs=False) if molblock else None
        return {
            "available": molblock is not None,
            "is_3d": bool(mol is not None and _is_conformer_3d(Chem, mol)),
            "molblock": molblock,
            "method": ensemble.get("method"),
            "force_field": ensemble.get("force_field"),
            # Zero by construction (the lowest conformer is the reference). The
            # useful number is the span across the ensemble.
            "relative_energy_kcal_per_mol": ensemble.get("lowest_energy_relative_kcal_per_mol"),
            "atom_count": mol.GetNumAtoms() if mol is not None else None,
            "conformer_count": ensemble.get("conformer_count"),
            "note": None,
        }

    # 3D is unavailable: give a flat 2D depiction rather than nothing at all.
    note = ensemble.get("note") or "3D structure generation is unavailable"
    shape = smiles_to_2d_structure(smiles)
    return {
        "available": shape.get("available", False),
        "is_3d": False,
        "molblock": shape.get("molblock"),
        "method": shape.get("method"),
        "force_field": None,
        "relative_energy_kcal_per_mol": None,
        "atom_count": shape.get("atom_count"),
        "conformer_count": 0,
        "note": note,
    }


def smiles_to_2d_structure(smiles: str) -> dict[str, Any]:
    """Flat 2D depiction built from RDKit's 2D coordinate generation."""
    mods = rdkit_compat()
    Chem = mods.get("Chem")
    rd_depictor = mods.get("rdDepictor")
    rdmolfiles = mods.get("rdmolfiles")
    if Chem is None or rd_depictor is None or rdmolfiles is None:
        note = rdkit_unavailable_reason(mods) if Chem is None else "2D depiction is unavailable"
        return {"available": False, "molblock": None, "note": note}

    mol = _prepare_mol(Chem, smiles)
    try:
        rd_depictor.Compute2DCoords(mol)
    except Exception as exc:
        return {"available": False, "molblock": None, "note": f"2D coordinate generation failed ({type(exc).__name__})"}

    try:
        molblock = rdmolfiles.MolToMolBlock(mol, confId=0)
    except Exception:
        try:
            molblock = rdmolfiles.MolToMolBlock(mol)
        except Exception as exc:
            return {"available": False, "molblock": None, "note": f"Molblock export failed ({type(exc).__name__})"}

    return {
        "available": bool(molblock),
        "molblock": molblock,
        "method": "RDKit 2D coordinate generation (flat depiction)",
        "atom_count": mol.GetNumHeavyAtoms(),
        "note": None,
    }
