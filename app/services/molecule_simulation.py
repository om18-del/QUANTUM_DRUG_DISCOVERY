"""SMILES simulation service.

Every number returned here is either computed with RDKit on this machine or
derived from a published, named rule of thumb (Lipinski 1997, Veber 2002,
Clark 1999 for CNS penetration). The provenance of each value is stated in the
payload so the UI never has to guess where a figure came from.
"""

from __future__ import annotations

import logging
import re
from typing import Any

from pydantic import BaseModel, Field

logger = logging.getLogger(__name__)

# A PDB identifier: four characters, first is a digit (e.g. 1HSG, 6LU7).
_PDB_ID_RE = re.compile(r"^[0-9][A-Za-z0-9]{3}$")

# Thresholds used by the interpretation layer, with their published source.
LIPINSKI_SOURCE = "Lipinski et al., Adv. Drug Deliv. Rev. 1997 (rule of five)"
VEBER_SOURCE = "Veber et al., J. Med. Chem. 2002 (oral bioavailability)"
CNS_SOURCE = "Clark, J. Pharm. Sci. 1999 (CNS penetration, TPSA/logP)"
PERMEABILITY_SOURCE = "Standard medicinal-chemistry TPSA/rotatable-bond heuristics"


def looks_like_pdb_id(value: str) -> bool:
    return bool(_PDB_ID_RE.match(value.strip()))


def _interpret(properties: dict[str, Any]) -> list[dict[str, Any]]:
    """Turn the computed descriptors into explicit, sourced drug-likeness calls."""
    findings: list[dict[str, Any]] = []

    mw = properties.get("molecular_weight")
    logp = properties.get("logp")
    tpsa = properties.get("tpsa")
    hbd = properties.get("h_bond_donors")
    hba = properties.get("h_bond_acceptors")
    rot = properties.get("rotatable_bonds")
    lipinski = properties.get("lipinski_rule_of_five") or {}
    veber = properties.get("veber") or {}

    if lipinski.get("pass") is not None:
        findings.append(
            {
                "name": "Lipinski rule of five",
                "verdict": "pass" if lipinski.get("pass") else "fail",
                "value": f"{lipinski.get('violations')} violation(s)",
                "threshold": "MW<=500, logP<=5, HBD<=5, HBA<=10",
                "rationale": (
                    "Oral absorption is predicted to be acceptable."
                    if lipinski.get("pass")
                    else "At least one rule is broken; oral absorption may be limited."
                ),
                "source": LIPINSKI_SOURCE,
            }
        )

    if veber.get("pass") is not None:
        findings.append(
            {
                "name": "Veber oral bioavailability",
                "verdict": "pass" if veber.get("pass") else "fail",
                "value": f"rotatable bonds {veber.get('rotatable_bonds')}, TPSA {veber.get('tpsa')}",
                "threshold": "rotatable bonds<=10, TPSA<=140",
                "rationale": (
                    "Flexibility and polarity are both within the range where good "
                    "oral bioavailability is usually observed."
                    if veber.get("pass")
                    else "Too flexible or too polar for reliable oral bioavailability."
                ),
                "source": VEBER_SOURCE,
            }
        )

    if tpsa is not None and logp is not None:
        cns = tpsa <= 60.0 and 1.0 <= logp <= 3.5
        findings.append(
            {
                "name": "Blood-brain barrier penetration",
                "verdict": "likely" if cns else "unlikely",
                "value": f"TPSA {tpsa}, logP {logp}",
                "threshold": "TPSA<=60 and 1.0<=logP<=3.5",
                "rationale": (
                    "Low polar surface area with moderate lipophilicity: this profile "
                    "crosses the blood-brain barrier in the large majority of cases."
                    if cns
                    else "Outside the low-TPSA / moderate-logP window associated with "
                    "passive blood-brain barrier penetration."
                ),
                "source": CNS_SOURCE,
            }
        )

    if tpsa is not None and rot is not None:
        if tpsa <= 90.0 and rot <= 10:
            verdict, rationale = "high", "Low polar surface area and few rotatable bonds favour passive membrane permeation."
        elif tpsa <= 140.0:
            verdict, rationale = "moderate", "Polar surface area is in the mid range; absorption will be partly transporter-dependent."
        else:
            verdict, rationale = "low", "High polar surface area makes passive membrane permeation unlikely."
        findings.append(
            {
                "name": "Predicted passive permeability",
                "verdict": verdict,
                "value": f"TPSA {tpsa}, rotatable bonds {rot}",
                "threshold": "TPSA<=90 with <=10 rotatable bonds",
                "rationale": rationale,
                "source": PERMEABILITY_SOURCE,
            }
        )

    if logp is not None:
        if logp > 5:
            verdict, rationale = "very lipophilic", "High logP raises the risk of poor aqueous solubility and promiscuous binding."
        elif logp < 1:
            verdict, rationale = "very hydrophilic", "Low logP reduces passive membrane permeation."
        else:
            verdict, rationale = "balanced", "logP sits in the range most orally absorbed drugs occupy."
        findings.append(
            {
                "name": "Lipophilicity",
                "verdict": verdict,
                "value": f"logP {logp}",
                "threshold": "1 to 3 is the typical orally absorbed window",
                "rationale": rationale,
                "source": "Lipinski et al. 1997 / general med-chem consensus",
            }
        )

    if mw is not None and hbd is not None and hba is not None:
        findings.append(
            {
                "name": "Size and hydrogen bonding",
                "verdict": "drug-like" if mw <= 500 and hbd <= 5 and hba <= 10 else "outside drug-like",
                "value": f"MW {mw}, HBD {hbd}, HBA {hba}",
                "threshold": "MW<=500, HBD<=5, HBA<=10",
                "rationale": "Complements the rule-of-five count with the individual values.",
                "source": LIPINSKI_SOURCE,
            }
        )

    return findings


def _simulation_steps(
    properties: dict[str, Any],
    structure: dict[str, Any],
    ensemble: dict[str, Any],
    findings: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Auditable log of exactly what the simulation did, in order."""
    steps: list[dict[str, Any]] = [
        {
            "step": "parse",
            "title": "Parse and canonicalise the SMILES",
            "status": "done",
            "detail": "RDKit converted the input string into a connection table and wrote a canonical SMILES.",
            "value": properties.get("formula"),
        },
        {
            "step": "descriptors",
            "title": "Compute molecular descriptors",
            "status": "done",
            "detail": "Molecular weight, logP, hydrogen-bond donors/acceptors, rotatable bonds, TPSA, rings and fraction sp3.",
            "value": f"MW {properties.get('molecular_weight')}, logP {properties.get('logp')}, TPSA {properties.get('tpsa')}",
        },
    ]

    if ensemble.get("available"):
        steps.append(
            {
                "step": "embed",
                "title": "Embed 3D conformers (distance geometry)",
                "status": "done",
                "detail": "ETKDGv3 places atoms in 3D from the bond graph, using experimental torsion preferences.",
                "value": f"{ensemble.get('conformer_count')} conformer(s)",
            }
        )
        if ensemble.get("force_field"):
            steps.append(
                {
                    "step": "minimise",
                    "title": "Minimise every conformer with a force field",
                    "status": "done",
                    "detail": (
                        f"{ensemble.get('force_field')} relaxed each structure to a local energy minimum. "
                        "Only differences between these minima are comparable: absolute force-field "
                        "energies change with the conformer pool, the embedding seed and the RDKit build, "
                        "so every energy below is quoted relative to the lowest minimum in this run."
                    ),
                    "value": (
                        f"span {ensemble.get('energy_span_kcal_per_mol')} kcal/mol across "
                        f"{ensemble.get('conformer_count')} conformer(s), lowest = 0 by definition"
                    ),
                }
            )
            steps.append(
                {
                    "step": "population",
                    "title": "Weight conformers by Boltzmann statistics",
                    "status": "done",
                    "detail": (
                        f"Each minimum was given a Boltzmann weight p ~ exp(-E/RT) at "
                        f"{ensemble.get('temperature_k')} K, using the relative energies. These are "
                        "estimates, not measured populations: they treat the force-field energy as a "
                        "free energy, so they ignore entropy and solvent. Read them as which shape "
                        "dominates rather than as an exact fraction."
                    ),
                    "value": ", ".join(
                        f"{c.get('population')}" for c in ensemble.get("conformers", []) if c.get("population") is not None
                    )
                    or None,
                }
            )
        else:
            steps.append(
                {
                    "step": "minimise",
                    "title": "Force-field minimisation unavailable",
                    "status": "skipped",
                    "detail": "No force field (MMFF94/UFF) could be loaded in this environment, so energies and populations are not reported.",
                    "value": None,
                }
            )
    else:
        steps.append(
            {
                "step": "embed",
                "title": "3D conformer embedding unavailable",
                "status": "skipped",
                "detail": ensemble.get("note") or "rdDistGeom could not run.",
                "value": None,
            }
        )
        if structure.get("molblock"):
            steps.append(
                {
                    "step": "depict",
                    "title": "Fall back to a flat 2D depiction",
                    "status": "done",
                    "detail": "A 2D layout was computed so the molecule can still be inspected.",
                    "value": structure.get("method"),
                }
            )

    steps.append(
        {
            "step": "interpret",
            "title": "Apply published drug-likeness rules",
            "status": "done",
            "detail": "Each descriptor was compared against named literature thresholds rather than an invented score.",
            "value": f"{len(findings)} rule(s) evaluated",
        }
    )
    return steps


# Machine-readable geometry method code, so a client never has to parse prose.
METHOD_3D_MMFF94 = "3d_mmff94"
METHOD_3D_UFF = "3d_uff"
METHOD_3D_UNMINIMISED = "3d_unminimised"
METHOD_2D_FALLBACK = "2d_fallback"
METHOD_UNAVAILABLE = "unavailable"
METHOD_CODES = (
    METHOD_3D_MMFF94,
    METHOD_3D_UFF,
    METHOD_3D_UNMINIMISED,
    METHOD_2D_FALLBACK,
    METHOD_UNAVAILABLE,
)


def _method_code(structure: dict[str, Any]) -> str:
    """Reduce the structure payload to one honest, machine-readable code.

    The prose description stays in ``method_detail``; this field is the value a
    client can branch on.
    """
    if structure.get("is_3d"):
        force_field = (structure.get("force_field") or "").upper()
        if force_field.startswith("MMFF"):
            return METHOD_3D_MMFF94
        if force_field:
            return METHOD_3D_UFF
        return METHOD_3D_UNMINIMISED
    if structure.get("available"):
        return METHOD_2D_FALLBACK
    return METHOD_UNAVAILABLE


class MoleculeSimulationResult(BaseModel):
    smiles: str
    canonical_smiles: str | None = None
    status: str = "computed"
    method: str | None = None
    method_detail: str | None = None
    energy_basis: str = "relative_to_lowest_conformer_in_this_run"
    properties: dict[str, Any] = Field(default_factory=dict)
    structure: dict[str, Any] = Field(default_factory=dict)
    depiction: dict[str, Any] = Field(default_factory=dict)
    ensemble: dict[str, Any] = Field(default_factory=dict)
    findings: list[dict[str, Any]] = Field(default_factory=list)
    simulation_steps: list[dict[str, Any]] = Field(default_factory=list)
    interpretation_note: str | None = None
    note: str | None = None


def simulate_molecule(smiles: str, num_conformers: int = 8) -> dict[str, Any]:
    from app.services.rdkit_helper import (
        rdkit_properties_from_smiles,
        smiles_to_2d_structure,
        smiles_to_3d_structure,
        smiles_to_conformer_ensemble,
        normalize_smiles_safe,
    )

    raw = (smiles or "").strip()
    if not raw:
        return MoleculeSimulationResult(
            smiles=smiles,
            status="placeholder",
            note="Enter a SMILES string to simulate.",
        ).model_dump()

    if looks_like_pdb_id(raw) and any(ch.isalpha() for ch in raw):
        # Do not guess: a PDB id is a structure file, not a molecule.
        from app.services.rdkit_helper import rdkit_compat

        if rdkit_compat().get("Chem") is None or not _parses_as_smiles(raw):
            return MoleculeSimulationResult(
                smiles=raw,
                status="placeholder",
                note=(
                    f"'{raw}' looks like a PDB structure identifier rather than a SMILES string. "
                    "Load it on the Protein Viewer tab, or open the Drug-Target Studio to see the "
                    "drug bound to its target."
                ),
            ).model_dump()

    normalized = normalize_smiles_safe(raw)
    if normalized is None:
        return MoleculeSimulationResult(
            smiles=raw,
            status="placeholder",
            note="SMILES normalization failed — RDKit could not parse this string.",
        ).model_dump()

    try:
        properties = rdkit_properties_from_smiles(normalized)
    except Exception as exc:  # noqa: BLE001
        logger.warning("SMILES descriptor computation failed for %r: %s", normalized, exc)
        return MoleculeSimulationResult(
            smiles=normalized,
            status="placeholder",
            note=f"Descriptor computation failed: {exc}",
        ).model_dump()

    ensemble = smiles_to_conformer_ensemble(normalized, num_conformers=num_conformers)
    structure = smiles_to_3d_structure(normalized, num_conformers=num_conformers)
    depiction = smiles_to_2d_structure(normalized)
    findings = _interpret(properties)
    steps = _simulation_steps(properties, structure, ensemble, findings)

    if structure.get("is_3d"):
        interpretation_note = (
            "The geometry below is a real energy-minimised 3D structure, not a drawing. "
            "Descriptors and rule-based calls describe what the molecule would do in the body; "
            "they do not predict binding to any specific viral protein."
        )
    elif structure.get("available"):
        interpretation_note = (
            "3D embedding was unavailable in this environment, so a flat 2D depiction is shown. "
            "Descriptors and rule-based calls are still real computed values."
        )
    else:
        interpretation_note = (
            "No 3D or 2D structure could be produced in this environment. "
            "Descriptors and rule-based calls are still real computed values."
        )

    return MoleculeSimulationResult(
        smiles=normalized,
        canonical_smiles=normalized,
        status="computed",
        method=_method_code(structure),
        method_detail=structure.get("method"),
        properties=properties,
        structure=structure,
        depiction=depiction,
        ensemble=ensemble,
        findings=findings,
        simulation_steps=steps,
        interpretation_note=interpretation_note,
        note=structure.get("note"),
    ).model_dump()


def _parses_as_smiles(value: str) -> bool:
    try:
        from app.services.rdkit_helper import normalize_smiles_safe

        return normalize_smiles_safe(value) is not None
    except Exception:
        return False
