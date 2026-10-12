"""Structural analysis of real PDB complexes.

Everything here is computed from the deposited crystallographic coordinates:
no interaction distances are invented. The analysis is a purely geometric
contact classification, which is stated plainly in the returned ``method``
field so the UI never presents it as a docking score or a binding free energy.
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass
from typing import Any

logger = logging.getLogger(__name__)

# Cutoffs in angstrom (heavy-atom..heavy-atom).
HBOND_MAX = 3.5
# The deposited files carry no hydrogens, so the D-H...A angle cannot be measured
# directly. Instead the angle is taken at the donor between its nearest covalent
# heavy-atom neighbour and the acceptor: the D-H bond points away from that
# neighbour, so an angle near 180 deg means the hydrogen can face the acceptor.
HBOND_MIN_ANGLE_DEG = 90.0
SALT_BRIDGE_MAX = 4.0
HYDROPHOBIC_MAX = 4.0
CONTACT_MAX = 4.5

# Interaction types, most specific first. Anything within CONTACT_MAX that is
# none of the first three is reported as a plain "close_contact".
CONTACT_TYPES = ("salt_bridge", "hydrogen_bond", "hydrophobic", "close_contact")

# How each label should be described. Hydrogen bonds and salt bridges rest on
# chemistry assumptions (no hydrogens or protonation states are deposited), so
# they are labelled as heuristics rather than measurements. Hydrophobic and close
# contacts are direct distance measurements between heavy atoms.
CONTACT_LABELS = {
    "hydrogen_bond": "geometric heuristic",
    "salt_bridge": "geometric heuristic",
    "hydrophobic": "measured geometry",
    "close_contact": "measured geometry",
}

# A donor/acceptor pair beyond HBOND_MAX but within CONTACT_MAX with the right
# angle is reported as a near miss, so a zero hydrogen-bond count can be
# explained instead of just displayed.
NEAR_MISS_MAX = CONTACT_MAX
NEAR_MISS_LIMIT = 8

CONTACT_RULES = {
    "hydrogen_bond": (
        "Donor/acceptor heavy-atom pair (N or O on both sides) with D...A <= "
        f"{HBOND_MAX} A AND the antecedent-D...A angle >= {HBOND_MIN_ANGLE_DEG:.0f} deg. "
        "The donor must be a hydrogen-bearing N/O (backbone N, Lys NZ, Arg N, Asn/Gln N, "
        "Trp NE1, Ser/Thr OG, Tyr OH, His N) and the acceptor a lone-pair N/O "
        "(backbone O, Asp/Glu O, Asn/Gln O, His N)."
    ),
    "salt_bridge": (
        f"Oppositely charged groups within {SALT_BRIDGE_MAX} A. Protein side: "
        "Asp OD1/OD2 or Glu OE1/OE2 (anion), or Lys NZ / Arg NE,NH1,NH2 / His ND1,NE2 "
        "(cation). Ligand side: no protonation state is deposited, so a ligand O is "
        "treated as a possible anion and a ligand N as a possible cation."
    ),
    "hydrophobic": f"Carbon-carbon contact only: C...C <= {HYDROPHOBIC_MAX} A.",
    "close_contact": (
        f"Any other heavy-atom pair within {CONTACT_MAX} A. Not a hydrogen bond, "
        "salt bridge or hydrophobic contact."
    ),
}

# Donor / acceptor character by PDB atom name. This is a chemistry-based
# assignment, not a measurement: the deposited coordinates contain no hydrogens
# and no protonation states.
DONOR_ATOM_NAMES = {
    "N", "NZ", "NE", "NH1", "NH2", "NE1", "ND1", "NE2", "OG", "OG1", "OH",
}
ACCEPTOR_ATOM_NAMES = {
    "O", "OXT", "OD1", "OD2", "OE1", "OE2", "SD", "OG", "OG1", "OH", "ND1", "NE2",
}
ANIONIC_ATOM_NAMES = {"OD1", "OD2", "OE1", "OE2"}
CATIONIC_ATOM_NAMES = {"NZ", "NE", "NH1", "NH2", "ND1", "NE2"}

WATER_CODES = {"HOH", "DOD", "WAT"}

# HETATM groups are separated by role so the UI can flag them instead of
# silently hiding them. A crystal entry routinely contains ions, sugars,
# cryoprotectants and cofactors alongside the actual drug.
ION_CODES = {
    "NA", "K", "CL", "BR", "IOD", "MG", "CA", "ZN", "MN", "FE", "FE2", "CU",
    "CU1", "NI", "CO", "CD", "HG", "CS", "RB", "SR", "BA", "AL", "YB", "GD",
    "SM", "AU", "PT", "AG", "MO", "W", "V", "CR", "LI", "CE", "PR", "TB",
}
SUGAR_CODES = {
    "NAG", "NDG", "BMA", "MAN", "FUC", "GAL", "GLA", "GLC", "BGC", "A2G",
    "XYS", "XYP", "SIA", "SLB", "TRE", "SUC", "FRU", "RIB", "RIP", "GCU",
    "IDR", "BDP", "ADA", "FBP", "G6P",
}
COFACTOR_CODES = {
    "NAD", "NAP", "NAI", "FAD", "FMN", "ADP", "ATP", "GDP", "GTP", "AMP",
    "ANP", "SAM", "SAH", "HEM", "HEA", "PLP", "TPP", "COA", "ACO", "UDP",
    "GSH", "BIO", "PQQ", "MTH", "THF", "FES", "SF4", "CLF", "BCL", "MQ7",
}
# Cryoprotectants, precipitants, buffer components and crystallisation additives.
BUFFER_CODES = {
    "SO4", "PO4", "NO3", "ACT", "FMT", "CIT", "EDO", "GOL", "PEG", "PG4",
    "MPD", "DMS", "TRS", "MES", "BME", "DTT", "TLA", "IMD", "EPE", "IPA",
    "ACE", "NH4", "FLC", "ACY", "AZI", "SCN", "PER", "BCT", "CO3", "PGE",
    "1PE", "2PE", "P6G", "SIN", "MLI", "CAC", "CIT", "TAR", "MAL", "LAC",
    "BOG", "DDQ", "HED", "PIN", "OXL", "TAM", "MRD", "SER", "BES", "CXS",
}

# Solvent, ions, sugars, cofactors and buffer components that are never the drug.
NON_LIGAND_CODES = (
    WATER_CODES | ION_CODES | SUGAR_CODES | COFACTOR_CODES | BUFFER_CODES
)

# Residues whose side-chain oxygens can carry a negative charge.
ACIDIC_RESIDUES = {"ASP", "GLU"}
# Residues whose side-chain nitrogens can carry a positive charge.
BASIC_RESIDUES = {"LYS", "ARG", "HIS"}

STANDARD_AMINO_ACIDS = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "MSE", "SEC", "PYL",
}


@dataclass
class PdbAtom:
    record: str
    atom_name: str
    element: str
    res_name: str
    chain: str
    res_seq: str
    x: float
    y: float
    z: float
    raw: str

    @property
    def residue_key(self) -> tuple[str, str, str]:
        return (self.chain, self.res_name, self.res_seq)

    @property
    def residue_label(self) -> str:
        return f"{self.res_name} {self.chain}{self.res_seq}".strip()

    @property
    def is_polar(self) -> bool:
        return self.element in ("N", "O")

    @property
    def is_carbon(self) -> bool:
        return self.element == "C"


def _element_from_name(atom_name: str) -> str:
    name = atom_name.strip()
    if not name:
        return "?"
    # PDB atom names right-pad single-letter elements (e.g. " CA " is carbon
    # alpha, while "CA  " is a calcium ion).
    if len(name) >= 2 and name[0].isalpha() and name[1].islower():
        return name[0].upper()
    if name[0].isdigit():
        return name[1].upper() if len(name) > 1 else "?"
    return name[0].upper()


def _parse_atom_line(line: str) -> PdbAtom | None:
    if len(line) < 54:
        return None
    try:
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
    except ValueError:
        return None

    atom_name = line[12:16]
    element = line[76:78].strip().upper() if len(line) >= 78 else ""
    if not element:
        element = _element_from_name(atom_name)

    return PdbAtom(
        record=line[:6].strip(),
        atom_name=atom_name.strip(),
        element=element,
        res_name=line[17:20].strip().upper(),
        chain=line[21:22].strip(),
        res_seq=line[22:26].strip(),
        x=x,
        y=y,
        z=z,
        raw=line,
    )


def parse_pdb(text: str) -> list[PdbAtom]:
    atoms: list[PdbAtom] = []
    for line in text.splitlines():
        record = line[:6]
        if record == "ATOM  " or record == "HETATM":
            atom = _parse_atom_line(line)
            if atom is not None:
                atoms.append(atom)
    return atoms


def ligand_category(res_name: str) -> str:
    """Classify a HETATM group so a caller can flag it rather than just drop it."""
    code = res_name.strip().upper()
    if code in WATER_CODES:
        return "water"
    if code in ION_CODES:
        return "ion"
    if code in SUGAR_CODES:
        return "sugar"
    if code in COFACTOR_CODES:
        return "cofactor"
    if code in BUFFER_CODES:
        return "buffer"
    if code in STANDARD_AMINO_ACIDS:
        return "amino_acid"
    return "drug_like"


def ligand_candidates(pdb_text: str) -> list[dict[str, Any]]:
    """List every HETATM group in the entry, in a deterministic order.

    Water is omitted (it is solvent, and listing it would swamp the list).
    Everything else is returned with its category and atom count so the caller
    can decide which group is the drug instead of guessing.
    """
    atoms = parse_pdb(pdb_text)
    groups: dict[tuple[str, str, str], list[PdbAtom]] = {}
    for atom in atoms:
        if atom.record != "HETATM":
            continue
        if atom.res_name in WATER_CODES:
            continue
        groups.setdefault(atom.residue_key, []).append(atom)

    candidates: list[dict[str, Any]] = []
    for key, group in groups.items():
        chain, res_name, res_seq = key
        heavy = [a for a in group if a.element != "H"]
        candidates.append(
            {
                "ligand_code": res_name,
                "chain": chain or None,
                "residue_number": res_seq,
                "atom_count": len(group),
                "heavy_atom_count": len(heavy),
                "category": ligand_category(res_name),
                "element_counts": _element_counts(heavy),
                "is_common_ion": res_name in ION_CODES,
                "is_sugar": res_name in SUGAR_CODES,
                "is_cofactor": res_name in COFACTOR_CODES,
                "is_buffer_or_additive": res_name in BUFFER_CODES,
                "is_standard_amino_acid": res_name in STANDARD_AMINO_ACIDS,
            }
        )

    # Most drug-like first (largest heavy-atom count), then alphabetically.
    candidates.sort(
        key=lambda c: (
            c["category"] != "drug_like",
            -c["heavy_atom_count"],
            c["ligand_code"],
            c["chain"] or "",
            c["residue_number"],
        )
    )
    return candidates


def _element_counts(atoms: list[PdbAtom]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for atom in atoms:
        counts[atom.element] = counts.get(atom.element, 0) + 1
    return dict(sorted(counts.items()))


def select_ligand(atoms: list[PdbAtom], preferred_code: str | None = None) -> list[PdbAtom]:
    """Pick the drug-like HETATM residue: the largest non-solvent, non-ionic,
    non-amino-acid HETATM group. ``preferred_code`` wins when present."""
    groups: dict[tuple[str, str, str], list[PdbAtom]] = {}
    for atom in atoms:
        if atom.record != "HETATM":
            continue
        if atom.res_name in WATER_CODES or atom.res_name in NON_LIGAND_CODES:
            continue
        if atom.res_name in STANDARD_AMINO_ACIDS:
            continue
        groups.setdefault(atom.residue_key, []).append(atom)

    if not groups:
        return []

    if preferred_code:
        code = preferred_code.strip().upper()
        for key, group in groups.items():
            if key[1] == code:
                return group

    # Largest heavy-atom count wins (ties broken deterministically).
    best_key = max(sorted(groups), key=lambda k: len(groups[k]))
    return groups[best_key]


def _distance(a: PdbAtom, b: PdbAtom) -> float:
    return math.sqrt((a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2)


def _angle_deg(a: PdbAtom, vertex: PdbAtom, c: PdbAtom) -> float:
    """Angle a-vertex-c in degrees."""
    v1 = (a.x - vertex.x, a.y - vertex.y, a.z - vertex.z)
    v2 = (c.x - vertex.x, c.y - vertex.y, c.z - vertex.z)
    n1 = math.sqrt(v1[0] ** 2 + v1[1] ** 2 + v1[2] ** 2)
    n2 = math.sqrt(v2[0] ** 2 + v2[1] ** 2 + v2[2] ** 2)
    if n1 == 0.0 or n2 == 0.0:
        return 0.0
    cosang = (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]) / (n1 * n2)
    cosang = max(-1.0, min(1.0, cosang))
    return math.degrees(math.acos(cosang))


def _antecedent(donor: PdbAtom, siblings: list[PdbAtom], max_bond: float = 1.9) -> PdbAtom | None:
    """Nearest covalent heavy-atom neighbour of ``donor`` within its own residue.

    Used as the stand-in for the missing hydrogen: the D-H bond points away from
    this atom, so the X-D...A angle carries the direction information.
    """
    best: PdbAtom | None = None
    best_d = max_bond
    for other in siblings:
        if other is donor:
            continue
        d = _distance(donor, other)
        if d <= best_d:
            best_d = d
            best = other
    return best


def _is_donor(atom: PdbAtom, is_ligand: bool) -> bool:
    if atom.element not in ("N", "O"):
        return False
    if is_ligand:
        # No hydrogen or protonation information exists for the ligand either.
        return True
    return atom.atom_name in DONOR_ATOM_NAMES


def _is_acceptor(atom: PdbAtom, is_ligand: bool) -> bool:
    if atom.element not in ("N", "O"):
        return False
    if is_ligand:
        return True
    return atom.atom_name in ACCEPTOR_ATOM_NAMES


def _donor_geometry_ok(
    donor: PdbAtom,
    acceptor: PdbAtom,
    siblings: list[PdbAtom],
) -> bool:
    """Apply the angle criterion; fall back to distance-only if it cannot run."""
    antecedent = _antecedent(donor, siblings)
    if antecedent is None:
        # Single-atom residue, or no neighbour inside bonding distance. The
        # angle is then unknown, so say so via the rule text rather than
        # silently inventing a direction.
        return True
    return _angle_deg(antecedent, donor, acceptor) >= HBOND_MIN_ANGLE_DEG


def _charge_role(atom: PdbAtom, is_ligand: bool) -> str | None:
    """Likely charge of a heavy atom: ``anion``, ``cation`` or ``None``."""
    if is_ligand:
        if atom.element == "O":
            return "anion"
        if atom.element == "N":
            return "cation"
        return None
    if atom.res_name in ACIDIC_RESIDUES and atom.atom_name in ANIONIC_ATOM_NAMES:
        return "anion"
    if atom.res_name in BASIC_RESIDUES and atom.atom_name in CATIONIC_ATOM_NAMES:
        return "cation"
    return None


def _near_miss_hydrogen_bond(
    protein_atom: PdbAtom,
    ligand_atom: PdbAtom,
    dist: float,
    protein_siblings: list[PdbAtom],
    ligand_siblings: list[PdbAtom],
) -> tuple[str, float | None] | None:
    """Right donor/acceptor geometry but too long for the cutoff.

    Returns ``(donor_side, angle_deg)`` or ``None``. This exists so the API can
    explain a zero hydrogen-bond count: e.g. nevirapine N14...O Lys101 in 1VRT is
    4.05 A at 164 degrees, which is ideal geometry outside the 3.5 A cutoff.
    """
    if dist <= HBOND_MAX or dist > NEAR_MISS_MAX:
        return None

    if _is_donor(protein_atom, False) and _is_acceptor(ligand_atom, True):
        antecedent = _antecedent(protein_atom, protein_siblings)
        angle = (
            _angle_deg(antecedent, protein_atom, ligand_atom)
            if antecedent is not None
            else None
        )
        if angle is None or angle >= HBOND_MIN_ANGLE_DEG:
            return "protein", angle

    if _is_donor(ligand_atom, True) and _is_acceptor(protein_atom, False):
        antecedent = _antecedent(ligand_atom, ligand_siblings)
        angle = (
            _angle_deg(antecedent, ligand_atom, protein_atom)
            if antecedent is not None
            else None
        )
        if angle is None or angle >= HBOND_MIN_ANGLE_DEG:
            return "ligand", angle

    return None


def _classify_contact(
    protein_atom: PdbAtom,
    ligand_atom: PdbAtom,
    dist: float,
    protein_siblings: list[PdbAtom],
    ligand_siblings: list[PdbAtom],
) -> str | None:
    """Strongest interaction type for one protein/ligand atom pair.

    Order matters: the most specific rule that matches wins, so a charged pair
    is a salt bridge rather than merely a close contact.
    """
    if dist <= SALT_BRIDGE_MAX:
        protein_role = _charge_role(protein_atom, is_ligand=False)
        ligand_role = _charge_role(ligand_atom, is_ligand=True)
        if protein_role and ligand_role and protein_role != ligand_role:
            return "salt_bridge"

    if dist <= HBOND_MAX:
        if _is_donor(protein_atom, False) and _is_acceptor(ligand_atom, True):
            if _donor_geometry_ok(protein_atom, ligand_atom, protein_siblings):
                return "hydrogen_bond"
        if _is_donor(ligand_atom, True) and _is_acceptor(protein_atom, False):
            if _donor_geometry_ok(ligand_atom, protein_atom, ligand_siblings):
                return "hydrogen_bond"

    if dist <= HYDROPHOBIC_MAX and protein_atom.is_carbon and ligand_atom.is_carbon:
        return "hydrophobic"

    if dist <= CONTACT_MAX:
        return "close_contact"

    return None


# Response-size cap on the per-residue atom_pairs list. atom_pair_counts keeps the
# full tally, so a truncated list is always flagged by atom_pairs_truncated.
ATOM_PAIRS_PER_RESIDUE_CAP = 12

# Ranked so a residue's headline interaction is its most specific one.
_CONTACT_RANK = {name: rank for rank, name in enumerate(CONTACT_TYPES)}


def analyse_complex(
    pdb_text: str,
    ligand_code: str | None = None,
) -> dict[str, Any]:
    """Geometric contact analysis of a protein/ligand complex."""
    atoms = parse_pdb(pdb_text)
    if not atoms:
        return {
            "status": "placeholder",
            "note": "No ATOM/HETATM records found in the structure file",
        }

    protein_atoms = [a for a in atoms if a.record == "ATOM"]
    if not protein_atoms:
        return {
            "status": "placeholder",
            "note": "No protein (ATOM) records found in the structure file",
        }

    ligand_atoms = select_ligand(atoms, ligand_code)
    if not ligand_atoms:
        return {
            "status": "placeholder",
            "protein_atom_count": len(protein_atoms),
            "note": "No drug-like HETATM ligand found in the structure file",
        }

    residues: dict[str, dict[str, Any]] = {}
    # Atom-pair tallies, kept for transparency.
    counts = {name: 0 for name in CONTACT_TYPES}
    # Residue-level tallies: the headline numbers count unique residues.
    residue_hits: dict[str, set[str]] = {name: set() for name in CONTACT_TYPES}
    # Donor/acceptor pairs rejected only because of the distance cutoff.
    near_misses: dict[tuple[str, str], dict[str, Any]] = {}

    protein_by_residue: dict[tuple[str, str, str], list[PdbAtom]] = {}
    for atom in protein_atoms:
        protein_by_residue.setdefault(atom.residue_key, []).append(atom)

    for lig_atom in ligand_atoms:
        for prot_atom in protein_atoms:
            dx = lig_atom.x - prot_atom.x
            dy = lig_atom.y - prot_atom.y
            dz = lig_atom.z - prot_atom.z
            # Cheap squared-distance reject before the sqrt.
            d2 = dx * dx + dy * dy + dz * dz
            if d2 > CONTACT_MAX * CONTACT_MAX:
                continue
            dist = math.sqrt(d2)
            siblings = protein_by_residue[prot_atom.residue_key]
            kind = _classify_contact(prot_atom, lig_atom, dist, siblings, ligand_atoms)

            if kind != "hydrogen_bond":
                miss = _near_miss_hydrogen_bond(
                    prot_atom, lig_atom, dist, siblings, ligand_atoms
                )
                if miss is not None:
                    donor_side, angle = miss
                    miss_key = (prot_atom.residue_label, lig_atom.atom_name)
                    previous = near_misses.get(miss_key)
                    if previous is None or dist < previous["distance_angstrom"]:
                        near_misses[miss_key] = {
                            "protein_residue": prot_atom.residue_label,
                            "protein_atom": prot_atom.atom_name,
                            "ligand_atom": lig_atom.atom_name,
                            "donor_side": donor_side,
                            "distance_angstrom": round(dist, 2),
                            "angle_deg": round(angle, 1) if angle is not None else None,
                            "cutoff_angstrom": HBOND_MAX,
                            "reason": (
                                "Donor/acceptor geometry is acceptable, but the distance is "
                                f"over the {HBOND_MAX} A hydrogen-bond cutoff, so this is not "
                                "counted as a hydrogen bond."
                            ),
                        }

            if kind is None:
                continue
            counts[kind] += 1
            residue_hits[kind].add(prot_atom.residue_label)
            key = prot_atom.residue_label
            entry = residues.get(key)
            pair = {
                "protein_atom": prot_atom.atom_name,
                "ligand_atom": lig_atom.atom_name,
                "distance_angstrom": round(dist, 2),
                "type": kind,
            }
            if entry is None:
                residues[key] = {
                    "residue": prot_atom.res_name,
                    "chain": prot_atom.chain,
                    "residue_number": prot_atom.res_seq,
                    "residue_label": key,
                    "types": [kind],
                    "contact_count": 1,
                    "min_distance_angstrom": round(dist, 2),
                    "atom_pairs": [pair],
                    "atom_pairs_truncated": False,
                }
            else:
                entry["contact_count"] += 1
                if kind not in entry["types"]:
                    entry["types"].append(kind)
                entry["min_distance_angstrom"] = min(
                    entry["min_distance_angstrom"], round(dist, 2)
                )
                if len(entry["atom_pairs"]) < ATOM_PAIRS_PER_RESIDUE_CAP:
                    entry["atom_pairs"].append(pair)
                else:
                    entry["atom_pairs_truncated"] = True

    ordered = sorted(
        residues.values(),
        key=lambda r: (r["min_distance_angstrom"], -r["contact_count"], r["residue_label"]),
    )
    for entry in ordered:
        entry["types"].sort(key=lambda t: _CONTACT_RANK.get(t, 9))

    binding_site = [r["residue_label"] for r in ordered]

    contact_summary = {name: len(residue_hits[name]) for name in CONTACT_TYPES}
    near_miss_list = sorted(
        near_misses.values(),
        key=lambda m: (m["distance_angstrom"], m["protein_residue"], m["ligand_atom"]),
    )[:NEAR_MISS_LIMIT]

    # The per-residue atom_pairs list is capped for response size. Say so, so a
    # truncated list is never mistaken for the complete set of contacts.
    atom_pairs_truncated = any(r.get("atom_pairs_truncated") for r in ordered)
    atom_pairs_truncated_note = (
        f"Per-residue atom_pairs are capped at {ATOM_PAIRS_PER_RESIDUE_CAP} per residue "
        "for response size; atom_pair_counts is the complete tally and is not truncated."
    )

    return {
        "status": "computed",
        "method": (
            "Geometric contact analysis of the deposited crystal coordinates. "
            "Hydrogen bonds are a GEOMETRIC HEURISTIC: a donor/acceptor N/O pair within "
            f"{HBOND_MAX} A with an antecedent-D...A angle >= {HBOND_MIN_ANGLE_DEG:.0f} deg. "
            "The deposited file contains no hydrogens and no protonation states, so this is an "
            "inference from distance, element and atom name, not an observed hydrogen bond. "
            f"Salt bridges are also a heuristic (oppositely charged groups within {SALT_BRIDGE_MAX} A); "
            f"hydrophobic contacts require C...C within {HYDROPHOBIC_MAX} A and are measured directly; "
            f"every other heavy-atom pair within {CONTACT_MAX} A is a plain close contact. "
            "Not a docking score or a binding free energy."
        ),
        "contact_rules": CONTACT_RULES,
        "contact_labels": CONTACT_LABELS,
        "near_miss_hydrogen_bonds": near_miss_list,
        "contact_summary": contact_summary,
        "contact_summary_basis": "unique_residues",
        "atom_pair_counts": counts,
        "atom_pairs_truncated": atom_pairs_truncated,
        "atom_pairs_truncated_note": atom_pairs_truncated_note,
        "ligand_code": ligand_atoms[0].res_name,
        "ligand_atom_count": len(ligand_atoms),
        "protein_atom_count": len(protein_atoms),
        "binding_site_residue_count": len(ordered),
        "interacting_residues": ordered,
        "binding_site_residues": binding_site,
        # State the cap when it actually affected the payload; otherwise leave the
        # note free for status commentary.
        "note": atom_pairs_truncated_note if atom_pairs_truncated else None,
    }


def build_viewer_pdb(pdb_text: str, ligand_code: str | None = None) -> str:
    """Trim the deposited PDB down to protein + the target ligand.

    Waters, ions and buffer components are dropped so the 3D viewer payload is
    small and the ligand is unambiguous.
    """
    atoms = parse_pdb(pdb_text)
    if not atoms:
        return pdb_text

    ligand_atoms = select_ligand(atoms, ligand_code)
    keep_keys = {a.residue_key for a in ligand_atoms}

    out: list[str] = []
    for atom in atoms:
        if atom.record == "ATOM":
            out.append(atom.raw)
        elif atom.record == "HETATM" and atom.residue_key in keep_keys:
            out.append(atom.raw)

    if not out:
        return pdb_text
    out.append("END")
    return "\n".join(out) + "\n"


def _residues_by_kind(analysis: dict[str, Any], kind: str, limit: int = 10) -> list[dict[str, Any]]:
    hits = [r for r in analysis.get("interacting_residues", []) if kind in r.get("types", [])]
    hits.sort(key=lambda r: (r.get("min_distance_angstrom", 999.0), r.get("residue_label", "")))
    return hits[:limit]


def _pair_text(residue: dict[str, Any], kind: str) -> str:
    pairs = [p for p in residue.get("atom_pairs", []) if p.get("type") == kind]
    if pairs:
        best = min(pairs, key=lambda p: p.get("distance_angstrom", 999.0))
        return f"{residue['residue_label']} ({best['protein_atom']}..{best['ligand_atom']})"
    return str(residue.get("residue_label", ""))


def build_mechanism(
    analysis: dict[str, Any],
    meta: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Turn the computed contacts into an ordered, auditable mechanism of action.

    Contact counts and distances come from ``analysis`` (computed from the
    deposited coordinates). Anything drawn from the literature — the target
    class, the drug class, which residues are catalytic — is carried in ``meta``
    and labelled ``source: literature`` so the two are never conflated.
    """
    meta = meta or {}
    if analysis.get("status") != "computed":
        return {
            "status": analysis.get("status", "placeholder"),
            "steps": [],
            "catalytic_residues": [],
            "note": analysis.get("note") or "No computed contacts available for a mechanism.",
        }

    counts = analysis.get("contact_summary", {})
    hbonds = _residues_by_kind(analysis, "hydrogen_bond")
    salts = _residues_by_kind(analysis, "salt_bridge")
    hydrophobics = _residues_by_kind(analysis, "hydrophobic")

    # Which curated catalytic residues does the computed contact set actually hit?
    catalytic: list[dict[str, Any]] = []
    for entry in meta.get("catalytic_residues", []) or []:
        match = str(entry.get("match", "")).strip()
        hit = next(
            (r for r in analysis.get("interacting_residues", []) if r.get("residue_label") == match),
            None,
        )
        catalytic.append(
            {
                "label": entry.get("label", match),
                "role": entry.get("role"),
                "match": match,
                "in_contact": hit is not None,
                "min_distance_angstrom": hit.get("min_distance_angstrom") if hit else None,
                "types": hit.get("types", []) if hit else [],
            }
        )

    occluded = [c for c in catalytic if c["in_contact"]]
    site_type = str(meta.get("binding_site_type") or "catalytic").lower()
    steps: list[dict[str, Any]] = []

    steps.append(
        {
            "step": "target",
            "title": "Identify the viral target and the drug",
            "status": "literature",
            "detail": meta.get("mechanism_summary")
            or "Target and drug identity taken from the curated complex entry.",
            "value": meta.get("target_class"),
        }
    )

    steps.append(
        {
            "step": "occupy",
            "title": "The drug settles into the binding pocket",
            "status": "computed",
            "detail": (
                f"The ligand ({analysis.get('ligand_code')}, {analysis.get('ligand_atom_count')} atoms) makes "
                f"close contact with {analysis.get('binding_site_residue_count')} protein residues. "
                "All distances below are measured from the deposited crystal coordinates, and every "
                "count is a number of residues, not of atom pairs."
            ),
            "value": f"{analysis.get('binding_site_residue_count')} residues in contact",
        }
    )

    if hbonds:
        steps.append(
            {
                "step": "anchor",
                "title": "Hydrogen bonds anchor the drug",
                "status": "computed",
                "detail": (
                    f"{counts.get('hydrogen_bond', 0)} residue(s) form a hydrogen bond. Rule: a "
                    "donor/acceptor N or O pair with D...A within 3.5 A and an antecedent-D...A "
                    "angle check (the deposited file contains no hydrogens). "
                    + "; ".join(_pair_text(r, "hydrogen_bond") for r in hbonds[:5])
                ),
                "value": f"{counts.get('hydrogen_bond', 0)} residues",
            }
        )

    if salts:
        steps.append(
            {
                "step": "charge",
                "title": "A charged pair locks the drug down",
                "status": "computed",
                "detail": (
                    f"{counts.get('salt_bridge', 0)} residue(s) form a salt bridge. Rule: oppositely "
                    "charged groups within 4.0 A (Asp/Glu carboxylate against a ligand N, or Lys/Arg/His "
                    "N against a ligand O). These are the strongest interactions in the complex: "
                    + "; ".join(_pair_text(r, "salt_bridge") for r in salts[:4])
                ),
                "value": f"{counts.get('salt_bridge', 0)} residues",
            }
        )

    if hydrophobics:
        steps.append(
            {
                "step": "hydrophobic",
                "title": "Hydrophobic packing excludes water",
                "status": "computed",
                "detail": (
                    f"{counts.get('hydrophobic', 0)} residue(s) make carbon-carbon contacts within "
                    "4.0 A, which buries non-polar surface and is what makes the association "
                    "favourable in water."
                ),
                "value": ", ".join(r["residue_label"] for r in hydrophobics[:6]),
            }
        )

    if site_type == "allosteric":
        # An allosteric inhibitor is *supposed* to miss the catalytic residues.
        # Reporting that as a failure to block would be wrong, so say what the
        # geometry actually shows: the curated pocket is occupied and the
        # catalytic residues are untouched.
        pocket = [
            {
                "label": entry.get("label"),
                "match": entry.get("match"),
                "in_contact": any(
                    r.get("residue_label") == str(entry.get("match", "")).strip()
                    for r in analysis.get("interacting_residues", [])
                ),
            }
            for entry in (meta.get("binding_pocket_residues") or [])
        ]
        pocket_hits = [p for p in pocket if p["in_contact"]]
        missing = ", ".join(c["label"] for c in catalytic if not c["in_contact"])
        detail = (
            "This inhibitor is allosteric, so it is expected to miss the catalytic residues. "
            f"Occupied pocket residues: {', '.join(p['label'] for p in pocket_hits) or 'none detected'}. "
            "Catalytic residues NOT in contact: "
            f"{missing or 'none (unexpected for an allosteric entry)'}. "
            "The blockade below is conformational, not a direct occupation of the active site."
        )
        steps.append(
            {
                "step": "blockade",
                "title": "The allosteric site is occupied, not the active site",
                "status": "computed",
                "detail": detail,
                "value": "allosteric",
            }
        )
    elif catalytic:
        if occluded:
            detail = (
                "The literature-defined catalytic residues below are all present in the computed contact set, "
                "so the drug is physically sitting on the catalytic machinery: "
                + ", ".join(
                    f"{c['label']} ({c['min_distance_angstrom']} A)" for c in occluded
                )
                + "."
            )
            verdict = "occluded"
        else:
            missing = ", ".join(c["label"] for c in catalytic if not c["in_contact"])
            detail = (
                "No curated catalytic residue appears within the 4.5 A contact cutoff in this structure. "
                f"Not in contact: {missing}. The drug binds nearby but this analysis does not show direct "
                "occupation of the catalytic centre."
            )
            verdict = "adjacent"
        steps.append(
            {
                "step": "blockade",
                "title": "The active site is blocked",
                "status": "computed",
                "detail": detail,
                "value": verdict,
            }
        )

    steps.append(
        {
            "step": "consequence",
            "title": "Consequence for the virus",
            "status": "literature",
            "detail": meta.get("consequence")
            or (
                "With the active site occupied the enzyme cannot process its substrate, so the "
                "downstream step of the viral life cycle stalls. This follows from the blockade above, "
                "not from a free-energy calculation."
            ),
            "value": None,
        }
    )

    return {
        "status": "computed",
        "target_class": meta.get("target_class"),
        "drug_class": meta.get("drug_class"),
        "target_origin": meta.get("target_origin"),
        "binding_site_type": site_type,
        "mechanism_summary": meta.get("mechanism_summary"),
        "catalytic_residues": catalytic,
        "catalytic_residues_occluded": len(occluded),
        "catalytic_residues_curated": len(catalytic),
        "steps": steps,
        "sources": meta.get("sources", []),
        "note": (
            "Computed steps come from the deposited crystal coordinates. Literature steps are labelled "
            "and are not derived from this structure."
        ),
    }
