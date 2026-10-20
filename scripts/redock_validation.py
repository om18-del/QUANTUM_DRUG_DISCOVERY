#!/usr/bin/env python
"""Redock known crystal complexes with AutoDock Vina and measure pose RMSD.

Run this on Linux or in Google Colab. The AutoDock Vina Python bindings and
meeko publish wheels for Linux/macOS only, so there is deliberately no Windows
path here -- exactly like ``precompute_hamiltonians.py`` and PySCF.

Install (the notebook installs the same pins):

    pip install "requests>=2.31" "numpy>=1.26" "rdkit" "vina==1.2.7" "meeko==0.8.0"

Usage:

    python scripts/redock_validation.py OUT.json
    python scripts/redock_validation.py OUT.json --targets 1VRT,2HU4,1HSG
    python scripts/redock_validation.py OUT.json --exhaustiveness 16 --seed 7

What it does, per PDB ID:

1. Download the entry from RCSB (files.rcsb.org) and its metadata
   (data.rcsb.org). The ligand is the single non-polymer entity with at least
   ``MIN_LIGAND_HEAVY_ATOMS`` heavy atoms according to the wwPDB chemical
   component definition; if that rule leaves zero or more than one candidate the
   target is recorded as **failed** with the candidate list. No residue code is
   ever guessed.
2. Keep the protein chains that have an atom within
   ``CHAIN_INTERFACE_CUTOFF_ANGSTROM`` of the crystal ligand, drop waters and
   every other heteroatom, and convert the receptor to PDBQT with meeko. The
   ligand is built from the deposited crystal coordinates with bond orders taken
   from the RCSB chemical component definition (never inferred from geometry).
3. Dock with a box centred on the crystal ligand centroid, 22 x 22 x 22 A, one
   fixed ``--exhaustiveness`` and ``--seed`` for every target, ``num_modes = 9``.
4. Measure the heavy-atom, symmetry-aware RMSD of every docked mode against the
   deposited crystal pose, with **no re-alignment**, via
   ``rdkit.Chem.rdMolAlign.CalcRMS`` (it searches symmetry-equivalent atom
   correspondences and reports the smallest RMSD; it does not move the probe).

Honest scope, repeated in the output file itself:

* A Vina score is a **prediction, not measured binding**. It is not a Kd, Ki or
  IC50 and it is not comparable with one.
* Redocking with a box around the **true** crystal site is an **easy** test. It
  says the search can find a pose that is already inside the box; it says nothing
  about affinity prediction, ranking, blind docking or virtual screening.
* Parameters are never retuned per target and no run is repeated until it
  passes. A target that cannot be prepared or docked is written out as
  ``status: "failed"`` with a ``failure_reason`` -- it is never omitted or tuned
  into a pass.

Exit status: 0 when the file was written (even with failed targets -- the
failures are the result), 2 for a usage error, 3 when the file was written but
every target failed, 1 when the file could not be written.

Standard library plus ``requests``, ``numpy``, ``rdkit``, ``vina`` and ``meeko``.
Every heavy import happens inside a function, so this module imports on a host
where none of them are installed (the import-safety test relies on that).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import re
import sys
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Sequence

SCHEMA_VERSION = "1.0"

# The three complexes from SPEC.md P3-3. This is only the default target list --
# nothing about the expected ligand codes is hardcoded anywhere in this script.
DEFAULT_TARGETS: tuple[str, ...] = ("1VRT", "2HU4", "1HSG")
EXPECTED_TARGETS: tuple[str, ...] = DEFAULT_TARGETS

PASS_THRESHOLD_ANGSTROM = 2.0
PASS_RULE = "top-ranked pose RMSD < 2.0 A"
BOX_SIZE_ANGSTROM: tuple[float, float, float] = (22.0, 22.0, 22.0)
NUM_MODES = 9
ENERGY_RANGE_KCAL_PER_MOL = 3.0
DEFAULT_EXHAUSTIVENESS = 16
DEFAULT_SEED = 7
# Vina's own default: use every available core. Recorded rather than left
# implicit, because bit-for-bit reproducibility needs the same CPU count as well
# as the same seed (AutoDock Vina FAQ).
VINA_CPU = 0
VINA_CPU_NOTE = (
    "cpu=0 means every available core, Vina's default. The seed makes the search "
    "reproducible given the same inputs and parameters, but bit-for-bit reproduction "
    "also needs the same CPU count on the same Vina version (AutoDock Vina FAQ)."
)

# One fixed setting per run, never per target (SPEC.md P3-3, and the validator
# checks that every target in a result file carries the same parameters).
RMSD_DECIMALS = 3
SCORE_DECIMALS = 3

# Receptor chain selection and ligand identification, both deterministic.
CHAIN_INTERFACE_CUTOFF_ANGSTROM = 5.0
MIN_LIGAND_HEAVY_ATOMS = 10
WATER_COMPONENT_CODES = ("HOH", "DOD", "WAT")

# A non-polymer entity is not automatically the bound drug: 1VRT also deposits a
# magnesium ion, and buffers/ions/cryoprotectants are common. The candidates are
# therefore filtered by the wwPDB chemical component definition's own heavy-atom
# count (Mg 1, NVP 20, MK1 45, G39 20), and a run that leaves 0 or >1 candidate
# fails rather than picking one.
LIGAND_IDENTIFICATION_RULE = (
    f"the single non-polymer entity of the entry whose wwPDB chemical component "
    f"definition (rcsb_chem_comp_info.atom_count_heavy) reports at least "
    f"{MIN_LIGAND_HEAVY_ATOMS} heavy atoms, excluding waters "
    f"({', '.join(WATER_COMPONENT_CODES)}); zero or several candidates fail the target"
)

RECEPTOR_PREPARATION = (
    f"keep the protein chains (ATOM records) with an atom within "
    f"{CHAIN_INTERFACE_CUTOFF_ANGSTROM:g} A of any crystal ligand atom; drop every "
    "HETATM record (waters, ions, the ligand itself and any modified residue "
    "deposited as HETATM), so the receptor is the deposited protein without the "
    "ligand; convert to PDBQT with meeko (Polymer + PDBQTWriterLegacy), which adds "
    "the missing hydrogens. No pH-based protonation assignment is performed: the "
    "deposited protonation states are used as they are."
)

LIGAND_PREPARATION = (
    "heavy atoms and their crystal coordinates come from the entry's HETATM records "
    "for the first deposited instance of the ligand residue; bond orders, formal "
    "charges and stereochemistry come from the RCSB chemical component definition "
    "(rcsb_chem_comp_descriptor.SMILES_stereo) applied with "
    "rdkit.Chem.AllChem.AssignBondOrdersFromTemplate. Bond orders are never guessed "
    "from geometry and the ligand is never re-embedded."
)

RMSD_METHOD = (
    "heavy-atom, symmetry-aware, no re-alignment: rdkit.Chem.rdMolAlign.CalcRMS(pose, "
    "crystal) with symmetrizeConjugatedTerminalGroups=True over the heavy-atom "
    "molecules. CalcRMS searches the symmetry-equivalent atom correspondences and "
    "reports the smallest root-mean-square distance for the coordinates as they "
    "stand -- it does not fit a rigid-body transform and does not move either "
    "molecule (that is what GetBestRMS/AlignMol would do). Both molecules are in "
    "the receptor frame, so no alignment is wanted."
)

BOX_NOTE = (
    "Standard but easy redocking setup: the box is centred on the centroid of the "
    "deposited crystal ligand (the true binding site) and is a 22 A cube, large "
    "enough to contain the known pose. Finding that pose here says the search works "
    "when the answer is already inside the box. It is not blind docking, it is not "
    "virtual screening, and it says nothing about finding or ranking new molecules."
)

ATOM_MAPPING_METHOD = (
    "the docked pose's heavy atoms are matched to the ligand molecule's heavy atoms "
    "in file order (Vina rewrites the atoms of the ligand it was given, and meeko "
    "writes the molecule's atoms in order); the correspondence is verified for every "
    "atom by comparing each PDBQT AutoDock type with the element of the molecule "
    "atom it is matched to, and the whole path is checked against the crystal pose "
    "itself (crystal_pose_self_check_rmsd_angstrom) before docking"
)

DISCLAIMERS: tuple[str, ...] = (
    "Scores are docking predictions, not measured binding: a Vina score is a "
    "force-field estimate for a pose in this rigid receptor, not an experimental "
    "Kd, Ki or IC50, and it must not be compared with one.",
    "Redocking is an easy test: the search box is centred on the true crystal "
    "ligand site and the ligand chemistry is taken from the crystal, so a low RMSD "
    "shows only that the search can rediscover a pose that is already inside the "
    "box. It is not evidence of virtual-screening performance, affinity prediction "
    "or lead-finding ability.",
    f"These results cover {len(EXPECTED_TARGETS)} complexes only "
    f"({', '.join(EXPECTED_TARGETS)}) and must not be generalised into a claim "
    "about the method, the force field or the tool.",
    "One fixed setting was used for every target (box, exhaustiveness, seed and "
    "num_modes as recorded). No parameter was retuned per target and no run was "
    "repeated until it passed; a target that failed is recorded as failed.",
    "RMSD is heavy-atom, symmetry-aware and computed without re-alignment against "
    "the deposited crystal pose; the crystal pose is a model of the experiment, not "
    "the experiment.",
)

PIP_INSTALL_LINE = (
    'pip install "requests>=2.31" "numpy>=1.26" "rdkit" "vina==1.2.7" "meeko==0.8.0"'
)

# Pins that are enforced (a mismatch is loud rather than silent), because the
# preparation and pose-reading code below is written against these APIs. The
# remaining tool versions are recorded but not enforced.
PINNED_VERSIONS = {"vina": "1.2.7", "meeko": "0.8.0"}

RCSB_DATA_API = "https://data.rcsb.org/rest/v1/core"
RCSB_FILES = "https://files.rcsb.org/download"
USER_AGENT = "qdde-offline-redock-validation/1.0"
REQUEST_TIMEOUT_SECONDS = 30.0
REQUEST_RETRIES = 3
REQUEST_BACKOFF_SECONDS = 2.0

MAX_SUMMARY_REASON_CHARACTERS = 60

# PDBQT atom types are an AutoDock vocabulary, not element symbols. This table
# maps the types that can appear in a ligand PDBQT onto elements and is used to
# verify that a pose atom really is the molecule atom it was matched to. An
# unknown type makes the target fail instead of producing an RMSD against an atom
# that may be the wrong one.
AUTODOCK_TYPE_TO_ELEMENT: dict[str, str] = {
    "H": "H",
    "HD": "H",
    "HS": "H",
    "HO": "H",
    "C": "C",
    "A": "C",
    "N": "N",
    "NA": "N",
    "NS": "N",
    "O": "O",
    "OA": "O",
    "OS": "O",
    "S": "S",
    "SA": "S",
    "P": "P",
    "F": "F",
    "Cl": "Cl",
    "CL": "Cl",
    "Br": "Br",
    "BR": "Br",
    "I": "I",
    # Metals/ions appear in receptor PDBQT files and in some ligands; they are
    # listed so a diff against a metal-containing ligand fails on chemistry
    # rather than on an unknown type. "NA" (nitrogen acceptor) and "Na" (sodium)
    # are case-sensitive AutoDock types.
    "Mg": "Mg",
    "MG": "Mg",
    "Mn": "Mn",
    "MN": "Mn",
    "Zn": "Zn",
    "ZN": "Zn",
    "Ca": "Ca",
    "CA": "Ca",
    "Fe": "Fe",
    "FE": "Fe",
    "K": "K",
    "Na": "Na",
}

PDB_ID_RE = re.compile(r"^[0-9][A-Za-z0-9]{3}$")


class TargetFailure(Exception):
    """A target that cannot be prepared, docked or measured.

    Raised for the expected failure modes (metadata cannot identify the ligand,
    a download fails, preparation fails, a pose cannot be matched to the ligand
    molecule). The caller records it as ``status: "failed"`` with the message as
    ``failure_reason`` -- it is never worked around and never hidden.
    """


# --------------------------------------------------------------------------- #
# small pure helpers (no third-party imports; directly unit-tested)
# --------------------------------------------------------------------------- #


def _rounded(value: float) -> float:
    """Round a score or an RMSD the way the file records it.

    Applied once and used for both the per-mode value and the target-level
    summary, so ``top1_rmsd`` and ``best_mode_rmsd`` cannot disagree with the
    modes they are supposed to summarise.
    """
    return round(float(value), RMSD_DECIMALS)


def _rounded_score(value: float) -> float:
    return round(float(value), SCORE_DECIMALS)


def parse_targets(raw: str) -> list[str]:
    """Parse ``--targets 1VRT,2HU4,1HSG`` into uppercase PDB IDs."""
    targets = [item.strip().upper() for item in raw.split(",") if item.strip()]
    if not targets:
        raise ValueError("--targets needs at least one PDB ID")
    for target in targets:
        if not PDB_ID_RE.match(target):
            raise ValueError(f"{target!r} is not a 4-character PDB ID")
    if len(set(targets)) != len(targets):
        raise ValueError(f"--targets repeats a PDB ID: {raw!r}")
    return targets


def sha256_of_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def ligand_centroid(atoms: Sequence[dict[str, float]]) -> list[float]:
    """Centroid of the ligand heavy atoms -- the centre of the docking box."""
    if not atoms:
        raise TargetFailure("the ligand has no atoms, so no box can be centred")
    count = len(atoms)
    return [sum(atom[axis] for atom in atoms) / count for axis in ("x", "y", "z")]


def _distance(first: dict[str, float], second: dict[str, float]) -> float:
    return math.sqrt(
        (first["x"] - second["x"]) ** 2
        + (first["y"] - second["y"]) ** 2
        + (first["z"] - second["z"]) ** 2
    )


def _pdb_atom(line: str) -> dict[str, Any]:
    """One ATOM/HETATM record as a dict. Raises on unreadable coordinates."""
    try:
        return {
            "record": line[:6].strip(),
            "serial": line[6:11].strip(),
            "name": line[12:16].strip(),
            "altloc": line[16],
            "residue": line[17:20].strip(),
            "chain": line[21],
            "resseq": line[22:27].strip(),
            "x": float(line[30:38]),
            "y": float(line[38:46]),
            "z": float(line[46:54]),
            "element": line[76:78].strip(),
            "line": line,
        }
    except ValueError as exc:
        raise TargetFailure(f"unreadable coordinate line: {line!r} ({exc})") from exc


def pdb_atoms(text: str, *, record: str | None = None) -> list[dict[str, Any]]:
    """ATOM and/or HETATM records of a PDB file, in file order.

    ``record="HETATM"`` or ``record="ATOM"`` restricts the result. Model 1 only:
    records after the first ``ENDMDL`` are ignored, because a docked-pose or
    multi-model entry is not what this reads.
    """
    atoms: list[dict[str, Any]] = []
    for line in text.splitlines():
        if line.startswith("ENDMDL"):
            break
        if line.startswith("HETATM") and record in (None, "HETATM"):
            atoms.append(_pdb_atom(line))
        elif line.startswith("ATOM") and record in (None, "ATOM"):
            atoms.append(_pdb_atom(line))
    return atoms


def ligand_instances(text: str, ligand_code: str) -> list[dict[str, Any]]:
    """Every deposited instance of ``ligand_code`` as HETATM records.

    Alternate locations: a residue is represented by its atoms with a blank
    altLoc plus those of the alphabetically first named altLoc it has (usually
    'A'). Deposited residues routinely leave shared atoms blank, so blank atoms
    are always kept; a residue with no named altLoc is represented by its blank
    atoms alone. The altLoc actually used is recorded on the instance rather than
    being chosen silently.
    """
    by_instance: dict[tuple[str, str], list[dict[str, Any]]] = {}
    order: list[tuple[str, str]] = []
    for atom in pdb_atoms(text, record="HETATM"):
        if atom["residue"] != ligand_code:
            continue
        key = (atom["chain"], atom["resseq"])
        if key not in by_instance:
            by_instance[key] = []
            order.append(key)
        by_instance[key].append(atom)

    instances: list[dict[str, Any]] = []
    for key in order:
        atoms = by_instance[key]
        named = sorted({atom["altloc"] for atom in atoms if atom["altloc"] != " "})
        altloc = named[0] if named else " "
        keep = {" ", altloc}
        selected = [atom for atom in atoms if atom["altloc"] in keep]
        if not selected:  # cannot happen with keep containing " ", but stay honest
            raise TargetFailure(
                f"no atoms could be selected for ligand residue {key[0]}:{key[1]}"
            )
        instances.append(
            {
                "chain": key[0],
                "resseq": key[1],
                "altloc": altloc,
                "atoms": selected,
                "hydrogen_atoms_in_deposit": sum(
                    1 for a in selected if a["element"] == "H"
                ),
            }
        )
    return instances


def chains_near_ligand(
    text: str,
    ligand_atoms: Sequence[dict[str, float]],
    cutoff: float = CHAIN_INTERFACE_CUTOFF_ANGSTROM,
) -> list[str]:
    """Protein chains with an atom within ``cutoff`` of any ligand heavy atom.

    Returned in the order the chains first appear in the file, so the choice is
    reproducible. A chain is kept whole, not trimmed to the residues inside the
    cutoff: a truncated chain would change the chemistry around the site for no
    documented reason.
    """
    protein = pdb_atoms(text, record="ATOM")
    chains: list[str] = []
    for atom in protein:
        if atom["chain"] not in chains:
            chains.append(atom["chain"])

    keep: list[str] = []
    for chain in chains:
        atoms = [a for a in protein if a["chain"] == chain]
        if any(_distance(atom, ligand) <= cutoff for atom in atoms for ligand in ligand_atoms):
            keep.append(chain)
    return keep


def receptor_pdb_text(text: str, chains: Sequence[str]) -> str:
    """The deposited protein chains as a PDB block: ATOM records only.

    Every HETATM record is dropped -- waters, ions, the ligand and any modified
    residue deposited as HETATM (1VRT chain A has one, CSD). Modified residues
    are therefore absent from the receptor, which is stated in the docs and in
    the per-target ``notes``.
    """
    wanted = set(chains)
    lines = ["REMARK   receptor prepared by scripts/redock_validation.py"]
    for atom in pdb_atoms(text, record="ATOM"):
        if atom["chain"] in wanted:
            lines.append(atom["line"].rstrip("\n"))
    lines.append("TER")
    lines.append("END")
    return "\n".join(lines) + "\n"


def parse_pdbqt_atoms(text: str) -> list[dict[str, Any]]:
    """ATOM/HETATM records of a PDBQT body, in file order, with their AD4 type.

    The AutoDock type is read from the trailing type column, with a fallback to
    the last whitespace-separated field of the line (both forms occur in the
    wild). A line whose type cannot be read keeps ``type: ""``, which fails the
    element check later rather than being replaced by a guess.
    """
    return [
        _pdbqt_atom(line)
        for line in text.splitlines()
        if line.startswith(("ATOM", "HETATM"))
    ]


def _pdbqt_atom(line: str) -> dict[str, Any]:
    """One PDBQT atom record with its AutoDock type attached."""
    atom = _pdb_atom(line)
    atom_type = line[77:80].strip()
    if not atom_type.isalpha():
        fields = [field for field in line[66:].split() if field.isalpha()]
        atom_type = fields[-1] if fields else ""
    atom["type"] = atom_type
    return atom


def parse_pdbqt_models(text: str) -> list[list[dict[str, Any]]]:
    """Split a multi-model PDBQT (Vina's ``--out`` file) into its models.

    Vina writes ``MODEL``/``ENDMDL`` per mode, best first, matching the order of
    ``energies()``. A file with no ``MODEL`` records is returned as one model.
    """
    models: list[list[dict[str, Any]]] = []
    current: list[dict[str, Any]] | None = None
    for line in text.splitlines():
        if line.startswith("MODEL"):
            if current:
                models.append(current)
            current = []
        elif line.startswith("ENDMDL"):
            if current is not None:
                models.append(current)
            current = None
        elif line.startswith(("ATOM", "HETATM")):
            if current is None:
                current = []
            current.append(_pdbqt_atom(line))
    if current:
        models.append(current)
    return models


def count_pdbqt_types(atoms: Iterable[dict[str, Any]]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for atom in atoms:
        key = atom.get("type", "") or "unreadable"
        counts[key] = counts.get(key, 0) + 1
    return counts


def element_of_pdbqt_type(atom_type: str) -> str | None:
    return AUTODOCK_TYPE_TO_ELEMENT.get(atom_type)


def format_summary_table(targets: Sequence[dict[str, Any]]) -> str:
    """The end-of-run table: one row per target, failures included."""
    header = (
        f"{'pdb_id':<8} {'ligand':<7} {'chains':<8} {'top1_rmsd':>10} "
        f"{'best_rmsd':>10} {'status':<8} {'passed':<7} note"
    )
    rows = [header, "-" * len(header)]
    for target in targets:
        top1 = target.get("top1_rmsd")
        best = target.get("best_mode_rmsd")
        reason = target.get("failure_reason") or ""
        if len(reason) > MAX_SUMMARY_REASON_CHARACTERS:
            reason = reason[: MAX_SUMMARY_REASON_CHARACTERS - 3] + "..."
        rows.append(
            f"{target['pdb_id']:<8} {str(target.get('ligand_code') or '-'):<7} "
            f"{','.join(target.get('receptor_chains_used') or []) or '-':<8} "
            f"{'-' if top1 is None else f'{top1:.3f}':>10} "
            f"{'-' if best is None else f'{best:.3f}':>10} "
            f"{target['status']:<8} {str(bool(target.get('passed'))):<7} {reason}"
        )
    return "\n".join(rows)


# --------------------------------------------------------------------------- #
# HTTP (requests, imported inside the function)
# --------------------------------------------------------------------------- #


def _require(module_name: str):
    """Import a third-party module or stop with the install line."""
    try:
        return __import__(module_name)
    except Exception as exc:  # noqa: BLE001 - report exactly what failed
        raise SystemExit(
            f"Cannot import {module_name!r} ({type(exc).__name__}: {exc}).\n"
            f"  platform : {platform.platform()}\n"
            f"  python   : {platform.python_version()}\n"
            "This script needs Vina and meeko, which have no Windows wheel.\n"
            f"Install them where they work (Linux or Google Colab):\n    {PIP_INSTALL_LINE}"
        ) from exc


def http_get_text(url: str, *, retries: int = REQUEST_RETRIES) -> str:
    """Download ``url`` with a timeout and retries.

    Retries only on transport errors, 429 and 5xx: a 404 or 403 will not fix
    itself by being asked again. Every attempt is bounded by
    ``REQUEST_TIMEOUT_SECONDS`` so a stalled host cannot hang the run.
    """
    requests = _require("requests")
    last_error = "no attempt was made"
    for attempt in range(1, retries + 1):
        try:
            response = requests.get(
                url,
                timeout=REQUEST_TIMEOUT_SECONDS,
                headers={"User-Agent": USER_AGENT},
            )
            if response.status_code == 200:
                return response.text
            last_error = f"HTTP {response.status_code}"
            if response.status_code < 500 and response.status_code != 429:
                break
        except Exception as exc:  # noqa: BLE001 - transport errors of any kind
            last_error = f"{type(exc).__name__}: {exc}"
        if attempt < retries:
            time.sleep(REQUEST_BACKOFF_SECONDS * attempt)
    raise TargetFailure(
        f"downloading {url} failed after {retries} attempt(s) ({last_error})"
    )


def http_get_json(url: str, *, retries: int = REQUEST_RETRIES) -> dict[str, Any]:
    text = http_get_text(url, retries=retries)
    try:
        payload = json.loads(text)
    except json.JSONDecodeError as exc:
        raise TargetFailure(f"{url} did not return JSON ({exc})") from exc
    if not isinstance(payload, dict):
        raise TargetFailure(f"{url} returned {type(payload).__name__}, expected an object")
    return payload


# --------------------------------------------------------------------------- #
# RCSB metadata: which residue code is the ligand
# --------------------------------------------------------------------------- #


def entry_metadata(pdb_id: str) -> dict[str, Any]:
    return http_get_json(f"{RCSB_DATA_API}/entry/{pdb_id}")


def nonpolymer_entity(pdb_id: str, entity_id: str) -> dict[str, Any]:
    return http_get_json(f"{RCSB_DATA_API}/nonpolymer_entity/{pdb_id}/{entity_id}")


def chem_comp(component_code: str) -> dict[str, Any]:
    return http_get_json(f"{RCSB_DATA_API}/chemcomp/{component_code}")


def identify_ligand(pdb_id: str, entry: dict[str, Any]) -> dict[str, Any]:
    """The ligand's HETATM residue code, from the RCSB entry metadata.

    Candidates are the entry's non-polymer entities, minus waters, minus
    entities whose chemical component definition has fewer than
    ``MIN_LIGAND_HEAVY_ATOMS`` heavy atoms. Exactly one candidate must remain:
    with none or several the target fails with the candidate list, because
    guessing a code is exactly what this must not do.
    """
    identifiers = entry.get("rcsb_entry_container_identifiers") or {}
    entity_ids = identifiers.get("non_polymer_entity_ids") or []
    if not entity_ids:
        raise TargetFailure(
            f"{pdb_id} lists no non-polymer entity, so no ligand can be identified "
            "from the metadata"
        )

    considered: list[dict[str, Any]] = []
    candidates: list[dict[str, Any]] = []
    for entity_id in entity_ids:
        entity = nonpolymer_entity(pdb_id, str(entity_id))
        nonpoly = entity.get("pdbx_entity_nonpoly") or {}
        code = str(nonpoly.get("comp_id") or "").strip()
        name = str(nonpoly.get("name") or "").strip()
        if not code:
            considered.append(
                {
                    "entity_id": str(entity_id),
                    "comp_id": None,
                    "excluded_because": "the entity metadata reports no comp_id",
                }
            )
            continue

        component = chem_comp(code)
        info = component.get("rcsb_chem_comp_info") or {}
        heavy = info.get("atom_count_heavy")
        entry_record = {
            "entity_id": str(entity_id),
            "comp_id": code,
            "name": name,
            "heavy_atom_count": heavy,
            "formula": (component.get("chem_comp") or {}).get("formula"),
        }

        if code in WATER_COMPONENT_CODES:
            entry_record["excluded_because"] = "water"
            considered.append(entry_record)
            continue
        if heavy is None:
            entry_record["excluded_because"] = (
                "the chemical component definition reports no atom_count_heavy"
            )
            considered.append(entry_record)
            continue
        if int(heavy) < MIN_LIGAND_HEAVY_ATOMS:
            entry_record["excluded_because"] = (
                f"only {heavy} heavy atoms (< {MIN_LIGAND_HEAVY_ATOMS}): an ion, "
                "buffer or other small molecule, not the ligand"
            )
            considered.append(entry_record)
            continue

        considered.append(entry_record)
        candidates.append(entry_record)

    if not candidates:
        raise TargetFailure(
            f"no non-polymer entity of {pdb_id} passes the ligand rule "
            f"(>= {MIN_LIGAND_HEAVY_ATOMS} heavy atoms, water excluded); "
            f"entities considered: {_brief(considered)}"
        )
    if len(candidates) > 1:
        raise TargetFailure(
            f"{pdb_id} has {len(candidates)} non-polymer entities that pass the ligand "
            f"rule ({_brief(candidates)}); the ligand cannot be identified with "
            "certainty and is not guessed"
        )

    selected = dict(candidates[0])
    selected["selection_rule"] = LIGAND_IDENTIFICATION_RULE
    selected["rcsb_entry_url"] = f"{RCSB_DATA_API}/entry/{pdb_id}"
    return {"selected": selected, "entities_considered": considered}


def chem_comp_definition(component_code: str) -> dict[str, Any]:
    """The RCSB chemical component definition used for the ligand's bonds."""
    component = chem_comp(component_code)
    descriptor = component.get("rcsb_chem_comp_descriptor") or {}
    smiles = descriptor.get("SMILES_stereo") or descriptor.get("SMILES")
    if not smiles:
        raise TargetFailure(
            f"the RCSB chemical component definition for {component_code} carries no "
            "SMILES, so the ligand's bond orders cannot be taken from it"
        )
    return {
        "comp_id": component_code,
        "smiles": str(smiles),
        "smiles_field": "SMILES_stereo" if descriptor.get("SMILES_stereo") else "SMILES",
        "inchi_key": descriptor.get("InChIKey"),
        "formula": (component.get("chem_comp") or {}).get("formula"),
        "atom_count_heavy": (component.get("rcsb_chem_comp_info") or {}).get(
            "atom_count_heavy"
        ),
        "source": f"{RCSB_DATA_API}/chemcomp/{component_code}",
    }


def _brief(records: Sequence[dict[str, Any]]) -> str:
    return "; ".join(
        f"{record.get('comp_id') or '?'} (heavy={record.get('heavy_atom_count')}"
        + (f", excluded: {record['excluded_because']}" if record.get("excluded_because") else "")
        + ")"
        for record in records
    )


# --------------------------------------------------------------------------- #
# structure preparation (rdkit + meeko, imported inside the functions)
# --------------------------------------------------------------------------- #


def prepare_crystal_ligand(
    ligand_code: str,
    instance: dict[str, Any],
    definition: dict[str, Any],
):
    """The deposited ligand as an RDKit molecule with the CCD's bond orders.

    Returns ``(heavy_mol, prep_mol, info)``: ``heavy_mol`` is the reference for
    the RMSD (heavy atoms, crystal coordinates), and ``prep_mol`` is the same
    molecule with hydrogens added for meeko and for the docked-pose mapping.
    Bond orders come from the chemical component definition through
    ``AssignBondOrdersFromTemplate``; the count of heavy atoms must agree, and
    that agreement is the condition for trusting the assignment.
    """
    try:
        from rdkit import Chem, RDLogger
        from rdkit.Chem import AllChem
    except Exception as exc:  # noqa: BLE001
        raise SystemExit(
            f"Cannot import rdkit ({type(exc).__name__}: {exc}).\n    {PIP_INSTALL_LINE}"
        ) from exc

    RDLogger.DisableLog("rdApp.warning")

    block = "\n".join(atom["line"].rstrip("\n") for atom in instance["atoms"]) + "\n"
    # sanitize=False: the PDB block has no reliable bond orders or valences; the
    # definition below supplies them.
    deposited = Chem.MolFromPDBBlock(block, sanitize=False, removeHs=True)
    if deposited is None:
        raise TargetFailure(
            f"rdkit could not read the deposited coordinates of {ligand_code} "
            f"(chain {instance['chain']} residue {instance['resseq']})"
        )

    template = Chem.MolFromSmiles(definition["smiles"])
    if template is None:
        raise TargetFailure(
            f"rdkit could not parse the chemical component definition for {ligand_code}"
        )

    try:
        assigned = AllChem.AssignBondOrdersFromTemplate(template, deposited)
    except Exception as exc:  # noqa: BLE001
        raise TargetFailure(
            f"assigning the {ligand_code} bond orders from its chemical component "
            f"definition failed ({type(exc).__name__}: {exc})"
        ) from exc
    if assigned is None:
        raise TargetFailure(
            f"the {ligand_code} bond-order assignment from its chemical component "
            "definition returned nothing; the deposited coordinates and the "
            "definition do not agree"
        )

    heavy_mol = Chem.Mol(assigned)
    if heavy_mol.GetNumConformers() == 0:
        raise TargetFailure(f"the {ligand_code} molecule has no 3D conformer")
    Chem.AssignStereochemistryFrom3D(heavy_mol)

    template_heavy = template.GetNumHeavyAtoms()
    if template_heavy != heavy_mol.GetNumHeavyAtoms():
        raise TargetFailure(
            f"the chemical component definition has {template_heavy} heavy atoms but the "
            f"deposited {ligand_code} residue has {heavy_mol.GetNumHeavyAtoms()}; the "
            "definition did not cover the deposited ligand, so its bond orders are not "
            "trustworthy"
        )
    if definition.get("atom_count_heavy") not in (None, template_heavy):
        raise TargetFailure(
            f"the chemical component definition reports "
            f"{definition['atom_count_heavy']} heavy atoms but parses to {template_heavy}; "
            "the definition is internally inconsistent"
        )

    prep_mol = Chem.AddHs(heavy_mol, addCoords=True)
    smiles_matches = Chem.MolToSmiles(Chem.RemoveHs(Chem.Mol(prep_mol)), isomericSmiles=False) == (
        Chem.MolToSmiles(Chem.RemoveHs(Chem.Mol(template)), isomericSmiles=False)
    )

    info = {
        "heavy_atom_count": heavy_mol.GetNumHeavyAtoms(),
        "atoms_with_hydrogens": prep_mol.GetNumAtoms(),
        "canonical_smiles_of_deposited_ligand": Chem.MolToSmiles(
            Chem.RemoveHs(Chem.Mol(prep_mol)), isomericSmiles=False
        ),
        "canonical_smiles_of_definition": Chem.MolToSmiles(
            Chem.RemoveHs(Chem.Mol(template)), isomericSmiles=False
        ),
        "deposited_matches_definition_smiles": bool(smiles_matches),
        "stereochemistry_source": "rdkit AssignStereochemistryFrom3D on the deposited coordinates",
    }
    return heavy_mol, prep_mol, info


def ligand_pdbqt(prep_mol) -> tuple[str, dict[str, Any]]:
    """Prepare the ligand PDBQT with meeko (default preparation)."""
    try:
        from meeko import MoleculePreparation, PDBQTWriterLegacy
    except Exception as exc:  # noqa: BLE001
        raise SystemExit(
            f"Cannot import meeko ({type(exc).__name__}: {exc}).\n    {PIP_INSTALL_LINE}"
        ) from exc

    preparation = MoleculePreparation()
    prepare = getattr(preparation, "prepare", None) or preparation
    try:
        setups = prepare(prep_mol)
    except Exception as exc:  # noqa: BLE001
        raise TargetFailure(
            f"meeko could not prepare the ligand ({type(exc).__name__}: {exc})"
        ) from exc
    # Both documented forms return the setup(s) created for the molecule; make a
    # single setup a list of one rather than assuming either shape.
    if not isinstance(setups, (list, tuple)):
        setups = [setups]
    if not setups or setups[0] is None:
        raise TargetFailure("meeko returned no molecule setup for the ligand")

    written = PDBQTWriterLegacy.write_string(setups[0])
    # write_string returns (pdbqt, success, error_message).
    if not isinstance(written, tuple) or len(written) < 2:
        raise TargetFailure(
            f"meeko's PDBQT writer returned {written!r}, which is not the documented "
            "(pdbqt, success, message) tuple"
        )
    pdbqt, ok = written[0], written[1]
    message = written[2] if len(written) > 2 else ""
    if not ok:
        raise TargetFailure(f"meeko refused to write the ligand PDBQT: {message}")

    atoms = parse_pdbqt_atoms(pdbqt)
    if not atoms:
        raise TargetFailure("the ligand PDBQT contains no atom records")
    return pdbqt, {
        "atom_count": len(atoms),
        "type_counts": count_pdbqt_types(atoms),
        "hydrogen_atom_count": sum(
            1 for atom in atoms if element_of_pdbqt_type(atom["type"]) == "H"
        ),
    }


def receptor_pdbqt(receptor_pdb_text: str) -> tuple[str, dict[str, Any]]:
    """Convert the prepared receptor PDB into PDBQT with meeko."""
    try:
        from meeko import PDBQTWriterLegacy, Polymer
    except Exception as exc:  # noqa: BLE001
        raise SystemExit(
            f"Cannot import meeko ({type(exc).__name__}: {exc}).\n    {PIP_INSTALL_LINE}"
        ) from exc

    try:
        polymer = Polymer.from_pdb_string(receptor_pdb_text)
        rigid, flexible = PDBQTWriterLegacy.write_string_from_polymer(polymer)
    except Exception as exc:  # noqa: BLE001
        raise TargetFailure(
            f"meeko could not prepare the receptor ({type(exc).__name__}: {exc})"
        ) from exc

    if flexible:
        raise TargetFailure(
            "meeko produced a flexible-residue receptor PDBQT, but this pipeline docks "
            "against a rigid receptor; the run would not match the recorded setup"
        )
    if not rigid:
        raise TargetFailure("meeko produced an empty receptor PDBQT")

    atoms = parse_pdbqt_atoms(rigid)
    if not atoms:
        raise TargetFailure("the receptor PDBQT contains no atom records")
    hydrogens = sum(1 for atom in atoms if element_of_pdbqt_type(atom["type"]) == "H")
    if hydrogens == 0:
        raise TargetFailure(
            "the prepared receptor PDBQT contains no hydrogens, so meeko did not add "
            "them as documented and the receptor is not what this pipeline claims to dock"
        )
    return rigid, {
        "atom_count": len(atoms),
        "hydrogen_atom_count": hydrogens,
        "type_counts": count_pdbqt_types(atoms),
    }


# --------------------------------------------------------------------------- #
# RMSD: heavy-atom, symmetry-aware, no re-alignment
# --------------------------------------------------------------------------- #


def map_pose_atoms_to_molecule(
    prep_mol, pose_atoms: Sequence[dict[str, Any]], *, where: str
) -> list[int]:
    """Indices of ``prep_mol``'s heavy atoms matched to the pose's heavy atoms.

    The correspondence is by file order, and it is **verified** rather than
    assumed: the count of heavy atoms must agree, and every pose atom's AutoDock
    type must map to the element of the molecule atom it is matched to. An
    unknown type or a single mismatch fails the target, because an RMSD computed
    against the wrong atom is worse than no RMSD.
    """
    heavy_indices = [
        atom.GetIdx()
        for atom in prep_mol.GetAtoms()
        if atom.GetAtomicNum() > 1  # implicit/explicit hydrogens excluded
    ]
    pose_heavy = [atom for atom in pose_atoms if element_of_pdbqt_type(atom["type"]) != "H"]

    unreadable = [atom["type"] for atom in pose_atoms if element_of_pdbqt_type(atom["type"]) is None]
    if unreadable:
        raise TargetFailure(
            f"{where}: AutoDock type(s) {sorted(set(unreadable))} are not in the "
            "documented type-to-element table, so the pose cannot be matched to the "
            "ligand molecule"
        )
    if len(pose_heavy) != len(heavy_indices):
        raise TargetFailure(
            f"{where}: the pose has {len(pose_heavy)} heavy atoms but the ligand molecule "
            f"has {len(heavy_indices)}; the pose cannot be mapped onto the molecule"
        )

    for atom, molecule_index in zip(pose_heavy, heavy_indices):
        pose_element = element_of_pdbqt_type(atom["type"])
        molecule_element = prep_mol.GetAtomWithIdx(molecule_index).GetSymbol()
        if pose_element != molecule_element:
            raise TargetFailure(
                f"{where}: pose heavy atom {atom.get('serial') or '?'} is typed "
                f"{atom['type']!r} ({pose_element}) but the molecule atom it would map to "
                f"is {molecule_element}; the atom correspondence is not trustworthy"
            )
    return heavy_indices


def pose_rmsd_angstrom(prep_mol, heavy_mol, pose_atoms: Sequence[dict[str, Any]], *, where: str):
    """Symmetry-aware heavy-atom RMSD of one pose against the crystal ligand.

    No re-alignment: ``rdMolAlign.CalcRMS`` is given the coordinates as they are.
    """
    heavy_indices = map_pose_atoms_to_molecule(prep_mol, pose_atoms, where=where)

    try:
        from rdkit import Chem
        from rdkit.Chem import rdMolAlign
        from rdkit.Geometry import Point3D
    except Exception as exc:  # noqa: BLE001
        raise SystemExit(
            f"Cannot import rdkit ({type(exc).__name__}: {exc}).\n    {PIP_INSTALL_LINE}"
        ) from exc

    pose_mol = Chem.Mol(heavy_mol)
    # The reference molecule is heavy-atom only, and AddHs appends hydrogens, so
    # the heavy atoms keep indices 0..n-1. The mapping must therefore be the
    # identity: if it is not, the coordinates would be written onto the wrong
    # atoms and the RMSD would be meaningless.
    if heavy_indices != list(range(pose_mol.GetNumAtoms())):
        raise TargetFailure(
            f"{where}: the ligand molecule's heavy atoms are not its first "
            f"{pose_mol.GetNumAtoms()} atoms ({heavy_indices[:8]}...), so the pose "
            "coordinates cannot be placed on the reference molecule"
        )
    conformer = pose_mol.GetConformer()
    pose_heavy = [a for a in pose_atoms if element_of_pdbqt_type(a["type"]) != "H"]
    for molecule_index, atom in zip(heavy_indices, pose_heavy):
        conformer.SetAtomPosition(
            molecule_index, Point3D(float(atom["x"]), float(atom["y"]), float(atom["z"]))
        )

    try:
        value = rdMolAlign.CalcRMS(
            pose_mol, heavy_mol, symmetrizeConjugatedTerminalGroups=True
        )
    except Exception as exc:  # noqa: BLE001
        raise TargetFailure(f"{where}: rdkit could not compute the RMSD ({exc})") from exc
    value = float(value)
    if not math.isfinite(value) or value < 0.0:
        raise TargetFailure(f"{where}: rdkit returned an unusable RMSD ({value!r})")
    return value


def crystal_pose_self_check(prep_mol, heavy_mol, ligand_pdbqt_text: str) -> float:
    """RMSD of the prepared ligand's own coordinates against the crystal pose.

    This is the internal consistency check of the whole path -- crystal block ->
    CCD bond orders -> meeko PDBQT -> pose reading -> RMSD. The answer must be
    ~0 A: anything else means the atom correspondence, not the docking, is wrong.
    It is not a docking result and is recorded as such.
    """
    own_atoms = parse_pdbqt_atoms(ligand_pdbqt_text)
    if not own_atoms:
        raise TargetFailure("the ligand PDBQT has no atom records to self-check")
    return pose_rmsd_angstrom(
        prep_mol, heavy_mol, own_atoms, where="crystal-pose self-check"
    )


# --------------------------------------------------------------------------- #
# Vina
# --------------------------------------------------------------------------- #


def dock(
    receptor_pdbqt_text: str,
    ligand_pdbqt_text: str,
    center: Sequence[float],
    *,
    exhaustiveness: int,
    seed: int,
    work_dir: Path,
) -> dict[str, Any]:
    """Run Vina once with the given fixed settings and read back its poses."""
    try:
        import numpy as np
        from vina import Vina
    except Exception as exc:  # noqa: BLE001
        raise SystemExit(
            f"Cannot import vina ({type(exc).__name__}: {exc}).\n    {PIP_INSTALL_LINE}"
        ) from exc

    receptor_path = work_dir / "receptor.pdbqt"
    poses_path = work_dir / "docked_poses.pdbqt"
    receptor_path.write_text(receptor_pdbqt_text, encoding="utf-8")

    engine = Vina(sf_name="vina", cpu=VINA_CPU, seed=seed, verbosity=0)
    engine.set_receptor(str(receptor_path))
    engine.set_ligand_from_string(ligand_pdbqt_text)
    engine.compute_vina_maps(center=list(center), box_size=list(BOX_SIZE_ANGSTROM))
    engine.dock(exhaustiveness=exhaustiveness, n_poses=NUM_MODES)

    scores = np.asarray(
        engine.energies(n_poses=NUM_MODES, energy_range=ENERGY_RANGE_KCAL_PER_MOL)
    )
    if scores.ndim != 2 or scores.shape[0] == 0:
        raise TargetFailure(f"Vina returned no scores for any pose (shape {scores.shape})")
    # Column 0 is the total score in kcal/mol, rows best-first.
    totals = [float(row[0]) for row in scores]
    if not all(math.isfinite(value) for value in totals):
        raise TargetFailure(f"Vina returned a non-finite score: {totals}")

    engine.write_poses(
        str(poses_path), n_poses=NUM_MODES, energy_range=ENERGY_RANGE_KCAL_PER_MOL, overwrite=True
    )
    poses_text = poses_path.read_text(encoding="utf-8")

    seed_used = None
    try:
        seed_used = int(engine.info()["seed"])
    except Exception:  # noqa: BLE001 - the actual seed is a bonus, not a requirement
        seed_used = None

    return {
        "scores_kcal_per_mol": totals,
        "poses_pdbqt": poses_text,
        "seed_used": seed_used,
        "receptor_pdbqt_sha256": sha256_of_text(receptor_pdbqt_text),
        "ligand_pdbqt_sha256": sha256_of_text(ligand_pdbqt_text),
        "docked_poses_pdbqt_sha256": sha256_of_text(poses_text),
        "receptor_pdbqt_path": str(receptor_path),
        "docked_poses_pdbqt_path": str(poses_path),
    }


# --------------------------------------------------------------------------- #
# one target, end to end
# --------------------------------------------------------------------------- #


def empty_target_record(pdb_id: str, *, exhaustiveness: int, seed: int) -> dict[str, Any]:
    """The record of a target: everything null/empty until it is measured.

    Written before any work so a failure has somewhere to be recorded, and every
    key is always present, so the validator reads the same shape for an ok and a
    failed target.
    """
    return {
        "pdb_id": pdb_id,
        "ligand_code": None,
        "ligand_identification": None,
        "ligand_definition": None,
        "ligand_instance": None,
        "ligand_preparation": None,
        "receptor_chains_used": [],
        "receptor_preparation": None,
        "box_center": None,
        "box_size": list(BOX_SIZE_ANGSTROM),
        "vina_params": {
            "exhaustiveness": exhaustiveness,
            "seed": seed,
            "num_modes": NUM_MODES,
            "energy_range_kcal_per_mol": ENERGY_RANGE_KCAL_PER_MOL,
            "cpu": VINA_CPU,
        },
        "modes": [],
        "top1_rmsd": None,
        "best_mode_rmsd": None,
        "passed": False,
        "status": "failed",
        "failure_reason": None,
        "failure_type": None,
        "inputs": None,
        "notes": [],
        "crystal_pose_self_check_rmsd_angstrom": None,
    }


def run_target(pdb_id: str, *, exhaustiveness: int, seed: int, work_dir: Path) -> dict[str, Any]:
    """Everything for one complex, with every failure recorded rather than raised.

    ``TargetFailure`` and any unexpected exception are both written into the
    record as ``status: "failed"`` with a reason. Nothing is retried with other
    parameters and nothing is skipped.
    """
    record = empty_target_record(pdb_id, exhaustiveness=exhaustiveness, seed=seed)
    try:
        _run_target(record, pdb_id, exhaustiveness=exhaustiveness, seed=seed, work_dir=work_dir)
    except TargetFailure as failure:
        record["failure_reason"] = str(failure)
        record["failure_type"] = "TargetFailure"
    except Exception as exc:  # noqa: BLE001 - an unexpected failure is still a failure
        record["failure_reason"] = f"unexpected {type(exc).__name__}: {exc}"
        record["failure_type"] = type(exc).__name__
    return record


def _run_target(
    record: dict[str, Any],
    pdb_id: str,
    *,
    exhaustiveness: int,
    seed: int,
    work_dir: Path,
) -> None:
    """The body of :func:`run_target`; raises on anything that cannot be done."""
    entry = entry_metadata(pdb_id)
    identification = identify_ligand(pdb_id, entry)
    record["ligand_identification"] = identification
    ligand_code = str(identification["selected"]["comp_id"])
    record["ligand_code"] = ligand_code
    print(f"  ligand: {ligand_code} ({identification['selected']['name']})")

    definition = chem_comp_definition(ligand_code)
    record["ligand_definition"] = definition

    entry_text = http_get_text(f"{RCSB_FILES}/{pdb_id}.pdb")
    instances = ligand_instances(entry_text, ligand_code)
    if not instances:
        raise TargetFailure(
            f"the entry's coordinates contain no HETATM record with residue code "
            f"{ligand_code}, although the metadata names it as the ligand"
        )
    instance = instances[0]
    record["ligand_instance"] = {
        "chain": instance["chain"],
        "resseq": instance["resseq"],
        "altloc": instance["altloc"],
        "heavy_atom_count": len(instance["atoms"]),
        "deposited_instances": len(instances),
        "selection_rule": (
            "the first deposited instance in the coordinate file; every instance "
            "carries the same residue code"
        ),
    }
    if instance["hydrogen_atoms_in_deposit"]:
        record["notes"].append(
            f"the deposited ligand instance includes "
            f"{instance['hydrogen_atoms_in_deposit']} hydrogen atom(s)"
        )

    heavy_mol, prep_mol, ligand_info = prepare_crystal_ligand(ligand_code, instance, definition)
    if not ligand_info["deposited_matches_definition_smiles"]:
        record["notes"].append(
            "the deposited ligand's canonical SMILES (without stereochemistry) differs "
            "from the chemical component definition's; the assigned bond orders come "
            "from the definition, and the difference is recorded rather than hidden"
        )
    record["ligand_preparation"] = ligand_info

    chains = chains_near_ligand(entry_text, instance["atoms"])
    if not chains:
        raise TargetFailure(
            f"no protein chain has an atom within {CHAIN_INTERFACE_CUTOFF_ANGSTROM:g} A of "
            f"the {ligand_code} residue, so there is no receptor site to dock into"
        )
    record["receptor_chains_used"] = chains
    record["notes"].append(
        f"receptor chains selected by a {CHAIN_INTERFACE_CUTOFF_ANGSTROM:g} A contact "
        "cutoff to the crystal ligand; all HETATM records are removed, so a modified "
        "residue deposited as HETATM is absent from the receptor"
    )

    receptor_text = receptor_pdb_text(entry_text, chains)
    receptor_text_pdbqt, receptor_info = receptor_pdbqt(receptor_text)
    record["receptor_preparation"] = {
        **receptor_info,
        "chains": chains,
        "pdb_atom_count": len(pdb_atoms(receptor_text, record="ATOM")),
        "pdb_sha256": sha256_of_text(receptor_text),
    }

    ligand_pdbqt_text, ligand_pdbqt_info = ligand_pdbqt(prep_mol)
    record["ligand_preparation"] = {**record["ligand_preparation"], **ligand_pdbqt_info}

    self_check = crystal_pose_self_check(prep_mol, heavy_mol, ligand_pdbqt_text)
    record["crystal_pose_self_check_rmsd_angstrom"] = _rounded(self_check)
    if self_check > 0.01:
        raise TargetFailure(
            f"the prepared ligand's own coordinates give an RMSD of {self_check:.3f} A "
            "against the crystal pose instead of ~0; the atom correspondence or the "
            "coordinate handling is wrong, so no docked-pose RMSD can be trusted"
        )

    center = ligand_centroid(instance["atoms"])
    record["box_center"] = [round(value, 3) for value in center]

    result = dock(
        receptor_text_pdbqt,
        ligand_pdbqt_text,
        center,
        exhaustiveness=exhaustiveness,
        seed=seed,
        work_dir=work_dir,
    )
    record["inputs"] = {
        "rcsb_entry_pdb_url": f"{RCSB_FILES}/{pdb_id}.pdb",
        "rcsb_entry_metadata_url": f"{RCSB_DATA_API}/entry/{pdb_id}",
        "rcsb_chem_comp_url": definition["source"],
        "receptor_pdbqt_sha256": result["receptor_pdbqt_sha256"],
        "ligand_pdbqt_sha256": result["ligand_pdbqt_sha256"],
        "docked_poses_pdbqt_sha256": result["docked_poses_pdbqt_sha256"],
    }
    record["vina_params"]["seed_used"] = result["seed_used"]

    models = parse_pdbqt_models(result["poses_pdbqt"])
    if len(models) != len(result["scores_kcal_per_mol"]):
        raise TargetFailure(
            f"Vina wrote {len(models)} pose model(s) but returned "
            f"{len(result['scores_kcal_per_mol'])} score(s); a score cannot be attributed "
            "to a pose, so no mode can be reported"
        )
    if not models:
        raise TargetFailure(
            "Vina returned no docked pose; the search evaluated nothing, which is a "
            "failed run rather than a result"
        )

    modes: list[dict[str, Any]] = []
    for rank, (model, score) in enumerate(zip(models, result["scores_kcal_per_mol"]), start=1):
        rmsd = pose_rmsd_angstrom(
            prep_mol, heavy_mol, model, where=f"docked mode {rank}"
        )
        modes.append(
            {
                "rank": rank,
                "score_kcal_per_mol": _rounded_score(score),
                "rmsd_to_crystal_angstrom": _rounded(rmsd),
            }
        )

    record["modes"] = modes
    record["top1_rmsd"] = modes[0]["rmsd_to_crystal_angstrom"]
    record["best_mode_rmsd"] = min(mode["rmsd_to_crystal_angstrom"] for mode in modes)
    record["passed"] = bool(record["top1_rmsd"] < PASS_THRESHOLD_ANGSTROM)
    record["status"] = "ok"
    print(
        f"  docked: {len(modes)} mode(s), top1 RMSD {record['top1_rmsd']:.3f} A, "
        f"best {record['best_mode_rmsd']:.3f} A, passed={record['passed']}"
    )


# --------------------------------------------------------------------------- #
# payload and CLI
# --------------------------------------------------------------------------- #


def collect_tool_versions() -> dict[str, str]:
    """The versions the result is attributable to; "unavailable" if not importable."""
    versions = {"python": platform.python_version()}
    for module_name in ("numpy", "rdkit", "vina", "meeko"):
        versions[module_name] = _version_of(module_name)
    return versions


def _version_of(module_name: str) -> str:
    try:
        module = __import__(module_name)
    except Exception:  # noqa: BLE001
        return "unavailable"
    value = getattr(module, "__version__", None)
    if value:
        return str(value)
    try:
        from importlib.metadata import version

        return str(version(module_name))
    except Exception:  # noqa: BLE001
        return "unknown"


def check_pinned_versions(versions: dict[str, str]) -> None:
    """Stop when an enforced pin is missing, naming what to install."""
    problems = [
        f"{name}: expected {expected}, found {versions.get(name, 'unavailable')}"
        for name, expected in PINNED_VERSIONS.items()
        if versions.get(name) != expected
    ]
    if problems:
        raise SystemExit(
            "Dependency versions do not match the pins this pipeline is written "
            "against.\n  "
            + "\n  ".join(problems)
            + f"\nInstall the pinned set, or change the pin here and in the Colab "
            f"notebook together:\n    {PIP_INSTALL_LINE}"
        )


def build_payload(
    targets: Sequence[dict[str, Any]],
    *,
    versions: dict[str, str],
    exhaustiveness: int,
    seed: int,
) -> dict[str, Any]:
    """Assemble the result file. Nothing here computes or adjusts a number."""
    return {
        "schema_version": SCHEMA_VERSION,
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds").replace(
            "+00:00", "Z"
        ),
        "generated_by": "scripts/redock_validation.py",
        "tool_versions": versions,
        "pass_threshold_angstrom": PASS_THRESHOLD_ANGSTROM,
        "pass_rule": PASS_RULE,
        "box_size_angstrom": list(BOX_SIZE_ANGSTROM),
        "setup": {
            "kind": "redocking of the deposited crystal pose",
            "box_center": "centroid of the deposited crystal ligand's heavy atoms (the true site)",
            "box_size_angstrom": list(BOX_SIZE_ANGSTROM),
            "note": BOX_NOTE,
        },
        "vina_params": {
            "exhaustiveness": exhaustiveness,
            "seed": seed,
            "num_modes": NUM_MODES,
            "energy_range_kcal_per_mol": ENERGY_RANGE_KCAL_PER_MOL,
            "cpu": VINA_CPU,
            "one_setting_for_all_targets": True,
            "reproducibility": VINA_CPU_NOTE,
        },
        "ligand_identification_rule": LIGAND_IDENTIFICATION_RULE,
        "ligand_preparation": LIGAND_PREPARATION,
        "receptor_preparation": RECEPTOR_PREPARATION,
        "rmsd_method": RMSD_METHOD,
        "atom_mapping_method": ATOM_MAPPING_METHOD,
        "crystal_pose_self_check": (
            "before docking, the ligand PDBQT's own coordinates are measured against the "
            "crystal ligand with the same RMSD code; the value must be ~0 A "
            "(crystal_pose_self_check_rmsd_angstrom) or the target fails. It checks the "
            "preparation and the pose reading, and is not a docking result."
        ),
        "disclaimers": list(DISCLAIMERS),
        "targets": list(targets),
    }


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Redock known crystal complexes with AutoDock Vina and measure the "
            "heavy-atom, symmetry-aware RMSD of every docked pose against the "
            "deposited crystal pose."
        )
    )
    parser.add_argument("out_path", help="where to write the JSON result file")
    parser.add_argument(
        "--targets",
        default=",".join(DEFAULT_TARGETS),
        help=f"comma-separated PDB IDs (default: {','.join(DEFAULT_TARGETS)})",
    )
    parser.add_argument(
        "--exhaustiveness",
        type=int,
        default=DEFAULT_EXHAUSTIVENESS,
        help=f"Vina exhaustiveness, one value for every target (default: {DEFAULT_EXHAUSTIVENESS})",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=DEFAULT_SEED,
        help=f"Vina seed, one value for every target (default: {DEFAULT_SEED})",
    )
    parser.add_argument(
        "--skip-version-check",
        action="store_true",
        help=(
            "record the tool versions without enforcing the vina/meeko pins "
            "(for a host whose packages cannot be changed)"
        ),
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        targets = parse_targets(args.targets)
    except ValueError as exc:
        print(f"usage error: {exc}", file=sys.stderr)
        return 2
    if args.exhaustiveness < 1:
        print("usage error: --exhaustiveness must be at least 1", file=sys.stderr)
        return 2

    versions = collect_tool_versions()
    if not args.skip_version_check:
        check_pinned_versions(versions)

    print("versions:", ", ".join(f"{k}={v}" for k, v in versions.items()))
    print(f"platform: {platform.platform()}  python {platform.python_version()}")
    print(
        f"targets: {', '.join(targets)}  exhaustiveness={args.exhaustiveness}  "
        f"seed={args.seed}  num_modes={NUM_MODES}  box={list(BOX_SIZE_ANGSTROM)}"
    )
    print(f"one fixed setting for every target; pass rule: {PASS_RULE}")
    print()

    records: list[dict[str, Any]] = []
    with tempfile.TemporaryDirectory(prefix="qdde_redock_") as temporary:
        for pdb_id in targets:
            print(f"{pdb_id}:")
            record = run_target(
                pdb_id,
                exhaustiveness=args.exhaustiveness,
                seed=args.seed,
                work_dir=Path(temporary),
            )
            records.append(record)
            if record["status"] != "ok":
                print(f"  FAILED ({record['failure_type']}): {record['failure_reason']}")
            print()

    payload = build_payload(
        records, versions=versions, exhaustiveness=args.exhaustiveness, seed=args.seed
    )

    out_path = Path(args.out_path)
    try:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    except OSError as exc:
        print(f"could not write {out_path}: {exc}", file=sys.stderr)
        return 1

    print(format_summary_table(records))
    print()
    ok_targets = sum(1 for record in records if record["status"] == "ok")
    passed = sum(1 for record in records if record["passed"])
    failed = len(records) - ok_targets
    print(f"wrote {out_path}")
    print(
        f"{ok_targets}/{len(records)} target(s) docked, {passed} passed the "
        f"{PASS_THRESHOLD_ANGSTROM:g} A rule, {failed} recorded as failed."
    )
    if failed:
        print("  failed targets are part of the result, not hidden; see failure_reason.")
    print("  Reminder: these are docking predictions for an easy redocking test, not")
    print("  measured binding, and not a virtual-screening result.")
    return 0 if ok_targets else 3


if __name__ == "__main__":
    raise SystemExit(main())
