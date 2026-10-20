from __future__ import annotations

import json
import logging
import re
from pathlib import Path
from typing import Any

from app.config import get_settings
from app.services.external_api import http_get_bytes_cached, http_get_json_cached
from app.services.structure_analysis import (
    analyse_complex,
    build_mechanism,
    build_viewer_pdb,
    drug_like_ligands,
    parse_pdb,
)

logger = logging.getLogger(__name__)

RSCB_PDB_DOWNLOAD_BASE = "https://files.rcsb.org/download"
RCSB_DATA_ENTRY_URL = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"

# A deposited PDB identifier: four characters, the first a digit (1-9), then
# letters or digits. Checked before any request, so a typo costs no network call.
PDB_ID_PATTERN = re.compile(r"^[0-9][A-Za-z0-9]{3}$")

# Caps for "explore any PDB entry". Both are project choices, stated in the API
# contract: the contact scan is O(ligand atoms x protein atoms) and a whole
# ribosome would keep a request busy for a very long time without telling the
# caller anything about its drug. The byte cap also bounds the viewer payload.
MAX_ENTRY_BYTES = 25 * 1024 * 1024
MAX_PROTEIN_ATOMS = 60_000

CURATION_NONE_NOTE = "Not curated: measured contacts only"
NO_MECHANISM_NOTE = (
    "No curated mechanism for this entry; values below are measured geometry."
)

EXPLORE_LIGAND_RULE = (
    "Deposited HETATM groups only; water, ions, sugars, cofactors, buffers/"
    "crystallisation additives and standard amino acids are excluded, and at least "
    "10 heavy atoms must be present in the deposited coordinates. The count is "
    "measured from the file, not taken from metadata."
)


class ComplexRequestError(RuntimeError):
    """A request the caller got wrong, carrying the standard error shape."""

    def __init__(self, code: str, message: str, status_code: int) -> None:
        self.code = code
        self.message = message
        self.status_code = status_code
        super().__init__(message)


def _load_complexes_manifest(path: Path) -> dict:
    if not path.exists():
        raise FileNotFoundError("Complexes manifest missing")
    return json.loads(path.read_text(encoding="utf-8"))


def list_complexes(manifest_path: Path) -> dict[str, Any]:
    manifest = _load_complexes_manifest(manifest_path)
    complexes = []
    for item in manifest.get("complexes", []):
        complexes.append({
            "id": item.get("id"),
            "display_name": item.get("display_name"),
            "protein_pdb": item.get("protein_pdb"),
            "ligand_name": item.get("ligand_name"),
            "ligand_code": item.get("ligand_code"),
            "target_class": item.get("target_class"),
            # Distinguishes viral examples from human ones, so the UI never
            # implies that a human protein is a virus.
            "target_origin": item.get("target_origin"),
            "category": item.get("category"),
            # "curated" carries a literature mechanism; "contacts_only" carries
            # measured geometry and nothing else.
            "tier": item.get("tier", "curated"),
            "resolution_angstrom": item.get("resolution_angstrom"),
            "binding_site_type": item.get("binding_site_type"),
            "status": item.get("status", "placeholder"),
        })
    return {"complexes": complexes, "count": len(complexes)}


def _find_complex(manifest_path: Path, complex_id: str) -> dict:
    manifest = _load_complexes_manifest(manifest_path)
    for item in manifest.get("complexes", []):
        if item.get("id") == complex_id:
            return dict(item)
    raise KeyError(complex_id)


def _fetch_pdb_text(pdb_id: str) -> str | None:
    """Download the deposited PDB-format coordinates for a complex.

    The PDB format is used (not mmCIF) because parsing ATOM/HETATM records and
    feeding 3Dmol a trimmed protein+ligand file is straightforward in it.
    """
    url = f"{RSCB_PDB_DOWNLOAD_BASE}/{pdb_id}.pdb"
    try:
        raw = http_get_bytes_cached(url, source="pdb", backoff_seconds=2.0)
    except Exception as exc:  # noqa: BLE001
        logger.warning("Complex structure download failed for %s: %s", pdb_id, exc)
        return None
    if not raw:
        return None
    return raw.decode("utf-8", "replace")


def fetch_complex_structure(complex_id: str, manifest_path: Path) -> dict[str, Any]:
    """Fetch the crystal structure and compute the drug's real contacts with it."""
    try:
        entry = _find_complex(manifest_path, complex_id)
    except KeyError:
        return {
            "complex_id": complex_id,
            "status": "placeholder",
            "note": f"No curated entry for complex {complex_id!r}",
        }

    pdb_id = entry.get("protein_pdb")
    ligand_code = entry.get("ligand_code")

    tier = entry.get("tier", "curated")

    base: dict[str, Any] = {
        "complex_id": complex_id,
        "display_name": entry.get("display_name"),
        "protein_pdb": pdb_id,
        "protein_title": entry.get("protein_title"),
        "ligand_name": entry.get("ligand_name"),
        "ligand_code": ligand_code,
        "complex_type": entry.get("complex_type"),
        "target_origin": entry.get("target_origin"),
        "category": entry.get("category"),
        "tier": tier,
        "curation": "curated" if tier == "curated" else "none",
        "curation_note": None if tier == "curated" else CURATION_NONE_NOTE,
        "resolution_angstrom": entry.get("resolution_angstrom"),
        "source_organism": entry.get("source_organism"),
        "binding_site_type": entry.get("binding_site_type"),
    }

    if not pdb_id:
        base.update({"status": "placeholder", "note": "Curated entry has no PDB identifier"})
        return base

    text = _fetch_pdb_text(str(pdb_id))
    if text is None:
        base.update({
            "status": "placeholder",
            "note": (
                f"Could not download the deposited coordinates for {pdb_id} "
                "(RCSB PDB unreachable or the entry has no PDB-format file)."
            ),
        })
        return base

    analysis = analyse_complex(text, ligand_code)
    if tier == "curated":
        mechanism = build_mechanism(analysis, entry)
    else:
        # A contacts-only entry has no literature narrative, so the mechanism
        # builder is not run at all -- there is nothing to populate it with, and
        # an empty-but-shaped step list would look like a curated one.
        mechanism = {
            "status": "not_curated",
            "steps": [],
            "catalytic_residues": [],
            "sources": [],
            "target_class": None,
            "drug_class": None,
            "note": NO_MECHANISM_NOTE,
        }
    viewer_pdb = build_viewer_pdb(text, ligand_code) if analysis.get("status") == "computed" else None

    base.update({
        "status": analysis.get("status", "placeholder"),
        "structure_file_reference": f"{pdb_id}.pdb",
        "structure_file_bytes": viewer_pdb,
        "protein_atom_count": analysis.get("protein_atom_count"),
        "ligand_atom_count": analysis.get("ligand_atom_count"),
        "binding_site_residue_count": analysis.get("binding_site_residue_count"),
        "contact_summary": analysis.get("contact_summary"),
        "contact_summary_basis": analysis.get("contact_summary_basis"),
        "contact_rules": analysis.get("contact_rules"),
        "contact_labels": analysis.get("contact_labels"),
        "near_miss_hydrogen_bonds": analysis.get("near_miss_hydrogen_bonds") or [],
        "atom_pair_counts": analysis.get("atom_pair_counts"),
        "atom_pairs_truncated": analysis.get("atom_pairs_truncated", False),
        "atom_pairs_truncated_note": analysis.get("atom_pairs_truncated_note"),
        "interacting_residues": analysis.get("interacting_residues", []),
        "binding_site_residues": analysis.get("binding_site_residues", []),
        "method": analysis.get("method"),
        "mechanism": mechanism,
        "citations": entry.get("citations") or [],
        # Each narrative statement carries its own LITERATURE tag and the source
        # it came from, so a reader can check one sentence without trusting the
        # whole entry.
        "literature_statements": entry.get("literature_statements") or [],
        "sources": entry.get("sources") or [],
        # Say it explicitly when a curated entry has no literature-verified
        # residue list, so a blank panel is never read as "no residue is in the
        # site" and no residue list is written from memory.
        "catalytic_residues_note": (
            None
            if tier != "curated" or entry.get("catalytic_residues")
            else "No literature-verified catalytic-residue list is curated for this entry."
        ),
        "note": analysis.get("note"),
    })
    return base


def fetch_complex_interactions(complex_id: str, manifest_path: Path) -> dict[str, Any]:
    structure = fetch_complex_structure(complex_id, manifest_path)
    if structure.get("status") != "computed":
        return {
            "complex_id": complex_id,
            "hydrogen_bonds": None,
            "hydrophobic_contacts": None,
            "salt_bridges": None,
            "interacting_residues": structure.get("interacting_residues") or [],
            "status": structure.get("status", "placeholder"),
            "note": structure.get("note"),
        }

    residues = structure.get("interacting_residues") or []
    return {
        "complex_id": complex_id,
        "ligand_code": structure.get("ligand_code"),
        # Each list holds unique residues, not atom pairs.
        "hydrogen_bonds": [r for r in residues if "hydrogen_bond" in r.get("types", [])],
        "hydrophobic_contacts": [r for r in residues if "hydrophobic" in r.get("types", [])],
        "salt_bridges": [r for r in residues if "salt_bridge" in r.get("types", [])],
        "close_contacts": [r for r in residues if "close_contact" in r.get("types", [])],
        "contact_summary": structure.get("contact_summary") or {},
        "contact_summary_basis": "unique_residues",
        "contact_rules": structure.get("contact_rules") or {},
        "contact_labels": structure.get("contact_labels") or {},
        "near_miss_hydrogen_bonds": structure.get("near_miss_hydrogen_bonds") or [],
        "atom_pair_counts": structure.get("atom_pair_counts") or {},
        "atom_pairs_truncated": structure.get("atom_pairs_truncated", False),
        "atom_pairs_truncated_note": structure.get("atom_pairs_truncated_note"),
        "interacting_residues": residues,
        "binding_site_residues": structure.get("binding_site_residues", []),
        "method": structure.get("method"),
        "status": "computed",
        "note": structure.get("note"),
    }


def fetch_entry_metadata(pdb_id: str) -> dict[str, Any] | None:
    """RCSB metadata for an entry: title, resolution and the experimental method.

    Cached like every other upstream call. ``None`` means the metadata could not
    be fetched; the coordinates may still be usable, so the caller decides rather
    than failing the whole request.
    """
    url = RCSB_DATA_ENTRY_URL.format(pdb_id=pdb_id)
    try:
        entry = http_get_json_cached(url, source="pdb", backoff_seconds=2.0)
    except Exception as exc:  # noqa: BLE001
        logger.warning("RCSB metadata fetch failed for %s: %s", pdb_id, exc)
        return None
    if not isinstance(entry, dict):
        return None
    info = entry.get("rcsb_entry_info") or {}
    resolution = info.get("resolution_combined") or []
    methods = entry.get("exptl") or []
    return {
        "title": (entry.get("struct") or {}).get("title"),
        "resolution_angstrom": resolution[0] if resolution else None,
        "method": methods[0].get("method") if methods else None,
        "deposited_atom_count": info.get("deposited_atom_count"),
        "nonpolymer_entity_ids": (
            entry.get("rcsb_entry_container_identifiers") or {}
        ).get("non_polymer_entity_ids"),
    }


def _explore_payload(
    text: str,
    pdb_id: str,
    ligand_code: str,
    metadata: dict[str, Any] | None,
    candidates: list[dict[str, Any]],
) -> dict[str, Any]:
    """The measured-contact payload for an uncurated entry, same shape as curated."""
    analysis = analyse_complex(text, ligand_code)
    viewer_pdb = build_viewer_pdb(text, ligand_code) if analysis.get("status") == "computed" else None
    title = (metadata or {}).get("title")
    return {
        "complex_id": f"explore:{pdb_id}:{ligand_code}",
        "display_name": title or pdb_id,
        "protein_pdb": pdb_id,
        "protein_title": title,
        "ligand_name": None,
        "ligand_code": ligand_code,
        "complex_type": None,
        "target_origin": None,
        "category": None,
        "tier": "contacts_only",
        "curation": "none",
        "curation_note": CURATION_NONE_NOTE,
        "resolution_angstrom": (metadata or {}).get("resolution_angstrom"),
        "source_organism": None,
        "binding_site_type": None,
        "ligand_candidates": candidates,
        "ligand_selection_required": False,
        "ligand_rule": EXPLORE_LIGAND_RULE,
        "status": analysis.get("status", "placeholder"),
        "structure_file_reference": f"{pdb_id}.pdb",
        "structure_file_bytes": viewer_pdb,
        "protein_atom_count": analysis.get("protein_atom_count"),
        "ligand_atom_count": analysis.get("ligand_atom_count"),
        "binding_site_residue_count": analysis.get("binding_site_residue_count"),
        "contact_summary": analysis.get("contact_summary"),
        "contact_summary_basis": analysis.get("contact_summary_basis"),
        "contact_rules": analysis.get("contact_rules"),
        "contact_labels": analysis.get("contact_labels"),
        "near_miss_hydrogen_bonds": analysis.get("near_miss_hydrogen_bonds") or [],
        "atom_pair_counts": analysis.get("atom_pair_counts"),
        "atom_pairs_truncated": analysis.get("atom_pairs_truncated", False),
        "atom_pairs_truncated_note": analysis.get("atom_pairs_truncated_note"),
        "interacting_residues": analysis.get("interacting_residues", []),
        "binding_site_residues": analysis.get("binding_site_residues", []),
        "method": analysis.get("method"),
        "mechanism": {
            "status": "not_curated",
            "steps": [],
            "catalytic_residues": [],
            "sources": [],
            "note": NO_MECHANISM_NOTE,
        },
        "catalytic_residues_note": None,
        "citations": [],
        "sources": [],
        "note": analysis.get("note"),
    }


def explore_entry(pdb_id: str, ligand_code: str | None = None) -> dict[str, Any]:
    """Measure contacts for any deposited PDB entry, with no curation at all.

    Validation order is deliberate: the identifier format is checked before any
    network call, the size caps before any analysis, and the ligand rule before
    the contact scan. A request that cannot be answered honestly raises
    ``ComplexRequestError`` with the standard error shape; nothing is guessed and
    no placeholder number is produced.
    """
    cleaned = (pdb_id or "").strip().upper()
    if not PDB_ID_PATTERN.match(cleaned):
        raise ComplexRequestError(
            "invalid_pdb_id",
            f"{pdb_id!r} is not a PDB identifier: expected four characters, the first a "
            "digit, then letters or digits (for example 1M17).",
            422,
        )

    text = _fetch_pdb_text(cleaned)
    if text is None:
        raise ComplexRequestError(
            "structure_unavailable",
            f"Could not download the deposited coordinates for {cleaned} (RCSB PDB "
            "unreachable, or the entry has no PDB-format file). Nothing was invented.",
            503,
        )

    if len(text.encode("utf-8")) > MAX_ENTRY_BYTES:
        raise ComplexRequestError(
            "entry_too_large",
            f"{cleaned} is larger than the {MAX_ENTRY_BYTES // (1024 * 1024)} MB cap for "
            "this endpoint.",
            422,
        )

    protein_atoms = sum(1 for atom in parse_pdb(text) if atom.record == "ATOM")
    if protein_atoms > MAX_PROTEIN_ATOMS:
        raise ComplexRequestError(
            "entry_too_large",
            f"{cleaned} has {protein_atoms} protein atoms, above the "
            f"{MAX_PROTEIN_ATOMS}-atom cap for an interactive contact scan.",
            422,
        )

    candidates = drug_like_ligands(text)
    if not candidates:
        raise ComplexRequestError(
            "no_drug_like_ligand",
            f"No drug-like ligand was found in {cleaned}: after excluding water, ions, "
            "sugars, cofactors, buffers/additives and standard amino acids, no "
            "non-polymer group has at least 10 heavy atoms.",
            404,
        )

    metadata = fetch_entry_metadata(cleaned)
    available = {candidate["ligand_code"]: candidate for candidate in candidates}
    requested = (ligand_code or "").strip().upper() or None

    if requested is not None and requested not in available:
        raise ComplexRequestError(
            "ligand_not_found",
            f"{requested} is not a drug-like ligand of {cleaned}. Candidates: "
            + ", ".join(sorted(available))
            + ".",
            404,
        )

    if requested is None and len(candidates) > 1:
        return {
            "complex_id": f"explore:{cleaned}",
            "display_name": (metadata or {}).get("title") or cleaned,
            "protein_pdb": cleaned,
            "resolution_angstrom": (metadata or {}).get("resolution_angstrom"),
            "status": "ligand_selection_required",
            "tier": "contacts_only",
            "curation": "none",
            "curation_note": CURATION_NONE_NOTE,
            "ligand_selection_required": True,
            "ligand_candidates": candidates,
            "ligand_rule": EXPLORE_LIGAND_RULE,
            "note": (
                f"{cleaned} contains {len(candidates)} drug-like non-polymer groups. "
                "Choose one; no ligand is assumed."
            ),
        }

    chosen = requested or candidates[0]["ligand_code"]
    return _explore_payload(text, cleaned, chosen, metadata, candidates)


def fetch_complex_docking_validation(complex_id: str, manifest_path: Path) -> dict[str, Any]:
    try:
        entry = _find_complex(manifest_path, complex_id)
    except KeyError:
        return {
            "complex_id": complex_id,
            "status": "placeholder",
            "note": f"No curated entry for complex {complex_id!r}",
        }

    return {
        "complex_id": complex_id,
        "crystal_pose_file": None,
        "docked_pose_file": None,
        "rmsd_angstrom": None,
        "status": "placeholder",
        "note": (
            "No docking run has been performed for this complex, so there is no RMSD to report. "
            "The Drug-Target Studio instead reports measured contacts from the deposited crystal "
            "coordinates, which require no docking."
        ),
    }
