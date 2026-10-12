from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any

from app.config import get_settings
from app.services.external_api import http_get_bytes_cached
from app.services.structure_analysis import analyse_complex, build_mechanism, build_viewer_pdb

logger = logging.getLogger(__name__)

RSCB_PDB_DOWNLOAD_BASE = "https://files.rcsb.org/download"


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
            # Distinguishes the viral examples from the human-enzyme example, so
            # the UI never implies that a human enzyme is a virus.
            "target_origin": item.get("target_origin"),
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

    base: dict[str, Any] = {
        "complex_id": complex_id,
        "display_name": entry.get("display_name"),
        "protein_pdb": pdb_id,
        "protein_title": entry.get("protein_title"),
        "ligand_name": entry.get("ligand_name"),
        "ligand_code": ligand_code,
        "complex_type": entry.get("complex_type"),
        "target_origin": entry.get("target_origin"),
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
    mechanism = build_mechanism(analysis, entry)
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
