from __future__ import annotations

from pathlib import Path

from app.services.rscb_service import fetch_protein as rcsb_fetch_protein

LOCAL_PDB_DIR = Path(__file__).resolve().parents[1] / "local_pdb"


def fetch_protein(pdb_id: str, timeout_seconds: float) -> dict:
    """Fetch protein metadata + optional local structure.

    Tries RCSB first. When a matching local ``.cif``/``.pdb`` file exists under
    ``local_pdb/``, its bytes and the resolved reference name are merged into the
    result so the frontend can render a local structure when RCSB returns a
    placeholder or is unavailable.
    """
    result = rcsb_fetch_protein(pdb_id, timeout_seconds)

    local_path = _resolve_local_pdb(pdb_id)
    if local_path is not None:
        try:
            result["structure_file_reference"] = local_path.name
            result["structure_file_bytes"] = local_path.read_bytes()
            if result.get("status") == "placeholder":
                result["status"] = "computed"
                result["note"] = result.get("note") or "RCSB unavailable; using local structure"
        except Exception as exc:
            pass

    return result


def _resolve_local_pdb(pdb_id: str) -> Path | None:
    """Return the first matching local structure file for *pdb_id*, if any.

    Checks ``local_pdb/{pdb_id}.cif`` then ``local_pdb/{pdb_id}.pdb`` (case-insensitive
    on the directory listing when the exact name is not found).
    """
    if not LOCAL_PDB_DIR.is_dir():
        return None

    candidates = [
        LOCAL_PDB_DIR / f"{pdb_id}.cif",
        LOCAL_PDB_DIR / f"{pdb_id}.pdb",
    ]
    for candidate in candidates:
        if candidate.is_file():
            return candidate

    lowered = pdb_id.lower()
    try:
        for entry in LOCAL_PDB_DIR.iterdir():
            if entry.is_file() and entry.name.lower().startswith(lowered):
                stem = entry.name.rsplit(".", 1)[0]
                if stem.lower() == lowered:
                    return entry
    except OSError:
        return None

    return None
