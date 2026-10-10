from __future__ import annotations

import logging
from typing import Any

from pydantic import BaseModel

from app.config import get_settings
from app.services.external_api import http_get_json_cached, http_get_bytes_cached

logger = logging.getLogger(__name__)

# New RCSB PDB REST API (data.rcsb.org) replaces the deprecated files.rcsb.org/json/{id} endpoint.
# The old https://files.rcsb.org/download/{pdb_id} also 404s; the working PDB download URL is
# https://files.rcsb.org/download/{pdb_id}.cif
RSCB_PDB_API_BASE = "https://data.rcsb.org/rest/v1"
RSCB_PDB_DOWNLOAD_BASE = "https://files.rcsb.org/download"


class ProteinResult(BaseModel):
    pdb_id: str
    title: str | None = None
    resolution: float | None = None
    organism: str | None = None
    chains: list[str] | None = None
    ligands: list[str] | None = None
    structure_file_reference: str = "placeholder.pdb"
    structure_file_bytes: bytes | None = None
    source: str = "rcsb"
    status: str = "computed"
    note: str | None = None


def pdb_entry_json(pdb_id: str) -> dict[str, Any]:
    """Fetch the core entry record from the new RCSB PDB REST API.

    The response IS the entry record itself (not wrapped in {"entry": ...}).
    """
    url = f"{RSCB_PDB_API_BASE}/core/entry/{pdb_id}"
    data = http_get_json_cached(url, source="pdb", backoff_seconds=2.0)
    return data


def _fetch_polymer_entity(pdb_id: str, entity_id: str) -> dict[str, Any]:
    """Fetch a single polymer entity detail record."""
    url = f"{RSCB_PDB_API_BASE}/core/polymer_entity/{pdb_id}/{entity_id}"
    return http_get_json_cached(url, source="pdb", backoff_seconds=2.0)


def pdb_entry_meta(entry: dict[str, Any], pdb_id: str) -> dict[str, Any]:
    """Extract title, resolution, and organism from the new API entry shape.

    - title: struct.title (fallback to rcsb_entry_info.struct_type_description when missing)
    - resolution: rcsb_entry_info.resolution_combined[0]
    - organism: rcsb_entity_source_organism[0].scientific_name from the first
      polymer entity (the protein's source organism). Falls back to the host
      organism when no source organism is available.
    """
    title = entry.get("struct", {}).get("title")
    if not title:
        title = entry.get("rcsb_entry_info", {}).get("struct_type_description")

    resolution = None
    res_combined = entry.get("rcsb_entry_info", {}).get("resolution_combined")
    if isinstance(res_combined, list) and res_combined:
        try:
            resolution = float(res_combined[0])
        except (TypeError, ValueError):
            resolution = None

    organism: str | None = None
    # Determine the source organism via the first polymer entity's taxonomy.
    polymer_entity_ids = entry.get("rcsb_entry_container_identifiers", {}).get(
        "polymer_entity_ids", []
    )
    if polymer_entity_ids:
        try:
            pe = _fetch_polymer_entity(pdb_id, polymer_entity_ids[0])
            src_orgs = pe.get("rcsb_entity_source_organism", [])
            if src_orgs:
                organism = (
                    src_orgs[0].get("common_name")
                    or src_orgs[0].get("scientific_name")
                )
            if not organism:
                host_orgs = pe.get("rcsb_entity_host_organism", [])
                if host_orgs:
                    organism = (
                        host_orgs[0].get("common_name")
                        or host_orgs[0].get("scientific_name")
                    )
        except Exception as exc:
            logger.warning(
                "Failed to fetch organism taxonomy for %s (entity %s): %s",
                pdb_id,
                polymer_entity_ids[0],
                exc,
            )

    return {
        "title": title,
        "resolution": resolution,
        "organism": organism,
    }


def pdb_entity_chains(entry: dict[str, Any], pdb_id: str) -> list[str]:
    """Derive chain identifiers from the new API's polymer entity records.

    Each polymer entity has ``entity_poly.pdbx_strand_id`` which is the PDB
    chain ID (e.g. ``"A"``, ``"B"``). We fetch each polymer entity referenced
    in ``rcsb_entry_container_identifiers.polymer_entity_ids`` and collect their
    strand IDs.
    """
    polymer_entity_ids = entry.get("rcsb_entry_container_identifiers", {}).get(
        "polymer_entity_ids", []
    )
    chains: list[str] = []
    for entity_id in polymer_entity_ids:
        try:
            pe = _fetch_polymer_entity(pdb_id, entity_id)
            strand_id = pe.get("entity_poly", {}).get("pdbx_strand_id")
            if strand_id:
                chains.append(str(strand_id))
        except Exception as exc:
            logger.warning(
                "Failed to fetch polymer entity %s for chain lookup: %s", entity_id, exc
            )
    chains.sort()
    return chains


def pdb_ligands(
    entry: dict[str, Any],
    pdb_id: str,
    structure_file_bytes: bytes | None = None,
) -> list[str]:
    """Return non-polymer small-molecule ligand names.

    Tries the new RCSB API ``/core/entity/{pdb_id}/{entity_id}`` endpoint for
    each non-polymer entity referenced in the entry. When the entity endpoint
    is unavailable (404), falls back to extracting ligand names directly from
    the downloaded structure file (``.cif`` format) by scanning the
    ``_chem_comp.name`` category.
    """
    import re

    from app.services.external_api import http_get_json_cached as _get

    entity_ids = entry.get("rcsb_entry_container_identifiers", {}).get(
        "entity_ids", []
    )
    polymer_entity_ids = set(
        entry.get("rcsb_entry_container_identifiers", {}).get(
            "polymer_entity_ids", []
        )
    )
    ligands: list[str] = []

    # Phase 1: try the per-entity API for every non-polymer entity.
    for entity_id in entity_ids:
        if entity_id in polymer_entity_ids:
            continue
        try:
            entity_data = _get(
                f"{RSCB_PDB_API_BASE}/core/entity/{pdb_id}/{entity_id}",
                source="pdb",
                backoff_seconds=2.0,
            )
            chem_comp = (
                entity_data.get("rcsb_entity_nonpolymer_entity", {})
                .get("chem_comp", {})
            )
            name = chem_comp.get("name") or chem_comp.get("id")
            if name:
                ligands.append(str(name))
        except Exception as exc:
            logger.debug(
                "Entity API failed for %s/%s (will try CIF fallback): %s",
                pdb_id,
                entity_id,
                exc,
            )

    # Phase 2: if the API produced nothing, parse the CIF content directly.
    if not ligands and structure_file_bytes:
        ligands = _extract_ligands_from_cif(structure_file_bytes)

    return sorted(ligands)


def _extract_ligands_from_cif(content_bytes: bytes) -> list[str]:
    """Parse a PDB CIF file and return unique small-molecule ligand names.

    Parses the ``loop_`` ``_chem_comp`` category block and returns the
    ``_chem_comp.name`` values for non-standard compounds (excluding water
    and standard L-peptide amino acids). Handles multi-line semicolon-quoted
    values and the column-aligned CIF data row format.
    """
    try:
        text = content_bytes.decode("utf-8", errors="replace")
    except UnicodeDecodeError:
        return []

    ligands: list[str] = []
    lines = text.splitlines()
    n = len(lines)
    i = 0

    while i < n:
        # Detect a ``loop_`` block that declares ``_chem_comp`` columns.
        if lines[i].strip() == "loop_":
            # Collect column declarations until the first data row.
            j = i + 1
            columns: list[str] = []
            while j < n and lines[j].strip().startswith("_chem_comp."):
                columns.append(lines[j].strip())
                j += 1

            if not columns:
                i += 1
                continue

            # Locate the _chem_comp.name column index.
            name_idx = None
            for idx, col in enumerate(columns):
                if col == "_chem_comp.name":
                    name_idx = idx
                    break
            if name_idx is None:
                i = j
                continue

            # Also record the _chem_comp.type and _chem_comp.mon_nstd_flag
            # column indices for filtering.
            type_idx = None
            flag_idx = None
            for idx, col in enumerate(columns):
                if col == "_chem_comp.type":
                    type_idx = idx
                elif col == "_chem_comp.mon_nstd_flag":
                    flag_idx = idx

            # Parse data rows until a blank/comment line or new keyword.
            k = j
            while k < n:
                row = lines[k].strip()
                if not row or row.startswith("#") or row.startswith("loop_"):
                    break
                # Split respecting single-quoted values (CIF quoting).
                parts: list[str] = []
                buf = ""
                in_quote = False
                for ch in row:
                    if ch == "'" and not in_quote:
                        in_quote = True
                        buf += ch
                    elif ch == "'" and in_quote:
                        in_quote = False
                        buf += ch
                    elif ch.isspace() and not in_quote:
                        if buf:
                            parts.append(buf)
                            buf = ""
                    else:
                        buf += ch
                if buf:
                    parts.append(buf)

                if len(parts) <= name_idx:
                    k += 1
                    continue

                raw_name = parts[name_idx]
                # Strip surrounding quotes if present.
                name = raw_name
                if len(name) >= 2 and name[0] == "'" and name[-1] == "'":
                    name = name[1:-1]

                # Filter out water and standard amino acids.
                comp_type = parts[type_idx].strip() if type_idx is not None and len(parts) > type_idx else ""
                if comp_type == "non-polymer" and name.upper() == "WATER":
                    k += 1
                    continue
                if "peptide linking" in comp_type or "L-peptide linking" in comp_type:
                    k += 1
                    continue

                if name and name not in ligands:
                    ligands.append(name)

                k += 1

            i = k
            continue

        i += 1

    return sorted(ligands)


def fetch_protein(pdb_id: str, timeout_seconds: float) -> dict[str, Any]:
    try:
        entry = pdb_entry_json(pdb_id)
    except Exception as exc:
        logger.warning("RCSB entry fetch failed for %s: %s", pdb_id, exc)
        return ProteinResult(
            pdb_id=pdb_id,
            status="placeholder",
            note=f"RCSB lookup failed: {exc}",
        ).model_dump()

    meta = pdb_entry_meta(entry, pdb_id)
    chains = pdb_entity_chains(entry, pdb_id)

    # Download the structure file BEFORE calling pdb_ligands, so the CIF
    # content can be used as a fallback when the entity API 404s.
    structure_file_bytes: bytes | None = None
    structure_file_reference = f"{pdb_id}.cif"
    try:
        url = f"{RSCB_PDB_DOWNLOAD_BASE}/{pdb_id}.cif"
        structure_file_bytes = http_get_bytes_cached(url, source="pdb", backoff_seconds=2.0)
        if len(structure_file_bytes) == 0:
            structure_file_bytes = None
    except Exception as exc:
        logger.warning("RCSB PDB file fetch failed for %s: %s", pdb_id, exc)
        structure_file_reference = "placeholder.pdb"

    ligands = pdb_ligands(entry, pdb_id, structure_file_bytes)

    return ProteinResult(
        pdb_id=pdb_id,
        title=meta["title"],
        resolution=meta["resolution"],
        organism=meta["organism"],
        chains=chains if chains else None,
        ligands=ligands if ligands else None,
        structure_file_reference=structure_file_reference,
        structure_file_bytes=structure_file_bytes,
        source="rcsb",
        status="computed" if meta["title"] is not None else "placeholder",
        note=None,
    ).model_dump()
