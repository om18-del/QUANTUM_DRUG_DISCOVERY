from __future__ import annotations

import json
import logging
import re
from pathlib import Path
from typing import Any, Literal
from urllib.parse import quote

from pydantic import BaseModel, Field, field_validator

from app.config import get_settings
from app.services.external_api import NotFoundError, http_get_bytes_cached, http_get_json_cached

logger = logging.getLogger(__name__)

PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

# The molfile dimension code is the two trailing characters of the program line
# (line index 1), for example "  -OEChem-10122600532D". PubChem and other writers
# put it there, so we look there first and fall back to the rest of the header.
_SDF_DIMENSION_RE = re.compile(r"(3D|2D)\s*$")
# MDL V2000 atom line: x is columns 1-10, y 11-20 and z 21-30 (0-indexed 20:30).
_SDF_Z_COLUMNS = slice(20, 30)
_SDF_Z_SPAN_MIN = 1e-4  # Å; below this, the z spread is effectively zero (flat record).


def _sdf_header_dimension(text: str) -> str | None:
    """The ``2D``/``3D`` marker from the molfile header, or ``None`` when absent."""
    lines = text.splitlines()
    for idx in (1, 0, 2, 3, 4, 5):
        if idx >= len(lines):
            continue
        match = _SDF_DIMENSION_RE.search(lines[idx].rstrip())
        if match:
            return match.group(1).upper()
    return None


def _sdf_atom_count(text: str) -> int:
    """Number of atom lines in the first SDF block."""
    lines = text.splitlines()
    if len(lines) < 5:
        return 0
    try:
        return max(0, int(lines[3].strip()[:3]))
    except Exception:
        return 0


def _sdf_z_extent(text: str, n_atoms: int) -> tuple[float | None, float | None]:
    """Min/max z coordinate over the first SDF block's atom lines."""
    if n_atoms <= 0:
        return None, None
    lines = text.splitlines()
    zs: list[float] = []
    for i in range(n_atoms):
        line = lines[4 + i] if 4 + i < len(lines) else ""
        try:
            z = float(line[_SDF_Z_COLUMNS].strip() or "0")
        except Exception:
            z = 0.0
        zs.append(z)
    if not zs:
        return None, None
    return min(zs), max(zs)



class CompoundNotFound(Exception):
    """PubChem has no record for this name (HTTP 404, or an empty CID list).

    Carries only the query, never the upstream URL, so the message built from it
    can be shown to a person without leaking where the answer came from.
    """

    def __init__(self, query: str):
        self.query = query
        super().__init__(query)


def compound_not_found_message(query: str) -> str:
    """The one user-facing wording for a name the registry does not know."""
    return f"No compound found for '{query}'. Check the spelling."


def normalize_query(query: str) -> str:
    return query.strip()


def pubchem_cid_from_name(name: str, timeout_seconds: float) -> int:
    # PubChem's PUG REST now requires an explicit output format in the URL;
    # the bare "/cid" form returns 400 "Output format unspecified or invalid".
    url = f"{PUBCHEM_BASE}/compound/name/{quote(name, safe='')}/cids/JSON"
    try:
        data = http_get_json_cached(url, source="pubchem", backoff_seconds=2.0)
    except NotFoundError as exc:
        # 404 means the name simply is not in the registry. Convert it into a
        # query-bearing error so no URL reaches the caller.
        raise CompoundNotFound(name) from exc
    cids = data.get("IdentifierList", {}).get("CID", [])
    if not cids:
        # Some 200 responses carry an empty identifier list; same meaning.
        raise CompoundNotFound(name)
    if len(cids) != 1:
        raise KeyError("ambiguous pubchem name")
    return int(cids[0])


def pubchem_properties(cid: int) -> dict[str, Any]:
    props = ("MolecularWeight", "MolecularFormula", "CanonicalSMILES", "IsomericSMILES", "IUPACName")
    url = f"{PUBCHEM_BASE}/compound/cid/{cid}/property/{','.join(props)}/JSON"
    data = http_get_json_cached(url, source="pubchem", backoff_seconds=2.0)
    props_list = data.get("PropertyTable", {}).get("Properties", [{}])
    if not props_list:
        raise KeyError("no pubchem property data")
    return dict(props_list[0])


def _detect_sdf_dimension(sdf_bytes: bytes) -> Literal["3D", "2D"]:
    """Detect 3D vs 2D from the SDF record itself (header marker + real z spread).

    PubChem's SDF header line carries a ``3D``/``2D`` marker, and we confirm it
    against the actual coordinate spread: a header that says ``2D``, or a ``3D``
    header whose atom block has an effectively zero z span, is reported as 2D so a
    flat record is never drawn as a 3D structure. Both signals must agree on 3D.
    """
    if not sdf_bytes:
        return "2D"
    try:
        text = sdf_bytes.decode("utf-8", errors="replace")
    except Exception:
        return "2D"

    marker = _sdf_header_dimension(text)
    n_atoms = _sdf_atom_count(text)
    min_z, max_z = _sdf_z_extent(text, n_atoms)
    z_spread = (max_z - min_z) if (min_z is not None and max_z is not None) else 0.0

    if marker == "3D" and z_spread > _SDF_Z_SPAN_MIN:
        return "3D"
    return "2D"


def fetch_pubchem_sdf(cid: int) -> dict[str, Any]:
    """Fetch PubChem SDF for *cid*, preferring a real 3D conformer.

    PubChem's default ``/compound/cid/{cid}/SDF`` endpoint frequently returns a flat
    2D record (all z == 0). The ``?record_type=3d`` variant returns a real 3D
    conformer when one is deposited. We request the 3D record first and fall back to
    the default record when PubChem has no 3D record for that CID.

    Returns a dict with ``sdf_bytes``, ``sdf_dimension`` ("3D"/"2D"), and
    ``note`` describing what was actually returned.
    """
    sdf_3d_url = f"{PUBCHEM_BASE}/compound/cid/{cid}/SDF?record_type=3d"
    sdf_default_url = f"{PUBCHEM_BASE}/compound/cid/{cid}/SDF"

    # 1) Ask PubChem for the deposited 3D conformer first.
    candidate: bytes | None = None
    try:
        candidate = http_get_bytes_cached(sdf_3d_url, source="pubchem", backoff_seconds=2.0)
    except NotFoundError:
        # 404/410: PubChem has no 3D record for this CID. That is an answer, so
        # fall back to the default record rather than treating it as a failure.
        logger.info("PubChem has no 3D record for CID %d (HTTP 404); using the default record", cid)
    except Exception as exc:
        logger.warning("PubChem 3D SDF request failed for CID %d: %s", cid, exc)

    used_3d_request = False
    sdf_bytes: bytes | None = None
    if candidate:
        # 2) Only keep it when the record itself is 3D. A record whose header says
        #    2D (or whose atom block is flat) is not a 3D conformer.
        if _detect_sdf_dimension(candidate) == "3D":
            sdf_bytes = candidate
            used_3d_request = True
        else:
            logger.info(
                "PubChem 3D request for CID %d returned a 2D record; using the default record",
                cid,
            )

    # 3) Fall back to the default record when the 3D one is missing or flat.
    if sdf_bytes is None:
        try:
            sdf_bytes = http_get_bytes_cached(
                sdf_default_url, source="pubchem", backoff_seconds=2.0
            )
        except Exception as exc:
            raise RuntimeError(f"PubChem SDF unavailable for CID {cid}") from exc

    if not sdf_bytes:
        raise RuntimeError(f"PubChem SDF returned an empty body for CID {cid}")

    sdf_dimension = _detect_sdf_dimension(sdf_bytes)
    if sdf_dimension == "3D":
        note = "3D conformer from PubChem"
    elif used_3d_request:
        note = "PubChem 3D endpoint returned a flat record; treated as 2D"
    else:
        note = "PubChem returned a 2D record; no 3D conformer available"

    return {
        "sdf_bytes": sdf_bytes,
        "sdf_dimension": sdf_dimension,
        "note": note,
    }


class MoleculeSearchResult(BaseModel):
    query: str
    pubchem_cid: int | None = None
    smiles: str | None = None
    canonical_smiles: str | None = None
    molecular_weight: float | None = None
    molecular_formula: str | None = None
    iupac_name: str | None = None
    sdf_bytes: bytes | None = None
    sdf_dimension: Literal["3D", "2D"] | None = None
    sdf_loaded: bool = False
    source: str = "pubchem"
    status: str = "computed"
    note: str | None = None

    @field_validator("sdf_bytes", mode="before")
    @classmethod
    def _coerce_sdf(cls, v: Any) -> bytes | None:
        if isinstance(v, bytes):
            return v
        return None


def search_molecule(query: str, timeout_seconds: float) -> dict[str, Any]:
    q = normalize_query(query)
    if not q:
        return {"query": q, "status": "placeholder", "note": "Empty query"}

    try:
        cid = pubchem_cid_from_name(q, timeout_seconds)
    except CompoundNotFound:
        # A real "not found" must reach the router intact so it can be reported
        # with its own code, instead of being flattened into a generic failure.
        raise
    except Exception as exc:
        logger.warning("PubChem name search failed for %r: %s", q, exc)
        return {
            "query": q,
            "pubchem_cid": None,
            "smiles": None,
            "canonical_smiles": None,
            "molecular_weight": None,
            "molecular_formula": None,
            "iupac_name": None,
            "sdf_bytes": None,
            "sdf_loaded": False,
            "status": "placeholder",
            "note": f"PubChem lookup failed: {exc}",
        }

    try:
        props = pubchem_properties(cid)
    except Exception as exc:
        logger.warning("PubChem properties failed for cid %d: %s", cid, exc)
        return {
            "query": q,
            "pubchem_cid": cid,
            "smiles": None,
            "canonical_smiles": None,
            "molecular_weight": None,
            "molecular_formula": None,
            "iupac_name": None,
            "sdf_bytes": None,
            "sdf_loaded": False,
            "status": "placeholder",
            "note": f"PubChem properties failed: {exc}",
        }

    # PubChem returns the SMILES under "SMILES" (with "ConnectivitySMILES");
    # "CanonicalSMILES"/"IsomericSMILES" are no longer returned.
    smiles = (
        props.get("SMILES")
        or props.get("ConnectivitySMILES")
        or props.get("CanonicalSMILES")
        or props.get("IsomericSMILES")
    )
    molecular_weight = props.get("MolecularWeight")
    molecular_formula = props.get("MolecularFormula")
    iupac_name = props.get("IUPACName")

    sdf_result: dict[str, Any] | None = None
    sdf_error: Exception | None = None
    try:
        sdf_result = fetch_pubchem_sdf(cid)
    except Exception as exc:
        logger.warning("PubChem SDF fetch failed for cid %d: %s", cid, exc)
        sdf_error = exc

    sdf_bytes = None
    sdf_dimension = "2D"
    sdf_note: str | None = None
    sdf_loaded = False
    if sdf_result is not None:
        sdf_bytes = sdf_result.get("sdf_bytes") or None
        sdf_dimension = str(sdf_result.get("sdf_dimension") or "2D")
        sdf_note = sdf_result.get("note") or None
        sdf_loaded = bool(sdf_bytes and len(sdf_bytes) > 0)

    return MoleculeSearchResult(
        query=q,
        pubchem_cid=cid,
        smiles=smiles,
        canonical_smiles=smiles,
        molecular_weight=molecular_weight,
        molecular_formula=molecular_formula,
        iupac_name=iupac_name,
        sdf_bytes=sdf_bytes,
        sdf_dimension=sdf_dimension,
        sdf_loaded=sdf_loaded,
        source="pubchem",
        status="computed",
        note=sdf_note,
    ).model_dump()


class MoleculePropertiesResult(BaseModel):
    smiles: str
    molecular_weight: float | None = None
    logp: float | None = None
    h_bond_donors: int | None = None
    h_bond_acceptors: int | None = None
    rotatable_bonds: int | None = None
    formula: str | None = None
    exact_mass: float | None = None
    tpsa: float | None = None
    lipinski_rule_of_five: dict[str, Any] = Field(default_factory=lambda: {"pass": None, "violations": None})
    status: str = "computed"
    note: str | None = None


def compute_molecule_properties(smiles: str, timeout_seconds: float) -> dict[str, Any]:
    from app.services.rdkit_helper import normalize_smiles_safe, rdkit_properties_from_smiles

    normalized = normalize_smiles_safe(smiles)
    if normalized is None:
        return MoleculePropertiesResult(
            smiles=smiles,
            status="placeholder",
            note="SMILES normalization failed",
        ).model_dump()

    try:
        props = rdkit_properties_from_smiles(normalized)
    except Exception as exc:  # noqa: BLE001
        logger.warning("RDKit properties failed for %r: %s", normalized, exc)
        return MoleculePropertiesResult(
            smiles=normalized,
            status="placeholder",
            note=f"RDKit property computation failed: {exc}",
        ).model_dump()

    return MoleculePropertiesResult(
        smiles=normalized,
        molecular_weight=props.get("molecular_weight"),
        logp=props.get("logp"),
        h_bond_donors=props.get("h_bond_donors"),
        h_bond_acceptors=props.get("h_bond_acceptors"),
        rotatable_bonds=props.get("rotatable_bonds"),
        formula=props.get("formula"),
        exact_mass=props.get("exact_mass"),
        tpsa=props.get("tpsa"),
        lipinski_rule_of_five=props.get("lipinski_rule_of_five", {"pass": None, "violations": None}),
        status="computed",
        note=None,
    ).model_dump()
