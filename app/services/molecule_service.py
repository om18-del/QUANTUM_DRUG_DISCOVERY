from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any
from urllib.parse import quote

from pydantic import BaseModel, Field, field_validator

from app.config import get_settings
from app.services.external_api import NotFoundError, http_get_bytes_cached, http_get_json_cached

logger = logging.getLogger(__name__)

PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"


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


def pubchem_sdf_bytes(cid: int) -> bytes:
    url = f"{PUBCHEM_BASE}/compound/cid/{cid}/SDF"
    return http_get_bytes_cached(url, source="pubchem", backoff_seconds=2.0)


class MoleculeSearchResult(BaseModel):
    query: str
    pubchem_cid: int | None = None
    smiles: str | None = None
    canonical_smiles: str | None = None
    molecular_weight: float | None = None
    molecular_formula: str | None = None
    iupac_name: str | None = None
    sdf_bytes: bytes | None = None
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

    sdf_bytes: bytes | None = None
    sdf_loaded = False
    try:
        sdf_bytes = pubchem_sdf_bytes(cid)
        sdf_loaded = len(sdf_bytes) > 0
    except Exception as exc:
        logger.warning("PubChem SDF fetch failed for cid %d: %s", cid, exc)

    return MoleculeSearchResult(
        query=q,
        pubchem_cid=cid,
        smiles=smiles,
        canonical_smiles=smiles,
        molecular_weight=molecular_weight,
        molecular_formula=molecular_formula,
        iupac_name=iupac_name,
        sdf_bytes=sdf_bytes,
        sdf_loaded=sdf_loaded,
        source="pubchem",
        status="computed",
        note=None,
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
