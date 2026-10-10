from __future__ import annotations

import logging
from typing import Union

from fastapi import APIRouter, HTTPException, Path, Query
from fastapi.requests import Request
from fastapi.responses import JSONResponse
from pydantic import ValidationError
from app.config import get_settings
from app.services.molecule_service import (
    CompoundNotFound,
    compound_not_found_message,
    compute_molecule_properties,
    search_molecule,
)
from app.services.molecule_simulation import simulate_molecule

router = APIRouter(prefix="/api/molecule", tags=["molecule"])
logger = logging.getLogger(__name__)

MoleculeResponse = Union[dict, JSONResponse]


def _error_response(code: str, message: str, status_code: int) -> JSONResponse:
    return JSONResponse(status_code=status_code, content={"error": {"code": code, "message": message}})


@router.get("/search", response_model=None)
def molecule_search(q: str = Query(..., min_length=1, max_length=200, description="Search query")) -> MoleculeResponse:
    settings = get_settings()
    try:
        return search_molecule(q, settings.external_timeout_seconds)
    except CompoundNotFound:
        # Deliberately worded once, in molecule_service, and deliberately free of
        # the upstream URL: the person searched for a compound, not for a link.
        return _error_response("compound_not_found", compound_not_found_message(q), 404)
    except Exception:
        # The raw exception can embed an upstream URL, so it stays in the log.
        logger.warning("Molecule search failed for %r", q, exc_info=True)
        return _error_response(
            "search_error",
            f"Molecule search for '{q}' is unavailable right now. Try again shortly.",
            502,
        )


@router.post("/properties", response_model=None)
async def molecule_properties(request: Request) -> MoleculeResponse:
    settings = get_settings()
    try:
        body = await request.json()
    except Exception:
        return _error_response("invalid_body", "Request body must be valid JSON", 400)

    if not isinstance(body, dict):
        return _error_response("invalid_body", "Request body must be a JSON object", 400)

    smiles = body.get("smiles")
    if not isinstance(smiles, str) or not smiles.strip():
        return _error_response("invalid_smiles", "Field 'smiles' is required and must be a non-empty string", 422)

    if len(smiles) > settings.max_smiles_length:
        return _error_response("invalid_smiles", f"SMILES too long (max {settings.max_smiles_length} chars)", 422)

    try:
        return compute_molecule_properties(smiles, settings.external_timeout_seconds)
    except Exception as exc:
        return _error_response("properties_error", f"Molecule property computation failed: {exc}", 500)


@router.post("/simulate", response_model=None)
async def molecule_simulate(request: Request) -> MoleculeResponse:
    """Real 3D conformer search plus a sourced drug-likeness interpretation."""
    settings = get_settings()
    try:
        body = await request.json()
    except Exception:
        return _error_response("invalid_body", "Request body must be valid JSON", 400)

    if not isinstance(body, dict):
        return _error_response("invalid_body", "Request body must be a JSON object", 400)

    smiles = body.get("smiles")
    if not isinstance(smiles, str) or not smiles.strip():
        return _error_response("invalid_smiles", "Field 'smiles' is required and must be a non-empty string", 422)

    if len(smiles) > settings.max_smiles_length:
        return _error_response("invalid_smiles", f"SMILES too long (max {settings.max_smiles_length} chars)", 422)

    try:
        num_conformers = int(body.get("num_conformers", 8))
    except (TypeError, ValueError):
        num_conformers = 8
    num_conformers = max(1, min(32, num_conformers))

    try:
        return simulate_molecule(smiles, num_conformers=num_conformers)
    except Exception as exc:
        return _error_response("simulation_error", f"Molecule simulation failed: {exc}", 500)
