from __future__ import annotations

from http import HTTPStatus
from typing import Union

from fastapi import APIRouter, HTTPException, Query
from fastapi.responses import JSONResponse

from app.config import get_settings
from app.services.complexes_service import (
    ComplexRequestError,
    explore_entry,
    fetch_complex_interactions,
    fetch_complex_docking_validation,
    fetch_complex_structure,
    list_complexes,
)

router = APIRouter(prefix="/api/complexes", tags=["complexes"])

ComplexesResponse = Union[dict, JSONResponse]


def _error_response(code: str, message: str, status_code: int) -> JSONResponse:
    return JSONResponse(status_code=status_code, content={"error": {"code": code, "message": message}})


@router.get("/list", response_model=None)
def get_complexes_list() -> ComplexesResponse:
    settings = get_settings()
    try:
        return list_complexes(settings.complexes_manifest_path)
    except FileNotFoundError:
        return _error_response(
            "manifest_unavailable",
            "Complexes manifest is unavailable",
            int(HTTPStatus.SERVICE_UNAVAILABLE),
        )


@router.get("/explore/{pdb_id}", response_model=None)
def get_explore_entry(
    pdb_id: str,
    ligand: str | None = Query(None, max_length=8),
) -> ComplexesResponse:
    """Measure contacts for any deposited PDB entry, with no curation.

    Errors use the standard envelope: ``invalid_pdb_id`` and ``entry_too_large``
    are 422s, ``no_drug_like_ligand`` and ``ligand_not_found`` are 404s, and an
    unreachable RCSB is a 503. Several drug-like ligands without ``?ligand=``
    returns the candidate list instead of choosing one.
    """
    try:
        return explore_entry(pdb_id, ligand)
    except ComplexRequestError as exc:
        return _error_response(exc.code, exc.message, exc.status_code)


@router.get("/{complex_id}/interactions", response_model=None)
def get_complex_interactions(complex_id: str) -> ComplexesResponse:
    settings = get_settings()
    try:
        return fetch_complex_interactions(complex_id, settings.complexes_manifest_path)
    except FileNotFoundError:
        return _error_response(
            "manifest_unavailable",
            "Complexes manifest is unavailable",
            int(HTTPStatus.SERVICE_UNAVAILABLE),
        )


@router.get("/{complex_id}/structure", response_model=None)
def get_complex_structure(complex_id: str) -> ComplexesResponse:
    """Deposited coordinates trimmed to protein + target ligand, plus the
    computed interactions and the resulting mechanism of action."""
    settings = get_settings()
    try:
        return fetch_complex_structure(complex_id, settings.complexes_manifest_path)
    except FileNotFoundError:
        return _error_response(
            "manifest_unavailable",
            "Complexes manifest is unavailable",
            int(HTTPStatus.SERVICE_UNAVAILABLE),
        )


@router.get("/{complex_id}/docking-validation", response_model=None)
def get_complex_docking_validation(complex_id: str) -> ComplexesResponse:
    settings = get_settings()
    try:
        return fetch_complex_docking_validation(complex_id, settings.complexes_manifest_path)
    except FileNotFoundError:
        return _error_response(
            "manifest_unavailable",
            "Complexes manifest is unavailable",
            int(HTTPStatus.SERVICE_UNAVAILABLE),
        )
