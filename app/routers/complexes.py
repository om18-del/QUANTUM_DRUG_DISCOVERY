from __future__ import annotations

from http import HTTPStatus
from typing import Union

from fastapi import APIRouter, HTTPException
from fastapi.responses import JSONResponse

from app.config import get_settings
from app.services.complexes_service import (
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
            int(HTTPStatus.service_unavailable),
        )


@router.get("/{complex_id}/interactions", response_model=None)
def get_complex_interactions(complex_id: str) -> ComplexesResponse:
    settings = get_settings()
    try:
        return fetch_complex_interactions(complex_id, settings.complexes_manifest_path)
    except FileNotFoundError:
        return _error_response(
            "manifest_unavailable",
            "Complexes manifest is unavailable",
            int(HTTPStatus.service_unavailable),
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
            int(HTTPStatus.service_unavailable),
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
            int(HTTPStatus.service_unavailable),
        )
