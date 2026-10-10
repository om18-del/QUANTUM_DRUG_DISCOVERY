from __future__ import annotations

import re
from http import HTTPStatus
from typing import Union

from fastapi import APIRouter, HTTPException, Query, Request
from fastapi.responses import JSONResponse
from app.config import get_settings
from app.services.protein_service import fetch_protein

router = APIRouter(prefix="/api/protein", tags=["protein"])

ProteinResponse = Union[dict, JSONResponse]
_PDB_ID_RE = re.compile(r"^[0-9][A-Za-z0-9]{3}$")


def _error_response(code: str, message: str, status_code: int) -> JSONResponse:
    return JSONResponse(status_code=status_code, content={"error": {"code": code, "message": message}})


def _normalize_pdb_id(raw: str) -> Union[str, JSONResponse]:
    value = raw.strip().upper()
    if not _PDB_ID_RE.match(value):
        return _error_response(
            "invalid_pdb_id",
            f"Invalid PDB id: {raw!r}. Expected 4 characters, starting with a digit.",
            400,
        )
    return value


@router.get("/{pdb_id}", response_model=None)
def get_protein(pdb_id: str, request: Request) -> ProteinResponse:
    normalized = _normalize_pdb_id(pdb_id)
    if isinstance(normalized, JSONResponse):
        return normalized

    pdb_id = normalized
    settings = get_settings()
    try:
        return fetch_protein(pdb_id, timeout_seconds=settings.external_timeout_seconds)
    except Exception as exc:
        return _error_response(
            "protein_lookup_error",
            f"Protein lookup failed: {exc}",
            502,
        )
