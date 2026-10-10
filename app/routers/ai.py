from __future__ import annotations

from typing import Union

from fastapi import APIRouter, HTTPException, Query
from fastapi.requests import Request
from fastapi.responses import JSONResponse

from app.config import get_settings
from app.services.ai_service import predict_properties, resolve_query
from app.services.rdkit_helper import normalize_smiles_safe

router = APIRouter(prefix="/api/ai", tags=["ai"])

AiResponse = Union[dict, JSONResponse]


def _error_response(code: str, message: str, status_code: int) -> JSONResponse:
    return JSONResponse(status_code=status_code, content={"error": {"code": code, "message": message}})


async def _extract_query_body(request: Request) -> dict | JSONResponse:
    try:
        body = await request.json()
    except Exception:
        return _error_response("invalid_body", "Request body must be valid JSON", 400)

    if not isinstance(body, dict):
        return _error_response("invalid_body", "Request body must be a JSON object", 400)

    query = body.get("query")
    if not isinstance(query, str) or not query.strip():
        return _error_response("invalid_query", "Field 'query' is required and must be a non-empty string", 422)

    return {"query": query}


@router.post("/resolve", response_model=None)
async def post_resolve(request: Request) -> AiResponse:
    try:
        body = await request.json()
    except Exception:
        return _error_response("invalid_body", "Request body must be valid JSON", 400)

    if not isinstance(body, dict):
        return _error_response("invalid_body", "Request body must be a JSON object", 400)

    body_or_error = await _extract_query_body(request)
    if isinstance(body_or_error, JSONResponse):
        return body_or_error

    return resolve_query(body_or_error["query"])


@router.post("/predict-properties", response_model=None)
async def post_predict_properties(request: Request) -> AiResponse:
    try:
        body = await request.json()
    except Exception:
        return _error_response("invalid_body", "Request body must be valid JSON", 400)

    if not isinstance(body, dict):
        return _error_response("invalid_body", "Request body must be a JSON object", 400)

    smiles = body.get("smiles")
    if not isinstance(smiles, str) or not smiles.strip():
        return _error_response("invalid_smiles", "Field 'smiles' is required and must be a non-empty string", 422)

    normalized = normalize_smiles_safe(smiles)
    if normalized is None:
        return _error_response("invalid_smiles", "SMILES could not be normalized", 422)

    settings = get_settings()
    return predict_properties(normalized, settings.models_manifest_path)
