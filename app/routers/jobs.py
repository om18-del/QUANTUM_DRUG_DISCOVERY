from __future__ import annotations

from http import HTTPStatus

from typing import Union

from fastapi import APIRouter, HTTPException
from fastapi.responses import JSONResponse

from app.config import get_settings
from app.services.jobs_service import get_job

router = APIRouter(prefix="/api/jobs", tags=["jobs"])

JobResponse = Union[dict, JSONResponse]


def _error_response(code: str, message: str, status_code: int) -> JSONResponse:
    return JSONResponse(status_code=status_code, content={"error": {"code": code, "message": message}})


@router.get("/{job_id}", response_model=None)
def get_job_status(job_id: str) -> JobResponse:
    settings = get_settings()
    try:
        return get_job(job_id, settings.job_store_path)
    except KeyError:
        return _error_response(
            "job_not_found",
            f"No job with id {job_id!r}",
            int(HTTPStatus.NOT_FOUND),
        )
