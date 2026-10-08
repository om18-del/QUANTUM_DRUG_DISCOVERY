from __future__ import annotations

import logging

from fastapi import FastAPI, HTTPException, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from fastapi.exceptions import RequestValidationError
from starlette.exceptions import HTTPException as StarletteHTTPException

from app.config import get_settings

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def error_body(code: str, message: str) -> dict[str, dict[str, str]]:
    return {"error": {"code": code, "message": message}}


def _normalize_http_exception(exc: HTTPException) -> JSONResponse:
    if isinstance(exc.detail, str):
        detail = exc.detail
    else:
        detail = str(exc.detail)

    if exc.status_code == 404:
        return JSONResponse(status_code=404, content=error_body("not_found", detail))
    if exc.status_code == 405:
        return JSONResponse(status_code=405, content=error_body("method_not_allowed", detail))
    return JSONResponse(status_code=exc.status_code, content=error_body("http_error", detail))


def create_app() -> FastAPI:
    settings = get_settings()

    app = FastAPI(
        title=settings.app_name,
        openapi_url=f"{settings.api_prefix}/openapi.json",
        docs_url=f"{settings.api_prefix}/docs",
        redoc_url=f"{settings.api_prefix}/redoc",
    )

    app.add_middleware(
        CORSMiddleware,
        allow_origins=settings.cors_origin_list,
        allow_credentials=settings.cors_allow_credentials,
        allow_methods=["*"],
        allow_headers=["*"],
    )

    @app.exception_handler(RequestValidationError)
    async def request_validation_exception_handler(
        request: Request, exc: RequestValidationError
    ) -> JSONResponse:
        logger.warning("Validation error: %s", exc)
        message = "; ".join(
            f"{err.get('loc', [])} {err.get('msg', '')}" for err in exc.errors()
        ) or "Invalid request"
        return JSONResponse(status_code=422, content=error_body("validation_error", message))

    @app.exception_handler(StarletteHTTPException)
    async def starlette_http_exception_handler(
        request: Request, exc: StarletteHTTPException
    ) -> JSONResponse:
        return _normalize_http_exception(exc)

    @app.exception_handler(Exception)
    async def generic_exception_handler(request: Request, exc: Exception) -> JSONResponse:
        logger.exception("Unexpected error: %s", exc)
        return JSONResponse(
            status_code=500, content=error_body("internal_error", "An unexpected error occurred")
        )

    @app.get("/health", tags=["system"])
    def health() -> dict[str, str]:
        return {"status": "ok", "app": settings.app_name}

    return app


app = create_app()
