from __future__ import annotations

import logging
from pathlib import Path

from fastapi import FastAPI, HTTPException, Request
from fastapi.exceptions import RequestValidationError
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import HTMLResponse, JSONResponse, RedirectResponse
from starlette.exceptions import HTTPException as StarletteHTTPException
from starlette.responses import Response
from starlette.staticfiles import StaticFiles

from app.config import get_settings


class NoCacheStaticFiles(StaticFiles):
    """Serve frontend static files with Cache-Control: no-cache so browsers
    never keep a stale asset around between edits.

    This must cover CSS and HTML, not just JavaScript: without an explicit
    directive Chromium applies heuristic freshness to stylesheets and will keep
    serving an old styles.css after it has been rewritten on disk.
    """

    async def get_response(self, path: str, request: Request) -> Response:
        response = await super().get_response(path, request)
        if response.status_code == 200:
            response.headers["Cache-Control"] = "no-cache, must-revalidate, max-age=0"
            response.headers["Expires"] = "0"
        return response


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

    from app.routers import ai
    from app.routers import complexes
    from app.routers import jobs
    from app.routers import molecules
    from app.routers import proteins
    from app.routers import quantum

    app.include_router(quantum.router)
    app.include_router(molecules.router)
    app.include_router(proteins.router)
    app.include_router(complexes.router)
    app.include_router(jobs.router)
    app.include_router(ai.router)

    frontend_dir = Path(__file__).resolve().parents[1] / "frontend" / "static"

    if frontend_dir.is_dir():
        stylesheet_version = "10"

        async def _index(_request: Request) -> HTMLResponse:
            index = frontend_dir / "index.html"
            if not index.is_file():
                raise HTTPException(status_code=404, detail="Frontend index.html not found")
            content = index.read_text(encoding="utf-8")
            # Version the stylesheet URL so a browser holding an old cached copy
            # fetches the current file instead of reusing it.
            html = content.replace(
                "</head>",
                f'<link rel="stylesheet" href="/app/styles.css?v={stylesheet_version}" /></head>',
                1,
            )
            return HTMLResponse(content=html, status_code=200)

        @app.get("/", include_in_schema=False)
        async def _root(request: Request) -> HTMLResponse:
            return await _index(request)

        @app.get("/index.html", include_in_schema=False)
        async def _index_html(request: Request) -> HTMLResponse:
            return await _index(request)

        @app.get("/docs", include_in_schema=False)
        async def _redirect_docs() -> RedirectResponse:
            return RedirectResponse(url="/api/docs")

        app.mount(
            "/app",
            NoCacheStaticFiles(directory=str(frontend_dir), html=False),
            name="frontend_static",
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
