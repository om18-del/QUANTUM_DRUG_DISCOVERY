from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
from typing import Any

_JOB_STORE: dict[str, dict[str, Any]] = {}


def _now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def get_job(job_id: str, job_store_path: Path) -> dict[str, Any]:
    """Return job status.

    In this build the job store is in-memory only. The ``job_store_path`` is
    reserved for a SQLite store wired later.
    """
    job = _JOB_STORE.get(job_id)
    if job is None:
        raise KeyError(job_id)

    return {
        "job_id": job_id,
        "state": job["state"],
        "created_at": job["created_at"],
        "updated_at": job["updated_at"],
        "result": job.get("result"),
        "error": job.get("error"),
    }


def enqueue_job(job_id: str, state: str = "queued", result: Any = None, error: Any = None) -> dict[str, Any]:
    """Record a job for the demo. Not part of the public router surface in this build."""
    job = {
        "state": state,
        "created_at": _now_iso(),
        "updated_at": _now_iso(),
        "result": result,
        "error": error,
    }
    _JOB_STORE[job_id] = job
    return {"job_id": job_id, "status": state}
