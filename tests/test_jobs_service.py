"""Tests for the internal job store's state validation."""

from __future__ import annotations

from pathlib import Path

import pytest

from app.services.jobs_service import VALID_JOB_STATES, enqueue_job, get_job


def test_documented_states_match_spec_section_5() -> None:
    # SPEC section 5: state is ``queued | running | done | failed``.
    assert VALID_JOB_STATES == frozenset({"queued", "running", "done", "failed"})


@pytest.mark.parametrize("state", ["queued", "running", "done", "failed"])
def test_documented_states_are_accepted(state: str) -> None:
    job_id = f"job-{state}"
    record = enqueue_job(job_id, state=state)

    assert record == {"job_id": job_id, "status": state}
    assert get_job(job_id, Path("."))["state"] == state


def test_unknown_state_is_rejected() -> None:
    with pytest.raises(ValueError) as excinfo:
        enqueue_job("job-bad", state="completed")

    message = str(excinfo.value)
    assert "Unknown job state" in message
    assert "completed" in message
    assert "queued" in message


def test_rejected_state_stores_nothing() -> None:
    with pytest.raises(ValueError):
        enqueue_job("job-nothing", state="pending")

    with pytest.raises(KeyError):
        get_job("job-nothing", Path("."))
