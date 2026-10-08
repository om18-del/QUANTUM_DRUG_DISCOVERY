from __future__ import annotations

import json
import re
from pathlib import Path

DOC_PATH = Path(__file__).resolve().parents[1] / "docs" / "API_CONTRACT.md"

SPEC_SECTION_5_ENDPOINTS = [
    "/health",
    "/api/quantum/h2-curve",
    "/api/quantum/simulate",
    "/api/quantum/hadamard-demo",
    "/api/quantum/resources",
    "/api/molecule/search",
    "/api/molecule/properties",
    "/api/protein/{pdb_id}",
    "/api/complexes/list",
    "/api/complexes/{id}/interactions",
    "/api/complexes/{id}/docking-validation",
    "/api/jobs/{id}",
    "/api/ai/resolve",
    "/api/ai/predict-properties",
]


def _json_blocks(text: str) -> list[str]:
    return re.findall(r"```json\s*\n(.*?)```", text, flags=re.DOTALL)


def test_all_json_blocks_parse():
    text = DOC_PATH.read_text(encoding="utf-8")
    blocks = _json_blocks(text)
    assert blocks, "Expected at least one json block in API_CONTRACT.md"
    for block in blocks:
        parsed = json.loads(block)
        assert isinstance(parsed, dict)


def test_all_spec_section_5_endpoints_mentioned():
    text = DOC_PATH.read_text(encoding="utf-8")
    for endpoint in SPEC_SECTION_5_ENDPOINTS:
        assert endpoint in text, f"Missing endpoint {endpoint} in API_CONTRACT.md"


def test_jobs_endpoint_is_correct():
    text = DOC_PATH.read_text(encoding="utf-8")
    assert "/api/jobs/{id}" in text
    assert "#### GET /api/jobs/{id}" in text
