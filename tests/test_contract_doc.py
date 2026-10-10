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
    "/api/molecule/simulate",
    "/api/protein/{pdb_id}",
    "/api/complexes/list",
    "/api/complexes/{id}/interactions",
    "/api/complexes/{id}/structure",
    "/api/complexes/{id}/docking-validation",
    "/api/jobs/{id}",
    "/api/ai/resolve",
    "/api/ai/predict-properties",
]

# Facts the contract must state, not just endpoint names.
REQUIRED_CONTRACT_STATEMENTS = [
    # Unique-residue contact counting and the rule behind each label.
    "contact_summary_basis",
    "contact_rules",
    # Every label must say whether it is inferred or measured (honesty pass), and a
    # rejected donor/acceptor pair must be reported rather than silently dropped.
    "contact_labels",
    "geometric heuristic",
    "near_miss_hydrogen_bonds",
    # The permissive ligand N/O typing must be declared an upper bound, not a count.
    "upper bound",
    # Energy is relative and populations are estimates (A5).
    "relative",
    "estimate",
    # The geometry method is a machine-readable code (A2).
    "3d_mmff94",
    # The frontend is served by FastAPI in this deployment.
    "frontend/static/index.html",
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


def test_contract_states_the_non_obvious_facts() -> None:
    text = DOC_PATH.read_text(encoding="utf-8")
    for statement in REQUIRED_CONTRACT_STATEMENTS:
        assert statement in text, f"API_CONTRACT.md must mention {statement!r}"
