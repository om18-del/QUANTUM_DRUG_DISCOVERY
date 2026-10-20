"""Opt-in live-network checks for the Drug-Target Studio.

These tests download real deposited coordinates and real RCSB metadata, so they
are skipped unless ``QDDE_LIVE_NETWORK=1`` is set. The rest of the suite must
stay offline: a skipped test here is not a pass, it means the live path was not
exercised in that run.

What they check: each curated entry's pinned ligand code is really in the entry,
the entry analyses against the live file, and the explore endpoint finds at
least one drug-like ligand for a known structure.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import pytest

from app.services.complexes_service import ComplexRequestError, explore_entry, fetch_complex_structure

LIVE = os.environ.get("QDDE_LIVE_NETWORK") == "1"

pytestmark = pytest.mark.skipif(
    not LIVE,
    reason=(
        "live RCSB download; set QDDE_LIVE_NETWORK=1 to run. A skip here is not a "
        "pass: the live path was not exercised."
    ),
)

MANIFEST = Path(__file__).resolve().parents[1] / "app" / "data" / "curated_complexes" / "manifest.json"
ENTRIES = json.loads(MANIFEST.read_text(encoding="utf-8"))["complexes"]


@pytest.mark.parametrize("entry", ENTRIES, ids=[entry["id"] for entry in ENTRIES])
def test_curated_entry_analyses_against_the_live_deposit(entry: dict) -> None:
    body = fetch_complex_structure(entry["id"], MANIFEST)

    assert body["status"] == "computed", body.get("note")
    assert body["ligand_code"] == entry["ligand_code"]
    assert body["binding_site_residue_count"] >= 1
    assert body["contact_summary"]["close_contact"] >= 1
    if entry["tier"] == "curated":
        assert body["mechanism"]["status"] == "computed"
    else:
        assert body["mechanism"]["status"] == "not_curated"


@pytest.mark.parametrize("pdb_id,code", [("4EY7", "E20"), ("1HWK", "117"), ("2H42", "VIA")])
def test_explore_finds_the_expected_drug_like_ligand(pdb_id: str, code: str) -> None:
    body = explore_entry(pdb_id)

    candidates = {candidate["ligand_code"] for candidate in body.get("ligand_candidates", [])}
    assert code in candidates or body.get("ligand_code") == code


def test_explore_rejects_an_identifier_the_registry_cannot_have() -> None:
    with pytest.raises(ComplexRequestError) as excinfo:
        explore_entry("ZZZZ")
    assert excinfo.value.code == "invalid_pdb_id"
