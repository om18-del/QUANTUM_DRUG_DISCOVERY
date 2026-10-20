"""Tests for the expanded Drug-Target Studio: tiers, categories and explore.

Everything here runs offline. The deposited-coordinate download and the RCSB
metadata call are monkeypatched with small hand-written fixtures, so no test
touches the network; the live-network checks live in
``tests/test_complexes_live_network.py`` and are opt-in via ``QDDE_LIVE_NETWORK=1``.

The fixtures are geometries invented for the test, not structures of real
complexes. Nothing in this file is a measured interaction of the proteins named
in the manifest.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest
from fastapi.testclient import TestClient

from app.main import app
from app.services import complexes_service
from app.services.structure_analysis import MIN_DRUG_HEAVY_ATOMS, drug_like_ligands

client = TestClient(app)
REPO_ROOT = Path(__file__).resolve().parents[1]
MANIFEST = json.loads(
    (REPO_ROOT / "app" / "data" / "curated_complexes" / "manifest.json").read_text(encoding="utf-8")
)
ENTRIES = {entry["id"]: entry for entry in MANIFEST["complexes"]}
NEW_PDB_IDS = {"1M17", "2HYY", "3ERT", "4EY7", "2H42", "3HS4", "1HWK", "6W63"}


def _pdb_atom(
    record: str,
    serial: int,
    name: str,
    res: str,
    chain: str,
    resseq: int,
    x: float,
    y: float,
    z: float,
    element: str,
) -> str:
    return (
        f"{record:<6}{serial:>5} {name:<4} {res:>3} {chain}{resseq:>4}    "
        f"{x:>8.3f}{y:>8.3f}{z:>8.3f}{1.0:>6.2f}{20.0:>6.2f}{'':>10}{element:>2}"
    )


def _protein_residue(serial_start: int, resname: str, resseq: int, *, x: float = 3.0) -> list[str]:
    """One protein residue, deliberately many atoms and near the ligand."""
    names = ["N", "CA", "C", "O", "CB", "CG", "CD", "CE", "CZ", "OH", "NE", "NH1", "NH2", "OD1", "OD2"]
    elements = ["N", "C", "C", "O", "C", "C", "C", "C", "C", "O", "N", "N", "N", "O", "O"]
    return [
        _pdb_atom("ATOM", serial_start + i, name, resname, "A", resseq, x, 0.05 * i, 0.0, element)
        for i, (name, element) in enumerate(zip(names, elements))
    ]


def _ligand(serial_start: int, code: str, atoms: int = 10, resseq: int = 400) -> list[str]:
    """A drug-like HETATM group with ``atoms`` heavy atoms."""
    return [
        _pdb_atom(
            "HETATM", serial_start + i, f"C{i + 1}", code, "A", resseq, 0.0, 0.2 * i, 0.0, "C"
        )
        for i in range(atoms)
    ]


def _build(*extra: list[str]) -> str:
    lines: list[str] = []
    lines += _protein_residue(1, "SER", 10)
    lines += _ligand(200, "MK1")
    for block in extra:
        lines += block
    return "\n".join(lines + ["END"]) + "\n"


FIXTURE_SINGLE = _build(
    [
        _pdb_atom("HETATM", 400, "ZN", "ZN", "A", 500, 10.0, 10.0, 10.0, "ZN"),
        _pdb_atom("HETATM", 401, "O", "HOH", "A", 600, 20.0, 20.0, 20.0, "O"),
        _pdb_atom("HETATM", 402, "C1", "EDO", "A", 601, 30.0, 30.0, 30.0, "C"),
    ]
)

# Two chemically distinct drug-like groups plus an ion: the caller must choose.
FIXTURE_MULTIPLE = _build(_ligand(500, "AB1", atoms=12, resseq=500))

# Only solvent, ions, buffers and a too-small fragment: no drug-like ligand.
FIXTURE_NO_LIGAND = "\n".join(
    _protein_residue(1, "SER", 10)
    + [
        _pdb_atom("HETATM", 300, "ZN", "ZN", "A", 500, 10.0, 10.0, 10.0, "ZN"),
        _pdb_atom("HETATM", 301, "O", "HOH", "A", 600, 20.0, 20.0, 20.0, "O"),
        _pdb_atom("HETATM", 302, "C1", "EDO", "A", 601, 30.0, 30.0, 30.0, "C"),
        _pdb_atom("HETATM", 303, "C1", "AC2", "A", 602, 1.0, 1.0, 1.0, "C"),
    ]
    + ["END"]
) + "\n"


@pytest.fixture
def patched(monkeypatch: pytest.MonkeyPatch):
    """Serve the fixture coordinates and metadata without any network call."""

    def fake_text(pdb_id: str) -> str:
        return FIXTURE_SINGLE

    def fake_metadata(pdb_id: str) -> dict[str, Any]:
        return {
            "title": f"Fixture entry {pdb_id}",
            "resolution_angstrom": 1.85,
            "method": "X-RAY DIFFRACTION",
            "deposited_atom_count": 26,
            "nonpolymer_entity_ids": [],
        }

    monkeypatch.setattr(complexes_service, "_fetch_pdb_text", fake_text)
    monkeypatch.setattr(complexes_service, "fetch_entry_metadata", fake_metadata)
    return monkeypatch


# --------------------------------------------------------------------------- #
# The curated manifest
# --------------------------------------------------------------------------- #


def test_manifest_has_the_expected_entries_and_fields() -> None:
    assert len(ENTRIES) == 12
    assert {entry["protein_pdb"] for entry in ENTRIES.values()} >= NEW_PDB_IDS

    allowed_categories = {"viral", "human_enzyme", "human_receptor"}
    allowed_tiers = {"curated", "contacts_only"}
    for entry in ENTRIES.values():
        assert entry["category"] in allowed_categories, entry["id"]
        assert entry["tier"] in allowed_tiers, entry["id"]
        assert str(entry["ligand_code"]).strip(), entry["id"]
        # Every entry must say where it is from; the filter and cards rely on it.
        assert entry["target_origin"] in {"viral", "human"}, entry["id"]


def test_every_new_entry_has_a_resolution_and_organism_from_rcsb() -> None:
    for pdb_id in NEW_PDB_IDS:
        entry = next(e for e in ENTRIES.values() if e["protein_pdb"] == pdb_id)
        assert isinstance(entry["resolution_angstrom"], (int, float)), pdb_id
        assert entry["source_organism"], pdb_id


def test_curated_entries_carry_citations_with_a_doi_or_pmid() -> None:
    curated = [e for e in ENTRIES.values() if e["tier"] == "curated"]
    assert len(curated) >= 10
    for entry in curated:
        assert entry.get("mechanism_summary"), entry["id"]
        assert entry.get("consequence"), entry["id"]
        assert entry.get("sources"), entry["id"]
        # The new entries carry structured citations because each was checked
        # against a fetched record; the four original entries keep the source
        # strings they shipped with, untouched.
        if entry["protein_pdb"] in NEW_PDB_IDS:
            citations = entry.get("citations") or []
            assert citations, entry["id"]
            for citation in citations:
                assert citation.get("title"), entry["id"]
                assert citation.get("doi") or citation.get("pmid"), entry["id"]
                if citation.get("doi"):
                    assert citation["doi"].startswith("10."), (entry["id"], citation["doi"])


def test_every_new_curated_statement_is_tagged_literature_with_its_source() -> None:
    for pdb_id in NEW_PDB_IDS:
        entry = next(e for e in ENTRIES.values() if e["protein_pdb"] == pdb_id)
        if entry["tier"] != "curated":
            continue
        statements = entry.get("literature_statements") or []
        assert statements, pdb_id
        dois = {c.get("doi") for c in entry["citations"]}
        for item in statements:
            assert item["tag"] == "LITERATURE", pdb_id
            assert item["statement"].strip(), pdb_id
            assert item["source"]["doi"] in dois, pdb_id
            assert item["source"].get("pmid"), pdb_id


def test_payload_exposes_tagged_literature_statements(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(complexes_service, "_fetch_pdb_text", lambda _pdb_id: FIXTURE_SINGLE)
    body = client.get("/api/complexes/abl-imatinib/structure").json()

    assert body["literature_statements"]
    for item in body["literature_statements"]:
        assert item["tag"] == "LITERATURE"
        assert item["source"]["doi"].startswith("10.")


def test_contacts_only_entry_has_no_narrative() -> None:
    contacts_only = [e for e in ENTRIES.values() if e["tier"] == "contacts_only"]
    assert contacts_only, "the tier must be demonstrated by at least one entry"
    for entry in contacts_only:
        assert not entry.get("mechanism_summary")
        assert not entry.get("consequence")
        assert not entry.get("catalytic_residues")
        assert not entry.get("citations")
        assert not entry.get("literature_statements")
        assert entry.get("curation_note")


def test_existing_entries_keep_their_narrative() -> None:
    """Adding the tier fields must not have rewritten the four original entries."""
    assert ENTRIES["hiv-protease-indinavir"]["mechanism_summary"].startswith("HIV-1 protease")
    assert ENTRIES["hiv-rt-nevirapine"]["binding_site_type"] == "allosteric"
    assert ENTRIES["influenza-neuraminidase-oseltamivir"]["ligand_code"] == "G39"
    assert ENTRIES["ace-lisinopril"]["target_origin"] == "human"
    for entry_id in (
        "hiv-protease-indinavir",
        "hiv-rt-nevirapine",
        "influenza-neuraminidase-oseltamivir",
        "ace-lisinopril",
    ):
        assert ENTRIES[entry_id]["tier"] == "curated"


def test_list_endpoint_exposes_category_tier_and_resolution() -> None:
    body = client.get("/api/complexes/list").json()
    by_id = {item["id"]: item for item in body["complexes"]}

    assert body["count"] == 12
    assert by_id["abl-imatinib"]["category"] == "human_enzyme"
    assert by_id["egfr-erlotinib"]["category"] == "human_receptor"
    assert by_id["sars-cov2-mpro-x77"]["tier"] == "contacts_only"
    assert by_id["sars-cov2-mpro-x77"]["resolution_angstrom"] == 2.1
    assert by_id["hiv-protease-indinavir"]["category"] == "viral"
    assert by_id["ca2-acetazolamide"]["resolution_angstrom"] == 1.1


# --------------------------------------------------------------------------- #
# The ligand rule
# --------------------------------------------------------------------------- #


def test_ligand_rule_excludes_ions_water_and_additives() -> None:
    ligands = drug_like_ligands(FIXTURE_SINGLE)
    codes = [ligand["ligand_code"] for ligand in ligands]

    assert codes == ["MK1"]
    assert "ZN" not in codes and "HOH" not in codes and "EDO" not in codes
    assert ligands[0]["heavy_atom_count"] >= MIN_DRUG_HEAVY_ATOMS


def test_ligand_rule_rejects_a_fragment_below_the_heavy_atom_threshold() -> None:
    assert drug_like_ligands(FIXTURE_NO_LIGAND) == []


def test_ligand_rule_collapses_repeated_copies_by_code() -> None:
    """1HWK deposits four atorvastatin copies; the caller chooses a chemical ligand."""
    text = _build(
        _ligand(500, "117", atoms=12, resseq=500),
        _ligand(600, "117", atoms=12, resseq=600),
    )
    ligands = drug_like_ligands(text)
    by_code = {ligand["ligand_code"]: ligand for ligand in ligands}

    assert by_code["117"]["instance_count"] == 2
    # The default MK1 copy in the shared fixture is still its own chemical ligand.
    assert by_code["MK1"]["instance_count"] == 1


# --------------------------------------------------------------------------- #
# Explore endpoint
# --------------------------------------------------------------------------- #


def test_explore_single_ligand_returns_measured_contacts(patched: pytest.MonkeyPatch) -> None:
    response = client.get("/api/complexes/explore/1ABC")
    assert response.status_code == 200
    body = response.json()

    assert body["status"] == "computed"
    assert body["ligand_code"] == "MK1"
    assert body["curation"] == "none"
    assert body["tier"] == "contacts_only"
    assert "measured contacts only" in body["curation_note"]
    assert body["mechanism"]["status"] == "not_curated"
    assert body["mechanism"]["steps"] == []
    assert body["mechanism"]["catalytic_residues"] == []
    assert body["protein_pdb"] == "1ABC"
    assert body["resolution_angstrom"] == 1.85
    # Same payload as a curated complex, including the truncation flags.
    for field in (
        "contact_summary",
        "contact_rules",
        "contact_labels",
        "interacting_residues",
        "binding_site_residues",
        "atom_pair_counts",
        "atom_pairs_truncated",
        "atom_pairs_truncated_note",
        "method",
        "structure_file_bytes",
    ):
        assert field in body, field
    assert isinstance(body["atom_pairs_truncated"], bool)


def test_explore_does_not_pick_the_ion(patched: pytest.MonkeyPatch) -> None:
    body = client.get("/api/complexes/explore/1abc").json()

    assert body["ligand_code"] != "ZN"
    assert [c["ligand_code"] for c in body["ligand_candidates"]] == ["MK1"]


def test_explore_returns_candidates_instead_of_choosing(patched: pytest.MonkeyPatch) -> None:
    patched.setattr(complexes_service, "_fetch_pdb_text", lambda _pdb_id: FIXTURE_MULTIPLE)
    body = client.get("/api/complexes/explore/1ABC").json()

    assert body["status"] == "ligand_selection_required"
    assert body["ligand_selection_required"] is True
    assert {c["ligand_code"] for c in body["ligand_candidates"]} == {"MK1", "AB1"}
    assert body["status"] != "computed"

    chosen = client.get("/api/complexes/explore/1ABC?ligand=ab1").json()
    assert chosen["status"] == "computed"
    assert chosen["ligand_code"] == "AB1"


def test_explore_rejects_an_unknown_ligand_code(patched: pytest.MonkeyPatch) -> None:
    response = client.get("/api/complexes/explore/1ABC?ligand=ZZZ")
    assert response.status_code == 404
    assert response.json()["error"]["code"] == "ligand_not_found"
    assert "MK1" in response.json()["error"]["message"]


@pytest.mark.parametrize("bad", ["ABC", "1AB", "1ABC2", "ABCD", "12"])
def test_explore_invalid_pdb_id(bad: str) -> None:
    response = client.get(f"/api/complexes/explore/{bad}")
    assert response.status_code == 422
    assert response.json()["error"]["code"] == "invalid_pdb_id"


def test_explore_empty_id_is_rejected_without_a_request() -> None:
    with pytest.raises(complexes_service.ComplexRequestError) as excinfo:
        complexes_service.explore_entry("")
    assert excinfo.value.code == "invalid_pdb_id"
    assert excinfo.value.status_code == 422


def test_explore_entry_without_a_drug_like_ligand(patched: pytest.MonkeyPatch) -> None:
    patched.setattr(complexes_service, "_fetch_pdb_text", lambda _pdb_id: FIXTURE_NO_LIGAND)
    response = client.get("/api/complexes/explore/1ABC")

    assert response.status_code == 404
    assert response.json()["error"]["code"] == "no_drug_like_ligand"


def test_explore_reports_unreachable_rcsb(patched: pytest.MonkeyPatch) -> None:
    patched.setattr(complexes_service, "_fetch_pdb_text", lambda _pdb_id: None)
    response = client.get("/api/complexes/explore/1ABC")

    assert response.status_code == 503
    error = response.json()["error"]
    assert error["code"] == "structure_unavailable"
    assert "Nothing was invented" in error["message"]


def test_explore_enforces_the_atom_cap(patched: pytest.MonkeyPatch) -> None:
    patched.setattr(complexes_service, "MAX_PROTEIN_ATOMS", 2)
    response = client.get("/api/complexes/explore/1ABC")

    assert response.status_code == 422
    assert response.json()["error"]["code"] == "entry_too_large"


def test_explore_enforces_the_byte_cap(patched: pytest.MonkeyPatch) -> None:
    patched.setattr(complexes_service, "MAX_ENTRY_BYTES", 10)
    response = client.get("/api/complexes/explore/1ABC")

    assert response.status_code == 422
    assert response.json()["error"]["code"] == "entry_too_large"


def test_explore_flags_truncated_atom_pairs(patched: pytest.MonkeyPatch) -> None:
    """One residue with many atoms against a 10-atom ligand exceeds the 12-pair cap."""
    body = client.get("/api/complexes/explore/1ABC").json()

    assert body["atom_pairs_truncated"] is True
    assert body["atom_pairs_truncated_note"]
    assert body["contact_summary"]["hydrophobic"] == 1


# --------------------------------------------------------------------------- #
# Curated payload tiers
# --------------------------------------------------------------------------- #


def test_contacts_only_entry_produces_no_mechanism(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(complexes_service, "_fetch_pdb_text", lambda _pdb_id: FIXTURE_SINGLE)
    body = client.get("/api/complexes/sars-cov2-mpro-x77/structure").json()

    assert body["status"] == "computed"
    assert body["tier"] == "contacts_only"
    assert body["curation"] == "none"
    assert body["mechanism"]["status"] == "not_curated"
    assert body["mechanism"]["steps"] == []
    assert body["mechanism"]["catalytic_residues"] == []
    assert "measured geometry" in body["mechanism"]["note"]


def test_curated_entry_without_a_residue_list_says_so(monkeypatch: pytest.MonkeyPatch) -> None:
    """No residue list is written from memory, so the payload must admit the gap."""
    monkeypatch.setattr(complexes_service, "_fetch_pdb_text", lambda _pdb_id: FIXTURE_SINGLE)
    body = client.get("/api/complexes/egfr-erlotinib/structure").json()

    assert body["tier"] == "curated"
    assert body["curation"] == "curated"
    assert body["mechanism"]["status"] == "computed"
    assert body["mechanism"]["catalytic_residues"] == []
    assert body["catalytic_residues_note"]
    assert body["citations"][0]["doi"].startswith("10.")


def test_curated_payload_exposes_citations_and_category(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(complexes_service, "_fetch_pdb_text", lambda _pdb_id: FIXTURE_SINGLE)
    body = client.get("/api/complexes/pde5-sildenafil/structure").json()

    assert body["category"] == "human_enzyme"
    assert body["resolution_angstrom"] == 2.3
    assert body["source_organism"] == "Homo sapiens"
    assert body["citations"][0]["pmid"] == "16735511"
    # The one literature-verified residue for this entry, and its measured status.
    glycine = body["mechanism"]["catalytic_residues"][0]
    assert glycine["match"] == "GLY A659"
    assert glycine["in_contact"] is False
