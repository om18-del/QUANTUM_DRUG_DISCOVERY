from __future__ import annotations

import json
import time

import pytest

from app.config import get_settings
from app.services import external_api
from app.services.external_api import (
    NotFoundError,
    UpstreamError,
    _cache,
    _cache_key,
    http_get_bytes_cached,
    http_get_json_cached,
    invalidate_cached_url,
)
from app.services.molecule_service import (
    CompoundNotFound,
    compound_not_found_message,
    pubchem_cid_from_name,
    pubchem_properties,
    search_molecule,
)
from app.services.rdkit_helper import normalize_smiles_safe, rdkit_properties_from_smiles
from app.services.rscb_service import fetch_protein
from tests.availability import require_live_network, require_rdkit_chem, require_rdkit_descriptors

CACHE_PROBE_KEY = "probe:test:external_api"

ASPIRIN_SMILES = "CC(=O)Oc1ccccc1C(=O)O"
ASPIRIN_CID = 2244
PUBCHEM_CIDS_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/aspirin/cids/JSON"

# A large, real drug whose rule-of-five profile is deliberately poor (MW > 500).
ATORVASTATIN_SMILES = (
    "CC(C)C1=C(C(=C(N1CCC(CC(CC(=O)O)O)O)C2=CC=C(C=C2)F)C3=CC=CC=C3)C4=CC=CC=C4"
)


def _clear_probe() -> None:
    try:
        _cache().delete(CACHE_PROBE_KEY)
    except Exception:
        pass


def test_cache_stores_and_returns_value() -> None:
    cache = _cache()
    cache.set(CACHE_PROBE_KEY, {"n": 1}, expire=60)
    value = cache.get(CACHE_PROBE_KEY)
    assert value == {"n": 1}
    cache.delete(CACHE_PROBE_KEY)


def test_http_get_bytes_cached_returns_pubchem_sdf() -> None:
    require_live_network()
    url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/SDF"
    data = http_get_bytes_cached(url, source="pubchem", backoff_seconds=2.0)
    assert isinstance(data, bytes)
    assert len(data) > 0
    assert b"2244" in data
    assert b"OEChem" in data


def test_http_get_json_cached_returns_pubchem_properties() -> None:
    require_live_network()
    props = pubchem_properties(ASPIRIN_CID)
    assert props.get("MolecularFormula") == "C9H8O4"
    assert float(props.get("MolecularWeight")) == pytest.approx(180.16, abs=0.1)


def test_http_get_json_cached_respects_ttl() -> None:
    require_live_network()
    # The cache is persistent, so drop any entry left by an earlier run first:
    # otherwise the expiry we measure belongs to when it was originally cached.
    invalidate_cached_url(PUBCHEM_CIDS_URL, source="pubchem")
    http_get_json_cached(PUBCHEM_CIDS_URL, source="pubchem", backoff_seconds=2.0)

    value, expire_time = _cache().get(_cache_key("pubchem", PUBCHEM_CIDS_URL), expire_time=True)
    assert value is not None, "the response should be cached"

    ttl = get_settings().cache_ttl_pubchem_hours * 3600
    now = time.time()
    assert expire_time > 0, "pubchem entries must carry a finite TTL, not an infinite one"
    assert now + ttl - 120 <= expire_time <= now + ttl + 120


def test_invalidate_cached_url_removes_entry() -> None:
    # Pure cache behaviour: no network required.
    url = "https://example.invalid/invalidate-probe"
    cache = _cache()
    key = _cache_key("pubchem", url)
    cache.set(key, {"probe": True}, expire=300)
    assert cache.get(key) == {"probe": True}

    invalidate_cached_url(url, source="pubchem")

    assert cache.get(key) is None


def test_pubchem_search_returns_computed_for_aspirin() -> None:
    require_live_network()
    result = search_molecule("aspirin", timeout_seconds=20.0)
    assert result["status"] == "computed", result.get("note")
    assert result["pubchem_cid"] == ASPIRIN_CID
    assert result["molecular_formula"] == "C9H8O4"
    assert result["smiles"]
    assert result["canonical_smiles"]
    assert isinstance(result["sdf_bytes"], bytes)
    assert len(result["sdf_bytes"]) > 0
    assert result["sdf_loaded"] is True


def test_pubchem_search_raises_compound_not_found_for_unknown_name() -> None:
    """An unknown name is now a 404 answer, not a successful empty placeholder.

    Replaces ``test_pubchem_search_returns_placeholder_for_unknown_name``, whose
    expectation predated the distinction between "not found" and "not working".
    """
    require_live_network()
    query = "xyzzynonexistentsearchterm"

    with pytest.raises(CompoundNotFound) as excinfo:
        search_molecule(query, timeout_seconds=20.0)

    # The exception carries the query and nothing else: no URL, no HTTP detail.
    assert excinfo.value.query == query
    assert excinfo.value.args == (query,)


def test_pubchem_name_lookup_is_served_from_cache() -> None:
    """The name -> CID lookup must be cached, so the second call is a hit.

    This replaces a test named ``..._uses_cached_cid`` whose name described a
    fallback that ``search_molecule`` does not implement: when the name lookup
    fails it returns a placeholder rather than reusing a previously cached CID.
    """
    require_live_network()
    first = pubchem_cid_from_name("aspirin", timeout_seconds=20.0)
    assert first == ASPIRIN_CID
    assert _cache().get(_cache_key("pubchem", PUBCHEM_CIDS_URL)) is not None

    second = pubchem_cid_from_name("aspirin", timeout_seconds=20.0)
    assert second == first


def test_rdkit_normalize_smiles_safe_returns_canonical() -> None:
    require_rdkit_chem()
    assert normalize_smiles_safe("CCO") == "CCO"
    assert normalize_smiles_safe("OCC") == "CCO"
    assert normalize_smiles_safe("c1ccccc1") == "c1ccccc1"
    assert normalize_smiles_safe("this-is-not-a-smiles") is None


def test_rdkit_properties_from_smiles_returns_descriptors() -> None:
    require_rdkit_descriptors()
    props = rdkit_properties_from_smiles(ASPIRIN_SMILES)

    assert props["formula"] == "C9H8O4"
    # NOTE: the service fills ``molecular_weight`` from CalcExactMolWt, so this
    # is the monoisotopic mass (180.04), not the average MW (180.16).
    assert props["molecular_weight"] == pytest.approx(180.04, abs=0.03)
    assert props["exact_mass"] == pytest.approx(180.04, abs=0.03)
    assert props["h_bond_donors"] == 1
    assert props["h_bond_acceptors"] == 4
    assert props["rotatable_bonds"] == 3
    assert props["tpsa"] == pytest.approx(63.6, abs=0.5)
    # Crippen logP needs rdPartialCharges, which an OS policy may block.
    if props["logp"] is not None:
        assert props["logp"] == pytest.approx(1.31, abs=0.35)


def test_rdkit_properties_lipinski_pass_for_aspirin() -> None:
    require_rdkit_descriptors()
    props = rdkit_properties_from_smiles(ASPIRIN_SMILES)
    assert props["lipinski_rule_of_five"]["violations"] == 0
    assert props["lipinski_rule_of_five"]["pass"] is True


def test_rdkit_properties_lipinski_violations_for_large_molecule() -> None:
    require_rdkit_descriptors()
    props = rdkit_properties_from_smiles(ATORVASTATIN_SMILES)
    assert props["molecular_weight"] > 500
    assert props["lipinski_rule_of_five"]["violations"] >= 1
    assert props["lipinski_rule_of_five"]["pass"] is False


def test_rcsb_fetch_protein_6lu7_returns_computed() -> None:
    require_live_network()
    result = fetch_protein("6LU7", timeout_seconds=30.0)
    assert result["pdb_id"] == "6LU7"
    assert result["status"] == "computed", result.get("note")
    assert result["title"]
    assert isinstance(result["structure_file_bytes"], bytes)
    assert len(result["structure_file_bytes"]) > 0


def test_rcsb_fetch_protein_1hsg_returns_computed() -> None:
    require_live_network()
    result = fetch_protein("1HSG", timeout_seconds=30.0)
    assert result["status"] == "computed", result.get("note")
    assert result["title"]
    assert result["structure_file_reference"].endswith(".cif")


def test_rcsb_fetch_protein_invalid_id_returns_placeholder() -> None:
    require_live_network()
    # Well-formed but nonexistent: RCSB must be reported as unavailable, not invented.
    result = fetch_protein("9ZZZ", timeout_seconds=30.0)
    assert result["status"] == "placeholder"
    assert result["pdb_id"] == "9ZZZ"
    assert result["title"] is None
    assert result["note"]


# ---------------------------------------------------------------------------
# A1: what is worth retrying
#
# The rule under test: 404/410 are answers, 400 is our fault, and only a rate
# limit (429) or a server error (5xx) can change by asking again.
# ---------------------------------------------------------------------------
class _StubResponse:
    def __init__(self, status_code: int, payload: dict | None = None):
        self.status_code = status_code
        self._payload = payload
        self.content = json.dumps(payload).encode() if payload is not None else b""

    def json(self):
        if self._payload is None:
            raise ValueError("no JSON body")
        return self._payload


class _StubClient:
    """A stand-in httpx client that records exactly how many attempts it gets."""

    def __init__(self, *responses: _StubResponse):
        self._responses = responses
        self.calls = 0

    def get(self, *args, **kwargs):
        index = min(self.calls, len(self._responses) - 1)
        self.calls += 1
        return self._responses[index]


def _no_sleep(monkeypatch) -> list[float]:
    slept: list[float] = []
    monkeypatch.setattr(external_api, "time_sleep", lambda seconds: slept.append(seconds))
    return slept


def test_http_404_is_not_retried(monkeypatch) -> None:
    """A 404 means "not found"; asking five times would still mean "not found"."""
    slept = _no_sleep(monkeypatch)
    client = _StubClient(_StubResponse(404))

    with pytest.raises(NotFoundError):
        external_api._http_get_json(
            client,
            "https://example.invalid/compound/name/nope/cids/JSON",
            "pubchem",
            retries=5,
            backoff_seconds=0.25,
        )

    assert client.calls == 1, f"404 was retried {client.calls} times"
    assert slept == [], "a not-found must not sleep or back off"


def test_rate_limit_and_server_errors_are_retried(monkeypatch) -> None:
    """429 and 5xx are the statuses that can recover, so they are retried."""
    slept = _no_sleep(monkeypatch)
    client = _StubClient(
        _StubResponse(429),
        _StubResponse(500),
        _StubResponse(200, {"IdentifierList": {"CID": [2244]}}),
    )

    data = external_api._http_get_json(
        client,
        "https://example.invalid/compound/name/aspirin/cids/JSON",
        "pubchem",
        retries=5,
        backoff_seconds=0.01,
    )

    assert data == {"IdentifierList": {"CID": [2244]}}
    assert client.calls == 3, f"expected 3 attempts (429, 500, 200), got {client.calls}"
    assert len(slept) == 2, "each retry should back off once"


def test_non_retryable_statuses_stop_after_one_attempt(monkeypatch) -> None:
    """400 and 403 cannot be fixed by repeating the request."""
    _no_sleep(monkeypatch)

    for status in (400, 403, 404, 410):
        client = _StubClient(_StubResponse(status))
        with pytest.raises((NotFoundError, UpstreamError)):
            external_api._http_get_json(
                client,
                "https://example.invalid/x",
                "pubchem",
                retries=5,
                backoff_seconds=0.01,
            )
        assert client.calls == 1, f"HTTP {status} was retried {client.calls} times"


def test_retry_exhaustion_raises_upstream_error(monkeypatch) -> None:
    """A 500 that never recovers is an upstream failure, not a not-found."""
    _no_sleep(monkeypatch)
    client = _StubClient(_StubResponse(503))

    with pytest.raises(UpstreamError):
        external_api._http_get_json(
            client, "https://example.invalid/x", "pubchem", retries=3, backoff_seconds=0.01
        )
    assert client.calls == 3


def test_compound_not_found_returns_clean_error(monkeypatch) -> None:
    """The wire contract for a name PubChem does not know.

    404, one error envelope, the exact agreed wording, and no URL or library
    detail leaking into what the person searched for is told.
    """
    from fastapi.testclient import TestClient

    from app.main import app

    query = "xyzzynonexistentsearchterm"

    def _not_found(*args, **kwargs):
        raise CompoundNotFound(query)

    monkeypatch.setattr(
        "app.services.molecule_service.http_get_json_cached", _not_found
    )

    response = TestClient(app).get("/api/molecule/search", params={"q": query})

    assert response.status_code == 404
    body = response.json()
    assert set(body) == {"error"}, f"expected a single error envelope, got {set(body)}"
    assert body["error"]["code"] == "compound_not_found"

    expected = f"No compound found for '{query}'. Check the spelling."
    assert body["error"]["message"] == expected
    assert body["error"]["message"] == compound_not_found_message(query)

    message = body["error"]["message"]
    for leak in ("http://", "https://", "pubchem", "pugrest", "traceback", "exception"):
        assert leak not in message.lower(), f"user-facing message leaked {leak!r}"


def test_compound_not_found_is_not_reported_as_a_retryable_failure(monkeypatch) -> None:
    """A 404 from PubChem must surface as compound_not_found, not a generic 502."""
    from fastapi.testclient import TestClient

    from app.main import app

    monkeypatch.setattr(
        "app.services.molecule_service.http_get_json_cached",
        lambda *args, **kwargs: (_ for _ in ()).throw(NotFoundError("pubchem 404")),
    )

    response = TestClient(app).get(
        "/api/molecule/search", params={"q": "zzzznotarealcompoundzzzz"}
    )

    assert response.status_code == 404
    assert response.json()["error"]["code"] == "compound_not_found"
