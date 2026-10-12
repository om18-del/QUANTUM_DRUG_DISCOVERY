"""Offline tests for PubChem SDF fetching and 2D/3D detection.

The synthetic molfiles below mirror the real PubChem record shape: the dimension
code is the two trailing characters of the program line (line index 1), the atom
counts live on line index 3, and each atom line carries three 10-character
coordinate fields with z in columns 21-30. No network is used: the HTTP helper is
always replaced with a fake.
"""

from __future__ import annotations

from app.services.external_api import NotFoundError
from app.services.molecule_service import (
    PUBCHEM_BASE,
    _detect_sdf_dimension,
    fetch_pubchem_sdf,
    search_molecule,
)

URL_3D = f"{PUBCHEM_BASE}/compound/cid/2244/SDF?record_type=3d"
URL_DEFAULT = f"{PUBCHEM_BASE}/compound/cid/2244/SDF"


def _atom_line(x: float, y: float, z: float, symbol: str = "C") -> str:
    """One MDL V2000 atom line: three F10.4 coordinates, then the element symbol."""
    return f"{x:10.4f}{y:10.4f}{z:10.4f} {symbol:<2}"


def _sdf(dimension: str, coords: list[tuple[float, float, float]]) -> bytes:
    """A minimal but realistic MDL molfile.

    ``dimension`` is written as the trailing dimension code of the program line,
    exactly as PubChem writes it (for example ``  -OEChem-10122600532D``).
    """
    n = len(coords)
    lines = [
        "2244",
        f"  -OEChem-0101240000{dimension}",
        "",
        f"{n:3d}{n:3d}  0     0  0  0  0  0  0999 V2000",
    ]
    lines += [_atom_line(x, y, z) for (x, y, z) in coords]
    lines.append("M  END")
    return ("\n".join(lines) + "\n").encode("utf-8")


# z spread ~2.35 (the real CID 2244 3D record spans roughly -1.48 to 0.87).
_SPREAD_3D = [(0.0, 0.0, -1.4821), (1.0, 0.0, 0.8685)]
# A flat 2D layout: every z is exactly zero.
_FLAT_2D = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)]


class _FakeFetch:
    """Stand-in for ``http_get_bytes_cached`` that serves canned bytes by URL."""

    def __init__(self, responses: dict[str, object]) -> None:
        self.responses = responses
        self.urls: list[str] = []

    def __call__(self, url: str, **_kwargs: object) -> bytes:
        self.urls.append(url)
        value = self.responses.get(url)
        if isinstance(value, Exception):
            raise value
        if value is None:
            raise AssertionError(f"unexpected URL requested: {url}")
        assert isinstance(value, bytes)
        return value


class TestDetectSdfDimension:
    def test_3d_marker_with_real_z_spread_is_3d(self) -> None:
        assert _detect_sdf_dimension(_sdf("3D", _SPREAD_3D)) == "3D"

    def test_2d_marker_with_flat_atoms_is_2d(self) -> None:
        assert _detect_sdf_dimension(_sdf("2D", _FLAT_2D)) == "2D"

    def test_3d_marker_with_flat_atoms_is_2d(self) -> None:
        # A degraded record: the header says 3D but the atom block is flat, so it
        # must not be presented as a 3D structure.
        assert _detect_sdf_dimension(_sdf("3D", _FLAT_2D)) == "2D"

    def test_empty_bytes_is_2d(self) -> None:
        assert _detect_sdf_dimension(b"") == "2D"

    def test_header_without_atom_block_is_2d(self) -> None:
        # No atom block means the coordinates cannot confirm 3D, so it is 2D.
        text = "2244\n  -OEChem-01012400003D\n\n"
        assert _detect_sdf_dimension(text.encode("utf-8")) == "2D"


class TestFetchPubchemSdf:
    def test_requests_the_3d_record_first_and_keeps_it(self, monkeypatch) -> None:
        fetch = _FakeFetch({URL_3D: _sdf("3D", _SPREAD_3D), URL_DEFAULT: _sdf("2D", _FLAT_2D)})
        monkeypatch.setattr("app.services.molecule_service.http_get_bytes_cached", fetch)

        result = fetch_pubchem_sdf(2244)

        assert fetch.urls[0] == URL_3D
        assert result["sdf_dimension"] == "3D"
        assert result["note"] == "3D conformer from PubChem"

    def test_falls_back_to_the_default_record_on_404(self, monkeypatch) -> None:
        fetch = _FakeFetch(
            {URL_3D: NotFoundError("pubchem has no 3D record"), URL_DEFAULT: _sdf("2D", _FLAT_2D)}
        )
        monkeypatch.setattr("app.services.molecule_service.http_get_bytes_cached", fetch)

        result = fetch_pubchem_sdf(2244)

        assert result["sdf_dimension"] == "2D"
        assert URL_DEFAULT in fetch.urls
        assert "2D" in result["note"]

    def test_falls_back_when_the_3d_request_returns_a_2d_record(self, monkeypatch) -> None:
        # PubChem answered the 3D request, but the record is a flat 2D one, so the
        # default record is fetched instead of passing the flat one off as 3D.
        fetch = _FakeFetch({URL_3D: _sdf("2D", _FLAT_2D), URL_DEFAULT: _sdf("2D", _FLAT_2D)})
        monkeypatch.setattr("app.services.molecule_service.http_get_bytes_cached", fetch)

        result = fetch_pubchem_sdf(2244)

        assert result["sdf_dimension"] == "2D"
        assert URL_DEFAULT in fetch.urls


class TestSearchMoleculeReportsDimension:
    def _patch_pubchem(self, monkeypatch, fetch: _FakeFetch) -> None:
        monkeypatch.setattr(
            "app.services.molecule_service.pubchem_cid_from_name",
            lambda *_args, **_kwargs: 2244,
        )
        monkeypatch.setattr(
            "app.services.molecule_service.pubchem_properties",
            lambda _cid: {
                "SMILES": "CC(=O)Oc1ccccc1C(=O)O",
                "MolecularWeight": "180.16",
                "MolecularFormula": "C9H8O4",
                "IUPACName": "aspirin",
            },
        )
        monkeypatch.setattr("app.services.molecule_service.http_get_bytes_cached", fetch)

    def test_search_reports_3d_dimension(self, monkeypatch) -> None:
        fetch = _FakeFetch({URL_3D: _sdf("3D", _SPREAD_3D), URL_DEFAULT: _sdf("2D", _FLAT_2D)})
        self._patch_pubchem(monkeypatch, fetch)

        result = search_molecule("aspirin", timeout_seconds=5.0)

        assert result["status"] == "computed"
        assert result["sdf_dimension"] == "3D"
        assert result["sdf_bytes"]

    def test_search_reports_2d_dimension_when_only_2d_is_available(self, monkeypatch) -> None:
        fetch = _FakeFetch({URL_3D: NotFoundError("no 3D record"), URL_DEFAULT: _sdf("2D", _FLAT_2D)})
        self._patch_pubchem(monkeypatch, fetch)

        result = search_molecule("aspirin", timeout_seconds=5.0)

        assert result["status"] == "computed"
        assert result["sdf_dimension"] == "2D"
        assert result["sdf_bytes"]
