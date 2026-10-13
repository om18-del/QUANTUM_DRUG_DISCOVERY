"""Tests for the generalised molecule artifact layer.

The artifact layer used to be H2-only. These tests exercise the registry, the
readers taking a molecule key, the ``/api/quantum/molecules`` endpoint, the curve
for a non-H2 molecule (served from stored data plus an exact NumPy reference) and
the validator's molecule-key / grid-from-metadata behaviour.

Every non-H2 payload here is a **synthetic fixture**: it satisfies the documented
schema and the energy cross-check by construction (a single identity term ``c * I``
has ground energy ``c``, and ``c + nuclear_repulsion_ha`` is set equal to the stored
FCI), but it is not LiH chemistry and is never presented as such. The real LiH
artifact does not exist yet; the H2 tests skip if the H2 artifact is absent.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pytest
from fastapi.testclient import TestClient

from app.core.exact import (
    MAX_DENSE_QUBITS,
    MOLECULE_REGISTRY,
    TooManyQubitsError,
    UnknownMoleculeError,
    artifact_filename,
    ensure_dense_diagonalisable,
    h2_hamiltonian_from_artifact,
    molecule_label,
    normalize_molecule_key,
)
from app.main import app
from app.services.quantum_service import h2_artifact_path

client = TestClient(app)
REPO_ROOT = Path(__file__).resolve().parents[1]
VALIDATOR = REPO_ROOT / "scripts" / "validate_hamiltonian_file.py"
LI_H_ARTIFACT = REPO_ROOT / "app" / "data" / "hamiltonians" / "lih_sto3g_active.json"

ACTIVE_SPACE_TEXT = (
    "LiH STO-3G: 6 spatial orbitals in the full basis; froze the Li 1s core and removed the 2 "
    "highest virtual spatial orbitals (zero-based indices [4, 5]); active space = 3 spatial / 6 "
    "spin orbitals, 2 active electrons."
)


def _point_quantum_service_at(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(
        "app.services.quantum_service.get_settings",
        lambda: SimpleNamespace(h2_hamiltonian_dir=tmp_path),
    )


def _require_checked_in_h2() -> None:
    if not h2_artifact_path().is_file():
        pytest.skip("h2_sto3g.json is not present in this checkout")


def _lih_artifact(
    n_qubits: int = 4,
    *,
    lengths: tuple[float, ...] = (1.0, 1.1, 1.2),
    active_space: str | None = ACTIVE_SPACE_TEXT,
) -> dict[str, Any]:
    """A schema-valid synthetic non-H2 artifact (consistency, not chemistry)."""
    points: list[dict[str, Any]] = []
    for bond_length in lengths:
        fci = -7.5 + 0.02 * (bond_length - 1.595) ** 2
        nuclear = 1.0 / bond_length
        # A single identity term: its ground energy is the coefficient itself.
        coefficient = fci - nuclear
        points.append(
            {
                "bond_length_angstrom": bond_length,
                "nuclear_repulsion_ha": nuclear,
                "hf_energy_ha": fci + 0.02,
                "fci_energy_ha": fci,
                "qubit_energy_ha": coefficient,
                "qubit_energy_plus_nuclear_ha": fci,
                "deviation_from_fci_ha": 0.0,
                "n_qubits": n_qubits,
                "pauli_terms": [{"pauli": "I" * n_qubits, "coefficient": coefficient}],
                "nuclear_repulsion_in_qubit_operator": False,
                "energy_span_two_state_ha": None,
            }
        )

    payload: dict[str, Any] = {
        "molecule": "LiH",
        "basis": "sto-3g",
        "mapping": "parity",
        "two_qubit_reduction": True,
        "n_qubits": n_qubits,
        "units": "hartree",
        "generated_by": "tests/test_molecule_registry.py (synthetic)",
        "library_versions": {"qiskit": "1.4.6", "qiskit-nature": "0.8.0", "pyscf": "2.14.0"},
        "convention": "synthetic fixture; the identity term's ground energy is its coefficient.",
        "max_deviation_from_fci_ha": 0.0,
        "point_count": len(points),
        "bond_length_grid": [float(length) for length in lengths],
        "points": points,
    }
    if active_space is not None:
        payload["active_space"] = active_space
    return payload


# ---------------------------------------------------------------------------
# Registry and readers
# ---------------------------------------------------------------------------
def test_registry_lists_h2_and_lih() -> None:
    assert set(MOLECULE_REGISTRY) == {"h2", "lih"}
    assert MOLECULE_REGISTRY["h2"] == {"file": "h2_sto3g.json", "label": "H2"}
    assert MOLECULE_REGISTRY["lih"]["file"] == "lih_sto3g_active.json"
    # The default key is h2, so every existing caller keeps working.
    assert artifact_filename() == "h2_sto3g.json"
    assert artifact_filename("LIh") == "lih_sto3g_active.json"
    assert molecule_label("lih") == "LiH (active space)"


def test_unknown_molecule_key_is_a_clear_error() -> None:
    with pytest.raises(UnknownMoleculeError, match="unknown molecule"):
        normalize_molecule_key("caffeine")
    with pytest.raises(UnknownMoleculeError):
        artifact_filename("aspirin")
    with pytest.raises(UnknownMoleculeError):
        h2_hamiltonian_from_artifact({"points": []}, 0.74, "nope")


def test_readers_generalise_beyond_two_qubits_and_five_terms() -> None:
    payload = _lih_artifact(n_qubits=4)
    hamiltonian = h2_hamiltonian_from_artifact(payload, 1.0, "lih")

    assert hamiltonian.n_qubits == 4
    assert hamiltonian.matrix.shape == (16, 16)
    assert hamiltonian.total_energy_ha == pytest.approx(
        payload["points"][0]["fci_energy_ha"], abs=1e-12
    )
    assert hamiltonian.active_space == ACTIVE_SPACE_TEXT
    # A single Pauli term is fine: the reader never assumed five of them.
    assert set(hamiltonian.pauli_terms) == {"IIII"}


def test_reader_requires_active_space_for_a_non_h2_artifact() -> None:
    payload = _lih_artifact(active_space=None)
    with pytest.raises(ValueError, match="active_space"):
        h2_hamiltonian_from_artifact(payload, 1.0, "lih")


def test_dense_diagonalisation_guard() -> None:
    ensure_dense_diagonalisable(MAX_DENSE_QUBITS, "lih")
    with pytest.raises(TooManyQubitsError, match=f"{MAX_DENSE_QUBITS}-qubit"):
        ensure_dense_diagonalisable(MAX_DENSE_QUBITS + 1, "lih")


# ---------------------------------------------------------------------------
# GET /api/quantum/molecules
# ---------------------------------------------------------------------------
def test_molecules_endpoint_reports_availability(tmp_path: Path, monkeypatch) -> None:
    (tmp_path / "h2_sto3g.json").write_text("{}", encoding="utf-8")
    _point_quantum_service_at(tmp_path, monkeypatch)

    response = client.get("/api/quantum/molecules")

    assert response.status_code == 200
    body = response.json()
    by_key = {entry["molecule"]: entry for entry in body["molecules"]}
    assert set(by_key) == {"h2", "lih"}
    assert by_key["h2"]["available"] is True
    assert by_key["lih"]["available"] is False
    assert by_key["lih"]["file"] == "lih_sto3g_active.json"
    assert by_key["lih"]["label"] == "LiH (active space)"
    assert body["available_count"] == 1


def test_molecules_endpoint_empty_directory(tmp_path: Path, monkeypatch) -> None:
    _point_quantum_service_at(tmp_path, monkeypatch)
    body = client.get("/api/quantum/molecules").json()
    assert all(entry["available"] is False for entry in body["molecules"])
    assert body["available_count"] == 0


# ---------------------------------------------------------------------------
# Curve endpoint: unknown molecule, missing artifact, LiH data, guard
# ---------------------------------------------------------------------------
def test_curve_unknown_molecule_is_404() -> None:
    response = client.get("/api/quantum/h2-curve", params={"molecule": "caffeine"})

    assert response.status_code == 404
    data = response.json()
    assert set(data) == {"error"}
    assert data["error"]["code"] == "unknown_molecule"
    assert "unknown molecule" in data["error"]["message"]


def test_curve_lih_is_503_without_the_artifact(tmp_path: Path, monkeypatch) -> None:
    _point_quantum_service_at(tmp_path, monkeypatch)

    response = client.get("/api/quantum/h2-curve", params={"molecule": "lih"})

    assert response.status_code == 503
    data = response.json()
    assert data["error"]["code"] == "hamiltonian_unavailable"
    assert "lih_sto3g_active.json" in data["error"]["message"]


def test_curve_lih_serves_stored_data_and_an_exact_reference(
    tmp_path: Path, monkeypatch
) -> None:
    payload = _lih_artifact(lengths=(1.0, 1.1, 1.2))
    (tmp_path / "lih_sto3g_active.json").write_text(json.dumps(payload), encoding="utf-8")
    _point_quantum_service_at(tmp_path, monkeypatch)

    response = client.get("/api/quantum/h2-curve", params={"molecule": "lih"})

    assert response.status_code == 200, response.text
    data = response.json()
    assert data["status"] == "computed"
    assert data["molecule"] == "LiH"
    assert data["n_qubits"] == 4
    assert data["active_space"] == ACTIVE_SPACE_TEXT
    assert data["artifact"] == "lih_sto3g_active.json"
    assert data["point_count"] == len(data["points"]) == 3
    assert data["bond_lengths_angstrom"] == [1.0, 1.1, 1.2]
    assert len(data["hf_energies_ha"]) == len(data["fci_energies_ha"]) == 3
    assert len(data["nuclear_repulsion_energies_ha"]) == len(data["exact_energies_ha"]) == 3
    # No VQE for this molecule: only stored data and the exact reference.
    assert "vqe_energies_ha" not in data
    assert "vqe_error_mha" not in data

    for point, exact in zip(data["points"], data["exact_energies_ha"]):
        assert point["exact_energy_ha"] == pytest.approx(exact, abs=1e-12)
        assert exact == pytest.approx(point["fci_energy_ha"], abs=1e-8)
    assert "No VQE runs for this molecule" in data["note"]


def test_curve_lih_refuses_more_than_eight_qubits(tmp_path: Path, monkeypatch) -> None:
    payload = _lih_artifact(n_qubits=MAX_DENSE_QUBITS + 1, lengths=(1.0, 1.1))
    (tmp_path / "lih_sto3g_active.json").write_text(json.dumps(payload), encoding="utf-8")
    _point_quantum_service_at(tmp_path, monkeypatch)

    response = client.get("/api/quantum/h2-curve", params={"molecule": "lih"})

    assert response.status_code == 422
    data = response.json()
    assert set(data) == {"error"}
    assert data["error"]["code"] == "too_many_qubits"


def test_curve_lih_without_active_space_is_a_fault(tmp_path: Path, monkeypatch) -> None:
    payload = _lih_artifact(active_space=None)
    (tmp_path / "lih_sto3g_active.json").write_text(json.dumps(payload), encoding="utf-8")
    _point_quantum_service_at(tmp_path, monkeypatch)

    response = client.get("/api/quantum/h2-curve", params={"molecule": "lih"})

    assert response.status_code == 500
    assert response.json()["error"]["code"] == "simulation_error"


# ---------------------------------------------------------------------------
# Simulate endpoint: LiH is 422; the default molecule is unchanged
# ---------------------------------------------------------------------------
def test_simulate_lih_is_422_vqe_not_supported() -> None:
    response = client.post(
        "/api/quantum/simulate",
        content=json.dumps({"bond_length": 1.5, "molecule": "lih"}),
    )

    assert response.status_code == 422
    data = response.json()
    assert set(data) == {"error"}
    assert data["error"]["code"] == "vqe_not_supported_for_molecule"
    message = data["error"]["message"]
    assert "curve" in message
    assert "exact diagonalisation" in message


def test_simulate_default_molecule_is_unchanged() -> None:
    _require_checked_in_h2()
    response = client.post(
        "/api/quantum/simulate",
        content=json.dumps({"bond_length": 0.74, "seed": 7}),
    )
    assert response.status_code == 200, response.text
    data = response.json()
    assert data["status"] == "computed"
    assert data["molecule"] == "H2"
    assert data["within_chemical_accuracy"] is True


def test_curve_default_molecule_is_unchanged() -> None:
    _require_checked_in_h2()
    data = client.get("/api/quantum/h2-curve").json()
    assert data["status"] == "computed"
    assert data["molecule"] == "H2"
    assert len(data["vqe_energies_ha"]) == 58


# ---------------------------------------------------------------------------
# Validator: molecule key and grid read from the file's own metadata
# ---------------------------------------------------------------------------
def test_validator_accepts_the_molecule_key_h2() -> None:
    _require_checked_in_h2()
    result = subprocess.run(
        [sys.executable, str(VALIDATOR), "h2"],
        cwd=str(REPO_ROOT),
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    assert "OK: every check passed." in result.stdout


def test_validator_lih_key_names_its_own_artifact_when_absent() -> None:
    if LI_H_ARTIFACT.is_file():
        pytest.skip("lih_sto3g_active.json is present in this checkout")
    result = subprocess.run(
        [sys.executable, str(VALIDATOR), "lih"],
        cwd=str(REPO_ROOT),
        capture_output=True,
        text=True,
    )
    assert result.returncode == 1
    assert "lih_sto3g_active.json" in result.stderr


def test_validator_reads_the_grid_from_the_file_metadata(tmp_path: Path) -> None:
    """A declared grid of 3 points is checked against itself, not against H2's 58."""
    payload = _lih_artifact(lengths=(1.0, 1.1, 1.2))
    path = tmp_path / "lih_synthetic.json"
    path.write_text(json.dumps(payload), encoding="utf-8")

    result = subprocess.run(
        [sys.executable, str(VALIDATOR), str(path)],
        cwd=str(REPO_ROOT),
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stderr
    assert "points                 : 3 (1 .. 1.2 A)" in result.stdout
    assert "OK: every check passed." in result.stdout
