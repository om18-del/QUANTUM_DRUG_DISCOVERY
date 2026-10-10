"""Tests for the offline-precomputed H2 Hamiltonian artifact.

``app/data/hamiltonians/h2_sto3g.json`` is produced by
``scripts/precompute_hamiltonians.ipynb`` (Colab or WSL), because PySCF publishes
no Windows wheel. Until that file has been generated and saved into the repo,
these tests skip with an explicit reason instead of passing vacuously.

The central claim being checked is the one the file promises: for every bond
length, the smallest eigenvalue of the stored qubit operator plus the stored
nuclear repulsion reproduces the stored PySCF FCI energy to 1e-6 Ha.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pytest

ARTIFACT = (
    Path(__file__).resolve().parents[1]
    / "app"
    / "data"
    / "hamiltonians"
    / "h2_sto3g.json"
)

# The grid the artifact must contain: 0.20 .. 3.00 angstrom in 0.05 steps (57
# uniform points) plus 0.74 angstrom -- the UI default and H2's real equilibrium
# bond length, which the 0.05 step misses (it would need index 10.8). 58 points in
# total, ascending, with no duplicates.
UNIFORM_POINTS = 57
EXPECTED_POINTS = 58
EXPECTED_BOND_LENGTHS = sorted(
    {*[round(0.20 + 0.05 * i, 2) for i in range(UNIFORM_POINTS)], 0.74}
)
TOLERANCE_HA = 1e-6
PINNED_VERSIONS = {"qiskit": "1.4.6", "qiskit-nature": "0.8.0", "pyscf": "2.14.0"}

TOP_LEVEL_FIELDS = (
    "molecule",
    "basis",
    "mapping",
    "two_qubit_reduction",
    "n_qubits",
    "library_versions",
    "point_count",
    "points",
)

POINT_FIELDS = (
    "bond_length_angstrom",
    "nuclear_repulsion_ha",
    "hf_energy_ha",
    "fci_energy_ha",
    "qubit_energy_plus_nuclear_ha",
    "deviation_from_fci_ha",
    "n_qubits",
    "pauli_terms",
)

_PAULI = {
    "I": np.eye(2, dtype=complex),
    "X": np.array([[0, 1], [1, 0]], dtype=complex),
    "Y": np.array([[0, -1j], [1j, 0]], dtype=complex),
    "Z": np.array([[1, 0], [0, -1]], dtype=complex),
}


def _load() -> dict[str, Any]:
    if not ARTIFACT.is_file():
        pytest.skip(
            f"{ARTIFACT.name} is not present. Run scripts/precompute_hamiltonians.ipynb "
            "(or scripts/precompute_hamiltonians.py) in Colab/WSL, then save the result to "
            "app/data/hamiltonians/h2_sto3g.json"
        )
    return json.loads(ARTIFACT.read_text(encoding="utf-8"))


def _matrix_from_terms(terms: list[dict[str, Any]]) -> np.ndarray:
    """Dense matrix of a Pauli-sum operator, built without Qiskit."""
    n_qubits = len(terms[0]["pauli"])
    matrix = np.zeros((2**n_qubits, 2**n_qubits), dtype=complex)
    for term in terms:
        op = np.array([[1.0 + 0.0j]])
        for char in term["pauli"]:
            op = np.kron(op, _PAULI[char])
        matrix = matrix + term["coefficient"] * op
    return matrix


def _ground_energy(terms: list[dict[str, Any]]) -> float:
    matrix = _matrix_from_terms(terms)
    hermitian = (matrix + matrix.conj().T) / 2.0
    return float(np.linalg.eigvalsh(hermitian)[0].real)


def test_point_count_is_58() -> None:
    data = _load()
    assert data["point_count"] == EXPECTED_POINTS
    assert len(data["points"]) == EXPECTED_POINTS
    assert EXPECTED_POINTS == 58


def test_top_level_fields_present() -> None:
    data = _load()
    for field in TOP_LEVEL_FIELDS:
        assert field in data, f"missing top-level field {field!r}"


def test_declares_parity_mapping_and_two_qubits() -> None:
    data = _load()
    assert data["molecule"] == "H2"
    assert data["basis"].lower() == "sto-3g"
    assert data["mapping"] == "parity"
    assert data["two_qubit_reduction"] is True
    assert data["n_qubits"] == 2


def test_library_versions_match_the_pins() -> None:
    data = _load()
    for name, expected in PINNED_VERSIONS.items():
        assert data["library_versions"].get(name) == expected, (
            f"{name} should be pinned to {expected}, got {data['library_versions'].get(name)}"
        )


def test_bond_length_grid_is_57_uniform_points_plus_0_74() -> None:
    """The full grid, so an extra, missing or reordered point cannot slip through."""
    data = _load()
    lengths = [point["bond_length_angstrom"] for point in data["points"]]
    assert lengths == EXPECTED_BOND_LENGTHS


def test_grid_contains_the_ui_default_0_74() -> None:
    """0.74 A is the UI default and H2's equilibrium: it must be storable exactly."""
    data = _load()
    lengths = [point["bond_length_angstrom"] for point in data["points"]]
    assert 0.74 in lengths
    assert lengths == sorted(lengths)
    assert len(set(lengths)) == len(lengths)


def test_every_point_has_required_fields() -> None:
    data = _load()
    for point in data["points"]:
        for field in POINT_FIELDS:
            assert field in point, f"missing {field!r} at r={point.get('bond_length_angstrom')}"


def test_pauli_terms_are_two_qubit_with_real_coefficients() -> None:
    data = _load()
    for point in data["points"]:
        terms = point["pauli_terms"]
        assert terms, f"no Pauli terms at r={point['bond_length_angstrom']}"
        assert point["n_qubits"] == 2
        for term in terms:
            assert len(term["pauli"]) == 2
            assert set(term["pauli"]) <= set("IXYZ")
            # A complex coefficient would mean the mapping produced something
            # that cannot be read as a real qubit Hamiltonian.
            assert isinstance(term["coefficient"], float)


@pytest.mark.parametrize("index", range(EXPECTED_POINTS))
def test_qubit_energy_plus_nuclear_repulsion_matches_fci(index: int) -> None:
    """diag(qubit Hamiltonian) + nuclear repulsion == PySCF FCI, to 1e-6 Ha."""
    data = _load()
    point = data["points"][index]
    computed = _ground_energy(point["pauli_terms"]) + point["nuclear_repulsion_ha"]
    assert computed == pytest.approx(point["fci_energy_ha"], abs=TOLERANCE_HA), (
        f"r={point['bond_length_angstrom']} A: qubit+nuc={computed:.12f} "
        f"vs stored FCI={point['fci_energy_ha']:.12f}"
    )


def test_recorded_deviation_is_within_tolerance() -> None:
    data = _load()
    assert data["max_deviation_from_fci_ha"] <= TOLERANCE_HA
    for point in data["points"]:
        assert abs(point["qubit_energy_plus_nuclear_ha"] - point["fci_energy_ha"]) <= TOLERANCE_HA
        assert point["deviation_from_fci_ha"] <= TOLERANCE_HA


def test_h2_minimum_near_0_74_angstrom_is_about_minus_1_137_ha() -> None:
    """SPEC.md section 6 reference: H2 in STO-3G near 0.74 A is about -1.137 Ha."""
    data = _load()
    nearest = min(
        data["points"], key=lambda p: abs(p["bond_length_angstrom"] - 0.74)
    )
    assert nearest["bond_length_angstrom"] == pytest.approx(0.74, abs=1e-9)
    assert nearest["fci_energy_ha"] == pytest.approx(-1.137, abs=0.005)

    # The curve must also be a genuine dissociation curve: the minimum should sit
    # near 0.74 A, well below the energy at the shortest and longest separations.
    energies = {p["bond_length_angstrom"]: p["fci_energy_ha"] for p in data["points"]}
    argmin = min(energies, key=energies.get)
    assert 0.5 <= argmin <= 1.0, f"equilibrium bond length looks wrong: {argmin} A"
    assert energies[argmin] < energies[0.20]
    assert energies[argmin] < energies[3.00]


def test_hf_energy_is_above_fci_at_every_point() -> None:
    """A variational sanity check that does not depend on our conventions."""
    data = _load()
    for point in data["points"]:
        assert point["hf_energy_ha"] >= point["fci_energy_ha"] - 1e-9, (
            f"HF below FCI at r={point['bond_length_angstrom']}"
        )
