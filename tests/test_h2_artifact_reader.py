"""Tests for reading the artifact the writer actually produces.

``scripts/precompute_hamiltonians.py`` stores, per bond length, a list of Pauli
terms plus a separate nuclear repulsion. ``h2_sto3g_hamiltonian_from_json`` read a
different shape (``n_qubits`` plus a dense ``matrix``), so before these readers the
application could not read the artifact it was built to consume. Everything here is
a **synthetic fixture**: the real ``app/data/hamiltonians/h2_sto3g.json`` is
generated offline (Colab/WSL, because PySCF has no Windows wheel) and is not in the
repository, so no test here needs it.

The fixture operator is chosen so its ground energy is known analytically:

    H = a * ZI + c * XX,  c > 0,  a^2 = d^2 - c^2,  d = nuclear_repulsion - fci

Its 4x4 matrix has eigenvalues ``+sqrt(a^2 + c^2)`` and ``-sqrt(a^2 + c^2)`` (each
twice), so the ground energy is exactly ``-d`` and the stored ``fci`` is
``ground + nuclear`` by construction. No diagonaliser is used to build the fixture,
which is what makes the cross-check below a real cross-check: the reader's NumPy
assembly and the validator's own NumPy assembly are two implementations compared
against one analytic number, not one implementation compared with itself.
"""

from __future__ import annotations

import copy
import importlib.util
import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from app.core.exact import (
    ArtifactPointHamiltonian,
    h2_artifact_bond_lengths,
    h2_artifact_pauli_terms,
    h2_artifact_point,
    h2_hamiltonian_from_artifact,
    h2_sto3g_hamiltonian_from_json,
    nearest_stored_bond_length,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
VALIDATOR = REPO_ROOT / "scripts" / "validate_hamiltonian_file.py"

# 57 uniform points 0.20 .. 3.00 in 0.05 steps, plus the UI default 0.74.
UNIFORM_POINTS = 57
EXPECTED_BOND_LENGTHS = sorted(
    {*[round(0.20 + 0.05 * i, 2) for i in range(UNIFORM_POINTS)], 0.74}
)

# The toy grid used by the reader unit tests: three stored points, one of which is
# the UI default. 0.80 and 1.234 are deliberately not on it. Every mutation below
# targets TARGET_INDEX, which is the point the reader is asked for, so a mutation
# that silently lands on another point cannot make a rejection test pass vacuously.
TOY_LENGTHS = [0.70, 0.74, 0.75]
TARGET_INDEX = 1

XX_COEFFICIENT = 0.05
HF_MINUS_FCI_HA = 0.02


def _analytic_point(bond_length: float, *, nuclear_repulsion: float = 0.70) -> dict:
    """One stored point whose ground energy is ``-(nuclear - fci)`` by construction."""
    fci = -1.137 + 0.05 * (bond_length - 0.74) ** 2
    d = nuclear_repulsion - fci
    assert d > XX_COEFFICIENT, "the fixture needs d > c so a^2 stays positive"
    a = math.sqrt(d**2 - XX_COEFFICIENT**2)
    return {
        "bond_length_angstrom": bond_length,
        "nuclear_repulsion_ha": nuclear_repulsion,
        "hf_energy_ha": fci + HF_MINUS_FCI_HA,
        "fci_energy_ha": fci,
        # d = -(ground energy), so these record the same claim the reader makes.
        "qubit_energy_ha": -d,
        "qubit_energy_plus_nuclear_ha": fci,
        "deviation_from_fci_ha": 0.0,
        "n_qubits": 2,
        "pauli_terms": [
            {"pauli": "ZI", "coefficient": a},
            {"pauli": "XX", "coefficient": XX_COEFFICIENT},
        ],
        "nuclear_repulsion_in_qubit_operator": False,
        "energy_span_two_state_ha": None,
    }


def _artifact(lengths: list[float], *, nuclear_repulsion: float = 0.70) -> dict:
    points = [_analytic_point(length, nuclear_repulsion=nuclear_repulsion) for length in lengths]
    return {
        "molecule": "H2",
        "basis": "sto-3g",
        "mapping": "parity",
        "two_qubit_reduction": True,
        "n_qubits": 2,
        "units": "hartree",
        "generated_by": "tests/test_h2_artifact_reader.py (synthetic)",
        "library_versions": {"qiskit": "1.4.6", "qiskit-nature": "0.8.0", "pyscf": "2.14.0"},
        "convention": "synthetic fixture; not a real PySCF run",
        "max_deviation_from_fci_ha": 0.0,
        "point_count": len(points),
        "points": points,
    }


def _validator_module():
    """Import the validator as a module so its recomputation can be compared directly."""
    spec = importlib.util.spec_from_file_location("validate_hamiltonian_file", VALIDATOR)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write(tmp_path: Path, payload: dict) -> Path:
    path = tmp_path / "synthetic_artifact.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


# ---------------------------------------------------------------------------
# Reading one stored bond length
# ---------------------------------------------------------------------------
def test_reader_reproduces_the_analytic_fixture_energy() -> None:
    """Ground state + nuclear repulsion must equal the stored FCI, per point."""
    data = _artifact(TOY_LENGTHS)

    for point in data["points"]:
        hamiltonian = h2_hamiltonian_from_artifact(data, point["bond_length_angstrom"])
        assert hamiltonian.total_energy_ha == pytest.approx(point["fci_energy_ha"], abs=1e-9)
        assert hamiltonian.ground_energy_ha == pytest.approx(
            point["fci_energy_ha"] - point["nuclear_repulsion_ha"], abs=1e-9
        )


def test_reader_keeps_nuclear_repulsion_and_units_as_separate_fields() -> None:
    data = _artifact(TOY_LENGTHS)
    hamiltonian = h2_hamiltonian_from_artifact(data, 0.74)

    assert isinstance(hamiltonian, ArtifactPointHamiltonian)
    assert hamiltonian.nuclear_repulsion_ha == 0.70
    assert hamiltonian.units == "hartree"
    assert hamiltonian.total_energy_ha == pytest.approx(
        hamiltonian.ground_energy_ha + hamiltonian.nuclear_repulsion_ha, abs=1e-15
    )
    assert hamiltonian.n_qubits == 2
    assert hamiltonian.matrix.shape == (4, 4)
    assert np.allclose(hamiltonian.matrix, hamiltonian.matrix.conj().T)
    assert hamiltonian.pauli_terms["XX"] == XX_COEFFICIENT


def test_reader_selects_the_requested_bond_length_not_just_any() -> None:
    data = _artifact(TOY_LENGTHS)

    at_0_70 = h2_hamiltonian_from_artifact(data, 0.70)
    at_0_74 = h2_hamiltonian_from_artifact(data, 0.74)

    assert at_0_70.bond_length_angstrom == 0.70
    assert at_0_74.bond_length_angstrom == 0.74
    # A lookup that returned the first point regardless would not be detectable
    # from the length alone, so the coefficients must differ -- and the 0.74 point's
    # coefficient must be the one derived for 0.74, not for any other stored point.
    assert at_0_70.pauli_terms != at_0_74.pauli_terms
    expected_at_0_74 = math.sqrt((0.70 - (-1.137)) ** 2 - XX_COEFFICIENT**2)
    assert at_0_74.pauli_terms["ZI"] == pytest.approx(expected_at_0_74, abs=1e-12)
    expected_at_0_70 = math.sqrt((0.70 - (-1.137 + 0.05 * 0.04**2)) ** 2 - XX_COEFFICIENT**2)
    assert at_0_70.pauli_terms["ZI"] == pytest.approx(expected_at_0_70, abs=1e-12)


def test_reader_rejects_an_off_grid_bond_length_and_names_the_nearest() -> None:
    data = _artifact(TOY_LENGTHS)

    with pytest.raises(ValueError, match="not on the artifact's grid"):
        h2_hamiltonian_from_artifact(data, 0.80)
    with pytest.raises(ValueError, match="0.75 A"):
        h2_hamiltonian_from_artifact(data, 0.80)
    with pytest.raises(ValueError, match="never interpolates"):
        h2_artifact_point(data, 1.234)


def test_nearest_stored_bond_length_never_interpolates() -> None:
    """The documented snap rule returns a stored length, never a computed one."""
    data = _artifact(TOY_LENGTHS)
    stored = set(h2_artifact_bond_lengths(data))

    assert nearest_stored_bond_length(data, 0.80) == 0.75
    assert nearest_stored_bond_length(data, 1.234) == 0.75
    assert nearest_stored_bond_length(data, 0.20) == 0.70
    assert nearest_stored_bond_length(data, 0.74) == 0.74

    for request in (0.20, 0.71, 0.80, 2.9):
        assert nearest_stored_bond_length(data, request) in stored


# ---------------------------------------------------------------------------
# Shape validation: the reader must not accept something the writer cannot emit
# ---------------------------------------------------------------------------
def test_reader_reads_the_writers_shape_and_rejects_the_legacy_dense_matrix() -> None:
    """The two shapes are distinct, and each reader accepts only its own."""
    data = _artifact(TOY_LENGTHS)

    hamiltonian = h2_hamiltonian_from_artifact(data, 0.74)
    assert set(hamiltonian.pauli_terms) == {"ZI", "XX"}

    legacy = {"n_qubits": 2, "matrix": np.eye(4).tolist()}
    # The old reader still reads the dense-matrix shape it was written for...
    assert h2_sto3g_hamiltonian_from_json(legacy).shape == (4, 4)
    # ...and the artifact reader refuses it: it has no points at all.
    with pytest.raises(ValueError, match="points"):
        h2_hamiltonian_from_artifact(legacy, 0.74)

    # A legacy point *inside* an otherwise writer-shaped artifact is refused too,
    # and the message names the field the writer would have written.
    mixed = {
        "points": [
            {"bond_length_angstrom": 0.74, "nuclear_repulsion_ha": 0.7, "matrix": np.eye(4).tolist()}
        ]
    }
    with pytest.raises(ValueError, match="pauli_terms"):
        h2_hamiltonian_from_artifact(mixed, 0.74)


@pytest.mark.parametrize(
    ("mutate", "expected"),
    [
        (lambda p: p["points"][TARGET_INDEX].pop("pauli_terms"), "pauli_terms"),
        (lambda p: p["points"][TARGET_INDEX].update(pauli_terms=[]), "non-empty 'pauli_terms'"),
        (lambda p: p["points"][TARGET_INDEX]["pauli_terms"][0].pop("coefficient"), "coefficient"),
        (
            lambda p: p["points"][TARGET_INDEX]["pauli_terms"][0].update(pauli="ZQ"),
            "invalid Pauli string",
        ),
        (
            lambda p: p["points"][TARGET_INDEX]["pauli_terms"][0].update(coefficient="0.5"),
            "real number",
        ),
        (
            lambda p: p["points"][TARGET_INDEX]["pauli_terms"][0].update(
                coefficient=float("nan")
            ),
            "finite",
        ),
        (
            lambda p: p["points"][TARGET_INDEX]["pauli_terms"].append(
                {"pauli": "Z", "coefficient": 0.1}
            ),
            "same length",
        ),
        (lambda p: p["points"][TARGET_INDEX].update(n_qubits=3), "n_qubits"),
        (lambda p: p["points"][TARGET_INDEX].pop("nuclear_repulsion_ha"), "nuclear_repulsion_ha"),
        (
            lambda p: p["points"][TARGET_INDEX].update(nuclear_repulsion_ha=True),
            "nuclear_repulsion_ha",
        ),
        (lambda p: p.pop("points"), "non-empty 'points'"),
        (
            lambda p: p["points"][TARGET_INDEX].pop("bond_length_angstrom"),
            "bond_length_angstrom",
        ),
    ],
    ids=[
        "missing-terms",
        "empty-terms",
        "missing-coefficient",
        "invalid-pauli-character",
        "non-numeric-coefficient",
        "non-finite-coefficient",
        "mixed-widths",
        "n-qubits-mismatch",
        "missing-nuclear-repulsion",
        "boolean-nuclear-repulsion",
        "no-points",
        "point-without-bond-length",
    ],
)
def test_reader_rejects_a_payload_the_writer_could_not_produce(mutate, expected: str) -> None:
    payload = copy.deepcopy(_artifact(TOY_LENGTHS))
    # Guard: every mutation below targets the point the reader is asked for, so a
    # mutation that lands on a different point cannot make this test pass vacuously.
    assert payload["points"][TARGET_INDEX]["bond_length_angstrom"] == 0.74
    mutate(payload)

    with pytest.raises(ValueError, match=expected):
        h2_hamiltonian_from_artifact(payload, 0.74)


def test_reader_sums_duplicate_pauli_strings() -> None:
    """A repeated Pauli string is added, not silently overwritten."""
    point = _analytic_point(0.74)
    point["pauli_terms"].append({"pauli": "XX", "coefficient": 0.01})

    terms = h2_artifact_pauli_terms(point)

    assert terms["XX"] == pytest.approx(XX_COEFFICIENT + 0.01, abs=1e-15)


def test_reader_bond_lengths_are_the_58_point_grid() -> None:
    """The reader reports the grid the writer now produces, 0.74 included."""
    data = _artifact(EXPECTED_BOND_LENGTHS)

    lengths = h2_artifact_bond_lengths(data)

    assert lengths == EXPECTED_BOND_LENGTHS
    assert len(lengths) == 58
    assert 0.74 in lengths


# ---------------------------------------------------------------------------
# The cross-check item 2 asks for: both paths agree on the same fixture
# ---------------------------------------------------------------------------
def test_reader_and_validator_recompute_the_same_energies(tmp_path: Path) -> None:
    """The core path and the validator's independent NumPy path must agree.

    The fixture is a complete artifact (58 points, the writer's field set), so this
    also proves the validator accepts the shape the reader consumes -- the two ends
    of the mismatch fixed here.
    """
    data = _artifact(EXPECTED_BOND_LENGTHS)
    validator = _validator_module()

    for point in data["points"]:
        hamiltonian = h2_hamiltonian_from_artifact(data, point["bond_length_angstrom"])
        recomputed = validator.ground_energy_from_terms(point["pauli_terms"])

        assert hamiltonian.ground_energy_ha == pytest.approx(recomputed, abs=1e-12)
        assert hamiltonian.total_energy_ha == pytest.approx(
            recomputed + point["nuclear_repulsion_ha"], abs=1e-12
        )
        assert hamiltonian.total_energy_ha == pytest.approx(point["fci_energy_ha"], abs=1e-9)


def test_validator_accepts_the_synthetic_artifact_the_reader_consumes(tmp_path: Path) -> None:
    path = _write(tmp_path, _artifact(EXPECTED_BOND_LENGTHS))

    result = subprocess.run(
        [sys.executable, str(VALIDATOR), str(path)],
        cwd=str(REPO_ROOT),
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stderr
    assert "OK: every check passed." in result.stdout
    assert "points                 : 58" in result.stdout
    assert "r=0.74 A" in result.stdout
