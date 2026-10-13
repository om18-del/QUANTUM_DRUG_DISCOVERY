"""Tests for ``scripts/validate_hamiltonian_file.py``.

The real artifact (``app/data/hamiltonians/h2_sto3g.json``) is generated offline
in Colab or WSL and is committed; ``tests/test_h2_hamiltonian_artifact.py`` checks
that file. These tests instead exercise the validator against small synthetic
payloads written to a temporary directory, because most of what matters here is
what the validator rejects. The synthetic curve is a
parabola with its minimum at 0.74 A, which is now on the 58-point grid (0.20 to
3.00 in 0.05 steps, 57 uniform points, plus 0.74 A -- the UI default).

The validator is run as a subprocess rather than imported, because its contract
is a process exit code plus a message on stderr -- importing the functions would
skip exactly the part the caller depends on.
"""

from __future__ import annotations

import copy
import json
import subprocess
import sys
from pathlib import Path
from typing import Any

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
VALIDATOR = REPO_ROOT / "scripts" / "validate_hamiltonian_file.py"
# Same 58-point grid the writer and the validator use: 57 uniform points of 0.20 ..
# 3.00 angstrom in 0.05 steps, plus 0.74 angstrom.
UNIFORM_POINTS = 57
EXPECTED_BOND_LENGTHS = sorted(
    {*[round(0.20 + 0.05 * i, 2) for i in range(UNIFORM_POINTS)], 0.74}
)

NUCLEAR_REPULSION_HA = 0.70
HF_MINUS_FCI_HA = 0.02


def _payload() -> dict[str, Any]:
    """A physically consistent artifact: FCI below HF, qubit operator matching FCI.

    ``ZI`` is diagonal, so ``c * ZI`` has ground energy ``-|c|``; choosing
    ``c = nuclear - fci`` (positive, because nuclear repulsion exceeds the FCI
    energy everywhere here) makes ``ground + nuclear`` equal the stored FCI
    exactly, up to floating point.
    """
    points: list[dict[str, Any]] = []
    for bond_length in EXPECTED_BOND_LENGTHS:
        fci = -1.137 + 0.05 * (bond_length - 0.74) ** 2
        hf = fci + HF_MINUS_FCI_HA
        points.append(
            {
                "bond_length_angstrom": bond_length,
                "nuclear_repulsion_ha": NUCLEAR_REPULSION_HA,
                "hf_energy_ha": hf,
                "fci_energy_ha": fci,
                "qubit_energy_ha": fci - NUCLEAR_REPULSION_HA,
                "qubit_energy_plus_nuclear_ha": fci,
                "deviation_from_fci_ha": 0.0,
                "n_qubits": 2,
                "pauli_terms": [
                    {"pauli": "ZI", "coefficient": NUCLEAR_REPULSION_HA - fci}
                ],
                "nuclear_repulsion_in_qubit_operator": False,
                "energy_span_two_state_ha": None,
            }
        )

    return {
        "molecule": "H2",
        "basis": "sto-3g",
        "mapping": "parity",
        "two_qubit_reduction": True,
        "n_qubits": 2,
        "units": "hartree",
        "generated_by": "tests/test_hamiltonian_validator.py (synthetic)",
        "library_versions": {"qiskit": "1.4.6", "qiskit-nature": "0.8.0", "pyscf": "2.14.0"},
        "convention": "synthetic fixture; not a real PySCF run",
        "max_deviation_from_fci_ha": 0.0,
        "point_count": len(points),
        "points": points,
    }


def _write(tmp_path: Path, payload: Any, name: str = "synthetic.json") -> Path:
    path = tmp_path / name
    if isinstance(payload, str):
        path.write_text(payload, encoding="utf-8")
    else:
        path.write_text(json.dumps(payload), encoding="utf-8")
    return path


def _run(path: Path) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(VALIDATOR), str(path)],
        cwd=str(REPO_ROOT),
        capture_output=True,
        text=True,
    )


def test_valid_synthetic_file_passes(tmp_path: Path) -> None:
    """Positive control: without this, the failure cases below prove nothing."""
    result = _run(_write(tmp_path, _payload()))

    assert result.returncode == 0, result.stderr
    assert "OK: every check passed." in result.stdout
    # The equilibrium energy must be printed, not merely asserted.
    assert "nearest point to 0.74 A" in result.stdout
    assert "FCI=" in result.stdout


def test_printed_equilibrium_energy_is_the_parabola_minimum(tmp_path: Path) -> None:
    result = _run(_write(tmp_path, _payload()))

    assert result.returncode == 0, result.stderr
    assert "r=0.74 A" in result.stdout
    assert "-1.137000" in result.stdout


def test_missing_file_fails_and_is_not_created(tmp_path: Path) -> None:
    missing = tmp_path / "does_not_exist.json"

    result = _run(missing)

    assert result.returncode != 0
    assert "does not exist" in result.stderr
    assert not missing.exists(), "the validator must not create or invent the artifact"


def test_malformed_json_fails(tmp_path: Path) -> None:
    result = _run(_write(tmp_path, "{not json at all"))

    assert result.returncode != 0
    assert "not valid JSON" in result.stderr


def _drop_a_point(payload: dict[str, Any]) -> dict[str, Any]:
    payload["points"].pop()
    payload["point_count"] = len(payload["points"])
    return payload


def _duplicate_a_bond_length(payload: dict[str, Any]) -> dict[str, Any]:
    payload["points"][5] = copy.deepcopy(payload["points"][4])
    payload["points"][5]["bond_length_angstrom"] = payload["points"][4]["bond_length_angstrom"]
    return payload


def _off_grid_step(payload: dict[str, Any]) -> dict[str, Any]:
    payload["points"][3]["bond_length_angstrom"] = 0.33
    return payload


def _missing_point_field(payload: dict[str, Any]) -> dict[str, Any]:
    del payload["points"][7]["hf_energy_ha"]
    return payload


def _missing_top_level_field(payload: dict[str, Any]) -> dict[str, Any]:
    del payload["library_versions"]
    return payload


def _fci_above_hf(payload: dict[str, Any]) -> dict[str, Any]:
    point = payload["points"][20]
    point["hf_energy_ha"] = point["fci_energy_ha"] - 0.01
    return payload


def _qubit_operator_does_not_reproduce_fci(payload: dict[str, Any]) -> dict[str, Any]:
    payload["points"][20]["pauli_terms"][0]["coefficient"] += 0.05
    return payload


def _recorded_max_deviation_too_large(payload: dict[str, Any]) -> dict[str, Any]:
    payload["max_deviation_from_fci_ha"] = 0.5
    return payload


def _wrong_mapping(payload: dict[str, Any]) -> dict[str, Any]:
    payload["mapping"] = "jordan-wigner"
    return payload


@pytest.mark.parametrize(
    ("mutate", "expected"),
    [
        (_drop_a_point, "expected exactly 58 bond lengths"),
        (_duplicate_a_bond_length, "duplicate bond length"),
        (_off_grid_step, "not 0.20 to 3.00 A in 0.05 steps"),
        (_missing_point_field, "missing field"),
        (_missing_top_level_field, "missing top-level field"),
        (_fci_above_hf, "above HF"),
        (_qubit_operator_does_not_reproduce_fci, "qubit operator + nuclear repulsion"),
        (_recorded_max_deviation_too_large, "max_deviation_from_fci_ha"),
        (_wrong_mapping, "parity"),
    ],
    ids=[
        "wrong-point-count",
        "duplicate-bond-length",
        "off-grid-step",
        "missing-point-field",
        "missing-top-level-field",
        "fci-above-hf",
        "qubit-operator-mismatch",
        "recorded-deviation-too-large",
        "wrong-mapping",
    ],
)
def test_invalid_payloads_fail_with_a_clear_message(
    tmp_path: Path, mutate, expected: str
) -> None:
    payload = mutate(_payload())

    result = _run(_write(tmp_path, payload))

    assert result.returncode != 0, f"expected a non-zero exit for {expected!r}"
    assert "FAIL:" in result.stderr
    assert expected in result.stderr, f"expected {expected!r} in:\n{result.stderr}"


def test_usage_error_for_extra_arguments(tmp_path: Path) -> None:
    result = subprocess.run(
        [sys.executable, str(VALIDATOR), "a.json", "b.json"],
        cwd=str(REPO_ROOT),
        capture_output=True,
        text=True,
    )

    assert result.returncode == 2
    assert "usage:" in result.stderr
