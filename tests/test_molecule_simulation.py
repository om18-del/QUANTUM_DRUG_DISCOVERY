"""Unit tests for the SMILES simulation method field and PDB-id detection.

These run without RDKit, because the machine-readable ``method`` code is derived
from the structure payload alone. The RDKit-dependent paths are covered in
``tests/test_endpoints.py`` behind an availability guard.
"""

from __future__ import annotations

import pytest

from app.services.molecule_simulation import (
    METHOD_2D_FALLBACK,
    METHOD_3D_MMFF94,
    METHOD_3D_UFF,
    METHOD_3D_UNMINIMISED,
    METHOD_CODES,
    METHOD_UNAVAILABLE,
    _method_code,
    looks_like_pdb_id,
)
from app.services.rdkit_helper import (
    ABSOLUTE_ENERGY_NOTE,
    POPULATION_ESTIMATE_NOTE,
    boltzmann_populations,
    relative_energies,
)


@pytest.mark.parametrize(
    "structure,expected",
    [
        ({"is_3d": True, "available": True, "force_field": "MMFF94"}, METHOD_3D_MMFF94),
        ({"is_3d": True, "available": True, "force_field": "UFF"}, METHOD_3D_UFF),
        ({"is_3d": True, "available": True, "force_field": None}, METHOD_3D_UNMINIMISED),
        ({"is_3d": False, "available": True, "force_field": None}, METHOD_2D_FALLBACK),
        ({"is_3d": False, "available": False, "force_field": None}, METHOD_UNAVAILABLE),
        ({}, METHOD_UNAVAILABLE),
    ],
)
def test_method_code_is_machine_readable(structure: dict, expected: str) -> None:
    assert _method_code(structure) == expected
    assert _method_code(structure) in METHOD_CODES


def test_method_code_never_returns_prose() -> None:
    for structure in (
        {"is_3d": True, "available": True, "force_field": "MMFF94"},
        {"is_3d": False, "available": True},
        {},
    ):
        code = _method_code(structure)
        assert " " not in code
        assert code == code.lower()


# --- A5: energies are relative, populations are estimates -------------------


def test_relative_energies_shifts_lowest_to_zero() -> None:
    assert relative_energies([10.0, 11.0, 12.5]) == [0.0, 1.0, 2.5]
    assert relative_energies([5.0]) == [0.0]
    assert relative_energies([]) == []
    # Order of the input must not matter.
    assert relative_energies([12.5, 10.0, 11.0]) == [2.5, 0.0, 1.0]


def test_relative_energies_are_never_negative() -> None:
    assert all(e >= 0.0 for e in relative_energies([-1.0, 4.0, 0.5]))


def test_boltzmann_populations_sum_to_one() -> None:
    pops = boltzmann_populations([0.0, 1.0, 2.0])
    assert len(pops) == 3
    assert sum(pops) == pytest.approx(1.0, abs=0.002)
    # Lower energy must never get the smaller share.
    assert pops[0] > pops[1] > pops[2]


def test_boltzmann_populations_degenerate_conformers_are_equal() -> None:
    assert boltzmann_populations([3.0, 3.0, 3.0]) == [0.3333, 0.3333, 0.3333]


def test_boltzmann_populations_empty_and_large_gap() -> None:
    assert boltzmann_populations([]) == []
    pops = boltzmann_populations([0.0, 100.0])
    assert pops[0] == pytest.approx(1.0, abs=0.001)
    assert pops[1] == pytest.approx(0.0, abs=0.001)


def test_energy_and_population_notes_say_they_are_not_absolute() -> None:
    assert "not comparable" in ABSOLUTE_ENERGY_NOTE
    assert "relative" in ABSOLUTE_ENERGY_NOTE
    assert "estimate" in POPULATION_ESTIMATE_NOTE.lower()


@pytest.mark.parametrize("value", ["1HSG", "6LU7", "9xyz", "1hsg"])
def test_looks_like_pdb_id_accepts(value: str) -> None:
    assert looks_like_pdb_id(value) is True


@pytest.mark.parametrize("value", ["CCO", "aspirin", "CC(=O)Oc1ccccc1C(=O)O", "", "ABC", "12345"])
def test_looks_like_pdb_id_rejects(value: str) -> None:
    assert looks_like_pdb_id(value) is False
