"""STO-3G spatial-orbital counts derived from a molecular formula.

FORMULA-DERIVED, NOT YET CONFIRMED AGAINST PySCF. The expected totals here are
computed inside this file by independent routes (a separately written formula
parser, plus hand-read atom counts summed against a separately transcribed basis
table) rather than copied from the implementation. That guards against a bug in
``app.core.resources``; it does **not** verify the basis itself. Nothing here has
been checked against ``mol.nao``, because the offline PySCF notebook
(``scripts/precompute_hamiltonians.ipynb``) has never been run.
"""

from __future__ import annotations

import re

import pytest

from app.core.resources import (
    STO3G_BASIS_FUNCTIONS_PER_ATOM,
    parse_formula_counts,
    sto3g_spatial_orbitals,
)

# Independent transcription of the STO-3G basis: H/He contribute one 1s function,
# Li to Ne contribute 1s, 2s, 2px, 2py, 2pz = five, and Na to Ar add 3s and
# 3px, 3py, 3pz on top of that = nine.
BASIS_FUNCTIONS_PER_ATOM = {
    "H": 1,
    "He": 1,
    "Li": 5,
    "Be": 5,
    "B": 5,
    "C": 5,
    "N": 5,
    "O": 5,
    "F": 5,
    "Ne": 5,
    "Na": 9,
    "Mg": 9,
    "Al": 9,
    "Si": 9,
    "P": 9,
    "S": 9,
    "Cl": 9,
    "Ar": 9,
}

# Atom counts read off each formula by hand, independent of any parser.
FORMULA_ATOM_COUNTS = {
    "H2": {"H": 2},
    "LiH": {"Li": 1, "H": 1},
    "H2O": {"H": 2, "O": 1},
    "CH4": {"C": 1, "H": 4},
    "C9H8O4": {"C": 9, "H": 8, "O": 4},
    "NaCl": {"Na": 1, "Cl": 1},
    "CH3Cl": {"C": 1, "H": 3, "Cl": 1},
}

_FORMULA_TOKEN = re.compile(r"([A-Z][a-z]?)(\d*)")


def _independent_parse(formula: str) -> dict[str, int]:
    """Parser written here, separate from the one under test."""
    counts: dict[str, int] = {}
    for symbol, digits in _FORMULA_TOKEN.findall(formula):
        counts[symbol] = counts.get(symbol, 0) + (int(digits) if digits else 1)
    return counts


def _independent_total(counts: dict[str, int]) -> int:
    return sum(count * BASIS_FUNCTIONS_PER_ATOM[element] for element, count in counts.items())


@pytest.mark.parametrize("formula", sorted(FORMULA_ATOM_COUNTS))
def test_total_matches_two_independent_routes(formula: str) -> None:
    # Route 1: a parser written separately in this file.
    parsed = _independent_parse(formula)
    assert parsed == FORMULA_ATOM_COUNTS[formula], "the two formula readers disagree"
    # Route 2: sum the hand-read atom counts, touching no parser at all.
    expected = _independent_total(FORMULA_ATOM_COUNTS[formula])
    assert _independent_total(parsed) == expected

    result = sto3g_spatial_orbitals(formula)
    print(
        f"STO-3G {formula}: {result['atom_counts']} -> "
        f"{result['basis_functions_per_element']} = {result['spatial_orbitals']} "
        f"spatial orbitals (unverified against PySCF)"
    )

    assert result["atom_counts"] == FORMULA_ATOM_COUNTS[formula]
    assert result["spatial_orbitals"] == expected


def test_basis_table_matches_the_transcribed_definition() -> None:
    assert STO3G_BASIS_FUNCTIONS_PER_ATOM == BASIS_FUNCTIONS_PER_ATOM


def test_formula_parser_handles_implicit_and_explicit_counts() -> None:
    assert parse_formula_counts("H2") == {"H": 2}
    assert parse_formula_counts("LiH") == {"Li": 1, "H": 1}
    assert parse_formula_counts("H2O") == {"H": 2, "O": 1}
    assert parse_formula_counts("C9H8O4") == {"C": 9, "H": 8, "O": 4}
    assert parse_formula_counts("C6H6") == {"C": 6, "H": 6}
    assert parse_formula_counts("NaCl") == {"Na": 1, "Cl": 1}
    assert parse_formula_counts("CH3Cl") == {"C": 1, "H": 3, "Cl": 1}


def test_result_is_flagged_as_unverified() -> None:
    result = sto3g_spatial_orbitals("H2O")
    assert result["verified_against_pyscf"] is False
    assert "PySCF" in result["note"]


@pytest.mark.parametrize("formula", ["", "   ", "H2O)", "c6h6", "Xx2", "Fe2O3", "H0"])
def test_rejects_unparseable_or_out_of_scope_formulas(formula: str) -> None:
    with pytest.raises(ValueError):
        sto3g_spatial_orbitals(formula)
