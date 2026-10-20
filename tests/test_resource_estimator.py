"""Resource estimates counted from a molecular formula (STO-3G).

Every expected value below is written out by hand from two independent things:
the minimal-basis rule (H/He contribute 1 spatial orbital, Li-Ne 5, Na-Ar 9) and
the neutral atomic numbers. Nothing is copied from the implementation, so an
arithmetic bug in ``app.core.resources`` cannot define its own expected answer.

What these tests do **not** do is confirm the basis itself: the counts are
counted by rule. Five formulas (H2, LiH, H2O, CH4, C9H8O4) were separately
cross-checked against PySCF ``mol.nao`` in Colab and are listed in
``PYSCF_VERIFIED_COUNTS``; those payloads say ``verified_against_pyscf: True``.
Every other formula stays ``False``, and tests below assert both cases.
"""

from __future__ import annotations

import json
import re

import pytest

from app.core.resources import (
    ATOMIC_NUMBERS,
    INTERACTIVE_CAP_QUBITS,
    LAPTOP_STATEVECTOR_CEILING_QUBITS,
    PYSCF_VERIFIED_COUNTS,
    STO3G_BASIS_FUNCTIONS_PER_ATOM,
    ActiveOrbitalsError,
    InvalidFormulaError,
    ResourceEstimateError,
    UnsupportedElementError,
    electron_count,
    estimate_resources,
    format_byte_count,
    normalise_formula,
    sto3g_spatial_orbitals,
)

# The project's worked examples: formula, spatial orbitals, electrons, full
# qubits, reduced qubits, verdict, human-readable statevector memory.
# (orbitals = hand-summed basis functions; electrons = hand-summed atomic
# numbers; full = 2 x orbitals; reduced = full - 2 for the parity mapping with
# two-qubit reduction; memory = 2**reduced x 16 bytes.)
EXAMPLES = [
    ("H2", 2, 2, 4, 2, "simulable_in_app", "64 bytes"),
    ("LiH", 6, 4, 12, 10, "precompute_only", "16.0 KB"),
    ("H2O", 7, 10, 14, 12, "precompute_only", "64.0 KB"),
    ("CH4", 9, 10, 18, 16, "precompute_only", "1.0 MB"),
    ("C9H8O4", 73, 94, 146, 144, "beyond_classical_simulation", "more than 10^15 bytes"),
]


@pytest.mark.parametrize(
    "formula,orbitals,electrons,full,reduced,verdict,memory",
    EXAMPLES,
)
def test_worked_examples(
    formula: str,
    orbitals: int,
    electrons: int,
    full: int,
    reduced: int,
    verdict: str,
    memory: str,
) -> None:
    payload = estimate_resources(formula)

    assert payload["status"] == "computed"
    assert payload["basis"] == "STO-3G"
    assert payload["n_spatial_orbitals"] == orbitals
    assert payload["n_electrons"] == electrons
    assert payload["n_qubits_full"] == full
    assert payload["n_qubits_after_parity_reduction"] == reduced
    assert payload["parity_reduction_removed_qubits"] == 2
    assert payload["feasibility"]["verdict"] == verdict
    assert payload["statevector_memory_human"] == memory


@pytest.mark.parametrize(
    "formula,orbitals,electrons,full,reduced,verdict,memory",
    EXAMPLES,
)
def test_arithmetic_is_internally_consistent(
    formula: str,
    orbitals: int,
    electrons: int,
    full: int,
    reduced: int,
    verdict: str,
    memory: str,
) -> None:
    """The identities, checked independently of the literal expected values."""
    payload = estimate_resources(formula)

    assert payload["n_spatial_orbitals"] == sto3g_spatial_orbitals(formula)["spatial_orbitals"]
    assert payload["n_electrons"] == electron_count(formula)
    assert payload["n_qubits_full"] == 2 * payload["n_spatial_orbitals"]
    assert payload["n_qubits_after_parity_reduction"] == payload["n_qubits_full"] - 2
    assert payload["statevector_memory_bytes"] == 2 ** payload["n_qubits_after_parity_reduction"] * 16
    assert payload["hilbert_dimension"] == str(2 ** payload["n_qubits_after_parity_reduction"])
    assert payload["subject"]["formula"] == formula


def test_hilbert_dimension_is_a_string_and_memory_is_an_int() -> None:
    """A 2**144 dimension would lose precision as a JSON number in a browser."""
    payload = estimate_resources("C9H8O4")

    assert isinstance(payload["hilbert_dimension"], str)
    assert payload["hilbert_dimension"] == str(2**144)
    assert len(payload["hilbert_dimension"]) == 44
    assert isinstance(payload["statevector_memory_bytes"], int)
    assert payload["statevector_memory_bytes"] == 2**144 * 16


@pytest.mark.parametrize("formula", [example[0] for example in EXAMPLES])
def test_payload_is_json_safe_and_round_trips(formula: str) -> None:
    payload = estimate_resources(formula)

    text = json.dumps(payload)
    assert json.loads(text) == payload
    # The huge integer survives the round trip exactly (a browser could not hold
    # it as a float, which is why the dimension is also given as a string).
    assert json.loads(text)["statevector_memory_bytes"] == payload["statevector_memory_bytes"]


@pytest.mark.parametrize("formula", [example[0] for example in EXAMPLES])
def test_every_figure_carries_an_exact_or_estimated_tag(formula: str) -> None:
    payload = estimate_resources(formula)
    tags = payload["tags"]

    assert tags, "a payload with untagged figures is not allowed"
    for key, tag in tags.items():
        assert tag in ("exact", "estimated"), (key, tag)
        # The tag must name a field that is actually there, so a renamed field
        # cannot leave a tag pointing at nothing.
        assert key.split(".")[0] in payload, key

    for field in (
        "n_spatial_orbitals",
        "n_electrons",
        "n_qubits_full",
        "n_qubits_after_parity_reduction",
        "statevector_memory_bytes",
        "statevector_memory_human",
        "hilbert_dimension",
    ):
        assert tags[field] == "exact", field

    # The one deliberately rough number and the one judgement are estimated.
    assert tags["feasibility.laptop_statevector_ceiling_qubits"] == "estimated"
    assert tags["feasibility.verdict"] == "estimated"


def test_no_invented_pauli_or_circuit_numbers() -> None:
    """A formula cannot give a term count, so none is reported."""
    payload = estimate_resources("C9H8O4")

    for invented in ("hamiltonian_terms", "estimated_circuit_gates", "estimated_circuit_depth"):
        assert invented not in payload, invented
    assert "roughly N^4" in payload["notes"]["pauli_term_scaling"]
    assert "estimated" in payload["notes"]["pauli_term_scaling"]
    assert payload["tags"]["notes.pauli_term_scaling"] == "estimated"


def test_states_what_it_is_and_is_not() -> None:
    payload = estimate_resources("C2H6O")

    assert payload["basis"] == "STO-3G"
    assert payload["method"] == "computed_by_formula"
    assert payload["classical_simulation"] is True
    assert "no PySCF run" in payload["method_note"]
    assert "project limits and rough estimates" in payload["method_note"]
    assert "not modelled" in payload["notes"]["spin_and_charge"]
    # Ethanol is not one of the five cross-checked formulas.
    assert payload["verified_against_pyscf"] is False
    assert payload["feasibility"]["interactive_cap_qubits"] == INTERACTIVE_CAP_QUBITS
    assert (
        payload["feasibility"]["laptop_statevector_ceiling_qubits"]
        == LAPTOP_STATEVECTOR_CEILING_QUBITS
    )


def test_the_five_pyscf_verified_formulas_are_flagged_true() -> None:
    """Exactly the five recorded formulas report verification, with matching counts."""
    assert set(PYSCF_VERIFIED_COUNTS) == {"H2", "LiH", "H2O", "CH4", "C9H8O4"}

    for formula, recorded in PYSCF_VERIFIED_COUNTS.items():
        payload = estimate_resources(formula)
        assert payload["verified_against_pyscf"] is True, formula
        assert payload["n_spatial_orbitals"] == recorded["spatial_orbitals"], formula
        assert payload["n_electrons"] == recorded["electrons"], formula
        assert payload["subject"]["normalised_formula"] == formula
        # One canonicaliser: the payload field and the helper agree.
        assert normalise_formula(formula) == formula
        # The note names the check that was actually done, and invents no version.
        note = payload["notes"]["verification"]
        assert "PySCF mol.nao" in note
        assert "No PySCF version was recorded" in note
        assert not re.search(r"PySCF\s+v?\d", note)


def test_unverified_formulas_stay_false_with_the_explanation() -> None:
    for formula in ("C2H6O", "C6H6", "CH3Cl", "Ne", "NaCl"):
        payload = estimate_resources(formula)
        assert payload["verified_against_pyscf"] is False, formula
        assert "unverified" in payload["notes"]["verification"], formula


def test_verification_ignores_how_the_formula_is_spelled() -> None:
    """The lookup is on the normalised (Hill-order) formula, not the raw text."""
    assert normalise_formula("OH2") == "H2O"
    assert normalise_formula("H4C") == "CH4"
    assert normalise_formula("H2C9O4H6") == "C9H8O4"
    # A verified formula keeps its recorded spelling even when Hill order would
    # have spelled it differently.
    assert normalise_formula("HLi") == "LiH"
    assert normalise_formula("C6H6") == "C6H6"
    assert estimate_resources("C6H6")["subject"]["normalised_formula"] == "C6H6"

    renamed_water = estimate_resources("OH2")
    assert renamed_water["subject"]["formula"] == "OH2"
    assert renamed_water["subject"]["normalised_formula"] == "H2O"
    assert renamed_water["n_spatial_orbitals"] == 7
    assert renamed_water["verified_against_pyscf"] is True


def test_active_space_does_not_change_the_verification_flag() -> None:
    for formula in ("H2O", "C2H6O"):
        plain = estimate_resources(formula)
        with_active = estimate_resources(formula, active_orbitals=2)
        assert with_active["verified_against_pyscf"] is plain["verified_against_pyscf"], formula
        assert with_active["active_space"] is not None


def test_verdict_boundaries() -> None:
    """Exactly at each cap is inside it; one past is not."""
    # 2 spatial orbitals -> 2 reduced qubits (H2), well inside the app cap.
    assert estimate_resources("H2")["feasibility"]["verdict"] == "simulable_in_app"
    # 5 spatial orbitals -> 8 reduced qubits: exactly at the 8-qubit app cap.
    neon = estimate_resources("Ne")
    assert neon["n_qubits_after_parity_reduction"] == INTERACTIVE_CAP_QUBITS
    assert neon["feasibility"]["verdict"] == "simulable_in_app"
    # 6 spatial orbitals -> 10 reduced qubits: past the cap, inside the ceiling.
    assert estimate_resources("LiH")["feasibility"]["verdict"] == "precompute_only"
    # 21 spatial orbitals -> 40 reduced qubits: past the rough ceiling.
    ethanol = estimate_resources("C2H6O")
    assert ethanol["n_qubits_after_parity_reduction"] > LAPTOP_STATEVECTOR_CEILING_QUBITS
    assert ethanol["feasibility"]["verdict"] == "beyond_classical_simulation"
    assert ethanol["feasibility"]["reason"]


def test_unsupported_element_names_the_supported_ones() -> None:
    with pytest.raises(UnsupportedElementError) as excinfo:
        estimate_resources("Fe2O3")

    message = str(excinfo.value)
    assert "Fe" in message
    assert "Ar" in message and "H," in message
    assert excinfo.value.elements == ("Fe",)
    # Still a ValueError, so callers written against the earlier build keep working.
    assert isinstance(excinfo.value, ValueError)


@pytest.mark.parametrize("formula", ["", "   ", "c6h6", "h2o", "H2O)", "H2O(", "H0", "2H2", "H2.O"])
def test_rejects_unparseable_formulas(formula: str) -> None:
    with pytest.raises(InvalidFormulaError):
        estimate_resources(formula)


def test_formula_errors_share_one_base_class() -> None:
    assert issubclass(InvalidFormulaError, ResourceEstimateError)
    assert issubclass(UnsupportedElementError, ResourceEstimateError)
    assert issubclass(ResourceEstimateError, ValueError)
    with pytest.raises(ResourceEstimateError):
        estimate_resources("H2O)")
    with pytest.raises(ResourceEstimateError):
        estimate_resources("Xx2")


def test_active_space_block() -> None:
    payload = estimate_resources("H2O", active_orbitals=4)
    active = payload["active_space"]

    assert active["active_orbitals"] == 4
    assert active["n_spatial_orbitals"] == 7
    assert active["n_qubits"] == 8
    assert active["n_qubits_after_parity_reduction"] == 6
    assert active["statevector_memory_bytes"] == 2**6 * 16
    assert active["statevector_memory_human"] == "1.0 KB"
    assert active["feasibility"]["verdict"] == "simulable_in_app"
    assert "approximation" in active["note"]
    assert payload["tags"]["active_space.n_qubits"] == "estimated"
    assert payload["tags"]["active_space.feasibility.verdict"] == "estimated"
    # The full-molecule figures are unaffected by the active-space block.
    assert payload["n_qubits_after_parity_reduction"] == 12


def test_active_space_is_absent_unless_asked_for() -> None:
    assert estimate_resources("H2O")["active_space"] is None
    assert "active_space.n_qubits" not in estimate_resources("H2O")["tags"]


def test_active_space_size_lets_aspirin_fit_the_app() -> None:
    """Cutting aspirin down to k orbitals is exactly how it could ever be run."""
    payload = estimate_resources("C9H8O4", active_orbitals=4)

    assert payload["feasibility"]["verdict"] == "beyond_classical_simulation"
    assert payload["active_space"]["n_qubits_after_parity_reduction"] == 6
    assert payload["active_space"]["feasibility"]["verdict"] == "simulable_in_app"


@pytest.mark.parametrize("k", [0, -1, 8, 99])
def test_active_space_out_of_range(k: int) -> None:
    with pytest.raises(ActiveOrbitalsError) as excinfo:
        estimate_resources("H2O", active_orbitals=k)
    assert "between 1 and 7" in str(excinfo.value)


@pytest.mark.parametrize("k", [2.5, "4", True])
def test_active_space_must_be_an_integer(k) -> None:
    with pytest.raises(ActiveOrbitalsError):
        estimate_resources("H2O", active_orbitals=k)


def test_active_space_k_maximum_is_allowed() -> None:
    payload = estimate_resources("H2O", active_orbitals=7)
    assert payload["active_space"]["n_qubits"] == payload["n_qubits_full"]


def test_odd_electron_count_is_allowed_and_flagged() -> None:
    odd = estimate_resources("CH3")
    assert odd["n_electrons"] == 9
    assert "open_shell" in odd["notes"]
    assert "odd electron count (9)" in odd["notes"]["open_shell"]

    even = estimate_resources("CH4")
    assert even["n_electrons"] == 10
    assert "open_shell" not in even["notes"]


def test_charge_is_not_modeled() -> None:
    """The electron count takes no charge, so it is the neutral count."""
    assert estimate_resources("H2O")["n_electrons"] == 10
    assert estimate_resources("NaCl")["n_electrons"] == 11 + 17
    assert estimate_resources("NaCl")["n_spatial_orbitals"] == 9 + 9


def test_basis_table_and_atomic_numbers_agree() -> None:
    assert set(STO3G_BASIS_FUNCTIONS_PER_ATOM) == set(ATOMIC_NUMBERS)
    assert STO3G_BASIS_FUNCTIONS_PER_ATOM["H"] == 1
    assert STO3G_BASIS_FUNCTIONS_PER_ATOM["Ne"] == 5
    assert STO3G_BASIS_FUNCTIONS_PER_ATOM["Na"] == 9
    assert STO3G_BASIS_FUNCTIONS_PER_ATOM["Ar"] == 9
    assert ATOMIC_NUMBERS["Ar"] == 18


@pytest.mark.parametrize(
    "n_bytes,expected",
    [
        (0, "0 bytes"),
        (64, "64 bytes"),
        (1023, "1023 bytes"),
        (1024, "1.0 KB"),
        (16384, "16.0 KB"),
        (1024**2, "1.0 MB"),
        (1024**3, "1.0 GB"),
        (1024**4, "1.0 TB"),
        (1024**5 - 1, "1024.0 TB"),
        (1024**5, "more than 10^15 bytes"),
        (2**148, "more than 10^15 bytes"),
    ],
)
def test_memory_formatting(n_bytes: int, expected: str) -> None:
    assert format_byte_count(n_bytes) == expected
