"""Resource estimation for a molecule we never simulate.

Two different things live here and they must not be confused:

* ``estimate_h2_resources`` looks at an actual Hamiltonian **matrix** and reports
  its dimension and a rough gate/depth figure. It is an order-of-magnitude guide
  for the tiny H2 demo problem, not a compiled resource estimate.
* ``estimate_resources`` answers "what would this molecule need?" from its
  **molecular formula alone** for the STO-3G minimal basis. Nothing is built and
  no electronic-structure program runs: the orbital count is a rule-based count
  (see the table below), the electron count is the sum of the neutral atomic
  numbers, and the memory figure is the closed-form size of a complex128
  statevector. No Pauli-term count is invented anywhere.

Every field of an ``estimate_resources`` payload is tagged in its ``tags`` dict as
``exact`` (computed by a closed-form rule in code: the orbital and electron
counts, the qubit counts and the memory formula) or ``estimated`` (a judgement or
a rough limit: the feasibility verdict, the laptop statevector ceiling, and
anything about an active space the caller chose).
"""

from __future__ import annotations

import re
from typing import Any

import numpy as np
from numpy.typing import NDArray


def estimate_h2_resources(hamiltonian: NDArray) -> dict[str, Any]:
    """Rough resource estimate for a small Hamiltonian.

    This is not a compiled resource estimate. It is a rough order-of-magnitude
    guide for the demonstration.
    """
    if hamiltonian.ndim != 2 or hamiltonian.shape[0] != hamiltonian.shape[1]:
        raise ValueError("Hamiltonian must be a square matrix")

    dim = hamiltonian.shape[0]
    n_qubits = int(np.log2(dim))
    if 2**n_qubits != dim:
        raise ValueError("Hamiltonian dimension is not a power of two")

    n_terms = int(np.count_nonzero(np.abs(hamiltonian) > 1e-12))

    # Very rough circuit estimate for a tiny 2-qubit demo ansatz
    gates = 4 * n_qubits + n_terms
    depth = 2 * n_qubits + 2

    return {
        "n_qubits": n_qubits,
        "hamiltonian_dim": dim,
        "hamiltonian_terms": n_terms,
        "estimated_circuit_gates": gates,
        "estimated_circuit_depth": depth,
        "note": "Estimate only, not a full compilation",
    }


# ---------------------------------------------------------------------------
# STO-3G spatial orbitals and electrons from a molecular formula
# ---------------------------------------------------------------------------
#
# The standard minimal-basis rule: STO-3G has exactly one contracted basis
# function per occupied atomic shell, and a shell contributes one function per
# orbital (1 for s, 3 for p). So, element by element:
#
#   H, He       1s                        -> 1 function
#   Li ... Ne   1s + 2s + 2p(3)           -> 5 functions
#   Na ... Ar   1s + 2s + 2p(3) + 3s + 3p(3) -> 9 functions
#
# The number of spatial orbitals equals that number of basis functions. The
# electron count is the sum of the neutral atomic numbers of the atoms in the
# formula: no charge input is taken, so an ion or a radical is out of scope and
# the payload says so.
#
# !! COUNTED BY RULE, NOT YET CONFIRMED AGAINST PySCF !!
# These counts are derived from the basis definition above and are not verified
# against an electronic-structure program on this host. The offline notebook
# (scripts/precompute_hamiltonians.ipynb) cross-checks ``mol.nao`` for H2 and
# LiH against this table, but that notebook needs PySCF and has not been run
# here. The payload therefore carries ``verified_against_pyscf: False`` (see
# STO3G_UNVERIFIED_NOTE) so a reader is never told the table was machine-checked
# when it was not.

STO3G_BASIS_FUNCTIONS_PER_ATOM: dict[str, int] = {
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

# Neutral atomic numbers, for the same elements the basis table covers.
ATOMIC_NUMBERS: dict[str, int] = {
    "H": 1,
    "He": 2,
    "Li": 3,
    "Be": 4,
    "B": 5,
    "C": 6,
    "N": 7,
    "O": 8,
    "F": 9,
    "Ne": 10,
    "Na": 11,
    "Mg": 12,
    "Al": 13,
    "Si": 14,
    "P": 15,
    "S": 16,
    "Cl": 17,
    "Ar": 18,
}

# Lightest first, so an error message and the usage response both read as a
# periodic-table run rather than an arbitrary set order.
SUPPORTED_ELEMENTS: tuple[str, ...] = tuple(
    sorted(STO3G_BASIS_FUNCTIONS_PER_ATOM, key=lambda symbol: ATOMIC_NUMBERS[symbol])
)

STO3G_UNVERIFIED_NOTE = (
    "Formula-derived from the STO-3G basis definition, not yet confirmed against "
    "PySCF mol.nao. Treat as unverified until the offline notebook cross-check runs."
)

_FORMULA_TOKEN = re.compile(r"([A-Z][a-z]?)(\d*)")


class ResourceEstimateError(ValueError):
    """A resource estimate that cannot be answered from the formula alone.

    Subclasses ``ValueError`` so callers that only ever saw a ``ValueError`` (the
    earlier build) keep working unchanged.
    """


class InvalidFormulaError(ResourceEstimateError):
    """The text is not a plain element-symbol-plus-count formula."""


class UnsupportedElementError(ResourceEstimateError):
    """The formula names an element with no STO-3G count in this table."""

    def __init__(self, elements: list[str] | tuple[str, ...]) -> None:
        self.elements = tuple(elements)
        super().__init__(
            "No STO-3G basis-function count for element(s): "
            + ", ".join(self.elements)
            + ". Supported elements are "
            + ", ".join(SUPPORTED_ELEMENTS)
            + "."
        )


class ActiveOrbitalsError(ResourceEstimateError):
    """``active_orbitals`` is not an integer in ``1..n_spatial_orbitals``."""


def parse_formula_counts(formula: str) -> dict[str, int]:
    """Atom counts from a plain molecular formula such as ``C9H8O4``.

    Handles the usual element-symbol-plus-optional-count form only. Raises
    ``InvalidFormulaError`` (a ``ValueError``) on anything else rather than
    guessing at charge, hydrates or parentheses.
    """
    text = formula.strip()
    if not text:
        raise InvalidFormulaError("molecular formula is empty")

    counts: dict[str, int] = {}
    position = 0
    for match in _FORMULA_TOKEN.finditer(text):
        if match.start() != position:
            raise InvalidFormulaError(f"could not parse molecular formula {formula!r}")
        symbol = match.group(1)
        digits = match.group(2)
        if digits and int(digits) == 0:
            raise InvalidFormulaError(
                f"{symbol}0 in formula {formula!r} is not a valid atom count"
            )
        counts[symbol] = counts.get(symbol, 0) + (int(digits) if digits else 1)
        position = match.end()

    if position != len(text):
        raise InvalidFormulaError(f"could not parse molecular formula {formula!r}")
    return counts


def sto3g_spatial_orbitals(formula: str) -> dict[str, Any]:
    """STO-3G spatial-orbital count for a formula.

    Counted by the minimal-basis rule, not computed by PySCF: see the module
    note above. The count is deliberately returned alongside
    ``verified_against_pyscf: False`` so a caller cannot quietly present it as a
    checked number.
    """
    counts = parse_formula_counts(formula)
    unsupported = sorted(
        (e for e in counts if e not in STO3G_BASIS_FUNCTIONS_PER_ATOM),
        key=lambda symbol: (len(symbol), symbol),
    )
    if unsupported:
        raise UnsupportedElementError(unsupported)

    per_element = {
        element: counts[element] * STO3G_BASIS_FUNCTIONS_PER_ATOM[element]
        for element in sorted(counts)
    }
    return {
        "formula": formula,
        "atom_counts": {element: counts[element] for element in sorted(counts)},
        "basis_functions_per_element": per_element,
        "spatial_orbitals": sum(per_element.values()),
        "verified_against_pyscf": False,
        "note": STO3G_UNVERIFIED_NOTE,
    }


def electron_count(formula: str) -> int:
    """Neutral-atom electron count for a formula (sum of atomic numbers).

    No charge is taken as input, so this is the neutral count only.
    """
    counts = parse_formula_counts(formula)
    unsupported = sorted(
        (e for e in counts if e not in ATOMIC_NUMBERS), key=lambda symbol: (len(symbol), symbol)
    )
    if unsupported:
        raise UnsupportedElementError(unsupported)
    return sum(count * ATOMIC_NUMBERS[element] for element, count in counts.items())


# Our app's interactive ceiling, in qubits after the parity two-qubit reduction
# (app/services/quantum_service.py refuses a dense diagonalisation above 8).
# This is a PROJECT choice, not a physics result.
INTERACTIVE_CAP_QUBITS = 8

# A rough ceiling for a dense complex128 statevector on a laptop, in qubits.
# 2**25 x 16 bytes = 512 MiB, the point where a single dense vector stops being
# comfortable; it is an order-of-magnitude judgement, and it is tagged
# "estimated" wherever it appears.
LAPTOP_STATEVECTOR_CEILING_QUBITS = 25

# The statevector entry size the memory formula uses: one complex128 amplitude.
STATEVECTOR_BYTES_PER_AMPLITUDE = 16

# Above this many bytes there is no useful unit left to quote, so the string
# says only the honest floor (the threshold is 1024**5 = 1.126e15 bytes, so
# "more than 10^15 bytes" is literally true of every value that reaches it).
MAX_QUOTED_BYTES = 1024**5

_VERDICTS = ("simulable_in_app", "precompute_only", "beyond_classical_simulation")


def format_byte_count(n_bytes: int) -> str:
    """Human-readable size: KB/MB/GB/TB, or an honest floor for huge values."""
    if n_bytes >= MAX_QUOTED_BYTES:
        return "more than 10^15 bytes"
    for unit, factor in (("TB", 1024**4), ("GB", 1024**3), ("MB", 1024**2), ("KB", 1024)):
        if n_bytes >= factor:
            return f"{n_bytes / factor:.1f} {unit}"
    return f"{n_bytes} bytes"


def _feasibility(n_qubits_reduced: int) -> dict[str, Any]:
    """Project-limit verdict for a qubit count after the parity reduction."""
    if n_qubits_reduced <= INTERACTIVE_CAP_QUBITS:
        verdict = _VERDICTS[0]
        reason = (
            f"{n_qubits_reduced} qubits after the parity two-qubit reduction is at or below "
            f"this app's {INTERACTIVE_CAP_QUBITS}-qubit interactive limit, so a dense statevector "
            "calculation of this size fits the app's budget."
        )
    elif n_qubits_reduced <= LAPTOP_STATEVECTOR_CEILING_QUBITS:
        verdict = _VERDICTS[1]
        reason = (
            f"{n_qubits_reduced} qubits after the parity two-qubit reduction is above this app's "
            f"{INTERACTIVE_CAP_QUBITS}-qubit interactive limit but inside a rough "
            f"{LAPTOP_STATEVECTOR_CEILING_QUBITS}-qubit laptop statevector ceiling, so the numbers "
            "would have to be precomputed offline and served as an artifact."
        )
    else:
        verdict = _VERDICTS[2]
        reason = (
            f"{n_qubits_reduced} qubits after the parity two-qubit reduction is above a rough "
            f"{LAPTOP_STATEVECTOR_CEILING_QUBITS}-qubit laptop statevector ceiling, so a dense "
            "statevector calculation is out of reach for this approach."
        )

    return {
        "based_on_qubits": n_qubits_reduced,
        "interactive_cap_qubits": INTERACTIVE_CAP_QUBITS,
        "laptop_statevector_ceiling_qubits": LAPTOP_STATEVECTOR_CEILING_QUBITS,
        "verdict": verdict,
        "reason": reason,
        "note": (
            "These limits are project choices and rough estimates, not measured hardware or "
            "algorithm bounds."
        ),
    }


def _qubit_scale(n_qubits: int) -> dict[str, Any]:
    """Memory and Hilbert-space figures for one qubit count (closed form)."""
    n_bytes = 2**n_qubits * STATEVECTOR_BYTES_PER_AMPLITUDE
    return {
        "statevector_memory_bytes": n_bytes,
        "statevector_memory_human": format_byte_count(n_bytes),
        "hilbert_dimension": str(2**n_qubits),
        "statevector_basis": "n_qubits_after_parity_reduction",
        "bytes_per_amplitude": STATEVECTOR_BYTES_PER_AMPLITUDE,
    }


def estimate_resources(formula: str, active_orbitals: int | None = None) -> dict[str, Any]:
    """What would this molecule need? Counted from its formula alone.

    Every figure comes from a closed-form rule in this module: the STO-3G
    spatial-orbital count (the table at the top of the file), the neutral
    electron count, one qubit per spin orbital (``2 * n_spatial_orbitals``), the
    fixed two-qubit saving of the parity mapping with two-qubit reduction, and
    the ``2**n_qubits * 16`` byte size of a complex128 statevector. Nothing is
    built, diagonalised or looked up in a database, and no Pauli-term or
    circuit-size number is produced because none of those can be had from a
    formula.

    ``active_orbitals`` (1 <= k <= n_spatial_orbitals) adds an ``active_space``
    block describing the smaller problem a caller would pick by hand; that block
    is tagged ``estimated`` because choosing an active space is an approximation,
    not a measurement.

    Raises ``InvalidFormulaError``, ``UnsupportedElementError`` or
    ``ActiveOrbitalsError`` (all ``ValueError`` subclasses) on bad input.
    """
    basis = sto3g_spatial_orbitals(formula)
    n_spatial = basis["spatial_orbitals"]
    n_electrons = electron_count(formula)
    n_qubits_full = 2 * n_spatial
    n_qubits_reduced = n_qubits_full - 2

    if active_orbitals is not None:
        if (
            isinstance(active_orbitals, bool)
            or not isinstance(active_orbitals, int)
            or not (1 <= active_orbitals <= n_spatial)
        ):
            raise ActiveOrbitalsError(
                f"'active_orbitals' must be an integer between 1 and {n_spatial} "
                f"(the STO-3G spatial-orbital count for {basis['formula']}), "
                f"got {active_orbitals!r}"
            )

    notes: dict[str, str] = {
        "spin_and_charge": (
            "Spin and charge are not modelled: the electron count is the sum of the neutral "
            "atomic numbers in the formula, and no charge input is taken."
        ),
        "method": (
            "Computed by formula for the STO-3G basis, not by running PySCF: no Hamiltonian was "
            "built and no diagonalisation was performed for this molecule."
        ),
        "pauli_term_scaling": (
            "Pauli term growth is roughly N^4 in the number of orbitals (estimated, "
            "order-of-magnitude). No specific term count is computed here, because a term count "
            "depends on the mapped Hamiltonian, which this estimator does not build."
        ),
        "limits": (
            "The classical simulation limits in 'feasibility' are project limits and rough "
            "estimates; they are not measured hardware or algorithm bounds."
        ),
        "interactive_limit": (
            f"Above {INTERACTIVE_CAP_QUBITS} qubits this app refuses a dense diagonalisation "
            "(too_many_qubits), and H2 is the only molecule with a VQE ansatz, so "
            "'simulable_in_app' means the size fits the app's budget rather than that a "
            "simulation exists for this molecule."
        ),
        "verification": (
            "The STO-3G orbital counts are counted from the minimal-basis rule, not read from a "
            "basis file, and have not been cross-checked against PySCF on this host."
        ),
    }
    if n_electrons % 2:
        notes["open_shell"] = (
            f"This formula has an odd electron count ({n_electrons}), so it is an open-shell "
            "system; the counts below still apply, but spin and charge are not modelled."
        )

    payload: dict[str, Any] = {
        "status": "computed",
        "subject": {
            "input": formula,
            "formula": basis["formula"],
            "atom_counts": basis["atom_counts"],
        },
        "basis": "STO-3G",
        "method": "computed_by_formula",
        "method_note": (
            "Computed by formula for the STO-3G basis (no PySCF run, no simulation). The "
            "classical simulation limits below are project limits and rough estimates."
        ),
        "n_spatial_orbitals": n_spatial,
        "basis_functions_per_element": basis["basis_functions_per_element"],
        "n_electrons": n_electrons,
        "n_qubits_full": n_qubits_full,
        "n_qubits_after_parity_reduction": n_qubits_reduced,
        "parity_reduction_removed_qubits": 2,
    }
    payload.update(_qubit_scale(n_qubits_reduced))
    payload["feasibility"] = _feasibility(n_qubits_reduced)
    payload["active_space"] = None

    tags: dict[str, str] = {
        "n_spatial_orbitals": "exact",
        "n_electrons": "exact",
        "n_qubits_full": "exact",
        "n_qubits_after_parity_reduction": "exact",
        "parity_reduction_removed_qubits": "exact",
        "statevector_memory_bytes": "exact",
        "statevector_memory_human": "exact",
        "hilbert_dimension": "exact",
        "feasibility.interactive_cap_qubits": "exact",
        "feasibility.laptop_statevector_ceiling_qubits": "estimated",
        "feasibility.verdict": "estimated",
        "notes.pauli_term_scaling": "estimated",
    }

    if active_orbitals is not None:
        k = active_orbitals
        active_qubits = 2 * k
        active_reduced = active_qubits - 2
        active_block: dict[str, Any] = {
            "active_orbitals": k,
            "n_spatial_orbitals": n_spatial,
            "n_qubits": active_qubits,
            "n_qubits_after_parity_reduction": active_reduced,
            "feasibility": _feasibility(active_reduced),
            "note": (
                "An active space is an approximation chosen by the user: every figure in this "
                "block describes the k-orbital problem, not the full molecule, and the frozen "
                "orbitals it leaves out change the energy."
            ),
        }
        active_block.update(_qubit_scale(active_reduced))
        payload["active_space"] = active_block
        tags.update(
            {
                "active_space.n_qubits": "estimated",
                "active_space.n_qubits_after_parity_reduction": "estimated",
                "active_space.statevector_memory_bytes": "estimated",
                "active_space.statevector_memory_human": "estimated",
                "active_space.hilbert_dimension": "estimated",
                "active_space.feasibility.verdict": "estimated",
            }
        )

    payload["tags"] = tags
    payload["notes"] = notes
    payload["verified_against_pyscf"] = False
    payload["classical_simulation"] = True
    return payload
