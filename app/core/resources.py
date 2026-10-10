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
# STO-3G spatial-orbital counts from a molecular formula
# ---------------------------------------------------------------------------
#
# Not wired into the API or the UI yet: this is groundwork for the Phase 3
# resource estimator (SPEC.md section 11, P3-1).
#
# STO-3G is a minimal basis with one contracted function per occupied atomic
# shell, so H and He contribute a single 1s function and every first-row atom
# from Li to Ne contributes 1s + 2s + 2px + 2py + 2pz = five. The number of
# spatial orbitals equals the number of such basis functions.
#
# !! FORMULA-DERIVED, NOT YET CONFIRMED AGAINST PySCF !!
# These counts are hand-derived from the basis definition and are not verified
# against an electronic-structure program. Before anything reports them as
# exact, scripts/precompute_hamiltonians.ipynb must build each molecule and
# check ``mol.nao`` against ``spatial_orbitals`` below. Until then the result
# must be presented as unverified, never as an exact count.

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
}

STO3G_UNVERIFIED_NOTE = (
    "Formula-derived from the STO-3G basis definition, not yet confirmed against "
    "PySCF mol.nao. Treat as unverified until the offline notebook cross-check runs."
)

_FORMULA_TOKEN = re.compile(r"([A-Z][a-z]?)(\d*)")


def parse_formula_counts(formula: str) -> dict[str, int]:
    """Atom counts from a plain molecular formula such as ``C9H8O4``.

    Handles the usual element-symbol-plus-optional-count form only. Raises
    ``ValueError`` on anything else rather than guessing at charge, hydrates or
    parentheses.
    """
    text = formula.strip()
    if not text:
        raise ValueError("molecular formula is empty")

    counts: dict[str, int] = {}
    position = 0
    for match in _FORMULA_TOKEN.finditer(text):
        if match.start() != position:
            raise ValueError(f"could not parse molecular formula {formula!r}")
        symbol = match.group(1)
        digits = match.group(2)
        if digits and int(digits) == 0:
            raise ValueError(f"{symbol}0 in formula {formula!r} is not a valid atom count")
        counts[symbol] = counts.get(symbol, 0) + (int(digits) if digits else 1)
        position = match.end()

    if position != len(text):
        raise ValueError(f"could not parse molecular formula {formula!r}")
    return counts


def sto3g_spatial_orbitals(formula: str) -> dict[str, Any]:
    """STO-3G spatial-orbital count for a formula.

    Unverified: see the module note above. The count is deliberately returned
    alongside ``verified_against_pyscf: False`` so a caller cannot quietly
    present it as an exact, checked number.
    """
    counts = parse_formula_counts(formula)
    unsupported = sorted(e for e in counts if e not in STO3G_BASIS_FUNCTIONS_PER_ATOM)
    if unsupported:
        raise ValueError(
            "no STO-3G basis-function count for element(s): " + ", ".join(unsupported)
        )

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
