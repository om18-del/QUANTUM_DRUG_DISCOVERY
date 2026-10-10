"""Exact diagonalisation helpers for small Hamiltonians.

Two independent reference paths live here, and they are deliberately kept apart:

* ``diagonalize_ground_state`` takes an explicit square matrix.
* ``pauli_sum_ground_energy`` takes a Pauli-sum operator and assembles its matrix
  with NumPy only -- no Qiskit, no simulator.

The second path exists so the "exact" number a VQE run is graded against is
produced by completely different code from the one producing the VQE number.
If both sides shared a Qiskit operator constructor, a convention mistake would
cancel out and the comparison would prove nothing.

The artifact readers at the bottom take an already-parsed JSON object, never a
path, so nothing here touches the filesystem and the caller decides which file to
read. Nothing in this module is reachable from an API route
(``tests/test_quantum_core.py`` enforces that), and the readers are exercised only
on synthetic fixtures defined inside ``tests/test_h2_artifact_reader.py``.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from functools import reduce
from typing import Any, Sequence

import numpy as np
from numpy.typing import NDArray

# Single-qubit Pauli operators. ``I`` is included so a Pauli string can be built
# by a plain Kronecker product over its characters.
PAULI_MATRICES: dict[str, NDArray] = {
    "I": np.array([[1.0, 0.0], [0.0, 1.0]], dtype=np.complex128),
    "X": np.array([[0.0, 1.0], [1.0, 0.0]], dtype=np.complex128),
    "Y": np.array([[0.0, -1.0j], [1.0j, 0.0]], dtype=np.complex128),
    "Z": np.array([[1.0, 0.0], [0.0, -1.0]], dtype=np.complex128),
}

VALID_PAULI_CHARS = frozenset(PAULI_MATRICES)


def diagonalize_ground_state(hamiltonian: NDArray) -> float:
    """Return the lowest eigenvalue of a Hermitian Hamiltonian matrix.

    This is the exact ground energy for the given finite-basis Hamiltonian.
    """
    if hamiltonian.ndim != 2 or hamiltonian.shape[0] != hamiltonian.shape[1]:
        raise ValueError("Hamiltonian must be a square matrix")
    eigenvalues = np.linalg.eigvalsh(hamiltonian)
    return float(eigenvalues[0].real)


def pauli_string_matrix(pauli_string: str) -> NDArray:
    """Dense matrix for one Pauli string, e.g. ``"IX"`` or ``"ZZ"``.

    Qubit ordering matches Qiskit's ``SparsePauliOp.to_matrix()``: the string is
    read left to right and the **rightmost character acts on qubit 0**, so

        matrix = P[s[0]] (x) P[s[1]] (x) ... (x) P[s[-1]]

    which makes ``"IX"`` equal to ``kron(I, X)``. That agreement is asserted
    against Qiskit in the test suite rather than merely asserted here, because a
    silent mismatch would make every VQE-vs-exact comparison meaningless.
    """
    if not pauli_string:
        raise ValueError("Pauli string must be non-empty")

    matrices: list[NDArray] = []
    for index, char in enumerate(pauli_string):
        key = char.upper()
        if key not in VALID_PAULI_CHARS:
            raise ValueError(
                f"Invalid Pauli character {char!r} at position {index} of {pauli_string!r}; "
                f"expected one of {sorted(VALID_PAULI_CHARS)}"
            )
        matrices.append(PAULI_MATRICES[key])

    return reduce(np.kron, matrices)


def pauli_sum_matrix(
    terms: Mapping[str, float],
    n_qubits: int | None = None,
) -> NDArray:
    """Dense matrix of a Pauli-sum operator such as ``{"ZZ": 1.0, "IX": 1.0}``.

    ``terms`` maps a Pauli string to a real coefficient. All strings must have the
    same length; ``n_qubits`` is optional and is only used to reject a mismatched
    width early.
    """
    if not terms:
        raise ValueError("A Pauli sum needs at least one term")

    widths = {len(key) for key in terms}
    if len(widths) != 1:
        raise ValueError(
            f"All Pauli strings must have the same length, got widths {sorted(widths)}"
        )

    width = widths.pop()
    if n_qubits is not None and width != n_qubits:
        raise ValueError(
            f"Pauli strings are {width} characters wide but n_qubits={n_qubits}"
        )
    if width == 0:
        raise ValueError("Pauli strings must be non-empty")

    dimension = 2**width
    total = np.zeros((dimension, dimension), dtype=np.complex128)
    for pauli_string, coefficient in terms.items():
        if not np.isreal(coefficient):
            raise ValueError(
                f"Coefficient for {pauli_string!r} must be real, got {coefficient!r}"
            )
        total = total + float(np.real(coefficient)) * pauli_string_matrix(pauli_string)

    return total


def pauli_sum_ground_energy(
    terms: Mapping[str, float],
    n_qubits: int | None = None,
) -> float:
    """Exact lowest eigenvalue of a Pauli-sum operator."""
    return diagonalize_ground_state(pauli_sum_matrix(terms, n_qubits))


def h2_sto3g_hamiltonian_from_json(data: dict) -> NDArray:
    """Reconstruct the H2 STO-3G Hamiltonian matrix from a precomputed JSON blob.

    Expected keys:
        n_qubits: int
        matrix: list of lists representing the Hamiltonian in the computational basis
    """
    n_qubits = int(data["n_qubits"])
    matrix = np.array(data["matrix"], dtype=np.complex128)
    if matrix.shape != (2**n_qubits, 2**n_qubits):
        raise ValueError("Stored matrix shape does not match n_qubits")
    return matrix


def pauli_terms_from_sequence(
    terms: Sequence[tuple[str, float]],
) -> dict[str, float]:
    """Convert ``[("ZZ", 1.0), ("IX", 1.0)]`` into the mapping form used above.

    Duplicate strings are summed rather than silently overwriting each other.
    """
    merged: dict[str, float] = {}
    for pauli_string, coefficient in terms:
        merged[pauli_string] = merged.get(pauli_string, 0.0) + float(coefficient)
    return merged


# ---------------------------------------------------------------------------
# Offline artifact readers: the shape the writer actually produces
# ---------------------------------------------------------------------------
# ``scripts/precompute_hamiltonians.py`` stores one entry per bond length as
#
#     {"bond_length_angstrom": 0.74, "nuclear_repulsion_ha": 0.7075,
#      "pauli_terms": [{"pauli": "II", "coefficient": -1.05}, ...], ...}
#
# which is NOT the shape ``h2_sto3g_hamiltonian_from_json`` above reads
# (``n_qubits`` plus a dense ``matrix``). Before these readers existed the
# application could not read the artifact it was built to consume: the two shapes
# never met. The readers below speak the artifact's shape instead, and they keep
# ``nuclear_repulsion_ha`` as its own field because the stored Pauli terms are the
# *electronic* operator only -- the chemistry number is their ground state PLUS
# the nuclear repulsion. Every energy here is in Hartree, as ``units`` records.
#
# Nothing below interpolates. A bond length is either stored or it is not, so a
# request that is off the grid raises and names the nearest stored value; deciding
# what to do about that (SPEC.md section 6: snap to the nearest stored point and
# disclose it in the response, never invent an energy) belongs to the caller.

ARTIFACT_UNITS = "hartree"

# Bond lengths are floats parsed from JSON, so an exact-match lookup uses a
# tolerance far below any physically meaningful spacing rather than ``==``.
BOND_LENGTH_MATCH_TOLERANCE_ANGSTROM = 1e-9


@dataclass(frozen=True)
class ArtifactPointHamiltonian:
    """One stored bond length of the offline artifact, ready to use.

    Attributes:
        bond_length_angstrom: the stored bond length, in angstrom.
        n_qubits: width of the stored Pauli strings (2 for the H2 artifact).
        pauli_terms: ``{pauli_string: coefficient}``, the electronic operator only.
        nuclear_repulsion_ha: the stored nuclear repulsion, kept as its own field.
        matrix: dense matrix of ``pauli_terms`` (electronic part only).
        ground_energy_ha: lowest eigenvalue of ``matrix`` -- electronic only.
        units: ``"hartree"``, as the artifact records.

    ``total_energy_ha`` is the number to compare with a chemistry reference: the
    electronic ground state plus the nuclear repulsion.
    """

    bond_length_angstrom: float
    n_qubits: int
    pauli_terms: dict[str, float]
    nuclear_repulsion_ha: float
    matrix: NDArray
    ground_energy_ha: float
    units: str = ARTIFACT_UNITS

    @property
    def total_energy_ha(self) -> float:
        """Electronic ground state plus nuclear repulsion, in Hartree."""
        return self.ground_energy_ha + self.nuclear_repulsion_ha


def _artifact_points(data: Mapping[str, Any]) -> list[Mapping[str, Any]]:
    """The ``points`` list, with each entry checked enough to be looked up."""
    if not isinstance(data, Mapping):
        raise ValueError(f"An artifact must be a JSON object, got {type(data).__name__}")
    points = data.get("points")
    if not isinstance(points, Sequence) or isinstance(points, (str, bytes)) or not points:
        raise ValueError("An artifact needs a non-empty 'points' list")
    for index, point in enumerate(points):
        if not isinstance(point, Mapping) or "bond_length_angstrom" not in point:
            raise ValueError(f"points[{index}] must be an object with 'bond_length_angstrom'")
    return list(points)


def h2_artifact_bond_lengths(data: Mapping[str, Any]) -> list[float]:
    """Every stored bond length, in the artifact's own order."""
    return [float(point["bond_length_angstrom"]) for point in _artifact_points(data)]


def h2_artifact_pauli_terms(point: Mapping[str, Any]) -> dict[str, float]:
    """``{pauli_string: coefficient}`` for one stored point's ``pauli_terms`` list.

    Validates the writer's shape rather than trusting it: every entry must be an
    object with a real, finite ``coefficient`` and an ``IXYZ`` ``pauli`` string.
    Repeated Pauli strings are summed, so a duplicate cannot silently overwrite a
    coefficient.
    """
    if not isinstance(point, Mapping):
        raise ValueError(f"A stored point must be an object, got {type(point).__name__}")

    terms = point.get("pauli_terms")
    if not isinstance(terms, Sequence) or isinstance(terms, (str, bytes)) or not terms:
        raise ValueError(
            "A stored point needs a non-empty 'pauli_terms' list; "
            f"found {type(terms).__name__}"
        )

    merged: dict[str, float] = {}
    for index, term in enumerate(terms):
        if not isinstance(term, Mapping) or "pauli" not in term or "coefficient" not in term:
            raise ValueError(
                f"pauli_terms[{index}] must be an object with 'pauli' and 'coefficient', "
                f"got {term!r}"
            )

        pauli_string = term["pauli"]
        if (
            not isinstance(pauli_string, str)
            or not pauli_string
            or set(pauli_string.upper()) - VALID_PAULI_CHARS
        ):
            raise ValueError(f"pauli_terms[{index}] has an invalid Pauli string {pauli_string!r}")

        coefficient = term["coefficient"]
        if isinstance(coefficient, bool) or not isinstance(coefficient, (int, float)):
            raise ValueError(
                f"pauli_terms[{index}] coefficient must be a real number, got {coefficient!r}"
            )
        if not np.isfinite(float(coefficient)):
            raise ValueError(f"pauli_terms[{index}] coefficient must be finite, got {coefficient!r}")

        key = pauli_string.upper()
        merged[key] = merged.get(key, 0.0) + float(coefficient)

    return merged


def h2_artifact_point(
    data: Mapping[str, Any],
    bond_length_angstrom: float,
) -> Mapping[str, Any]:
    """The stored point for exactly this bond length. Never interpolates.

    Raises ``ValueError`` naming the nearest stored value when the grid does not
    contain the request. Applying SPEC.md section 6's "snap to the nearest stored
    point and disclose it" rule is the caller's job, not this function's.
    """
    requested = float(bond_length_angstrom)
    points = _artifact_points(data)

    for point in points:
        stored = float(point["bond_length_angstrom"])
        if abs(stored - requested) <= BOND_LENGTH_MATCH_TOLERANCE_ANGSTROM:
            return point

    nearest = min(points, key=lambda p: abs(float(p["bond_length_angstrom"]) - requested))
    raise ValueError(
        f"{requested:g} A is not on the artifact's grid; the nearest stored point is "
        f"{float(nearest['bond_length_angstrom']):g} A. This reader never interpolates."
    )


def nearest_stored_bond_length(data: Mapping[str, Any], bond_length_angstrom: float) -> float:
    """The stored bond length closest to the request. Lookup only, no interpolation.

    SPEC.md section 6 documents the simulate endpoint's handling of an off-grid
    request as "snap to the nearest stored point and disclose it", which needs a
    deterministic tie-break; ``min`` returns the first of equal distances, i.e. the
    shorter bond length. The endpoint itself is not wired to this module yet.
    """
    requested = float(bond_length_angstrom)
    points = _artifact_points(data)
    nearest = min(points, key=lambda p: abs(float(p["bond_length_angstrom"]) - requested))
    return float(nearest["bond_length_angstrom"])


def h2_hamiltonian_from_artifact(
    data: Mapping[str, Any],
    bond_length_angstrom: float,
) -> ArtifactPointHamiltonian:
    """Read one stored bond length of the offline artifact into a Hamiltonian.

    Builds the electronic operator from the stored ``pauli_terms`` (2-qubit Paulis
    with real coefficients for the H2 artifact), keeps ``nuclear_repulsion_ha`` as
    a separate field, and reports every energy in Hartree. The ground energy is
    the lowest eigenvalue of the electronic operator alone; add the nuclear
    repulsion (``total_energy_ha``) to compare with a chemistry reference.

    Raises ``ValueError`` for an off-grid bond length or a payload the writer
    could not have produced. No interpolation happens anywhere in this path.
    """
    point = h2_artifact_point(data, bond_length_angstrom)
    terms = h2_artifact_pauli_terms(point)

    widths = {len(key) for key in terms}
    if len(widths) != 1:
        raise ValueError(f"All Pauli strings must have the same length, got {sorted(widths)}")
    width = widths.pop()

    stored_n_qubits = point.get("n_qubits", data.get("n_qubits"))
    if stored_n_qubits is not None and int(stored_n_qubits) != width:
        raise ValueError(
            f"stored n_qubits={stored_n_qubits} but the Pauli terms are {width} characters wide"
        )

    nuclear = point.get("nuclear_repulsion_ha")
    if (
        isinstance(nuclear, bool)
        or not isinstance(nuclear, (int, float))
        or not np.isfinite(float(nuclear))
    ):
        raise ValueError(
            f"A stored point needs a finite 'nuclear_repulsion_ha', got {nuclear!r}"
        )

    matrix = pauli_sum_matrix(terms, width)
    return ArtifactPointHamiltonian(
        bond_length_angstrom=float(point["bond_length_angstrom"]),
        n_qubits=width,
        pauli_terms=terms,
        nuclear_repulsion_ha=float(nuclear),
        matrix=matrix,
        ground_energy_ha=diagonalize_ground_state(matrix),
    )
