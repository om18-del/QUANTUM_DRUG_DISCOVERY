#!/usr/bin/env python
"""Standalone validator for the offline H2 Hamiltonian artifact.

Usage:
    python scripts/validate_hamiltonian_file.py [path]

Defaults to ``app/data/hamiltonians/h2_sto3g.json``. Exits 0 only when every
check passes, and non-zero with a clear message otherwise. It never writes,
repairs or fabricates the file: a missing artifact is reported as missing.

Standard library plus numpy only, and it deliberately does not import ``app``.
This is an independent check of the file the application will read, so it
recomputes the operator's ground energy from the stored Pauli terms itself
instead of trusting any stored energy or reusing the writer's code path.

Two notes on the checks, because they are not all the same kind of claim:

* The literal request "the stored HF and FCI energies agree within 1e-6 Ha" is
  deliberately NOT implemented. For H2 in STO-3G the two differ by the
  correlation energy -- order 0.02 Ha near equilibrium, far above 1e-6 -- so such
  a check would reject a physically correct file, and it contradicts the
  variational ordering checked below. What is checked instead is the artifact's
  own documented promise (``convention`` in the payload, and the writer's
  post-check): the ground state of the stored qubit operator plus the stored
  nuclear repulsion reproduces the stored FCI energy to 1e-6 Ha. Where the file
  records two *sources* for the same quantity, that is the comparison worth
  making to 1e-6.
* The FCI/HF pair is checked for the variational relation instead: FCI must sit
  at or below HF at every bond length.

The required field lists below mirror the contract in
``tests/test_h2_hamiltonian_artifact.py`` and the payload written by
``scripts/precompute_hamiltonians.py``. ``app/core`` now reads that payload
(``h2_hamiltonian_from_artifact`` and friends in ``app/core/exact.py``), while the
older ``h2_sto3g_hamiltonian_from_json`` there still reads a different shape
(``n_qubits`` plus a dense ``matrix``). Routers and services may not import
``app.core`` (a test enforces that), so ``app/services/quantum_service.py`` reads
the payload itself and serves both the curve and the single point from it; only
the resource estimate still reports the artifact as "present but not yet wired
to its format".
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
HAMILTONIAN_DIR = REPO_ROOT / "app" / "data" / "hamiltonians"
DEFAULT_ARTIFACT = HAMILTONIAN_DIR / "h2_sto3g.json"

# The registry is mirrored here so a molecule key resolves to its file without
# importing the application (this validator stays standalone and independent).
MOLECULE_FILES = {"h2": "h2_sto3g.json", "lih": "lih_sto3g_active.json"}

# The H2 grid, used as the fallback only when a file does not declare its own grid
# (the committed H2 artifact predates the ``bond_length_grid`` field). New artifacts
# declare their grid, so nothing here assumes 58 points or a 0.05 step for them.
UNIFORM_POINTS = 57
EXPECTED_POINTS = 58
UNIFORM_BOND_LENGTHS = [round(0.20 + 0.05 * i, 2) for i in range(UNIFORM_POINTS)]
EXPECTED_BOND_LENGTHS = sorted({*UNIFORM_BOND_LENGTHS, 0.74})
assert len(EXPECTED_BOND_LENGTHS) == EXPECTED_POINTS
TOLERANCE_HA = 1e-6

# SPEC.md section 6: H2 in STO-3G near 0.74 A is about -1.137 Ha. The actual
# value is always printed; this tolerance only decides pass or fail.
EQUILIBRIUM_ANGSTROM = 0.74
SPEC_FCI_NEAR_EQUILIBRIUM_HA = -1.137
SPEC_FCI_TOLERANCE_HA = 0.005

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

PAULI_MATRICES = {
    "I": np.eye(2, dtype=complex),
    "X": np.array([[0, 1], [1, 0]], dtype=complex),
    "Y": np.array([[0, -1j], [1j, 0]], dtype=complex),
    "Z": np.array([[1, 0], [0, -1]], dtype=complex),
}


class Failure(Exception):
    """A single failed check, reported and accumulated rather than raised out."""


def ground_energy_from_terms(terms: list[dict[str, Any]]) -> float:
    """Lowest eigenvalue of a Pauli sum, built with numpy only.

    Independent of Qiskit and of the writer: the matrix is assembled here from
    the stored (pauli, coefficient) pairs and symmetrised before diagonalising,
    so a Hermiticity mistake shows up as a real difference rather than a
    complex eigenvalue.
    """
    width = len(terms[0]["pauli"])
    matrix = np.zeros((2**width, 2**width), dtype=complex)
    for term in terms:
        operator = np.array([[1.0 + 0.0j]])
        for char in term["pauli"]:
            try:
                operator = np.kron(operator, PAULI_MATRICES[char])
            except KeyError:
                raise Failure(f"invalid Pauli character {char!r} in term {term['pauli']!r}") from None
        matrix = matrix + complex(term["coefficient"]) * operator
    hermitian = (matrix + matrix.conj().T) / 2.0
    return float(np.linalg.eigvalsh(hermitian)[0].real)


def _number(value: Any) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise Failure(f"expected a number, got {type(value).__name__}: {value!r}")
    if not math.isfinite(float(value)):
        raise Failure(f"expected a finite number, got {value!r}")
    return float(value)


def expected_grid(data: dict) -> tuple[list[float], str]:
    """The grid to check against, read from the file's own metadata when declared.

    New artifacts declare ``bond_length_grid``, so no molecule's point count or
    step is hard-coded here. A file without it -- the committed H2 artifact
    predates the field -- falls back to the H2 grid, and only when it really is
    H2; another molecule must declare its own grid rather than be measured
    against H2's.
    """
    declared = data.get("bond_length_grid")
    if declared is not None:
        if not isinstance(declared, list) or not declared:
            raise Failure("bond_length_grid must be a non-empty list when it is declared")
        return [_number(value) for value in declared], "declared"

    molecule = str(data.get("molecule", "")).strip().upper()
    if molecule in ("", "H2"):
        return list(EXPECTED_BOND_LENGTHS), "h2-default"
    raise Failure(
        f"{data.get('molecule')!r} declares no bond_length_grid; a non-H2 artifact must "
        "declare its own grid rather than be checked against the H2 one"
    )


def check_top_level(data: Any) -> list[dict[str, Any]]:
    if not isinstance(data, dict):
        raise Failure(f"top level must be a JSON object, got {type(data).__name__}")

    missing = [field for field in TOP_LEVEL_FIELDS if field not in data]
    if missing:
        raise Failure("missing top-level field(s): " + ", ".join(sorted(missing)))

    if not isinstance(data["library_versions"], dict):
        raise Failure(f"library_versions must be an object, got {type(data['library_versions']).__name__}")
    if not isinstance(data["points"], list):
        raise Failure(f"points must be a list, got {type(data['points']).__name__}")
    return data["points"]


def check_grid(points: list[Any], expected: list[float], source: str) -> None:
    """The grid the file declares (or the H2 fallback), ascending, with no duplicates."""
    if len(points) != len(expected):
        raise Failure(f"expected exactly {len(expected)} bond lengths, found {len(points)}")

    lengths: list[float] = []
    for index, point in enumerate(points):
        if not isinstance(point, dict):
            raise Failure(f"points[{index}] must be an object, got {type(point).__name__}")
        if "bond_length_angstrom" not in point:
            raise Failure(f"points[{index}] is missing bond_length_angstrom")
        lengths.append(_number(point["bond_length_angstrom"]))

    duplicates = sorted({value for value in lengths if lengths.count(value) > 1})
    if duplicates:
        raise Failure("duplicate bond length(s): " + ", ".join(f"{value:g}" for value in duplicates))

    if lengths != sorted(lengths):
        raise Failure("bond lengths are not in ascending order")

    mismatches = [
        f"index {index}: expected {want:g}, found {found:g}"
        for index, (want, found) in enumerate(zip(expected, lengths))
        if want != found
    ]
    if mismatches:
        if source == "h2-default":
            raise Failure(
                "bond-length grid is not 0.20 to 3.00 A in 0.05 steps plus 0.74 A -- "
                + "; ".join(mismatches[:5])
            )
        raise Failure(
            "bond lengths do not match the file's declared bond_length_grid -- "
            + "; ".join(mismatches[:5])
        )

    first, last = lengths[0], lengths[-1]
    if (first, last) != (expected[0], expected[-1]):
        raise Failure(
            f"grid must start at {expected[0]:g} and end at {expected[-1]:g} A, got {first:g} .. {last:g}"
        )


def check_points(
    points: list[dict[str, Any]], n_qubits_expected: int
) -> dict[str, float]:
    """Field presence plus the physics claims. Returns measured summaries.

    The qubit count is the file's own top-level ``n_qubits``; no width is assumed.
    """
    problems: list[str] = []
    worst_qubit_deviation = 0.0
    # Tracked as a minimum, so it must start at +inf: starting at 0.0 made this
    # summary silently report 0.000000 whatever the file contained.
    smallest_hf_minus_fci = float("inf")

    for point in points:
        length = point["bond_length_angstrom"]

        missing = [field for field in POINT_FIELDS if field not in point]
        if missing:
            problems.append(f"r={length:g}: missing field(s) {', '.join(sorted(missing))}")
            continue

        terms = point["pauli_terms"]
        if not isinstance(terms, list) or not terms:
            problems.append(f"r={length:g}: pauli_terms must be a non-empty list")
            continue

        if point["n_qubits"] != n_qubits_expected:
            problems.append(
                f"r={length:g}: n_qubits must be {n_qubits_expected}, got {point['n_qubits']!r}"
            )
            continue

        widths = {len(str(term.get("pauli", ""))) for term in terms if isinstance(term, dict)}
        if len(widths) != 1:
            problems.append(
                f"r={length:g}: all Pauli terms must have the same length, got {sorted(widths)}"
            )
            continue
        width = widths.pop()
        if width != n_qubits_expected:
            problems.append(
                f"r={length:g}: Pauli terms are {width} characters wide but n_qubits="
                f"{n_qubits_expected}"
            )
            continue

        bad_term = next(
            (
                term
                for term in terms
                if not isinstance(term, dict)
                or "pauli" not in term
                or "coefficient" not in term
                or set(str(term.get("pauli", ""))) - set("IXYZ")
                or isinstance(term.get("coefficient"), bool)
                or not isinstance(term.get("coefficient"), (int, float))
            ),
            None,
        )
        if bad_term is not None:
            problems.append(f"r={length:g}: malformed Pauli term {bad_term!r}")
            continue

        nuclear = _number(point["nuclear_repulsion_ha"])
        hf = _number(point["hf_energy_ha"])
        fci = _number(point["fci_energy_ha"])

        # The artifact's documented promise, recomputed rather than trusted.
        try:
            computed = ground_energy_from_terms(terms) + nuclear
        except Failure as exc:
            problems.append(f"r={length:g}: {exc}")
            continue
        deviation = abs(computed - fci)
        worst_qubit_deviation = max(worst_qubit_deviation, deviation)

        if deviation > TOLERANCE_HA:
            problems.append(
                f"r={length:g}: qubit operator + nuclear repulsion = {computed:.12f} but stored "
                f"FCI = {fci:.12f} (deviation {deviation:.3e} Ha > {TOLERANCE_HA:g})"
            )

        recorded = _number(point["qubit_energy_plus_nuclear_ha"])
        if abs(recorded - computed) > TOLERANCE_HA:
            problems.append(
                f"r={length:g}: stored qubit_energy_plus_nuclear_ha {recorded:.12f} disagrees with "
                f"the recomputed {computed:.12f} by {abs(recorded - computed):.3e} Ha"
            )

        recorded_deviation = _number(point["deviation_from_fci_ha"])
        if abs(recorded_deviation - deviation) > TOLERANCE_HA:
            problems.append(
                f"r={length:g}: stored deviation_from_fci_ha {recorded_deviation:.3e} does not match "
                f"the recomputed {deviation:.3e}"
            )

        # NOT the check that "the stored HF and FCI energies agree within 1e-6 Ha".
        # That was never the right check: FCI lies *below* HF by the correlation
        # energy (order 0.02 Ha for H2 in STO-3G near equilibrium, not 1e-6), so it
        # would fail a physically correct file. The 1e-6 comparison belongs to the
        # two independent routes to one number -- the recomputed ground state of the
        # stored Pauli terms plus nuclear repulsion, versus the stored PySCF FCI
        # (checked above). What the HF/FCI pair is checked for here is the
        # variational ordering: a correlation method cannot sit above Hartree-Fock.
        hf_minus_fci = hf - fci
        smallest_hf_minus_fci = min(smallest_hf_minus_fci, hf_minus_fci)
        if hf_minus_fci < -1e-9:
            problems.append(
                f"r={length:g}: FCI ({fci:.12f}) is above HF ({hf:.12f}) -- a variational method "
                "cannot be higher than the mean-field energy"
            )

    if problems:
        raise Failure(f"{len(problems)} point-level problem(s):\n    " + "\n    ".join(problems[:12]))

    return {
        "worst_qubit_deviation": worst_qubit_deviation,
        "smallest_hf_minus_fci": smallest_hf_minus_fci,
    }


def check_equilibrium_energy(points: list[dict[str, Any]]) -> None:
    """The printed -1.137 Ha reference, checked but never hardcoded as a pass.

    0.74 A is on the grid by construction (see EXPECTED_BOND_LENGTHS above) and
    ``check_grid`` has already required an exact match to that grid, so the point
    nearest 0.74 A is 0.74 A itself. That is asserted rather than tolerated: the
    equilibrium entry is the one a reader is most likely to quote.
    """
    nearest = min(points, key=lambda point: abs(float(point["bond_length_angstrom"]) - EQUILIBRIUM_ANGSTROM))
    length = float(nearest["bond_length_angstrom"])
    actual = float(nearest["fci_energy_ha"])
    delta = abs(actual - SPEC_FCI_NEAR_EQUILIBRIUM_HA)

    print(
        f"  nearest point to {EQUILIBRIUM_ANGSTROM:g} A  : r={length:g} A  "
        f"FCI={actual:.6f} Ha  (SPEC.md reference {SPEC_FCI_NEAR_EQUILIBRIUM_HA:g} Ha, "
        f"difference {delta:.6f} Ha)"
    )
    if length != EQUILIBRIUM_ANGSTROM:
        raise Failure(
            f"the grid point nearest {EQUILIBRIUM_ANGSTROM:g} A is {length:g} A, but the grid "
            f"must contain {EQUILIBRIUM_ANGSTROM:g} A exactly: it is the UI default and H2's real "
            "equilibrium bond length, and the 0.05 step alone never produces it"
        )
    if delta > SPEC_FCI_TOLERANCE_HA:
        raise Failure(
            f"FCI near equilibrium is {actual:.6f} Ha, which is {delta:.6f} Ha away from the "
            f"SPEC.md reference of {SPEC_FCI_NEAR_EQUILIBRIUM_HA:g} Ha (tolerance "
            f"{SPEC_FCI_TOLERANCE_HA:g} Ha). The measured value above is what was stored."
        )


def check_recorded_claims(data: dict[str, Any], measured: dict[str, float]) -> None:
    """The payload's own summary fields, so they cannot drift from the points."""
    if data["point_count"] != len(data["points"]):
        raise Failure(f"point_count is {data['point_count']} but points holds {len(data['points'])}")

    n_qubits = data["n_qubits"]
    if isinstance(n_qubits, bool) or not isinstance(n_qubits, int) or n_qubits < 1:
        raise Failure(f"top-level n_qubits must be a positive integer, got {n_qubits!r}")
    point_widths = {
        len(str(term["pauli"]))
        for point in data["points"]
        for term in point.get("pauli_terms", [])
        if isinstance(term, dict) and "pauli" in term
    }
    if point_widths != {n_qubits}:
        raise Failure(
            f"top-level n_qubits is {n_qubits} but the stored Pauli terms are widths "
            f"{sorted(point_widths)}"
        )

    if str(data.get("molecule", "")).strip().upper() != "H2":
        active_space = data.get("active_space")
        if not isinstance(active_space, str) or not active_space.strip():
            raise Failure(
                f"top-level field 'active_space' is required for {data.get('molecule')!r}"
            )

    if str(data["mapping"]).lower() != "parity":
        raise Failure(f"expected the parity mapping, found {data['mapping']!r}")
    if data["two_qubit_reduction"] is not True:
        raise Failure(f"two_qubit_reduction must be true, found {data['two_qubit_reduction']!r}")

    recorded_max = data.get("max_deviation_from_fci_ha")
    if recorded_max is not None:
        recorded_max = _number(recorded_max)
        if recorded_max > TOLERANCE_HA:
            raise Failure(
                f"the file's own max_deviation_from_fci_ha {recorded_max:.3e} exceeds {TOLERANCE_HA:g}"
            )
        if abs(recorded_max - measured["worst_qubit_deviation"]) > TOLERANCE_HA:
            raise Failure(
                f"max_deviation_from_fci_ha is {recorded_max:.3e} but the worst recomputed "
                f"deviation is {measured['worst_qubit_deviation']:.3e}"
            )

    for point in data["points"]:
        if abs(float(point["deviation_from_fci_ha"])) > TOLERANCE_HA:
            raise Failure(
                f"r={point['bond_length_angstrom']:g}: recorded deviation_from_fci_ha "
                f"{point['deviation_from_fci_ha']:.3e} exceeds {TOLERANCE_HA:g}"
            )


def validate(path: Path) -> int:
    """Run every check. Returns a process exit code."""
    print(f"validating: {path}")

    if not path.is_file():
        print(
            f"FAIL: {path} does not exist.\n"
            "  Nothing to validate. Generate it offline with "
            "scripts/precompute_hamiltonians.ipynb (Colab or WSL, because PySCF has no\n"
            f"  Windows wheel) and save the result to app/data/hamiltonians/{path.name}.\n"
            "  This validator will not create or invent the file.",
            file=sys.stderr,
        )
        return 1

    try:
        raw = path.read_text(encoding="utf-8")
    except OSError as exc:
        print(f"FAIL: could not read {path}: {exc}", file=sys.stderr)
        return 1

    try:
        data = json.loads(raw)
    except json.JSONDecodeError as exc:
        print(f"FAIL: {path} is not valid JSON: {exc}", file=sys.stderr)
        return 1

    try:
        points = check_top_level(data)
        expected, source = expected_grid(data)
        check_grid(points, expected, source)
        measured = check_points(points, data["n_qubits"])
        check_recorded_claims(data, measured)
        # The 0.74 A / -1.137 Ha equilibrium reference is an H2 fact; another
        # molecule's grid has no reason to contain it.
        if str(data.get("molecule", "")).strip().upper() == "H2":
            check_equilibrium_energy(points)
    except Failure as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        return 1

    print(f"  points                 : {len(points)} ({expected[0]:g} .. {expected[-1]:g} A)")
    print(f"  worst qubit-vs-FCI gap : {measured['worst_qubit_deviation']:.3e} Ha")
    print(f"  smallest (HF - FCI)    : {measured['smallest_hf_minus_fci']:.6f} Ha  (must be >= 0)")
    print("OK: every check passed.")
    return 0


def resolve_target(target: str | None) -> Path:
    """A molecule key (h2, lih) or a file path, as the CLI accepts either."""
    if target is None:
        return DEFAULT_ARTIFACT
    key = target.strip().lower()
    if key in MOLECULE_FILES:
        return HAMILTONIAN_DIR / MOLECULE_FILES[key]
    return Path(target).resolve()


def main() -> int:
    parser = argparse.ArgumentParser(
        prog=Path(sys.argv[0]).name,
        description="Validate an offline Hamiltonian artifact, by molecule key or file path.",
    )
    parser.add_argument("target", nargs="?", default=None, help="molecule key (h2, lih) or a JSON path")
    parser.add_argument(
        "--molecule",
        choices=sorted(MOLECULE_FILES),
        default=None,
        help="molecule key; mutually exclusive with a path target",
    )
    args = parser.parse_args()

    if args.molecule is not None and args.target is not None:
        print("usage: pass a molecule key or a path, not both", file=sys.stderr)
        return 2
    if args.molecule is not None:
        return validate(HAMILTONIAN_DIR / MOLECULE_FILES[args.molecule])
    return validate(resolve_target(args.target))


if __name__ == "__main__":
    raise SystemExit(main())
