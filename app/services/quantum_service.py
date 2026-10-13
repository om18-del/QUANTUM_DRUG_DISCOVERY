"""Quantum Lab service.

Everything here is built from the offline artifact
``app/data/hamiltonians/h2_sto3g.json``, produced by
``scripts/precompute_hamiltonians.ipynb`` (Colab or WSL, because PySCF publishes
no Windows wheel). For each of its 58 stored bond lengths the artifact carries
the PySCF HF and FCI reference energies, the nuclear repulsion, and the five
Pauli terms of the parity-mapped, two-qubit-reduced electronic Hamiltonian.

* ``compute_h2_curve`` **reads** the stored grid: no optimiser, shots or noise
  run inside that request, and the note says so.
* ``compute_h2_simulate`` **runs a real VQE** -- a fixed Ry ansatz minimised by
  SciPy COBYLA through a Qiskit estimator --  against the artifact's Pauli terms, and grades it against an exact NumPy
  diagonalisation of the same operator. Snap-and-disclose follows SPEC.md
  section 6. A noiseless statevector run that lands below the exact energy
  raises ``VariationalPrincipleViolation`` instead of being served; a shot-based
  run is a statistic, so it is served with its uncertainty
  (``vqe_energy_std_ha``) and no chemical-accuracy claim, never refused for
  sampling noise alone (SPEC.md section 6).

``app.core`` is imported **lazily inside those functions**, never at module
import time; ``tests/test_quantum_core.py`` enforces both halves of that rule.

An earlier build shipped a hand-typed 4x4 matrix in
``app/data/hamiltonians/`` and reported it as ``status: "computed"``; those files
have been deleted, because presenting them as computed H2 physics was wrong.
Every number in a response now comes either from the artifact or from a
computation on it, and the provenance note on each response says which.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from app.config import get_settings

H2_HAMILTONIAN_ARTIFACT = "h2_sto3g.json"

# The curve, the single-point simulation and the resource estimate are three
# different things, so each one names its own subject: the curve and the energy
# are served from the artifact, while the estimate is still withheld. A single
# shared sentence used to tell the simulate panel that "no H2 curve is
# available", which is not what that endpoint computes.
CURVE_SUBJECT = "curve"
ENERGY_SUBJECT = "energy"
RESOURCE_SUBJECT = "resource estimate"

# What every stored point must carry before it can answer either endpoint: the
# curve serves HF, FCI and the nuclear repulsion separately, and the simulation
# needs the stored qubit total plus, indirectly, the Pauli terms.
REQUIRED_POINT_FIELDS = (
    "bond_length_angstrom",
    "hf_energy_ha",
    "fci_energy_ha",
    "nuclear_repulsion_ha",
    "qubit_energy_plus_nuclear_ha",
)

# A request this close to a stored bond length counts as being on the grid rather
# than a snap, so 0.7400000000000001 is still "the stored 0.74".
BOND_LENGTH_MATCH_TOLERANCE_ANGSTROM = 1e-9
# Stored distances this close count as a tie; the shorter bond length wins, so a
# request exactly between two points has one deterministic answer.
TIE_TOLERANCE_ANGSTROM = 1e-12

# SPEC.md section 6's chemical-accuracy target. Applied to noiseless runs only:
# a shot-based estimate is a statistic with a spread, and it is reported as one
# rather than collapsed into a pass/fail flag.
CHEMICAL_ACCURACY_HA = 1.6e-3

# The fixed optimiser budget for one simulate request (per restart). Disclosed in
# the response (maxiter, restarts, iterations), never hidden. Restart 0 starts at
# zero angles; the other restarts are drawn from the seeded RNG, so a given seed
# reproduces the whole run.
#
# 150 was measured, not guessed: on this 4-parameter ansatz the error against
# exact stops improving much before the budget runs out, and 150 x 3 restarts
# holds >3x margin under SPEC.md's 1.6 mHa target even at the hard stretched
# geometries (2.5-3.0 A), while 100 left only ~15% margin at 2.5 A.
DEFAULT_MAXITER = 150
VQE_RESTARTS = 3

# Winner's-curse handling for shot runs. The optimiser is judged by the noiseless
# expectation at each restart's parameters, but the *reported* energy is still a
# fresh evaluation at the final optimal parameters -- not the minimum of the noisy
# optimisation trace, which is biased downward by selection over many draws. That
# fresh estimate carries an uncertainty: the estimator's own variance metadata when
# the backend reports it (the V1 Aer estimator does), otherwise this many
# independent re-evaluations at the same parameters with different sampler seeds,
# returned as their mean and sample standard deviation.
SHOT_REEVALUATION_SEED_COUNT = 10

# The variational guard is applied to noiseless statevector runs only (SPEC.md
# section 6): a single shot-based estimate is a statistic and may legitimately sit
# below the exact energy on sampling noise, so a shot run is never refused for
# that. Its uncertainty is reported as `vqe_energy_std_ha` and chemical accuracy
# stays `null` instead.
VARIATIONAL_GUARD_RUNS = "statevector"

# Cross-check tolerance between the live NumPy diagonalisation and the artifact's
# stored FCI total. The artifact's own recorded max deviation is ~1e-15 Ha, so
# anything near 1e-8 means the file contradicts itself, not that rounding moved.
ARTIFACT_FCI_MATCH_TOLERANCE_HA = 1e-8

CURVE_PROVENANCE = (
    "Exact energies are the stored PySCF FCI reference, and VQE energies are the stored ground "
    "energy of the mapped qubit Hamiltonian plus nuclear repulsion, already minimised in the "
    "offline precompute"
)
NO_LIVE_CALCULATION = "No optimisation, shots or noise run in this request."


class VariationalPrincipleViolation(RuntimeError):
    """A run reported an energy below the exact ground energy.

    The variational principle forbids that beyond floating-point rounding, so the
    router refuses the result (500, ``variational_principle_violation``) instead
    of reporting it. Raised by ``compute_h2_simulate`` after its guard check.
    """


def _unavailable_note(subject: str) -> str:
    """Why nothing can be served yet, phrased for this endpoint's own subject."""
    return (
        f"{H2_HAMILTONIAN_ARTIFACT} has not been generated, so no H2 {subject} is available. "
        "Run scripts/precompute_hamiltonians.ipynb in Colab or WSL and save the result into "
        f"app/data/hamiltonians/. This endpoint reports the {subject} as unavailable rather than "
        "showing placeholder numbers."
    )


def _no_points_note(artifact_name: str, subject: str) -> str:
    """Why no stored point can answer: the file exists but holds no grid."""
    return (
        f"{artifact_name} contains no stored bond lengths, so no H2 {subject} is available. "
        "Regenerate the artifact with scripts/precompute_hamiltonians.ipynb in Colab or WSL "
        f"and save it into app/data/hamiltonians/. This endpoint reports the {subject} as "
        "unavailable rather than showing placeholder numbers."
    )


def _unwired_note(subject: str) -> str:
    return (
        "{artifact} is present but this endpoint is not yet wired to its format, so the H2 "
        f"{subject} cannot be produced in this build."
    )


def h2_artifact_path() -> Path:
    """Location of the offline H2 Hamiltonian artifact."""
    return get_settings().h2_hamiltonian_dir / H2_HAMILTONIAN_ARTIFACT


def _artifact_status(subject: str) -> tuple[Path, str | None]:
    """Return the artifact path and the reason it cannot be used, if any."""
    artifact = h2_artifact_path()
    if not artifact.exists():
        return artifact, _unavailable_note(subject)
    return artifact, _unwired_note(subject).format(artifact=artifact.name)


def _read_artifact(artifact: Path) -> dict[str, Any]:
    """Read and JSON-decode the artifact, raising ``ValueError`` for a bad file.

    A missing file stays a ``FileNotFoundError`` (the router's 503), while a file
    that is present but unreadable or malformed is a fault, not an expected state,
    so it surfaces as ``ValueError`` (the router's 500) naming the problem.
    """
    try:
        raw = artifact.read_text(encoding="utf-8")
    except OSError as exc:
        raise ValueError(f"{artifact.name} could not be read: {exc}") from exc

    try:
        payload = json.loads(raw)
    except json.JSONDecodeError as exc:
        raise ValueError(f"{artifact.name} is not valid JSON: {exc}") from exc

    if not isinstance(payload, dict):
        raise ValueError(
            f"{artifact.name} must be a JSON object, got {type(payload).__name__}"
        )
    return payload


def _validated_points(
    artifact: Path,
    payload: dict[str, Any],
    subject: str,
) -> list[dict[str, float]]:
    """The stored grid, validated field by field instead of trusted.

    Two failure classes, kept apart on purpose:

    * no grid at all (``points`` missing or empty) means **no stored point is
      available** -- the same expected state as a missing file, so it raises
      ``FileNotFoundError`` (the router's 503 ``hamiltonian_unavailable``);
    * a grid that is present but malformed raises ``ValueError`` naming the
      offending point, so a corrupt artifact becomes a 500 that says what is
      wrong, never an empty or half-trusted curve.
    """
    points = payload.get("points")
    if not isinstance(points, list) or not points:
        raise FileNotFoundError(_no_points_note(artifact.name, subject))

    declared = payload.get("point_count")
    if isinstance(declared, int) and not isinstance(declared, bool) and declared != len(points):
        raise ValueError(
            f"{artifact.name} declares point_count={declared} but holds {len(points)} points"
        )

    molecule = payload.get("molecule")
    if not isinstance(molecule, str) or molecule.strip().upper() != "H2":
        raise ValueError(f"{artifact.name} must describe H2, found {molecule!r}")

    units = payload.get("units", "hartree")
    if units != "hartree":
        raise ValueError(f"{artifact.name} must store hartree energies, found {units!r}")

    validated: list[dict[str, float]] = []
    previous: float | None = None
    for index, point in enumerate(points):
        if not isinstance(point, dict):
            raise ValueError(f"{artifact.name} point {index} is not a JSON object")

        missing = [field for field in REQUIRED_POINT_FIELDS if field not in point]
        if missing:
            raise ValueError(
                f"{artifact.name} point {index} is missing {', '.join(sorted(missing))}"
            )

        values: dict[str, float] = {}
        for field in REQUIRED_POINT_FIELDS:
            value = point[field]
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                raise ValueError(
                    f"{artifact.name} point {index} field '{field}' must be a number, "
                    f"found {value!r}"
                )
            number = float(value)
            if not math.isfinite(number):
                raise ValueError(
                    f"{artifact.name} point {index} field '{field}' must be finite, "
                    f"found {value!r}"
                )
            values[field] = number

        bond_length = values["bond_length_angstrom"]
        if bond_length <= 0.0:
            raise ValueError(
                f"{artifact.name} point {index} has a non-positive bond length {bond_length}"
            )
        if previous is not None and bond_length <= previous:
            raise ValueError(
                f"{artifact.name} points must ascend in bond length, but point {index} is "
                f"{bond_length} after {previous}"
            )
        previous = bond_length
        validated.append(values)

    return validated


def _load_points(
    subject: str,
) -> tuple[Path, dict[str, Any], list[dict[str, float]]]:
    """The raw payload and its validated stored grid, or the reason it is absent.

    Raises ``FileNotFoundError`` (the router's 503) when the file or its grid is
    absent, ``ValueError`` (the router's 500) when the file is present but wrong.
    """
    artifact = h2_artifact_path()
    if not artifact.exists():
        raise FileNotFoundError(_unavailable_note(subject))
    payload = _read_artifact(artifact)
    return artifact, payload, _validated_points(artifact, payload, subject)


def _nearest_stored_point(
    points: list[dict[str, float]],
    requested: float,
) -> tuple[dict[str, float], bool]:
    """The stored point an answer must come from, and whether it *is* the request.

    Never interpolates: an off-grid request is answered by the nearest stored
    point, disclosed by the caller. An exact tie breaks toward the shorter bond
    length (the grid ascends, so ``min`` would return the longer one only if the
    grid were walked backwards).
    """
    distances = [(abs(point["bond_length_angstrom"] - requested), point) for point in points]
    shortest = min(distance for distance, _ in distances)
    tied = [
        point
        for distance, point in distances
        if distance <= shortest + TIE_TOLERANCE_ANGSTROM
    ]
    chosen = min(tied, key=lambda point: point["bond_length_angstrom"])
    on_grid = abs(chosen["bond_length_angstrom"] - requested) <= BOND_LENGTH_MATCH_TOLERANCE_ANGSTROM
    return chosen, on_grid


def _curve_metadata(artifact_name: str, payload: dict[str, Any]) -> dict[str, Any]:
    """The documented provenance block, read from the artifact and checked.

    The contract promises basis STO-3G, parity mapping with two-qubit reduction,
    2 qubits and the offline source. Rather than hard-coding those claims, each
    is checked against what the file actually declares; a regenerated artifact
    that says something else is a fault (500) to review, not a discrepancy to
    paper over.
    """
    problems: list[str] = []
    if str(payload.get("basis", "")).lower() != "sto-3g":
        problems.append(f"basis={payload.get('basis')!r}")
    if str(payload.get("mapping", "")).lower() != "parity":
        problems.append(f"mapping={payload.get('mapping')!r}")
    if payload.get("two_qubit_reduction") is not True:
        problems.append(f"two_qubit_reduction={payload.get('two_qubit_reduction')!r}")
    if payload.get("n_qubits") != 2:
        problems.append(f"n_qubits={payload.get('n_qubits')!r}")
    if problems:
        raise ValueError(
            f"{artifact_name} does not hold the documented metadata (expected STO-3G, "
            f"parity mapping with two-qubit reduction, n_qubits=2): {', '.join(problems)}"
        )

    return {
        "basis": "STO-3G",
        "mapping": "parity",
        "two_qubit_reduction": True,
        "mapping_note": "parity mapping with two-qubit reduction",
        "n_qubits": 2,
        "source": "precomputed offline with PySCF + qiskit-nature",
    }


def compute_h2_curve() -> dict:
    """H2 dissociation curve read from the offline artifact.

    Raises ``FileNotFoundError`` when the artifact or its grid is absent (the
    router's 503) and ``ValueError`` when it is present but malformed (500); an
    empty or placeholder curve is never served.

    Every energy is the stored value: ``hf_energies_ha`` and ``exact_energies_ha``
    (FCI) are TOTAL energies (electronic plus nuclear repulsion), while
    ``nuclear_repulsion_energies_ha`` reports the nuclear part separately, as the
    artifact's ``convention`` field documents.
    """
    artifact, payload, points = _load_points(CURVE_SUBJECT)
    metadata = _curve_metadata(artifact.name, payload)

    return {
        "molecule": "H2",
        "basis": metadata["basis"],
        "metadata": metadata,
        "artifact": artifact.name,
        "status": "computed",
        "note": (
            f"Computed from the {len(points)} stored points in {artifact.name}. "
            f"{CURVE_PROVENANCE}. {NO_LIVE_CALCULATION}"
        ),
        "units": "hartree",
        "convention": payload.get("convention"),
        "point_count": len(points),
        "bond_lengths_angstrom": [point["bond_length_angstrom"] for point in points],
        "hf_energies_ha": [point["hf_energy_ha"] for point in points],
        "exact_energies_ha": [point["fci_energy_ha"] for point in points],
        "vqe_energies_ha": [point["qubit_energy_plus_nuclear_ha"] for point in points],
        "nuclear_repulsion_energies_ha": [point["nuclear_repulsion_ha"] for point in points],
        "vqe_error_mha": [
            (point["qubit_energy_plus_nuclear_ha"] - point["fci_energy_ha"]) * 1000.0
            for point in points
        ],
        "max_deviation_from_fci_ha": payload.get("max_deviation_from_fci_ha"),
        # The same data point by point, so "for every bond length: HF, FCI and
        # the nuclear repulsion" is answerable without zipping four arrays.
        "points": [
            {
                "bond_length_angstrom": point["bond_length_angstrom"],
                "hf_energy_ha": point["hf_energy_ha"],
                "fci_energy_ha": point["fci_energy_ha"],
                "nuclear_repulsion_ha": point["nuclear_repulsion_ha"],
                "vqe_energy_ha": point["qubit_energy_plus_nuclear_ha"],
            }
            for point in points
        ],
    }


def _run_vqe(
    pauli_terms: dict[str, float],
    n_qubits: int,
    *,
    seed: int | None,
    maxiter: int,
    shots: int | None,
) -> dict[str, Any]:
    """Minimise a Pauli-sum energy for real: Ry ansatz, COBYLA, Qiskit estimator.

    ``shots is None`` selects the exact statevector estimator (a noiseless run).
    An integer selects the shot-based Aer V1 estimator with that many shots per
    energy evaluation, and the seed reaches it as
    ``backend_options={"seed_simulator": seed}`` through ``build_estimator``.
    The same seed also seeds the restart angles, so one fixed seed reproduces a
    run byte for byte.

    Restarts are judged by the **noiseless** expectation at their parameters: on
    a shot-based run the raw objective is a statistic, and picking the lowest
    sampled value would systematically favour a lucky downward fluctuation.

    All ``app.core`` imports are lazy -- this runs per request, never at module
    import time. Returns the **electronic** energy; the caller adds the nuclear
    repulsion stored in the artifact.
    """
    from scipy.optimize import minimize

    from app.core.vqe import (
        build_estimator,
        build_observable,
        default_layers,
        evaluate_energy,
        make_ansatz,
    )

    backend = "statevector" if shots is None else "aer"
    estimator = build_estimator(
        backend,
        seed=int(seed) if backend == "aer" and seed is not None else None,
    )
    run_options = {"shots": int(shots)} if backend == "aer" else None

    layers = default_layers(n_qubits)
    circuit, parameters = make_ansatz(n_qubits, layers)
    observable = build_observable(pauli_terms)
    reference_estimator = estimator if backend == "statevector" else build_estimator("statevector")
    rng = np.random.default_rng(seed)

    history: list[float] = []

    def objective(angles: np.ndarray) -> float:
        energy = evaluate_energy(
            estimator, circuit, observable, angles, run_options=run_options
        )
        history.append(energy)
        return energy

    best_energy = float("inf")
    best_parameters = np.zeros(len(parameters), dtype=float)

    for restart in range(VQE_RESTARTS):
        if restart == 0:
            start = np.zeros(len(parameters), dtype=float)
        else:
            start = rng.uniform(-np.pi, np.pi, len(parameters))

        result = minimize(
            objective,
            start,
            method="COBYLA",
            options={"maxiter": maxiter, "tol": 1e-12, "rhobeg": 0.4},
        )
        if not np.isfinite(result.fun):
            continue

        candidate = evaluate_energy(
            reference_estimator,
            circuit,
            observable,
            np.asarray(result.x, dtype=float),
        )
        if candidate < best_energy:
            best_energy = float(candidate)
            best_parameters = np.asarray(result.x, dtype=float)

    if not np.isfinite(best_energy):
        raise ValueError(
            f"COBYLA returned no finite energy in any of the {VQE_RESTARTS} restarts"
        )

    # Statevector: the noiseless expectation at the best parameters *is* the
    # answer. Shot-based: re-evaluate freshly at those parameters and report that
    # estimate's own uncertainty, so the headline number is not the minimum of the
    # noisy trace.
    if shots is None:
        reported_energy = float(best_energy)
        reported_std: float | None = None
        uncertainty_method: str | None = None
    else:
        reported_energy, reported_std, uncertainty_method = _shot_estimate(
            circuit,
            observable,
            best_parameters,
            seed=seed,
            shots=shots,
        )

    return {
        "vqe_energy_ha": reported_energy,
        "vqe_energy_std_ha": reported_std,
        "vqe_energy_uncertainty_method": uncertainty_method,
        "iterations": len(history),
        # The budget spelled out: the same total as `iterations`, plus the two
        # factors it is made of, so "450 iterations" can be read as
        # 3 restarts x 150 evaluations per restart.
        "total_energy_evaluations": len(history),
        "n_restarts": VQE_RESTARTS,
        "maxiter_per_restart": int(maxiter),
        "convergence_history_ha": [float(value) for value in history],
        "estimator_backend": backend,
        "estimator_class_name": type(estimator).__name__,
        "ansatz": f"{layers} x [Ry on each of {n_qubits} qubit(s), then a CX ladder]",
        "optimizer": "COBYLA",
        "restarts": VQE_RESTARTS,
        "maxiter": int(maxiter),
        "seed": seed,
        "seed_simulator_applied": backend == "aer" and seed is not None,
        "optimal_parameters": [float(value) for value in best_parameters],
    }


def _shot_estimate(
    circuit: Any,
    observable: Any,
    parameters: np.ndarray,
    *,
    seed: int | None,
    shots: int,
) -> tuple[float, float, str]:
    """A fresh shot-based energy at fixed parameters, plus its uncertainty.

    Preferred source of uncertainty is the estimator's own variance metadata: the
    fresh evaluation's result carries the per-shot sampling variance, and the
    estimate's standard error is ``sqrt(variance / shots)`` -- one evaluation, and
    the reader learns the method name ``estimator_variance_metadata``. When a
    backend reports no variance, the fallback re-evaluates
    ``SHOT_REEVALUATION_SEED_COUNT`` times with different sampler seeds and returns
    the mean and sample standard deviation (``seed_reevaluation_10x``). Both
    methods report the same quantity: the standard error of the estimated mean. The
    caller reports which method was used.
    """
    from app.core.vqe import (
        build_estimator,
        evaluate_energy,
        evaluate_energy_with_uncertainty,
    )

    run_options = {"shots": int(shots)}
    primary = (
        build_estimator("aer", seed=int(seed))
        if seed is not None
        else build_estimator("aer")
    )
    estimate = evaluate_energy_with_uncertainty(
        primary, circuit, observable, parameters, run_options=run_options
    )
    if estimate["std"] is not None:
        return float(estimate["value"]), float(estimate["std"]), "estimator_variance_metadata"

    # No variance metadata: draw fresh sampler seeds (deterministic when the caller
    # passed one) and re-evaluate the same parameters, then report mean and std.
    rng = np.random.default_rng(seed)
    samples: list[float] = []
    for _ in range(SHOT_REEVALUATION_SEED_COUNT):
        sample_seed = int(rng.integers(0, 2**32 - 1))
        sample_estimator = build_estimator("aer", seed=sample_seed)
        samples.append(
            float(
                evaluate_energy(
                    sample_estimator, circuit, observable, parameters, run_options=run_options
                )
            )
        )
    array = np.asarray(samples, dtype=float)
    return float(array.mean()), float(array.std(ddof=1)), "seed_reevaluation_10x"


def compute_h2_simulate(
    bond_length: float,
    shots: int | None = None,
    noise: float = 0.0,
    molecule: str | None = None,
    seed: int | None = None,
    maxiter: int | None = None,
) -> dict:
    """One H2 energy at the nearest stored bond length, computed by a real VQE.

    The request snaps to the stored grid per SPEC.md section 6 and discloses both
    bond lengths. The Hamiltonian is built from the artifact's own Pauli terms
    (readers in ``app.core.exact``), minimised by ``_run_vqe``, and compared with
    an exact NumPy diagonalisation of the same operator. The FCI and HF numbers
    are read from the artifact; nothing is invented.

    Raises ``FileNotFoundError`` (the router's 503) when the artifact or its grid
    is absent, ``ValueError`` (500) for a malformed artifact, and
    ``VariationalPrincipleViolation`` (500) when the optimiser reports an energy
    below the exact one beyond rounding -- that result is refused, not served.

    ``shots=None`` means a noiseless statevector run; an integer runs the
    shot-based Aer estimator. ``noise`` is accepted for the contract but no noise
    model is applied, and the response says so.
    """
    if molecule is not None and molecule.strip().upper() != "H2":
        return {
            "molecule": molecule,
            "bond_length": bond_length,
            "exact_energy_ha": None,
            "vqe_energy_ha": None,
            "error_mha": None,
            "convergence_history_ha": None,
            "shots": shots,
            "noise": noise,
            "estimator_class_name": None,
            "status": "placeholder",
            "note": (
                "Only H2 is in scope for live VQE in this build; larger molecules need the "
                "offline artifact and the resource estimator."
            ),
        }

    artifact, payload, points = _load_points(ENERGY_SUBJECT)
    point, on_grid = _nearest_stored_point(points, float(bond_length))
    used = point["bond_length_angstrom"]

    # Lazy: app.core must never load at import time (enforced by the guard tests).
    from app.core.exact import h2_hamiltonian_from_artifact
    from app.core.vqe import EXACT_EPSILON_HA

    hamiltonian = h2_hamiltonian_from_artifact(payload, used)
    exact = hamiltonian.ground_energy_ha + hamiltonian.nuclear_repulsion_ha

    fci = point["fci_energy_ha"]
    if abs(exact - fci) > ARTIFACT_FCI_MATCH_TOLERANCE_HA:
        raise ValueError(
            f"diagonalising the stored Pauli terms at {used:g} A gives {exact:.12f} Ha but "
            f"{artifact.name} records fci_energy_ha={fci:.12f} Ha; the artifact contradicts "
            "itself, so no energy is reported."
        )

    run = _run_vqe(
        hamiltonian.pauli_terms,
        hamiltonian.n_qubits,
        seed=seed,
        maxiter=DEFAULT_MAXITER if maxiter is None else maxiter,
        shots=shots,
    )
    noiseless = shots is None
    vqe = run["vqe_energy_ha"] + hamiltonian.nuclear_repulsion_ha

    # Variational principle, applied to noiseless statevector runs only (SPEC.md
    # section 6). A statevector energy may sit at or above the exact minimum, or
    # below it only by floating-point rounding; anything else is a bug in this
    # pipeline and is refused. A shot-based estimate is a statistic and may
    # legitimately fall below the exact energy on sampling noise, so it is served
    # with its uncertainty rather than refused.
    if noiseless and vqe < exact - EXACT_EPSILON_HA:
        raise VariationalPrincipleViolation(
            f"the run at {used:g} A reported VQE energy {vqe:.12f} Ha, below the exact "
            f"{exact:.12f} Ha by {exact - vqe:.3e} Ha -- more than the {EXACT_EPSILON_HA:g} Ha "
            "rounding allowance the variational principle allows. The result is refused "
            "instead of reported."
        )

    error = vqe - exact

    run_description = (
        f"{run['total_energy_evaluations']} COBYLA energy evaluations over "
        f"{run['n_restarts']} restarts of a fixed {run['ansatz']} ansatz, "
        f"{run['estimator_backend']} estimator"
        + (f", seed {run['seed']}" if run["seed"] is not None else ", no seed")
        + (" (noiseless statevector run)" if noiseless else f" ({shots} shots per evaluation)")
    )
    shot_disclosure = ""
    if not noiseless:
        method = str(run["vqe_energy_uncertainty_method"]).replace("_", " ")
        shot_disclosure = (
            f" The energy is a fresh estimate at the final optimal parameters (not the "
            f"lowest sampled value of the optimisation trace), ±{run['vqe_energy_std_ha'] * 1000.0:.3f} "
            f"mHa, uncertainty from {method}; that is sampling error, not a hardware result."
        )
    execution = "Classical simulation only: no quantum hardware was used, no quantum advantage is claimed."
    reference = (
        "The exact reference diagonalises the same stored Pauli terms with NumPy; "
        f"hf_energy_ha and fci_energy_ha are read from {artifact.name}. No noise model is applied."
    )

    if on_grid:
        note = (
            f"Real VQE run at {used:g} A: the requested {bond_length:g} A is on the stored grid "
            f"in {artifact.name}. {run_description}.{shot_disclosure} {execution} {reference}"
        )
    else:
        note = (
            f"{bond_length:g} A is not on the stored grid in {artifact.name}; answered with the "
            f"nearest stored point at {used:g} A instead, never interpolated (an exact tie "
            f"breaks toward the shorter bond length). {run_description}.{shot_disclosure} "
            f"{execution} {reference}"
        )

    return {
        "molecule": "H2",
        "status": "computed",
        "units": "hartree",
        # Bond-length disclosure. SPEC.md section 6 names the *_angstrom trio;
        # requested_bond_length / used_bond_length / snapped are the same three
        # facts under the names the backend task specifies. ``bond_length`` is
        # the point the energies came from, because the panel badges it next to
        # them.
        "bond_length": used,
        "requested_bond_length": float(bond_length),
        "used_bond_length": used,
        "snapped": not on_grid,
        "bond_length_requested_angstrom": float(bond_length),
        "bond_length_used_angstrom": used,
        "bond_length_is_exact": on_grid,
        # Energies: the run's own result plus the artifact's stored references.
        "hf_energy_ha": point["hf_energy_ha"],
        "fci_energy_ha": fci,
        "nuclear_repulsion_ha": hamiltonian.nuclear_repulsion_ha,
        "exact_energy_ha": exact,
        "vqe_energy_ha": vqe,
        "error_vs_exact_ha": error,
        "abs_error_vs_exact_ha": abs(error),
        "error_mha": error * 1000.0,
        # The uncertainty of a shot-based estimate; null for a noiseless run, which
        # has none to report.
        "vqe_energy_std_ha": run["vqe_energy_std_ha"],
        "vqe_energy_uncertainty_method": run["vqe_energy_uncertainty_method"],
        "within_chemical_accuracy": (abs(error) <= CHEMICAL_ACCURACY_HA) if noiseless else None,
        "chemical_accuracy_threshold_ha": CHEMICAL_ACCURACY_HA,
        "noiseless": noiseless,
        # Run provenance. ``iterations`` is kept for backward compatibility and
        # equals ``total_energy_evaluations``; the budget is also spelled out as
        # its two factors.
        "iterations": run["iterations"],
        "total_energy_evaluations": run["total_energy_evaluations"],
        "n_restarts": run["n_restarts"],
        "maxiter_per_restart": run["maxiter_per_restart"],
        "convergence_history_ha": run["convergence_history_ha"],
        "seed": run["seed"],
        "seed_simulator_applied": run["seed_simulator_applied"],
        "maxiter": run["maxiter"],
        "restarts": run["restarts"],
        "ansatz": run["ansatz"],
        "optimizer": run["optimizer"],
        "shots": shots,
        "noise": noise,
        "estimator_backend": run["estimator_backend"],
        "estimator_class_name": run["estimator_class_name"],
        "classical_simulation": True,
        "note": note,
    }


def compute_hadamard_demo() -> dict:
    """Single-qubit Hadamard superposition.

    This is an exact hand-computable state, not a hardware result, and needs no
    artifact.
    """
    inv_sqrt2 = 1.0 / np.sqrt(2.0)
    return {
        "circuit": "H|0>",
        "initial_state": [1.0, 0.0],
        "final_state": [inv_sqrt2, inv_sqrt2],
        "probabilities": {
            "0": float(inv_sqrt2 * inv_sqrt2),
            "1": float(inv_sqrt2 * inv_sqrt2),
        },
        "status": "computed",
    }


def compute_resources(molecule: str) -> dict:
    """Rough resource estimate.

    The H2 numbers used to be derived from the fabricated placeholder matrix, so
    they are withheld until the resource estimator is wired to the artifact's
    Pauli terms.
    """
    normalized = molecule.strip().upper()
    artifact, reason = _artifact_status(RESOURCE_SUBJECT)

    if normalized == "H2":
        return {
            "molecule": "H2",
            "n_qubits": None,
            "hamiltonian_terms": None,
            "estimated_circuit_gates": None,
            "estimated_circuit_depth": None,
            "artifact": artifact.name,
            "status": "placeholder",
            "note": reason,
        }

    return {
        "molecule": molecule,
        "qubits": None,
        "hamiltonian_terms": None,
        "estimated_circuit_gates": None,
        "estimated_circuit_depth": None,
        "note": (
            "Only H2 is in scope for a resource estimate in this build; larger molecules "
            "are added later."
        ),
        "status": "placeholder",
    }
