"""Manual VQE: a simple ansatz, one estimator primitive, and SciPy COBYLA.

Scope, stated plainly because it bounds what this module proves:

* It is a **core algorithm module**. ``app.services.quantum_service`` reaches it
  through *lazy* imports inside functions, so starting the app does not load it;
  the test suite enforces both halves of that rule: no module-level import of
  ``app.core`` anywhere in ``app/routers`` or ``app/services`` (function-level
  imports are allowed), and ``app.core`` absent from ``sys.modules`` after
  ``import app.main`` in a fresh interpreter.
* There are **no H2 coefficients in this file**. The Hamiltonians the unit tests
  exercise are synthetic fixtures defined inside the test file, chosen because
  their ground energies are known analytically; the API path passes in Pauli
  terms read from the offline artifact instead.
* The estimator is a real Qiskit primitive, used noiselessly, and the exact
  reference comes from ``app.core.exact`` -- which assembles the Pauli sum with
  NumPy alone, so the comparison is not circular.

Estimator backends and their pinned interfaces in this environment (qiskit
1.4.6, qiskit-aer 0.17.2). These were probed, not assumed:

===============================  ==================  ====================================
import path                      primitive version   call form
===============================  ==================  ====================================
qiskit.primitives.              BaseEstimatorV2     run([(circuit, observable, params)])
  StatevectorEstimator
qiskit_aer.primitives.          BaseEstimatorV1     run([circuit], [obs], [[params]])
  Estimator
qiskit_aer.primitives.          BaseEstimatorV2     run([(circuit, observable, params)])
  EstimatorV2
===============================  ==================  ====================================

Note that ``qiskit.primitives`` exports no ``EstimatorV2`` symbol in this version:
the V2 statevector estimator is named ``StatevectorEstimator``. The V1
``qiskit.primitives.Estimator`` also exists but is deprecated, so it is not used
here.

Accuracy caveat, measured rather than assumed: ``qiskit_aer.primitives.Estimator``
(the V1 one) is **shot-based and therefore approximate**. It defaults to 1024 shots
and only warns, rather than raising, when ``shots=None`` is requested; asking for
more shots reduces but does not remove the deviation. It is kept here because the
adapter has to speak both the V1 and V2 calling conventions, but it must never be
used as the reference an accuracy target is judged against. ``StatevectorEstimator``
and ``qiskit_aer.primitives.EstimatorV2`` are exact and agree to ~1e-9; those are
the ones the fixture tests grade against.

Reproducibility, also measured: the V1 estimator takes a sampler seed
(``backend_options={"seed_simulator": n}``), which makes repeated evaluations
byte-identical and identical across processes. Seeding pins the sampling noise; it
does not make the estimator exact. ``build_estimator`` exposes it as ``seed=``, and
refuses a seed on the exact backends rather than accepting it as a silent no-op.
"""

from __future__ import annotations

import importlib
import logging
import math
from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np
from numpy.typing import NDArray

from app.core.exact import pauli_sum_ground_energy, pauli_sum_matrix

logger = logging.getLogger(__name__)

# Documented tolerance for the "a variational energy may sit infinitesimally below
# the exact one due to floating point" allowance. Anything below this is a bug,
# not rounding, and is reported as such rather than silently accepted.
EXACT_EPSILON_HA = 1e-9

# Tolerance the fixture tests require the optimiser to reach.
TARGET_ACCURACY_HA = 1e-6

DEFAULT_BACKEND = "statevector"

# backend -> (module, attribute). Kept as data so the import path can be reported
# and asserted without importing anything.
ESTIMATOR_BACKENDS: dict[str, tuple[str, str]] = {
    "statevector": ("qiskit.primitives", "StatevectorEstimator"),
    "aer": ("qiskit_aer.primitives", "Estimator"),
    "aer_v2": ("qiskit_aer.primitives", "EstimatorV2"),
}

# Backends whose result is a finite-sample estimate and therefore depends on a
# sampler seed -> the keyword argument that carries the seed.
#
# Probed in this pinned environment (qiskit-aer 0.17.2), not assumed:
#   qiskit_aer.primitives.Estimator(backend_options={"seed_simulator": 1234}) gives
#   byte-identical values across repeated runs and across separate processes, while
#   the unseeded estimator wanders (observed spread 2.7e-2 over 5 runs at a true
#   <Z> = 0.5). EstimatorV2 rejects that option ("Options.__init__() got an
#   unexpected keyword argument 'seed_simulator'") and is already exact, and
#   StatevectorEstimator is analytic, so neither takes a sampler seed here.
SEEDED_ESTIMATOR_OPTIONS: dict[str, str] = {"aer": "backend_options"}


def estimator_import_path(backend: str = DEFAULT_BACKEND) -> str:
    """Fully-qualified import path for a backend, e.g. ``qiskit.primitives.StatevectorEstimator``."""
    try:
        module, attribute = ESTIMATOR_BACKENDS[backend]
    except KeyError:
        raise ValueError(
            f"Unknown estimator backend {backend!r}; expected one of {sorted(ESTIMATOR_BACKENDS)}"
        ) from None
    return f"{module}.{attribute}"


def build_estimator(backend: str = DEFAULT_BACKEND, *, seed: int | None = None) -> Any:
    """Instantiate an estimator primitive.

    ``seed`` is a sampler seed for the shot-based backends (see
    ``SEEDED_ESTIMATOR_OPTIONS``). Passing one to an exact backend is an error
    rather than a silent no-op, so a caller cannot believe it has made an exact
    estimator reproducible when the seed could not have done anything.

    Raises ``ImportError`` if the backend is unavailable rather than falling back
    to a different one, so a missing dependency is never silently masked.
    """
    module_name, attribute = ESTIMATOR_BACKENDS.get(backend, (None, None))
    if module_name is None:
        raise ValueError(
            f"Unknown estimator backend {backend!r}; expected one of {sorted(ESTIMATOR_BACKENDS)}"
        )

    module = importlib.import_module(module_name)
    estimator_class = getattr(module, attribute)

    if seed is None:
        return estimator_class()

    option_name = SEEDED_ESTIMATOR_OPTIONS.get(backend)
    if option_name is None:
        raise ValueError(
            f"The {backend!r} estimator is exact and takes no sampler seed; seeding is "
            f"supported for {sorted(SEEDED_ESTIMATOR_OPTIONS)}"
        )
    return estimator_class(**{option_name: {"seed_simulator": int(seed)}})


def make_ansatz(n_qubits: int, layers: int = 1) -> tuple[Any, list[Any]]:
    """Hardware-efficient ansatz: ``layers`` x (Ry on every qubit, then a CX ladder).

    Deliberately simple and fully parameterised by ``Ry`` angles, which is enough
    to represent the ground state of the small Ising fixtures the tests use. It is
    not a chemistry-inspired ansatz and makes no claim to be one.
    """
    if n_qubits < 1:
        raise ValueError("n_qubits must be at least 1")
    if layers < 1:
        raise ValueError("layers must be at least 1")

    from qiskit import QuantumCircuit
    from qiskit.circuit import Parameter

    circuit = QuantumCircuit(n_qubits)
    parameters: list[Any] = []

    for layer in range(layers):
        for qubit in range(n_qubits):
            parameter = Parameter(f"ry_{layer}_{qubit}")
            parameters.append(parameter)
            circuit.ry(parameter, qubit)
        for qubit in range(n_qubits - 1):
            circuit.cx(qubit, qubit + 1)

    return circuit, parameters


def build_observable(terms: Mapping[str, float]) -> Any:
    """Qiskit ``SparsePauliOp`` for a Pauli sum, sorted for deterministic output."""
    if not terms:
        raise ValueError("A Pauli sum needs at least one term")

    from qiskit.quantum_info import SparsePauliOp

    return SparsePauliOp.from_list([(key, float(value)) for key, value in sorted(terms.items())])


def _extract_expectation(result: Any) -> float:
    """Pull a single expectation value out of a primitive result.

    The V2 result exposes ``.data.evs``, which is a **0-dimensional** array for a
    single parameter set, so it cannot be indexed directly; the V1 result exposes
    ``.values``. Both are normalised to one float here.
    """
    values = getattr(result, "values", None)
    if values is None:
        values = result.data.evs
    flat = np.atleast_1d(np.asarray(values)).reshape(-1)
    if flat.size != 1:
        raise RuntimeError(f"Expected exactly one expectation value, got {flat.size}")
    return float(flat[0])


def _extract_variance(result: Any) -> float | None:
    """**Per-shot** variance of one V1 estimate, or ``None`` when none is reported.

    The V1 Aer estimator reports ``variance`` per observable in ``result.metadata``
    (probed in this pinned environment: a 1024-shot ``<ZZ+IX+XI>`` estimate carried
    ``metadata=[{"shots": 1024, "variance": 1.875...}]``). That number is the
    variance of a single shot's value, NOT of the mean: 40 independent 1024-shot
    estimates of the same expectation had sample std 0.052, while
    ``sqrt(variance)`` was 1.39 and ``sqrt(variance/shots)`` was 0.043. The
    estimate's standard error is therefore ``sqrt(variance / shots)`` and the
    caller does that division (see ``evaluate_energy_with_uncertainty``). The
    analytic statevector estimator reports neither variance nor std, which is the
    signal that there is no sampling uncertainty to quote -- an exact value is not
    a statistic.
    """
    metadata = getattr(result, "metadata", None)
    if not isinstance(metadata, Sequence) or isinstance(metadata, (str, bytes)) or not metadata:
        return None
    entry = metadata[0]
    if not isinstance(entry, Mapping):
        return None
    raw = entry.get("variance")
    if isinstance(raw, bool) or not isinstance(raw, (int, float)):
        return None
    if not np.isfinite(float(raw)) or float(raw) < 0.0:
        return None
    return float(raw)


def evaluate_energy_with_uncertainty(
    estimator: Any,
    circuit: Any,
    observable: Any,
    parameters: Sequence[float] | NDArray,
    *,
    run_options: Mapping[str, Any] | None = None,
) -> dict[str, float | None]:
    """One energy evaluation plus its sampling uncertainty, when the estimator has one.

    Returns ``{"value": float, "std": float | None}``, where ``std`` is the standard
    error of the estimate (the standard deviation of the reported *mean*, not of a
    single shot). For the V1 Aer estimator that is ``sqrt(variance / shots)`` from
    its per-observable ``variance`` metadata and the ``shots`` in ``run_options``;
    a V2 shot-based estimator's ``stds`` already is the estimate's standard error.
    An exact estimator reports ``None`` because an analytic expectation value
    carries no sampling error.

    ``run_options`` are forwarded verbatim to the V1 primitive's ``run`` call, so a
    shot-based run can ask for a concrete sample count (for example
    ``{"shots": 1024}`` on the Aer estimator). The exact V2/statevector estimators
    accept no per-run options, so passing any is an error rather than a silently
    ignored keyword.
    """
    values = np.asarray(parameters, dtype=float).reshape(-1)

    from qiskit.primitives import BaseEstimatorV2

    if isinstance(estimator, BaseEstimatorV2):
        if run_options:
            raise ValueError(
                "The statevector estimator takes no per-run options, got "
                f"{sorted(run_options)}; shots belong on the shot-based aer backend"
            )
        job = estimator.run([(circuit, observable, values)])
        result = job.result()[0]
        std: float | None = None
        stds = getattr(getattr(result, "data", None), "stds", None)
        if stds is not None:
            flat = np.atleast_1d(np.asarray(stds, dtype=float)).reshape(-1)
            if flat.size == 1 and np.isfinite(flat[0]) and flat[0] >= 0.0:
                std = float(flat[0])
        evs = np.atleast_1d(np.asarray(result.data.evs)).reshape(-1)
        if evs.size != 1:
            raise RuntimeError(f"Expected exactly one expectation value, got {evs.size}")
        return {"value": float(evs[0]), "std": std}

    # V1 primitives take parallel lists of circuits, observables and parameter sets.
    options = dict(run_options or {})
    job = estimator.run([circuit], [observable], [values], **options)
    result = job.result()
    variance = _extract_variance(result)
    std: float | None = None
    if variance is not None:
        shots = options.get("shots")
        divisor = float(shots) if isinstance(shots, (int, float)) and shots else 1.0
        std = math.sqrt(variance / divisor)
    return {"value": _extract_expectation(result), "std": std}


def evaluate_energy(
    estimator: Any,
    circuit: Any,
    observable: Any,
    parameters: Sequence[float] | NDArray,
    *,
    run_options: Mapping[str, Any] | None = None,
) -> float:
    """One energy evaluation: <psi(theta)|H|psi(theta)> for a single parameter set.

    The point-estimate half of ``evaluate_energy_with_uncertainty``; use that one
    when the estimate's sampling uncertainty matters (a shot-based estimate) and
    this one when only the value matters.
    """
    return float(
        evaluate_energy_with_uncertainty(
            estimator, circuit, observable, parameters, run_options=run_options
        )["value"]
    )


def default_layers(n_qubits: int) -> int:
    """Layer count used when the caller does not choose one.

    Two qubits get two layers because the ``ZZ + X`` fixtures are entangled and a
    single Ry+CX layer does not span their ground state. One qubit needs one.
    """
    return 1 if n_qubits < 2 else 2


def vqe_ground_energy(
    terms: Mapping[str, float],
    n_qubits: int | None = None,
    *,
    layers: int | None = None,
    restarts: int = 4,
    maxiter: int = 1200,
    seed: int = 7,
    backend: str = DEFAULT_BACKEND,
    rhobeg: float = 0.4,
) -> dict[str, Any]:
    """Minimise a Pauli-sum energy with COBYLA over the simple ansatz.

    The optimiser is deterministic: restart 0 always starts from all-zero angles
    and any further restarts come from a seeded generator, so two runs of the same
    fixture give the same answer.

    Returns measured quantities plus the tolerances the caller should judge them
    against, so no acceptance threshold is hidden inside this function.
    """
    if not terms:
        raise ValueError("A Pauli sum needs at least one term")

    widths = {len(key) for key in terms}
    if len(widths) != 1:
        raise ValueError(f"All Pauli strings must have the same length, got {sorted(widths)}")
    inferred_qubits = widths.pop()
    if n_qubits is None:
        n_qubits = inferred_qubits
    elif n_qubits != inferred_qubits:
        raise ValueError(f"Terms are {inferred_qubits} qubits wide but n_qubits={n_qubits}")

    if restarts < 1:
        raise ValueError("restarts must be at least 1")

    exact_energy = pauli_sum_ground_energy(terms, n_qubits)

    resolved_layers = default_layers(n_qubits) if layers is None else layers
    circuit, parameters = make_ansatz(n_qubits, resolved_layers)
    observable = build_observable(terms)
    estimator = build_estimator(backend)

    rng = np.random.default_rng(seed)
    evaluation_count = 0
    best_energy = float("inf")
    best_angles = np.zeros(len(parameters), dtype=float)
    best_converged = False
    last_result: Any = None

    def objective(angles: np.ndarray) -> float:
        nonlocal evaluation_count
        energy = evaluate_energy(estimator, circuit, observable, angles)
        evaluation_count += 1
        return energy

    from scipy.optimize import minimize

    for restart in range(restarts):
        if restart == 0:
            start = np.zeros(len(parameters), dtype=float)
        else:
            start = rng.uniform(-np.pi, np.pi, len(parameters))

        result = minimize(
            objective,
            start,
            method="COBYLA",
            options={"maxiter": maxiter, "tol": 1e-12, "rhobeg": rhobeg},
        )
        last_result = result

        if np.isfinite(result.fun) and result.fun < best_energy:
            best_energy = float(result.fun)
            best_angles = np.asarray(result.x, dtype=float)
            best_converged = bool(result.success)

    abs_error = abs(best_energy - exact_energy)
    below_exact = best_energy < exact_energy - EXACT_EPSILON_HA
    if below_exact:
        logger.warning(
            "VQE energy %.12f is below the exact %.12f by more than the documented "
            "epsilon %g; this indicates a bug rather than rounding.",
            best_energy,
            exact_energy,
            EXACT_EPSILON_HA,
        )

    return {
        "exact_energy_ha": exact_energy,
        "vqe_energy_ha": best_energy,
        "abs_error_ha": abs_error,
        "error_mha": (best_energy - exact_energy) * 1000.0,
        # Tolerances the caller should assert against, surfaced rather than hidden.
        "target_accuracy_ha": TARGET_ACCURACY_HA,
        "converged_within_target": bool(abs_error <= TARGET_ACCURACY_HA),
        "exact_epsilon_ha": EXACT_EPSILON_HA,
        "below_exact_beyond_epsilon": bool(below_exact),
        # Provenance of the run.
        "n_qubits": n_qubits,
        "layers": resolved_layers,
        "n_parameters": len(parameters),
        "optimal_parameters": [float(value) for value in best_angles],
        "restarts": restarts,
        "maxiter": maxiter,
        "random_seed": seed,
        "energy_evaluations": evaluation_count,
        "optimizer_success": best_converged,
        "optimizer_message": getattr(last_result, "message", None),
        "estimator_backend": backend,
        "estimator_class_name": type(estimator).__name__,
        "estimator_module": type(estimator).__module__,
        "ansatz": (
            f"{resolved_layers} x [Ry on each of {n_qubits} qubit(s), then a CX ladder]"
        ),
        "hamiltonian_terms": {key: float(value) for key, value in sorted(terms.items())},
        "hamiltonian_matrix_hermitian": bool(
            np.allclose(
                pauli_sum_matrix(terms, n_qubits),
                pauli_sum_matrix(terms, n_qubits).conj().T,
            )
        ),
        "note": (
            "Variational upper bound on the ground energy of the given Pauli sum. "
            "This function is exercised on synthetic test fixtures; the API's VQE "
            "loop lives in app.services.quantum_service, which builds its own COBYLA "
            "run from these primitives and passes in Pauli terms read from the "
            "offline artifact."
        ),
    }
