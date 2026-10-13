"""Tests for the exact solver, the manual VQE and the resource estimator.

Every Hamiltonian in this file is a **synthetic test fixture** written out below.
They are not data, they are not chemistry, and they are not the H2 problem: each
one was chosen because its ground energy is known in closed form, so the exact
solver can be checked against arithmetic rather than against itself.

Qubit-ordering convention (matches Qiskit, asserted in
``test_pauli_string_matrix_matches_qiskit_ordering``): a Pauli string is read left
to right and the RIGHTMOST character acts on qubit 0, i.e. ``"IX"`` means X on
qubit 0 and ``"XI"`` means X on qubit 1.

These fixtures must never be reachable from the API *at import time*.
``app/core`` may only be pulled in by lazy imports inside functions (that is how
the quantum endpoints reach the exact solver and the VQE machinery);
``test_app_core_is_not_imported_eagerly_by_routers_or_services`` checks that
statically -- module-level imports forbidden, function-level imports allowed --
and ``test_app_core_is_absent_from_sys_modules_after_importing_the_app``
re-checks it dynamically in a fresh interpreter.

Every Hamiltonian below is a synthetic fixture written for this file; the real
H2 Hamiltonian lives in the offline artifact
(``app/data/hamiltonians/h2_sto3g.json``) and is not referenced here.
"""

from __future__ import annotations

import ast
import importlib
import math
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from app.core.exact import (
    PAULI_MATRICES,
    diagonalize_ground_state,
    h2_sto3g_hamiltonian_from_json,
    pauli_string_matrix,
    pauli_sum_ground_energy,
    pauli_sum_matrix,
    pauli_terms_from_sequence,
)
from app.core.resources import estimate_h2_resources
from app.core.vqe import (
    DEFAULT_BACKEND,
    ESTIMATOR_BACKENDS,
    EXACT_EPSILON_HA,
    TARGET_ACCURACY_HA,
    build_observable,
    build_estimator,
    default_layers,
    estimator_import_path,
    evaluate_energy,
    make_ansatz,
    vqe_ground_energy,
)

REPO_ROOT = Path(__file__).resolve().parents[1]

# The V1 Aer estimator is shot-based (1024 shots by default) and only warns when
# asked for shots=None, so its values carry sampling noise on top of the true
# expectation. This band is what that noise looks like at 1024 shots; it is
# explicitly NOT an accuracy claim about the VQE, which is graded against the
# exact estimator at 1e-6.
V1_SHOT_NOISE_BAND = 0.05

# The V1 estimator's sampler can be seeded, so the test below no longer depends on
# shot luck. Measured before this was added: the +/-0.05 band is about 1.85 sigma of
# the 1024-shot sampling noise, so the seam was crossed in roughly 1 run in 15
# (observed: 1 failure in 10 runs). Seeding was probed rather than assumed --
# ``qiskit_aer.primitives.Estimator(backend_options={"seed_simulator": 1234})``
# returns byte-identical values across repeated evaluations and across separate
# processes, where the unseeded estimator wandered by up to 2.7e-2. The band stays
# because it still describes the estimator's sampling error: the seed pins the noise
# so the test is reproducible, it does not make the estimator exact.
V1_SAMPLER_SEED = 1234


# ---------------------------------------------------------------------------
# Test fixtures: (label, Pauli sum, ground energy known in closed form)
#
# The closed forms are derived by hand, not read off a diagonaliser:
#   Z            -> eigenvalues +/-1                       -> -1
#   ZZ           -> eigenvalues +/-1                       -> -1
#   ZZ + X(0)    -> two 2x2 blocks [[+/-1, 1], [1, -1]]    -> -sqrt(2)
#   ZZ + X(0,1)  -> symmetric; ground state is the          -> -sqrt(5)
#                   uniform superposition, energy -sqrt(5)
# ---------------------------------------------------------------------------
FIXTURES: list[tuple[str, dict[str, float], float]] = [
    ("1q_Z", {"Z": 1.0}, -1.0),
    ("2q_ZZ", {"ZZ": 1.0}, -1.0),
    ("2q_ZZ_plus_X0", {"ZZ": 1.0, "IX": 1.0}, -math.sqrt(2.0)),
    ("2q_ZZ_plus_X0_plus_X1", {"ZZ": 1.0, "IX": 1.0, "XI": 1.0}, -math.sqrt(5.0)),
]

FIXTURE_IDS = [label for label, _, _ in FIXTURES]
FIXTURE_CASES = [(label, terms, expected) for label, terms, expected in FIXTURES]


# ---------------------------------------------------------------------------
# Exact solver: matrix path
# ---------------------------------------------------------------------------
def test_exact_solver_returns_lowest_eigenvalue() -> None:
    matrix = np.array([[1.0, 0.5], [0.5, 1.0]], dtype=np.float64)
    assert abs(diagonalize_ground_state(matrix) - 0.5) < 1e-9


def test_exact_solver_rejects_a_non_square_matrix() -> None:
    with pytest.raises(ValueError):
        diagonalize_ground_state(np.zeros((2, 3)))


def test_exact_solver_is_symmetric_real() -> None:
    rng = np.random.default_rng(17)
    raw = rng.normal(size=(4, 4))
    hermitian = (raw + raw.T) / 2.0
    assert diagonalize_ground_state(hermitian) == pytest.approx(
        float(np.linalg.eigvalsh(hermitian)[0].real), abs=1e-9
    )


# ---------------------------------------------------------------------------
# Exact solver: Pauli-sum path
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("pauli_string", ["I", "X", "Y", "Z", "IX", "XI", "ZZ", "XY", "ZY"])
def test_pauli_string_matrix_matches_qiskit_ordering(pauli_string: str) -> None:
    """The NumPy kron order must agree with Qiskit, or every VQE comparison is void."""
    from qiskit.quantum_info import SparsePauliOp

    qiskit_matrix = SparsePauliOp.from_list([(pauli_string, 1.0)]).to_matrix()
    assert np.allclose(pauli_string_matrix(pauli_string), qiskit_matrix), (
        f"Pauli ordering mismatch for {pauli_string!r}: the exact reference and the "
        f"VQE observable would disagree about which qubit each character acts on."
    )


def test_pauli_matrix_table_is_unitary_and_hermitian() -> None:
    for name, matrix in PAULI_MATRICES.items():
        assert np.allclose(matrix @ matrix.conj().T, np.eye(2)), f"{name} not unitary"
        assert np.allclose(matrix, matrix.conj().T), f"{name} not Hermitian"


@pytest.mark.parametrize("label,terms,expected", FIXTURE_CASES, ids=FIXTURE_IDS)
def test_fixture_ground_energies_match_the_closed_form(
    label: str, terms: dict[str, float], expected: float
) -> None:
    """The exact solver must reproduce the hand-derived ground energy."""
    assert pauli_sum_ground_energy(terms) == pytest.approx(expected, abs=1e-12)


@pytest.mark.parametrize("label,terms,expected", FIXTURE_CASES, ids=FIXTURE_IDS)
def test_fixture_pauli_sum_matrix_is_hermitian_and_dimensioned(
    label: str, terms: dict[str, float], expected: float
) -> None:
    width = len(next(iter(terms)))
    matrix = pauli_sum_matrix(terms)
    assert matrix.shape == (2**width, 2**width)
    assert np.allclose(matrix, matrix.conj().T)


def test_pauli_sum_matrix_rejects_mismatched_widths() -> None:
    with pytest.raises(ValueError, match="same length"):
        pauli_sum_matrix({"Z": 1.0, "ZZ": 1.0})


def test_pauli_sum_matrix_rejects_an_invalid_character() -> None:
    with pytest.raises(ValueError, match="Invalid Pauli character"):
        pauli_sum_matrix({"ZQ": 1.0})


def test_pauli_sum_matrix_rejects_empty_and_mismatched_qubit_count() -> None:
    with pytest.raises(ValueError):
        pauli_sum_matrix({})
    with pytest.raises(ValueError, match="n_qubits"):
        pauli_sum_matrix({"ZZ": 1.0}, n_qubits=3)


def test_pauli_terms_from_sequence_sums_duplicates() -> None:
    merged = pauli_terms_from_sequence([("ZZ", 1.0), ("ZZ", 0.5), ("IX", -1.0)])
    assert merged == {"ZZ": 1.5, "IX": -1.0}


# ---------------------------------------------------------------------------
# Estimator probes (the import paths and interfaces in THIS pinned environment)
# ---------------------------------------------------------------------------
def test_estimator_import_paths_are_the_probed_ones() -> None:
    """Pin the probed FQNs, so a dependency bump that moves them fails loudly."""
    assert estimator_import_path("statevector") == "qiskit.primitives.StatevectorEstimator"
    assert estimator_import_path("aer") == "qiskit_aer.primitives.Estimator"
    assert estimator_import_path("aer_v2") == "qiskit_aer.primitives.EstimatorV2"
    assert set(ESTIMATOR_BACKENDS) == {"statevector", "aer", "aer_v2"}
    assert DEFAULT_BACKEND == "statevector"

    with pytest.raises(ValueError):
        estimator_import_path("does-not-exist")


def test_estimator_classes_resolve_and_match_their_probed_names() -> None:
    """Every backend must import, and the class name must match the attribute name."""
    for backend, (module_name, attribute) in ESTIMATOR_BACKENDS.items():
        module = importlib.import_module(module_name)
        estimator_class = getattr(module, attribute)
        assert estimator_class.__name__ == attribute, f"{backend} resolved to a different class"
        instance = build_estimator(backend)
        assert type(instance).__name__ == attribute


def test_qiskit_primitives_has_no_estimatorv2_symbol() -> None:
    """Documents the trap: the V2 statevector estimator is StatevectorEstimator."""
    primitives = importlib.import_module("qiskit.primitives")
    assert hasattr(primitives, "StatevectorEstimator")
    assert not hasattr(primitives, "EstimatorV2")


def test_estimator_versions_are_v2_except_aer_estimator() -> None:
    """The adapter branches on BaseEstimatorV2, so the versions must be as probed."""
    from qiskit.primitives import BaseEstimatorV1, BaseEstimatorV2

    assert isinstance(build_estimator("statevector"), BaseEstimatorV2)
    assert isinstance(build_estimator("aer_v2"), BaseEstimatorV2)
    assert isinstance(build_estimator("aer"), BaseEstimatorV1)


def test_build_estimator_rejects_an_unknown_backend() -> None:
    with pytest.raises(ValueError):
        build_estimator("nope")


def test_make_ansatz_has_one_parameter_per_rotation() -> None:
    circuit, parameters = make_ansatz(2, layers=2)
    assert len(parameters) == 4
    assert circuit.num_qubits == 2
    with pytest.raises(ValueError):
        make_ansatz(0, layers=1)
    with pytest.raises(ValueError):
        make_ansatz(2, layers=0)


def test_default_layers_gives_entangled_fixtures_enough_depth() -> None:
    assert default_layers(1) == 1
    assert default_layers(2) == 2


def test_build_observable_rejects_an_empty_sum() -> None:
    with pytest.raises(ValueError):
        build_observable({})


def test_evaluate_energy_endpoints_are_exact() -> None:
    """<Z> must be +1 at theta=0 and -1 at theta=pi for Ry(theta)|0>."""
    estimator = build_estimator()
    circuit, parameters = make_ansatz(1, layers=1)
    observable = build_observable({"Z": 1.0})

    assert evaluate_energy(estimator, circuit, observable, [0.0]) == pytest.approx(1.0, abs=1e-12)
    assert evaluate_energy(estimator, circuit, observable, [math.pi]) == pytest.approx(
        -1.0, abs=1e-12
    )


def test_v2_estimators_agree_on_one_evaluation() -> None:
    """Both V2 backends are exact, so they must agree to floating-point precision."""
    circuit, _ = make_ansatz(2, layers=2)
    observable = build_observable({"ZZ": 1.0, "IX": 1.0, "XI": 1.0})
    angles = [0.3, -0.7, 1.1, 0.2]

    statevector_energy = evaluate_energy(
        build_estimator("statevector"), circuit, observable, angles
    )
    aer_v2_energy = evaluate_energy(build_estimator("aer_v2"), circuit, observable, angles)

    assert statevector_energy == pytest.approx(aer_v2_energy, abs=1e-9)


def test_v1_aer_adapter_uses_the_v1_call_form() -> None:
    """Exercises the V1 positional branch against an analytic expectation.

    ``<Z>`` after ``Ry(theta)|0>`` is exactly ``cos(theta)``, so this needs no
    diagonaliser. The loose band is the recorded shot noise of the V1 estimator,
    not slack in the algorithm: it checks the adapter wiring, and the 1e-6
    accuracy claim is made against the exact estimator elsewhere. The sampler is
    seeded so the band is crossed deterministically or not at all.
    """
    estimator = build_estimator("aer", seed=V1_SAMPLER_SEED)
    circuit, _ = make_ansatz(1, layers=1)
    observable = build_observable({"Z": 1.0})
    theta = math.pi / 3.0

    energy = evaluate_energy(estimator, circuit, observable, [theta])
    assert isinstance(energy, float)
    assert energy == pytest.approx(math.cos(theta), abs=V1_SHOT_NOISE_BAND)

    # Determinism, not luck: the same seeded sampler must reproduce the value on
    # every further evaluation. Without the seed this list is not constant, which is
    # exactly why this test used to fail occasionally.
    repeats = [evaluate_energy(estimator, circuit, observable, [theta]) for _ in range(3)]
    assert repeats == [energy] * 3, f"the seeded V1 sampler is not reproducible: {repeats}"

    # The V1 primitive rejects the V2 tuple form, which proves the adapter really
    # took the V1 branch instead of one call shape happening to satisfy both.
    with pytest.raises(TypeError):
        estimator.run([(circuit, observable, [0.0])])


def test_seeding_an_exact_estimator_is_refused_not_ignored() -> None:
    """A seed that cannot take effect must never look like it did."""
    for backend in ("statevector", "aer_v2"):
        with pytest.raises(ValueError, match="no sampler seed"):
            build_estimator(backend, seed=V1_SAMPLER_SEED)


# ---------------------------------------------------------------------------
# VQE: the two required properties
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("label,terms,expected", FIXTURE_CASES, ids=FIXTURE_IDS)
def test_vqe_reaches_exact_within_1e6(
    label: str, terms: dict[str, float], expected: float
) -> None:
    """VQE must land within 1e-6 of the exact ground energy on a noiseless fixture."""
    result = vqe_ground_energy(terms)
    assert result["abs_error_ha"] <= TARGET_ACCURACY_HA, (
        f"{label}: |VQE - exact| = {result['abs_error_ha']:.3e} exceeds "
        f"{TARGET_ACCURACY_HA:g} (VQE {result['vqe_energy_ha']}, exact {expected})"
    )
    assert result["converged_within_target"] is True
    assert result["vqe_energy_ha"] == pytest.approx(expected, abs=TARGET_ACCURACY_HA)


@pytest.mark.parametrize("label,terms,expected", FIXTURE_CASES, ids=FIXTURE_IDS)
def test_vqe_never_below_exact_beyond_epsilon(
    label: str, terms: dict[str, float], expected: float
) -> None:
    """A variational energy is an upper bound; going below it is a bug, not noise."""
    result = vqe_ground_energy(terms)
    assert result["below_exact_beyond_epsilon"] is False, (
        f"{label}: VQE {result['vqe_energy_ha']} is below exact {result['exact_energy_ha']} "
        f"by more than the documented epsilon {result['exact_epsilon_ha']:g}"
    )
    assert result["vqe_energy_ha"] >= result["exact_energy_ha"] - EXACT_EPSILON_HA
    assert result["exact_epsilon_ha"] == EXACT_EPSILON_HA


def test_vqe_reports_tolerances_and_provenance() -> None:
    terms = {"ZZ": 1.0, "IX": 1.0}
    result = vqe_ground_energy(terms)

    for key in (
        "exact_energy_ha",
        "vqe_energy_ha",
        "abs_error_ha",
        "error_mha",
        "target_accuracy_ha",
        "converged_within_target",
        "exact_epsilon_ha",
        "below_exact_beyond_epsilon",
        "n_qubits",
        "layers",
        "n_parameters",
        "optimal_parameters",
        "restarts",
        "energy_evaluations",
        "optimizer_success",
        "estimator_backend",
        "estimator_class_name",
        "estimator_module",
        "ansatz",
        "hamiltonian_terms",
        "note",
    ):
        assert key in result, f"missing {key}"

    assert result["n_qubits"] == 2
    assert result["n_parameters"] == len(result["optimal_parameters"])
    assert result["energy_evaluations"] > 0
    assert result["estimator_class_name"] == "StatevectorEstimator"
    assert result["estimator_module"] == "qiskit.primitives.statevector_estimator"
    assert result["hamiltonian_terms"] == terms
    assert isinstance(result["note"], str)


def test_vqe_is_deterministic_for_a_fixed_seed() -> None:
    terms = {"ZZ": 1.0, "IX": 1.0, "XI": 1.0}
    first = vqe_ground_energy(terms, seed=3)
    second = vqe_ground_energy(terms, seed=3)
    assert first["vqe_energy_ha"] == pytest.approx(second["vqe_energy_ha"], abs=1e-12)
    assert first["optimal_parameters"] == pytest.approx(second["optimal_parameters"], abs=1e-12)

    different = vqe_ground_energy(terms, seed=99)
    assert different["abs_error_ha"] <= TARGET_ACCURACY_HA


def test_vqe_bounds_and_arguments_are_validated() -> None:
    with pytest.raises(ValueError):
        vqe_ground_energy({})
    with pytest.raises(ValueError, match="same length"):
        vqe_ground_energy({"Z": 1.0, "ZZ": 2.0})
    with pytest.raises(ValueError, match="n_qubits"):
        vqe_ground_energy({"ZZ": 1.0}, n_qubits=3)
    with pytest.raises(ValueError, match="restarts"):
        vqe_ground_energy({"Z": 1.0}, restarts=0)


def test_vqe_optimum_is_a_lower_bound_on_the_fixture() -> None:
    """The variational energy must never sit below the true minimum."""
    for label, terms, expected in FIXTURES:
        result = vqe_ground_energy(terms)
        assert result["vqe_energy_ha"] >= expected - EXACT_EPSILON_HA, label


# -------------------------------------------------------------------------------
# app.core must never load at API import time. Lazy imports inside the service
# functions are allowed (that is how the endpoints reach the exact solver and the
# VQE loop); an import at module level -- including inside a top-level if/try --
# is what both tests below forbid.
# ---------------------------------------------------------------------------
def test_app_core_is_not_imported_eagerly_by_routers_or_services() -> None:
    """Static check: no API module may import app.core at module import time.

    Walks the AST instead of matching raw lines, so an import nested inside a
    function (a lazy, per-request import) passes while the same import at module
    level fails. Module-level ``if``/``try`` blocks still count as eager: their
    code runs when the module is imported.
    """
    offenders: list[str] = []
    for directory in ("app/routers", "app/services"):
        for path in sorted((REPO_ROOT / directory).glob("*.py")):
            tree = ast.parse(path.read_text(encoding="utf-8"))
            parents = {
                child: parent
                for parent in ast.walk(tree)
                for child in ast.iter_child_nodes(parent)
            }
            for node in ast.walk(tree):
                if not isinstance(node, (ast.Import, ast.ImportFrom)):
                    continue
                if isinstance(node, ast.Import):
                    modules = [alias.name for alias in node.names]
                else:
                    modules = [node.module or ""]
                if not any(m == "app.core" or m.startswith("app.core.") for m in modules):
                    continue
                # Find the nearest enclosing scope: a function body means the
                # import runs per call (lazy, allowed); anything else -- module,
                # class body, top-level if/try -- runs at import time (forbidden).
                ancestor = parents.get(node)
                while ancestor is not None and not isinstance(
                    ancestor, (ast.FunctionDef, ast.AsyncFunctionDef, ast.Lambda)
                ):
                    ancestor = parents.get(ancestor)
                if ancestor is None:
                    offenders.append(
                        f"{path.relative_to(REPO_ROOT)}:{node.lineno}: {ast.unparse(node)}"
                    )
    assert not offenders, (
        "API modules must not import app.core at module import time (lazy imports "
        f"inside functions are fine): {offenders}"
    )


def test_app_core_is_absent_from_sys_modules_after_importing_the_app() -> None:
    """Dynamic check in a fresh interpreter: importing the app does not load app.core."""
    code = (
        "import sys\n"
        "import app.main  # noqa: F401\n"
        "leaked = sorted(m for m in sys.modules if m.startswith('app.core'))\n"
        "assert not leaked, f'app.core loaded by the app: {leaked}'\n"
        "print('clean')\n"
    )
    completed = subprocess.run(
        [sys.executable, "-c", code],
        cwd=str(REPO_ROOT),
        capture_output=True,
        text=True,
        timeout=180,
    )
    assert completed.returncode == 0, (
        f"Fresh-interpreter check failed:\nstdout: {completed.stdout}\nstderr: {completed.stderr}"
    )
    assert "clean" in completed.stdout


# ---------------------------------------------------------------------------
# Pre-existing artifact loader and resource estimator (unchanged behaviour)
# ---------------------------------------------------------------------------
def _synthetic_two_qubit_matrix() -> np.ndarray:
    return np.array(
        [
            [-1.10, 0.00, 0.00, 0.08],
            [0.00, -0.50, 0.02, 0.00],
            [0.00, 0.02, -0.50, 0.00],
            [0.08, 0.00, 0.00, -0.20],
        ],
        dtype=np.float64,
    )


def test_matrix_loader_roundtrips_a_json_payload() -> None:
    payload = {"n_qubits": 2, "matrix": _synthetic_two_qubit_matrix().tolist()}
    matrix = h2_sto3g_hamiltonian_from_json(payload)
    assert matrix.shape == (4, 4)
    assert np.allclose(matrix.real, _synthetic_two_qubit_matrix())


def test_matrix_loader_rejects_a_shape_mismatch() -> None:
    with pytest.raises(ValueError):
        h2_sto3g_hamiltonian_from_json({"n_qubits": 2, "matrix": [[1.0, 0.0], [0.0, 1.0]]})


def test_h2_resources_estimate() -> None:
    resources = estimate_h2_resources(_synthetic_two_qubit_matrix())
    assert resources["n_qubits"] == 2
    assert resources["hamiltonian_dim"] == 4
    assert resources["estimated_circuit_gates"] >= 0
    assert resources["estimated_circuit_depth"] >= 0
