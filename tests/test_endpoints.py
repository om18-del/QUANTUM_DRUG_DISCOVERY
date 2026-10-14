from __future__ import annotations

import json
from types import SimpleNamespace

import pytest

from fastapi.testclient import TestClient

from app.main import app
from app.services.molecule_simulation import METHOD_CODES
from app.services.quantum_service import H2_HAMILTONIAN_ARTIFACT, h2_artifact_path
from tests.availability import (
    require_live_network,
    require_rdkit_descriptors,
)

client = TestClient(app)


def _assert_error_shape(response, expected_code: str, expected_status: int = 422) -> dict:
    assert response.status_code == expected_status, response.text
    data = response.json()
    assert "error" in data
    assert data["error"]["code"] == expected_code
    assert "message" in data["error"]
    return data


def test_health() -> None:
    response = client.get("/health")
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "ok"
    assert data["app"] == "Quantum Drug Discovery Explorer API"


def _point_quantum_service_at(tmp_path, monkeypatch) -> None:
    """Point the quantum service at a directory we control.

    Both branches below are then deterministic, whether or not a real
    h2_sto3g.json happens to exist in the working tree.
    """
    monkeypatch.setattr(
        "app.services.quantum_service.get_settings",
        lambda: SimpleNamespace(h2_hamiltonian_dir=tmp_path),
    )


def test_h2_curve_reports_unavailable_when_artifact_is_missing(tmp_path, monkeypatch) -> None:
    """No fabricated matrix may be served, so the curve is honestly unavailable."""
    _point_quantum_service_at(tmp_path, monkeypatch)

    response = client.get("/api/quantum/h2-curve")

    data = _assert_error_shape(response, "hamiltonian_unavailable", 503)
    message = data["error"]["message"]
    assert H2_HAMILTONIAN_ARTIFACT in message
    assert "not been generated" in message
    assert "Colab" in message
    # This endpoint really does draw a curve, so the word belongs here -- which
    # is exactly why the shared sentence used to leak it into the simulate panel.
    assert "no H2 curve is available" in message
    assert "curve" in message.lower()


def _require_checked_in_artifact() -> None:
    """The computed tests read the real artifact; skip, never fake, without it."""
    if not h2_artifact_path().is_file():
        pytest.skip(f"{H2_HAMILTONIAN_ARTIFACT} is not present in this checkout")


def test_h2_curve_is_computed_from_the_checked_in_artifact() -> None:
    """The curve must be the stored 58-point grid, not a placeholder."""
    _require_checked_in_artifact()

    response = client.get("/api/quantum/h2-curve")

    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "computed"
    assert data["artifact"] == H2_HAMILTONIAN_ARTIFACT
    assert data["units"] == "hartree"

    # Documented provenance, checked against what the artifact actually declares.
    metadata = data["metadata"]
    assert data["basis"] == metadata["basis"] == "STO-3G"
    assert metadata["mapping"] == "parity"
    assert metadata["two_qubit_reduction"] is True
    assert metadata["mapping_note"] == "parity mapping with two-qubit reduction"
    assert metadata["n_qubits"] == 2
    assert metadata["source"] == "precomputed offline with PySCF + qiskit-nature"

    lengths = data["bond_lengths_angstrom"]
    hf = data["hf_energies_ha"]
    exact = data["exact_energies_ha"]
    vqe = data["vqe_energies_ha"]
    nuclear = data["nuclear_repulsion_energies_ha"]
    errors = data["vqe_error_mha"]
    assert data["point_count"] == len(lengths) == 58
    assert len(hf) == len(exact) == len(vqe) == len(nuclear) == len(errors) == 58
    assert lengths == sorted(lengths), "the grid must ascend for the canvas x-axis"
    assert all(n > 0.0 for n in nuclear), "nuclear repulsion must be positive"
    assert data["convention"] and "nuclear_repulsion_ha" in data["convention"]

    # The same data point by point: for every bond length the endpoint answers
    # HF, FCI and the nuclear repulsion individually, not only as zipped arrays.
    points = data["points"]
    assert len(points) == 58
    for i, point in enumerate(points):
        assert point["bond_length_angstrom"] == lengths[i]
        assert point["hf_energy_ha"] == hf[i]
        assert point["fci_energy_ha"] == exact[i]
        assert point["nuclear_repulsion_ha"] == nuclear[i]
        assert point["vqe_energy_ha"] == vqe[i]

    # The UI default sits on the grid and holds the published equilibrium values.
    index = lengths.index(0.74)
    assert exact[index] == pytest.approx(-1.1372838, abs=1e-6)
    assert vqe[index] == pytest.approx(exact[index], abs=1e-6)
    # Stored HF and FCI are TOTAL energies (electronic plus nuclear repulsion),
    # so the known equilibrium references must match without any adjustment.
    assert hf[index] == pytest.approx(-1.1167593, abs=1e-6)
    assert nuclear[index] == pytest.approx(0.7151043, abs=1e-6)

    for bond, h, e, v, d in zip(lengths, hf, exact, vqe, errors):
        assert h >= e - 1e-9, f"HF reference below FCI at {bond} A"
        assert v >= e - 1e-9, f"variational energy below FCI at {bond} A"
        assert abs(d) <= 1.6, f"{abs(d):.6f} mHa exceeds chemical accuracy at {bond} A"

    assert H2_HAMILTONIAN_ARTIFACT in data["note"]
    assert "No optimisation, shots or noise" in data["note"]


def test_h2_curve_reports_a_malformed_point_as_a_fault(tmp_path, monkeypatch) -> None:
    """A point that is present but wrong is a fault (500), never an empty curve."""
    (tmp_path / H2_HAMILTONIAN_ARTIFACT).write_text(
        json.dumps({"molecule": "H2", "points": [{"bond_length_angstrom": 0.74}]}),
        encoding="utf-8",
    )
    _point_quantum_service_at(tmp_path, monkeypatch)

    response = client.get("/api/quantum/h2-curve")

    data = _assert_error_shape(response, "simulation_error", 500)
    assert "H2 curve computation failed" in data["error"]["message"]
    assert "point 0 is missing" in data["error"]["message"]


def test_h2_curve_with_an_empty_grid_reports_unavailable(tmp_path, monkeypatch) -> None:
    """A file with no stored points has no point to serve: the expected 503."""
    (tmp_path / H2_HAMILTONIAN_ARTIFACT).write_text(
        json.dumps({"molecule": "H2", "point_count": 0, "points": []}),
        encoding="utf-8",
    )
    _point_quantum_service_at(tmp_path, monkeypatch)

    response = client.get("/api/quantum/h2-curve")

    data = _assert_error_shape(response, "hamiltonian_unavailable", 503)
    message = data["error"]["message"]
    assert H2_HAMILTONIAN_ARTIFACT in message
    assert "no stored bond lengths" in message
    assert "no H2 curve is available" in message


def test_hadamard_demo() -> None:
    response = client.get("/api/quantum/hadamard-demo")
    assert response.status_code == 200
    data = response.json()
    assert data["circuit"] == "H|0>"
    assert data["probabilities"]["0"] == pytest.approx(0.5, abs=1e-9)
    assert data["probabilities"]["1"] == pytest.approx(0.5, abs=1e-9)


def test_simulate_missing_body() -> None:
    response = client.post("/api/quantum/simulate", content="not json")
    _assert_error_shape(response, "invalid_body", 400)


def test_simulate_invalid_field() -> None:
    response = client.post(
        "/api/quantum/simulate",
        content=json.dumps({"bond_length": "bad"}),
    )
    _assert_error_shape(response, "invalid_bond_length", 422)


def test_simulate_h2_reports_unavailable_without_the_artifact(tmp_path, monkeypatch) -> None:
    """Without the offline artifact there is no honest energy to report."""
    _point_quantum_service_at(tmp_path, monkeypatch)

    response = client.post(
        "/api/quantum/simulate",
        content=json.dumps({"bond_length": 0.74, "shots": 1024, "noise": 0.0}),
    )

    data = _assert_error_shape(response, "hamiltonian_unavailable", 503)
    message = data["error"]["message"]
    assert H2_HAMILTONIAN_ARTIFACT in message
    # The simulate panel reports one energy point. It must name that energy, not
    # the curve that a different endpoint draws.
    assert "no H2 energy is available" in message
    assert "curve" not in message.lower(), message


def _post_simulate(
    bond_length: float,
    shots: int | None = None,
    noise: float = 0.0,
    seed: int | None = 7,
    maxiter: int | None = None,
) -> dict:
    """POST /api/quantum/simulate; default is the noiseless statevector path."""
    payload = {"bond_length": bond_length, "shots": shots, "noise": noise, "seed": seed}
    if maxiter is not None:
        payload["maxiter"] = maxiter
    response = client.post("/api/quantum/simulate", content=json.dumps(payload))
    assert response.status_code == 200, response.text
    return response.json()


def test_simulate_h2_runs_a_real_vqe_at_the_stored_default_bond_length() -> None:
    """A fixed-seed noiseless run lands on -1.137 Ha within chemical accuracy."""
    _require_checked_in_artifact()

    data = _post_simulate(0.74)

    assert data["status"] == "computed"
    assert data["molecule"] == "H2"
    assert data["units"] == "hartree"
    # Snap-and-disclose: an on-grid answer says so under both field spellings
    # (SPEC.md section 6 names the *_angstrom trio; the task names the other).
    assert data["bond_length"] == pytest.approx(0.74, abs=1e-9)
    assert data["requested_bond_length"] == pytest.approx(0.74, abs=1e-12)
    assert data["used_bond_length"] == pytest.approx(0.74, abs=1e-9)
    assert data["snapped"] is False
    assert data["bond_length_requested_angstrom"] == pytest.approx(0.74, abs=1e-12)
    assert data["bond_length_used_angstrom"] == pytest.approx(0.74, abs=1e-9)
    assert data["bond_length_is_exact"] is True

    # The run's result, the exact reference, and the artifact's stored numbers.
    assert data["exact_energy_ha"] == pytest.approx(-1.1372838, abs=1e-6)
    assert data["vqe_energy_ha"] == pytest.approx(-1.1372838, abs=1.6e-3)
    assert data["vqe_energy_ha"] >= data["exact_energy_ha"] - 1e-9
    # The exact reference is a live NumPy diagonalisation of the same Pauli
    # terms; it must land on the artifact's stored FCI value.
    assert data["exact_energy_ha"] == pytest.approx(data["fci_energy_ha"], abs=1e-8)
    assert data["hf_energy_ha"] == pytest.approx(-1.1167593, abs=1e-6)
    assert data["nuclear_repulsion_ha"] == pytest.approx(0.7151043, abs=1e-6)
    assert abs(data["error_mha"]) <= 1.6
    assert data["error_vs_exact_ha"] == pytest.approx(
        data["vqe_energy_ha"] - data["exact_energy_ha"], abs=1e-15
    )
    assert data["abs_error_vs_exact_ha"] == pytest.approx(
        abs(data["error_vs_exact_ha"]), abs=1e-15
    )
    assert data["within_chemical_accuracy"] is True
    assert data["chemical_accuracy_threshold_ha"] == pytest.approx(1.6e-3)
    assert data["noiseless"] is True

    # A real optimiser run, not a stored number.
    assert data["iterations"] > 1
    assert len(data["convergence_history_ha"]) == data["iterations"]
    assert data["seed"] == 7
    assert data["maxiter"] >= 1
    assert data["restarts"] >= 1
    # The budget is spelled out, and the backward-compatible `iterations` equals
    # the total across all restarts.
    assert data["n_restarts"] == data["restarts"] == 3
    assert data["maxiter_per_restart"] == data["maxiter"]
    assert data["total_energy_evaluations"] == data["iterations"]
    # A noiseless run has no sampling uncertainty to report.
    assert data["vqe_energy_std_ha"] is None
    assert data["vqe_energy_uncertainty_method"] is None
    assert data["shots"] is None
    assert data["noise"] == 0.0
    assert data["estimator_backend"] == "statevector"
    assert data["estimator_class_name"] == "StatevectorEstimator"
    assert data["optimizer"] == "COBYLA"
    assert "Ry" in data["ansatz"]
    assert data["classical_simulation"] is True
    # The response must say what the run was, and claim no quantum advantage.
    assert "Classical simulation" in data["note"]
    assert "no quantum advantage is claimed" in data["note"]
    # The rule behind the shared sentence still holds: this endpoint reports an
    # energy, and must not describe itself in terms of the curve.
    assert "curve" not in data["note"].lower(), data["note"]


def test_simulate_h2_snaps_to_the_nearest_stored_point_and_discloses_it() -> None:
    """0.741 A is off the grid; the answer comes from 0.74 A and says so."""
    _require_checked_in_artifact()

    snapped = _post_simulate(0.741)
    on_grid = _post_simulate(0.74)

    assert snapped["status"] == "computed"
    # Both spellings of the disclosure, as specified by the task and SPEC sec. 6.
    assert snapped["requested_bond_length"] == pytest.approx(0.741, abs=1e-12)
    assert snapped["used_bond_length"] == pytest.approx(0.74, abs=1e-9)
    assert snapped["snapped"] is True
    assert snapped["bond_length_requested_angstrom"] == pytest.approx(0.741, abs=1e-12)
    assert snapped["bond_length_used_angstrom"] == pytest.approx(0.74, abs=1e-9)
    assert snapped["bond_length"] == pytest.approx(0.74, abs=1e-9)
    assert snapped["bond_length_is_exact"] is False
    assert "not on the stored grid" in snapped["note"]
    assert "0.74" in snapped["note"]
    # Same stored point and same seed, so the same run: nothing was interpolated.
    assert snapped["exact_energy_ha"] == pytest.approx(on_grid["exact_energy_ha"], abs=1e-12)
    assert snapped["vqe_energy_ha"] == pytest.approx(on_grid["vqe_energy_ha"], abs=1e-12)
    assert snapped["convergence_history_ha"] == on_grid["convergence_history_ha"]


def test_simulate_h2_tie_between_two_stored_points_prefers_the_shorter_bond_length() -> None:
    """0.745 A sits between 0.74 and 0.75 A; the documented rule picks 0.74."""
    _require_checked_in_artifact()

    data = _post_simulate(0.745)

    assert data["bond_length_used_angstrom"] == pytest.approx(0.74, abs=1e-9)
    assert data["bond_length_is_exact"] is False
    assert "shorter bond length" in data["note"]


def test_simulate_h2_reports_a_malformed_artifact(tmp_path, monkeypatch) -> None:
    """A present-but-unusable file is a fault (500), never placeholder energies."""
    (tmp_path / H2_HAMILTONIAN_ARTIFACT).write_text(
        json.dumps({"molecule": "H2", "points": [{"bond_length_angstrom": 0.74}]}),
        encoding="utf-8",
    )
    _point_quantum_service_at(tmp_path, monkeypatch)

    response = client.post(
        "/api/quantum/simulate",
        content=json.dumps({"bond_length": 0.74, "noise": 0.0}),
    )

    data = _assert_error_shape(response, "simulation_error", 500)
    assert "H2 simulation failed" in data["error"]["message"]
    assert "point 0 is missing" in data["error"]["message"]


def test_simulate_with_an_empty_grid_reports_unavailable(tmp_path, monkeypatch) -> None:
    """No stored points means no point can answer: the expected 503."""
    (tmp_path / H2_HAMILTONIAN_ARTIFACT).write_text(
        json.dumps({"molecule": "H2", "points": []}), encoding="utf-8"
    )
    _point_quantum_service_at(tmp_path, monkeypatch)

    response = client.post(
        "/api/quantum/simulate",
        content=json.dumps({"bond_length": 0.74}),
    )

    data = _assert_error_shape(response, "hamiltonian_unavailable", 503)
    message = data["error"]["message"]
    assert H2_HAMILTONIAN_ARTIFACT in message
    assert "no stored bond lengths" in message
    assert "no H2 energy is available" in message


def test_simulate_refuses_a_result_below_the_exact_energy(tmp_path, monkeypatch) -> None:
    """Variational guard: an energy below exact is a bug, so it is refused.

    ``_run_vqe`` is replaced with a runner that reports 50 mHa below the true
    minimum; the endpoint must answer an error envelope and ship no result field.
    """
    _require_checked_in_artifact()

    def broken_run(pauli_terms, n_qubits, *, seed, maxiter, shots):
        from app.core.exact import pauli_sum_ground_energy

        return {
            "vqe_energy_ha": pauli_sum_ground_energy(pauli_terms, n_qubits) - 0.05,
            "iterations": 1,
            "convergence_history_ha": [-999.0],
            "estimator_backend": "statevector",
            "estimator_class_name": "StatevectorEstimator",
            "ansatz": "1 x [Ry on each of 2 qubit(s), then a CX ladder]",
            "optimizer": "COBYLA",
            "restarts": 1,
            "maxiter": maxiter or 1,
            "seed": seed,
            "seed_simulator_applied": False,
            "optimal_parameters": [],
        }

    monkeypatch.setattr("app.services.quantum_service._run_vqe", broken_run)

    response = client.post(
        "/api/quantum/simulate",
        content=json.dumps({"bond_length": 0.74}),
    )

    data = _assert_error_shape(response, "variational_principle_violation", 500)
    assert "variational principle" in data["error"]["message"]
    assert set(data) == {"error"}, "a refused run must not ship any result fields"


def _shot_run_stub(below_exact_by: float):
    """A shot-based runner reporting an estimate this far below the exact energy."""

    def runner(pauli_terms, n_qubits, *, seed, maxiter, shots):
        from app.core.exact import pauli_sum_ground_energy

        return {
            "vqe_energy_ha": pauli_sum_ground_energy(pauli_terms, n_qubits) - below_exact_by,
            "vqe_energy_std_ha": 0.01,
            "vqe_energy_uncertainty_method": "estimator_variance_metadata",
            "iterations": 1,
            "total_energy_evaluations": 1,
            "n_restarts": 1,
            "maxiter_per_restart": 1,
            "convergence_history_ha": [-999.0],
            "estimator_backend": "aer",
            "estimator_class_name": "Estimator",
            "ansatz": "1 x [Ry on each of 2 qubit(s), then a CX ladder]",
            "optimizer": "COBYLA",
            "restarts": 1,
            "maxiter": 1,
            "seed": seed,
            "seed_simulator_applied": True,
            "optimal_parameters": [],
        }

    return runner


def test_simulate_shot_run_below_exact_is_not_refused(monkeypatch) -> None:
    """A shot estimate may sit below exact on sampling noise: served, not 500.

    The variational guard is scoped to noiseless statevector runs (SPEC.md section
    6), so a shot run that lands below exact reports its uncertainty instead of
    being refused.
    """
    _require_checked_in_artifact()
    monkeypatch.setattr(
        "app.services.quantum_service._run_vqe", _shot_run_stub(below_exact_by=0.05)
    )

    response = client.post(
        "/api/quantum/simulate",
        content=json.dumps({"bond_length": 0.74, "shots": 1024}),
    )

    assert response.status_code == 200, response.text
    data = response.json()
    assert data["status"] == "computed"
    assert data["vqe_energy_ha"] < data["exact_energy_ha"], "the stub is below exact"
    assert data["vqe_energy_std_ha"] == pytest.approx(0.01)
    assert data["vqe_energy_uncertainty_method"] == "estimator_variance_metadata"
    assert data["within_chemical_accuracy"] is None


def test_shot_run_reports_a_fresh_reevaluation_not_the_trace_minimum(monkeypatch) -> None:
    """The headline shot energy is a fresh evaluation, not the biased trace minimum."""
    import app.core.vqe as vqe_core

    from app.services.quantum_service import _run_vqe

    def biased_objective(estimator, circuit, observable, parameters, *, run_options=None):
        return -99.0

    def fresh_estimate(estimator, circuit, observable, parameters, *, run_options=None):
        return {"value": -1.5, "std": 0.02}

    monkeypatch.setattr(vqe_core, "evaluate_energy", biased_objective)
    monkeypatch.setattr(vqe_core, "evaluate_energy_with_uncertainty", fresh_estimate)

    run = _run_vqe({"ZZ": 1.0}, 2, seed=1, maxiter=3, shots=1024)

    assert min(run["convergence_history_ha"]) == pytest.approx(-99.0)
    assert run["vqe_energy_ha"] == pytest.approx(-1.5), "not a fresh re-evaluation"
    assert run["vqe_energy_std_ha"] == pytest.approx(0.02)
    assert run["vqe_energy_uncertainty_method"] == "estimator_variance_metadata"


def test_shot_run_falls_back_to_seed_reevaluation_without_metadata(monkeypatch) -> None:
    """With no variance metadata, 10 re-evaluations with fresh seeds give mean and std."""
    import app.core.vqe as vqe_core

    from app.services.quantum_service import _run_vqe

    real = vqe_core.evaluate_energy_with_uncertainty

    def strip_std(estimator, circuit, observable, parameters, *, run_options=None):
        result = dict(real(estimator, circuit, observable, parameters, run_options=run_options))
        result["std"] = None
        return result

    monkeypatch.setattr(vqe_core, "evaluate_energy_with_uncertainty", strip_std)

    run = _run_vqe({"ZZ": 1.0, "IX": 1.0}, 2, seed=3, maxiter=5, shots=512)

    assert run["vqe_energy_uncertainty_method"] == "seed_reevaluation_10x"
    assert run["vqe_energy_std_ha"] is not None
    assert run["vqe_energy_std_ha"] >= 0.0


def test_simulate_is_reproducible_with_a_fixed_seed() -> None:
    """Same seed, same run: byte-identical convergence history."""
    _require_checked_in_artifact()

    first = _post_simulate(0.74, seed=7)
    second = _post_simulate(0.74, seed=7)

    assert first["iterations"] == second["iterations"]
    assert first["convergence_history_ha"] == second["convergence_history_ha"]
    assert first["vqe_energy_ha"] == second["vqe_energy_ha"]


def test_simulate_runs_without_a_seed() -> None:
    """The seed is optional; a null seed still produces a complete run."""
    _require_checked_in_artifact()

    response = client.post(
        "/api/quantum/simulate",
        content=json.dumps({"bond_length": 0.74, "seed": None}),
    )

    assert response.status_code == 200, response.text
    data = response.json()
    assert data["status"] == "computed"
    assert data["seed"] is None
    assert data["noiseless"] is True
    assert data["iterations"] > 1
    assert data["vqe_energy_ha"] >= data["exact_energy_ha"] - 1e-9


def test_simulate_shot_based_run_discloses_sampling_not_a_noiseless_pass() -> None:
    """A shots run uses the Aer estimator and cannot claim a noiseless pass."""
    _require_checked_in_artifact()

    data = _post_simulate(0.74, shots=512, seed=7, maxiter=60)

    assert data["status"] == "computed"
    assert data["shots"] == 512
    assert data["noiseless"] is False
    assert data["estimator_backend"] == "aer"
    assert data["estimator_class_name"] == "Estimator"
    assert data["seed"] == 7
    assert data["seed_simulator_applied"] is True, "the seed must reach seed_simulator"
    # Chemical accuracy is a noiseless claim; a sampled run reports the number
    # and the threshold, but not a pass/fail flag.
    assert data["within_chemical_accuracy"] is None
    # A shot-based estimate is a statistic, so it is NOT held to the variational
    # floor; instead it reports its own sampling uncertainty and the method used.
    assert data["vqe_energy_std_ha"] is not None
    assert data["vqe_energy_std_ha"] >= 0.0
    assert data["vqe_energy_uncertainty_method"] in {
        "estimator_variance_metadata",
        "seed_reevaluation_10x",
    }
    assert data["n_restarts"] == data["restarts"]
    assert data["maxiter_per_restart"] == data["maxiter"]
    assert data["total_energy_evaluations"] == data["iterations"]
    assert len(data["convergence_history_ha"]) == data["iterations"]
    assert "shots per evaluation" in data["note"]
    assert "sampling error" in data["note"]


def test_simulate_rejects_bad_optional_fields() -> None:
    """seed, maxiter and shots are validated before any run starts."""
    cases = [
        ({"bond_length": 0.74, "seed": "abc"}, "invalid_seed"),
        ({"bond_length": 0.74, "seed": -1}, "invalid_seed"),
        ({"bond_length": 0.74, "maxiter": 0}, "invalid_maxiter"),
        ({"bond_length": 0.74, "maxiter": 2.5}, "invalid_maxiter"),
        ({"bond_length": 0.74, "shots": 0}, "invalid_shots"),
        ({"bond_length": 0.74, "shots": "many"}, "invalid_shots"),
    ]
    for payload, code in cases:
        response = client.post("/api/quantum/simulate", content=json.dumps(payload))
        _assert_error_shape(response, code, 422)


def test_simulate_accepts_bond_length_angstrom_alias() -> None:
    """The artifact-shaped input name is accepted alongside ``bond_length``."""
    _require_checked_in_artifact()

    response = client.post(
        "/api/quantum/simulate",
        content=json.dumps({"bond_length_angstrom": 0.74, "seed": 7}),
    )

    assert response.status_code == 200, response.text
    data = response.json()
    assert data["status"] == "computed"
    assert data["requested_bond_length"] == pytest.approx(0.74, abs=1e-12)
    assert data["used_bond_length"] == pytest.approx(0.74, abs=1e-9)


def test_simulate_unsupported_molecule_is_placeholder() -> None:
    response = client.post(
        "/api/quantum/simulate",
        content=json.dumps({"bond_length": 0.74, "shots": 1, "noise": 0.0, "molecule": "caffeine"}),
    )
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "placeholder"
    assert data["vqe_energy_ha"] is None


def test_resources_are_counted_from_the_formula_without_any_artifact(tmp_path, monkeypatch) -> None:
    """The estimate no longer depends on a Hamiltonian file being present.

    It used to answer H2 with ``status: "placeholder"`` and null qubit counts
    whenever the artifact was missing, because its numbers were derived from the
    file. It now counts from the formula, so pointing the service at an empty
    directory must not change the answer -- and it must never mention the curve.
    """
    _point_quantum_service_at(tmp_path, monkeypatch)

    response = client.get("/api/quantum/resources", params={"formula": "H2"})

    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "computed"
    assert data["n_spatial_orbitals"] == 2
    assert data["n_electrons"] == 2
    assert data["n_qubits_full"] == 4
    assert data["n_qubits_after_parity_reduction"] == 2
    assert data["statevector_memory_bytes"] == 64
    assert data["feasibility"]["verdict"] == "simulable_in_app"
    assert "curve" not in json.dumps(data).lower(), data


def test_resources_rejects_a_name_that_is_not_a_formula() -> None:
    """``caffeine`` is a name, not a formula, and the estimator takes formulas."""
    response = client.post(
        "/api/quantum/resources",
        content=json.dumps({"molecule": "caffeine"}),
    )
    data = _assert_error_shape(response, "invalid_formula", 422)
    assert "caffeine" in data["error"]["message"]
    assert "formula" in data["error"]["message"]


def test_resources_endpoint_returns_the_worked_examples() -> None:
    """The five documented examples, counted by the live endpoint.

    orbitals, electrons, full qubits, reduced qubits, verdict -- every value
    written out by hand here rather than read back from the response's own
    formulas.
    """
    expected = {
        "H2": (2, 2, 4, 2, "simulable_in_app"),
        "LiH": (6, 4, 12, 10, "precompute_only"),
        "H2O": (7, 10, 14, 12, "precompute_only"),
        "CH4": (9, 10, 18, 16, "precompute_only"),
        "C9H8O4": (73, 94, 146, 144, "beyond_classical_simulation"),
    }

    for formula, values in expected.items():
        response = client.get("/api/quantum/resources", params={"formula": formula})
        assert response.status_code == 200, response.text
        data = response.json()
        assert data["status"] == "computed"
        assert data["basis"] == "STO-3G"
        assert (
            data["n_spatial_orbitals"],
            data["n_electrons"],
            data["n_qubits_full"],
            data["n_qubits_after_parity_reduction"],
            data["feasibility"]["verdict"],
        ) == values, formula
        assert data["tags"]["n_qubits_full"] == "exact"


def test_resources_without_parameters_is_a_usage_response() -> None:
    """No parameters must not be answered with a guessed molecule."""
    response = client.get("/api/quantum/resources")

    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "usage"
    assert "formula" in data["message"]
    assert data["basis"] == "STO-3G"
    assert data["supported_elements"][:3] == ["H", "He", "Li"]
    assert set(data["supported_elements"]) == set(data["basis_functions_per_atom"])
    assert "n_qubits_full" not in data


def test_resources_post_alias_still_answers() -> None:
    """The body-shaped call an older client makes goes through the same path."""
    response = client.post("/api/quantum/resources", content=json.dumps({"molecule": "H2"}))

    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "computed"
    assert data["requested"]["molecule"] == "H2"
    assert data["n_qubits_full"] == 4

    empty = client.post("/api/quantum/resources", content=json.dumps({}))
    assert empty.status_code == 200
    assert empty.json()["status"] == "usage"


def test_resources_unsupported_element_is_a_422() -> None:
    response = client.get("/api/quantum/resources", params={"formula": "Fe2O3"})

    data = _assert_error_shape(response, "unsupported_element", 422)
    assert "Fe" in data["error"]["message"]
    assert "Ar" in data["error"]["message"]


def test_resources_invalid_formula_is_a_422() -> None:
    response = client.get("/api/quantum/resources", params={"formula": "c6h6"})
    _assert_error_shape(response, "invalid_formula", 422)


def test_resources_active_orbitals_bounds() -> None:
    ok = client.get("/api/quantum/resources", params={"formula": "H2O", "active_orbitals": 4})
    assert ok.status_code == 200
    assert ok.json()["active_space"]["n_qubits"] == 8
    assert ok.json()["tags"]["active_space.n_qubits"] == "estimated"

    too_big = client.get("/api/quantum/resources", params={"formula": "H2O", "active_orbitals": 99})
    _assert_error_shape(too_big, "invalid_active_orbitals", 422)

    orphan = client.get("/api/quantum/resources", params={"active_orbitals": 4})
    _assert_error_shape(orphan, "missing_formula", 422)


def test_resources_molecule_key_is_resolved_to_its_formula() -> None:
    response = client.get("/api/quantum/resources", params={"molecule": "lih"})

    assert response.status_code == 200
    data = response.json()
    assert data["requested"]["resolved_formula"] == "LiH"
    assert data["requested"]["resolved_from"] == "molecule"
    assert data["n_spatial_orbitals"] == 6


def test_resources_molecule_keys_match_the_registry() -> None:
    """A registry molecule must be nameable at the estimator without a formula."""
    from app.core.exact import MOLECULE_REGISTRY
    from app.services.quantum_service import MOLECULE_FORMULAS

    for key in MOLECULE_REGISTRY:
        assert key in MOLECULE_FORMULAS, key
        response = client.get("/api/quantum/resources", params={"molecule": key})
        assert response.status_code == 200, response.text
        assert response.json()["status"] == "computed"


def test_molecule_search_live_pubchem() -> None:
    require_live_network()


def test_molecule_search_unknown_name_returns_compound_not_found() -> None:
    """A name PubChem does not have is a 404 answer, not a retry and not a 200.

    Replaces ``test_molecule_search_returns_placeholder_when_pubchem_fails``:
    A1 separates "not found" from "not working", so the honest response is an
    error envelope with one clean, URL-free message.
    """
    query = "xyzzynonexistentsearchterm"
    response = client.get(f"/api/molecule/search?q={query}")

    assert response.status_code == 404
    data = response.json()
    assert set(data) == {"error"}, f"expected a single error envelope, got {set(data)}"
    assert data["error"]["code"] == "compound_not_found"
    assert data["error"]["message"] == (
        f"No compound found for '{query}'. Check the spelling."
    )
    message = data["error"]["message"].lower()
    for leak in ("http://", "https://", "pubchem", "traceback"):
        assert leak not in message


def test_molecule_properties_live_rdkit() -> None:
    """The endpoint must return real RDKit descriptors, not a placeholder."""
    require_rdkit_descriptors()
    response = client.post(
        "/api/molecule/properties",
        content=json.dumps({"smiles": "CC(=O)Oc1ccccc1C(=O)O"}),
    )
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "computed", data.get("note")
    assert data["formula"] == "C9H8O4"
    assert data["h_bond_donors"] == 1
    assert data["h_bond_acceptors"] == 4
    assert data["rotatable_bonds"] == 3
    assert data["lipinski_rule_of_five"]["violations"] == 0


def test_molecule_properties_returns_placeholder_when_rdkit_fails(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def _rdkit_unavailable(*_args: object, **_kwargs: object) -> None:
        raise RuntimeError("RDKit not available in this build")

    monkeypatch.setattr(
        "app.services.rdkit_helper.rdkit_properties_from_smiles", _rdkit_unavailable
    )
    response = client.post(
        "/api/molecule/properties",
        content=json.dumps({"smiles": "CC(=O)Oc1ccccc1C(=O)O"}),
    )
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "placeholder"
    assert data["smiles"] is not None
    assert data["note"] is not None
    assert "RDKit" in data["note"]


def test_molecule_properties_invalid_smiles_is_placeholder() -> None:
    response = client.post(
        "/api/molecule/properties",
        content=json.dumps({"smiles": "not a valid smiles"}),
    )
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "placeholder"
    assert data["note"]


def _pdb_atom(
    record: str,
    serial: int,
    name: str,
    res: str,
    chain: str,
    resseq: int,
    x: float,
    y: float,
    z: float,
    element: str,
) -> str:
    """Format one ATOM/HETATM line with the columns the PDB spec requires."""
    return (
        f"{record:<6}{serial:>5} {name:<4} {res:>3} {chain}{resseq:>4}    "
        f"{x:>8.3f}{y:>8.3f}{z:>8.3f}{1.0:>6.2f}{20.0:>6.2f}{'':>10}{element:>2}"
    )


# A physically plausible fragment: the backbone amide N of ASP A25 donates a
# hydrogen bond to a ligand oxygen placed in the direction the N-H bond points
# (the N -> CA direction is +x, so the acceptor sits on -x), and the ASP
# carboxylate oxygen is 3.11 A from a ligand nitrogen for the salt bridge. The
# geometry is what makes the classification unambiguous, and it is independent
# of any network round-trip or the machine's RDKit build.
PDB_FIXTURE = "\n".join([
    _pdb_atom("ATOM", 1, "N", "ASP", "A", 25, 0.000, 0.000, 0.000, "N"),
    _pdb_atom("ATOM", 2, "CA", "ASP", "A", 25, 1.460, 0.000, 0.000, "C"),
    _pdb_atom("ATOM", 3, "OD1", "ASP", "A", 25, 1.200, 1.400, 0.000, "O"),
    _pdb_atom("HETATM", 9, "O1", "MK1", "B", 400, -2.800, 0.200, 0.000, "O"),
    _pdb_atom("HETATM", 10, "N2", "MK1", "B", 400, 1.500, 4.500, 0.000, "N"),
    _pdb_atom("HETATM", 11, "C1", "MK1", "B", 400, 6.000, 6.000, 6.000, "C"),
    "END",
])


def test_complexes_list_reports_real_ligand_codes() -> None:
    response = client.get("/api/complexes/list")
    assert response.status_code == 200
    data = response.json()
    assert "complexes" in data
    assert data["count"] > 0
    codes = {item["id"]: item["ligand_code"] for item in data["complexes"]}
    # The manifest must carry the HETATM codes actually deposited in each entry.
    assert codes["hiv-protease-indinavir"] == "MK1"
    assert codes["ace-lisinopril"] == "LPR"
    # Codes verified against the deposited files for the viral additions.
    assert codes["hiv-rt-nevirapine"] == "NVP"
    assert codes["influenza-neuraminidase-oseltamivir"] == "G39"


def test_far_donor_acceptor_pair_is_a_near_miss_not_a_hydrogen_bond() -> None:
    """Ideal donor/acceptor geometry beyond 3.5 A must not be called a hydrogen bond.

    This is the nevirapine situation in 1VRT: N14...O Lys101 is 4.05 A at 164 deg.
    """
    from app.services.structure_analysis import HBOND_MAX, analyse_complex

    pdb = "\n".join(
        [
            _pdb_atom("ATOM", 1, "N", "ASP", "A", 25, 0.000, 0.000, 0.000, "N"),
            _pdb_atom("ATOM", 2, "CA", "ASP", "A", 25, 1.460, 0.000, 0.000, "C"),
            _pdb_atom("HETATM", 9, "O1", "MK1", "B", 400, -4.050, 0.100, 0.000, "O"),
            "END",
        ]
    )
    analysis = analyse_complex(pdb, "MK1")

    assert analysis["status"] == "computed"
    assert analysis["contact_summary"]["hydrogen_bond"] == 0

    misses = analysis["near_miss_hydrogen_bonds"]
    assert len(misses) == 1
    miss = misses[0]
    assert miss["protein_residue"] == "ASP A25"
    assert miss["ligand_atom"] == "O1"
    assert miss["distance_angstrom"] > HBOND_MAX
    assert miss["distance_angstrom"] == pytest.approx(4.05, abs=0.05)
    assert miss["angle_deg"] is not None and miss["angle_deg"] > 150.0
    assert "cutoff" in miss["reason"]


def test_contact_labels_mark_hydrogen_bonds_as_a_heuristic() -> None:
    from app.services.structure_analysis import CONTACT_LABELS

    assert CONTACT_LABELS["hydrogen_bond"] == "geometric heuristic"
    assert CONTACT_LABELS["salt_bridge"] == "geometric heuristic"
    assert CONTACT_LABELS["hydrophobic"] == "measured geometry"


def test_complexes_label_human_versus_viral_targets() -> None:
    """ACE must be presented as a human enzyme, never as a virus."""
    response = client.get("/api/complexes/list")
    data = response.json()
    by_id = {item["id"]: item for item in data["complexes"]}

    assert by_id["ace-lisinopril"]["target_origin"] == "human"
    assert "not a virus" in by_id["ace-lisinopril"]["display_name"].lower()

    for viral_id in (
        "hiv-protease-indinavir",
        "hiv-rt-nevirapine",
        "influenza-neuraminidase-oseltamivir",
    ):
        assert by_id[viral_id]["target_origin"] == "viral"

    assert by_id["hiv-rt-nevirapine"]["binding_site_type"] == "allosteric"
    assert by_id["hiv-protease-indinavir"]["binding_site_type"] == "catalytic"


def test_allosteric_mechanism_does_not_claim_active_site_blockade() -> None:
    """An allosteric inhibitor is expected to miss the catalytic residues."""
    from app.services.structure_analysis import build_mechanism

    analysis = {
        "status": "computed",
        "ligand_code": "NVP",
        "ligand_atom_count": 20,
        "binding_site_residue_count": 1,
        "contact_summary": {
            "hydrogen_bond": 0,
            "salt_bridge": 0,
            "hydrophobic": 5,
            "close_contact": 9,
        },
        "interacting_residues": [
            {
                "residue": "TYR",
                "chain": "A",
                "residue_number": "181",
                "residue_label": "TYR A181",
                "types": ["hydrophobic"],
                "contact_count": 3,
                "min_distance_angstrom": 3.04,
                "atom_pairs": [],
            }
        ],
        "binding_site_residues": ["TYR A181"],
    }
    meta = {
        "binding_site_type": "allosteric",
        "target_origin": "viral",
        "binding_pocket_residues": [{"label": "Tyr181", "match": "TYR A181"}],
        "catalytic_residues": [
            {"label": "Asp110 (catalytic)", "match": "ASP A110", "role": "catalytic"}
        ],
        "target_class": "HIV-1 reverse transcriptase",
        "drug_class": "NNRTI",
        "mechanism_summary": "summary",
        "consequence": "consequence",
    }

    mechanism = build_mechanism(analysis, meta)
    blockade = next(s for s in mechanism["steps"] if s["step"] == "blockade")

    assert blockade["value"] == "allosteric"
    assert "allosteric" in blockade["title"].lower()
    assert "not the active site" in blockade["title"].lower()
    assert "Tyr181" in blockade["detail"]
    # The catalytic residue is curated but must be reported as NOT in contact.
    assert mechanism["catalytic_residues_occluded"] == 0
    assert mechanism["catalytic_residues_curated"] == 1
    assert mechanism["binding_site_type"] == "allosteric"
    assert mechanism["catalytic_residues"][0]["in_contact"] is False


def test_catalytic_mechanism_still_reports_occlusion() -> None:
    """The catalytic branch must be unchanged by the allosteric branch."""
    from app.services.structure_analysis import build_mechanism

    analysis = {
        "status": "computed",
        "ligand_code": "MK1",
        "ligand_atom_count": 45,
        "binding_site_residue_count": 1,
        "contact_summary": {
            "hydrogen_bond": 1,
            "salt_bridge": 0,
            "hydrophobic": 3,
            "close_contact": 4,
        },
        "interacting_residues": [
            {
                "residue": "ASP",
                "chain": "A",
                "residue_number": "25",
                "residue_label": "ASP A25",
                "types": ["hydrogen_bond"],
                "contact_count": 1,
                "min_distance_angstrom": 2.77,
                "atom_pairs": [],
            }
        ],
        "binding_site_residues": ["ASP A25"],
    }
    meta = {
        "binding_site_type": "catalytic",
        "catalytic_residues": [
            {"label": "Asp25", "match": "ASP A25", "role": "catalytic aspartate"}
        ],
        "mechanism_summary": "summary",
        "consequence": "consequence",
    }

    mechanism = build_mechanism(analysis, meta)
    blockade = next(s for s in mechanism["steps"] if s["step"] == "blockade")
    assert blockade["value"] == "occluded"
    assert mechanism["catalytic_residues_occluded"] == 1


def test_complex_interactions_computed_from_coordinates(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        "app.services.complexes_service.http_get_bytes_cached",
        lambda *_args, **_kwargs: PDB_FIXTURE.encode("utf-8"),
    )
    response = client.get("/api/complexes/hiv-protease-indinavir/interactions")
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "computed", data.get("note")
    assert data["ligand_code"] == "MK1"
    assert data["contact_summary"]["hydrogen_bond"] >= 1
    assert data["contact_summary"]["salt_bridge"] >= 1
    assert "ASP A25" in data["binding_site_residues"]


def test_complex_interactions_placeholder_when_pdb_unavailable(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def _unreachable(*_args: object, **_kwargs: object) -> None:
        raise RuntimeError("RCSB unreachable")

    monkeypatch.setattr(
        "app.services.complexes_service.http_get_bytes_cached", _unreachable
    )
    response = client.get("/api/complexes/hiv-protease-indinavir/interactions")
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "placeholder"
    assert data["complex_id"] == "hiv-protease-indinavir"
    assert data["note"] is not None


def test_complex_structure_exposes_mechanism_and_coordinates(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        "app.services.complexes_service.http_get_bytes_cached",
        lambda *_args, **_kwargs: PDB_FIXTURE.encode("utf-8"),
    )
    response = client.get("/api/complexes/hiv-protease-indinavir/structure")
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "computed", data.get("note")
    # The viewer payload is the trimmed protein + ligand, not the raw deposit.
    assert "MK1" in data["structure_file_bytes"]
    assert data["structure_file_reference"] == "1HSG.pdb"
    mechanism = data["mechanism"]
    assert mechanism["status"] == "computed"
    steps = {step["step"] for step in mechanism["steps"]}
    assert {"target", "occupy", "anchor", "blockade", "consequence"} <= steps
    asp25 = next(c for c in mechanism["catalytic_residues"] if c["match"] == "ASP A25")
    assert asp25["in_contact"] is True
    assert asp25["min_distance_angstrom"] is not None


def test_complex_structure_unknown_id_is_placeholder() -> None:
    response = client.get("/api/complexes/not-a-real-complex/structure")
    assert response.status_code == 200
    assert response.json()["status"] == "placeholder"


def test_complex_docking_validation_placeholder() -> None:
    response = client.get("/api/complexes/hiv-protease-indinavir/docking-validation")
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "placeholder"
    assert data["complex_id"] == "hiv-protease-indinavir"
    assert data["rmsd_angstrom"] is None


def test_invalid_pdb_id() -> None:
    response = client.get("/api/protein/ABC")
    _assert_error_shape(response, "invalid_pdb_id", 400)


def test_protein_live_rcsb() -> None:
    require_live_network()


def test_protein_returns_placeholder_when_rcsb_fails(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def _rcsb_unreachable(*_args: object, **_kwargs: object) -> None:
        raise RuntimeError("RCSB unreachable")

    monkeypatch.setattr("app.services.rscb_service.pdb_entry_json", _rcsb_unreachable)
    response = client.get("/api/protein/6LU7")
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "placeholder"
    assert data["pdb_id"] == "6LU7"
    assert data["note"] is not None
    assert "RCSB" in data["note"]


def test_ai_resolve_no_match() -> None:
    response = client.post(
        "/api/ai/resolve",
        content=json.dumps({"query": "blahblah"}),
    )
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "no_match"
    assert data["resolved"] is None


def test_ai_resolve_scaffold_example() -> None:
    response = client.post(
        "/api/ai/resolve",
        content=json.dumps({"query": "covid protease"}),
    )
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "resolved"


def test_ai_predict_properties_not_trained() -> None:
    response = client.post(
        "/api/ai/predict-properties",
        content=json.dumps({"smiles": "CC(=O)Oc1ccccc1C(=O)O"}),
    )
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "not_trained"
    assert data["model"] == "not trained"
    assert data["method"] == "heuristic"
    assert data["predictions"] is None


def test_molecule_simulate_requires_smiles() -> None:
    response = client.post("/api/molecule/simulate", content=json.dumps({}))
    _assert_error_shape(response, "invalid_smiles", 422)


def test_molecule_simulate_unparseable_smiles_is_placeholder() -> None:
    response = client.post(
        "/api/molecule/simulate",
        content=json.dumps({"smiles": "not a valid smiles"}),
    )
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "placeholder"
    assert data["note"] is not None


def test_molecule_simulate_shape_with_mocked_chemistry(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Pin the simulation payload contract without depending on this machine's
    RDKit build, which varies with the OS application-control policy."""
    properties = {
        "molecular_weight": 180.16,
        "logp": 1.31,
        "h_bond_donors": 1,
        "h_bond_acceptors": 4,
        "rotatable_bonds": 3,
        "formula": "C9H8O4",
        "exact_mass": 180.04,
        "tpsa": 63.6,
        "heavy_atom_count": 13,
        "ring_count": 1,
        "aromatic_ring_count": 1,
        "fraction_csp3": 0.111,
        "veber": {"pass": True, "rotatable_bonds": 3, "tpsa": 63.6},
        "lipinski_rule_of_five": {"pass": True, "violations": 0},
    }
    conformers = [
        {
            "index": 0,
            "relative_energy_kcal_per_mol": 0.0,
            "minimised": True,
            "population": 0.7,
            "molblock": None,
        },
        {
            "index": 1,
            "relative_energy_kcal_per_mol": 1.0,
            "minimised": True,
            "population": 0.3,
            "molblock": None,
        },
    ]
    monkeypatch.setattr(
        "app.services.rdkit_helper.normalize_smiles_safe", lambda s: s
    )
    monkeypatch.setattr(
        "app.services.rdkit_helper.rdkit_properties_from_smiles", lambda s: properties
    )
    monkeypatch.setattr(
        "app.services.rdkit_helper.smiles_to_conformer_ensemble",
        lambda s, num_conformers=8, max_iters=400, temperature_k=298.15: {
            "available": True,
            "method": "mocked ensemble",
            "force_field": "MMFF94",
            "temperature_k": 298.15,
            "conformer_count": 2,
            "energy_span_kcal_per_mol": 1.0,
            "lowest_absolute_energy_kcal_per_mol": 10.0,
            "energy_reference": "relative_to_lowest_conformer_in_this_run",
            "absolute_energies_comparable": False,
            "absolute_energy_note": "absolute energies are not comparable between runs",
            "populations_are_estimates": True,
            "population_note": "populations are Boltzmann estimates",
            "conformers": conformers,
            "lowest_energy_index": 0,
            "lowest_energy_relative_kcal_per_mol": 0.0,
            "lowest_energy_molblock": "MOLBLOCK",
            "note": None,
        },
    )
    monkeypatch.setattr(
        "app.services.rdkit_helper.smiles_to_3d_structure",
        lambda s, num_conformers=8, max_iters=400: {
            "available": True,
            "is_3d": True,
            "molblock": "MOLBLOCK",
            "method": "mocked 3D",
            "force_field": "MMFF94",
            "relative_energy_kcal_per_mol": 0.0,
            "atom_count": 21,
            "conformer_count": 2,
            "note": None,
        },
    )
    monkeypatch.setattr(
        "app.services.rdkit_helper.smiles_to_2d_structure",
        lambda s: {"available": True, "molblock": "TWOD", "method": "mocked 2D", "atom_count": 13, "note": None},
    )

    response = client.post(
        "/api/molecule/simulate",
        content=json.dumps({"smiles": "CC(=O)Oc1ccccc1C(=O)O", "num_conformers": 2}),
    )
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "computed"
    assert data["structure"]["is_3d"] is True
    assert data["structure"]["molblock"] == "MOLBLOCK"
    assert data["ensemble"]["conformer_count"] == 2
    assert data["depiction"]["molblock"] == "TWOD"
    # Populations must be reported, and the sourced rules must be evaluated.
    assert [c["population"] for c in data["ensemble"]["conformers"]] == [0.7, 0.3]
    # A5: geometry method is a machine-readable code, energies are relative, and
    # populations are labelled as estimates.
    assert data["method"] in METHOD_CODES
    assert data["method_detail"]
    assert data["energy_basis"] == "relative_to_lowest_conformer_in_this_run"
    assert data["ensemble"]["absolute_energies_comparable"] is False
    assert data["ensemble"]["populations_are_estimates"] is True
    assert data["ensemble"]["population_note"]
    assert [c["relative_energy_kcal_per_mol"] for c in data["ensemble"]["conformers"]] == [0.0, 1.0]
    assert data["findings"], "drug-likeness findings should be present"
    assert all(f["source"] for f in data["findings"])
    steps = {step["step"] for step in data["simulation_steps"]}
    assert {"parse", "descriptors", "embed", "minimise", "population", "interpret"} <= steps
