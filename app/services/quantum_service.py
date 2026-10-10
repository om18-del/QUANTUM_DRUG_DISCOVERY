"""Quantum Lab service.

The H2 curve is meant to be built from the offline artifact
``app/data/hamiltonians/h2_sto3g.json``, produced by
``scripts/precompute_hamiltonians.ipynb`` (Colab or WSL, because PySCF publishes
no Windows wheel).

That artifact is not in the repository yet, so these endpoints say so instead of
serving numbers. An earlier build shipped a hand-typed 4x4 matrix in
``app/data/hamiltonians/`` and reported it as ``status: "computed"``; those files
have been deleted, because presenting them as computed H2 physics was wrong.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from app.config import get_settings

H2_HAMILTONIAN_ARTIFACT = "h2_sto3g.json"

# The curve, the single-point simulation and the resource estimate are three
# different things withheld for the same missing artifact, so each one names its
# own subject. A single shared sentence used to tell the simulate panel that "no
# H2 curve is available", which is not what that endpoint computes.
CURVE_SUBJECT = "curve"
ENERGY_SUBJECT = "energy"
RESOURCE_SUBJECT = "resource estimate"


def _unavailable_note(subject: str) -> str:
    """Why nothing can be served yet, phrased for this endpoint's own subject."""
    return (
        f"{H2_HAMILTONIAN_ARTIFACT} has not been generated, so no H2 {subject} is available. "
        "Run scripts/precompute_hamiltonians.ipynb in Colab or WSL and save the result into "
        f"app/data/hamiltonians/. This endpoint reports the {subject} as unavailable rather than "
        "showing placeholder numbers."
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


def compute_h2_curve() -> dict:
    """H2 dissociation curve, or an honest unavailable response.

    Raises ``FileNotFoundError`` when the offline artifact has not been generated,
    which the router turns into the standard error shape.
    """
    artifact = h2_artifact_path()
    if not artifact.exists():
        raise FileNotFoundError(_unavailable_note(CURVE_SUBJECT))

    return {
        "molecule": "H2",
        "basis": "sto-3g",
        "artifact": artifact.name,
        "status": "placeholder",
        "note": _unwired_note(CURVE_SUBJECT).format(artifact=artifact.name),
        "bond_lengths_angstrom": [],
        "exact_energies_ha": [],
        "vqe_energies_ha": [],
        "vqe_error_mha": [],
    }


def compute_h2_simulate(
    bond_length: float,
    shots: int,
    noise: float,
    molecule: str | None = None,
) -> dict:
    """One H2 VQE simulation, or an honest unavailable response.

    ``shots`` and ``noise`` are accepted for the contract and echoed back, but no
    shot or noise simulation is performed in this build.
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

    artifact = h2_artifact_path()
    if not artifact.exists():
        raise FileNotFoundError(_unavailable_note(ENERGY_SUBJECT))

    return {
        "molecule": "H2",
        "bond_length": bond_length,
        "exact_energy_ha": None,
        "vqe_energy_ha": None,
        "error_mha": None,
        "convergence_history_ha": None,
        "shots": shots,
        "noise": noise,
        "estimator_class_name": None,
        "status": "placeholder",
        "note": _unwired_note(ENERGY_SUBJECT).format(artifact=artifact.name),
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
    they are withheld until the real artifact is wired in.
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
