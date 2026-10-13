from __future__ import annotations

from typing import Union

from fastapi import APIRouter, HTTPException, Query
from fastapi.requests import Request
from fastapi.responses import JSONResponse

from app.services.quantum_service import (
    DEFAULT_MOLECULE,
    QuantumRequestError,
    VariationalPrincipleViolation,
    compute_curve,
    compute_h2_simulate,
    compute_hadamard_demo,
    compute_resources,
    list_molecules,
)

QuantumResponse = Union[dict, JSONResponse]

router = APIRouter(prefix="/api/quantum", tags=["quantum"])


def _error_response(code: str, message: str, status_code: int) -> JSONResponse:
    return JSONResponse(status_code=status_code, content={"error": {"code": code, "message": message}})


@router.get("/molecules", response_model=None)
def get_molecules() -> dict:
    """The molecule registry, with whether each artifact file is present."""
    return list_molecules()


@router.get("/h2-curve", response_model=None)
def get_h2_curve(molecule: str = Query(DEFAULT_MOLECULE)) -> QuantumResponse:
    """The stored curve for a molecule (default H2); the path is unchanged."""
    label = "H2" if str(molecule).strip().lower() == "h2" else str(molecule)
    try:
        return compute_curve(molecule)
    except QuantumRequestError as exc:
        return _error_response(exc.code, exc.message, exc.status_code)
    except FileNotFoundError as exc:
        return _error_response("hamiltonian_unavailable", str(exc), 503)
    except Exception as exc:
        return _error_response("simulation_error", f"{label} curve computation failed: {exc}", 500)


def _parse_optional_int(body: dict, field: str, code: str, low: int, high: int):
    """An optional integer in [low, high], ``None`` when absent, else an error."""
    value = body.get(field)
    if value is None:
        return None, None
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or int(value) != value
        or not (low <= int(value) <= high)
    ):
        return None, _error_response(
            code,
            f"Field '{field}' must be an integer between {low} and {high}",
            422,
        )
    return int(value), None


@router.post("/simulate", response_model=None)
async def post_h2_simulate(request: Request) -> QuantumResponse:
    try:
        body = await request.json()
    except Exception:
        return _error_response("invalid_body", "Request body must be valid JSON", 400)

    if not isinstance(body, dict):
        return _error_response("invalid_body", "Request body must be a JSON object", 400)

    # Both spellings are accepted; the artifact-shaped name wins when both appear.
    bond_length = (
        body["bond_length_angstrom"] if "bond_length_angstrom" in body else body.get("bond_length")
    )
    shots = body.get("shots")
    noise = body.get("noise")
    molecule = body.get("molecule")

    if molecule is not None and not isinstance(molecule, str):
        return _error_response("invalid_molecule", "Field 'molecule' must be a string", 422)
    if bond_length is None or not isinstance(bond_length, (int, float)) or isinstance(bond_length, bool):
        return _error_response(
            "invalid_bond_length",
            "Field 'bond_length' (or 'bond_length_angstrom') is required and must be a number",
            422,
        )
    bond_length = float(bond_length)
    if not (0.0 <= bond_length <= 5.0):
        return _error_response("invalid_bond_length", "bond_length must be between 0.0 and 5.0", 422)

    # shots: absent or null means the noiseless statevector run (SPEC/task);
    # a number is the per-evaluation sample count of the shot-based estimator.
    if shots is not None and (
        isinstance(shots, bool)
        or not isinstance(shots, (int, float))
        or int(shots) != shots
        or int(shots) < 1
    ):
        return _error_response(
            "invalid_shots",
            "Field 'shots' must be a positive integer, or null for a noiseless statevector run",
            422,
        )
    shots = int(shots) if shots is not None else None

    seed, error = _parse_optional_int(body, "seed", "invalid_seed", 0, 2**32 - 1)
    if error is not None:
        return error
    maxiter, error = _parse_optional_int(body, "maxiter", "invalid_maxiter", 1, 10000)
    if error is not None:
        return error

    noise = float(noise) if isinstance(noise, (int, float)) and not isinstance(noise, bool) else 0.0

    try:
        return compute_h2_simulate(
            bond_length=bond_length,
            shots=shots,
            noise=noise,
            molecule=molecule,
            seed=seed,
            maxiter=maxiter,
        )
    except QuantumRequestError as exc:
        # A request-level answer (unknown molecule, unsupported VQE, too many
        # qubits), never a fault: its own code and status, no result fields.
        return _error_response(exc.code, exc.message, exc.status_code)
    except VariationalPrincipleViolation as exc:
        # The optimiser broke the variational principle: refuse the result with
        # its own code instead of reporting a physically impossible energy.
        return _error_response("variational_principle_violation", str(exc), 500)
    except FileNotFoundError as exc:
        return _error_response("hamiltonian_unavailable", str(exc), 503)
    except Exception as exc:
        return _error_response("simulation_error", f"H2 simulation failed: {exc}", 500)


@router.get("/hadamard-demo")
def get_hadamard_demo() -> dict:
    return compute_hadamard_demo()


@router.post("/resources", response_model=None)
async def post_resources(request: Request) -> QuantumResponse:
    try:
        body = await request.json()
    except Exception:
        return _error_response("invalid_body", "Request body must be valid JSON", 400)

    if not isinstance(body, dict):
        return _error_response("invalid_body", "Request body must be a JSON object", 400)

    molecule = body.get("molecule")
    if not isinstance(molecule, str) or not molecule.strip():
        return _error_response("invalid_molecule", "Field 'molecule' is required and must be a non-empty string", 422)

    return compute_resources(molecule)
