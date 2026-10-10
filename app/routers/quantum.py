from __future__ import annotations

from typing import Union

from fastapi import APIRouter, HTTPException, Query
from fastapi.requests import Request
from fastapi.responses import JSONResponse

from app.services.quantum_service import (
    compute_h2_curve,
    compute_h2_simulate,
    compute_hadamard_demo,
    compute_resources,
)

QuantumResponse = Union[dict, JSONResponse]

router = APIRouter(prefix="/api/quantum", tags=["quantum"])


def _error_response(code: str, message: str, status_code: int) -> JSONResponse:
    return JSONResponse(status_code=status_code, content={"error": {"code": code, "message": message}})


@router.get("/h2-curve", response_model=None)
def get_h2_curve() -> QuantumResponse:
    try:
        return compute_h2_curve()
    except FileNotFoundError as exc:
        return _error_response("hamiltonian_unavailable", str(exc), 503)
    except Exception as exc:
        return _error_response("simulation_error", f"H2 curve computation failed: {exc}", 500)


@router.post("/simulate", response_model=None)
async def post_h2_simulate(request: Request) -> QuantumResponse:
    try:
        body = await request.json()
    except Exception:
        return _error_response("invalid_body", "Request body must be valid JSON", 400)

    if not isinstance(body, dict):
        return _error_response("invalid_body", "Request body must be a JSON object", 400)

    bond_length = body.get("bond_length")
    shots = body.get("shots")
    noise = body.get("noise")
    molecule = body.get("molecule")

    if molecule is not None and not isinstance(molecule, str):
        return _error_response("invalid_molecule", "Field 'molecule' must be a string", 422)
    if bond_length is None or not isinstance(bond_length, (int, float)):
        return _error_response("invalid_bond_length", "Field 'bond_length' is required and must be a number", 422)
    bond_length = float(bond_length)
    if not (0.0 <= bond_length <= 5.0):
        return _error_response("invalid_bond_length", "bond_length must be between 0.0 and 5.0", 422)

    shots = int(shots) if isinstance(shots, (int, float)) else 1024
    noise = float(noise) if isinstance(noise, (int, float)) else 0.0

    try:
        return compute_h2_simulate(bond_length=bond_length, shots=shots, noise=noise, molecule=molecule)
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
