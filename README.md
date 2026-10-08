# Quantum Drug Discovery Explorer

A hosted web app for studying quantum simulations, molecules, proteins, and drug-target
complexes. It is an exploration and education tool, not a drug discovery pipeline and not
medical advice.

## Project structure

- `app/` — FastAPI backend
- `tests/` — pytest suite
- `scripts/` — offline tooling such as Hamiltonian precomputation
- `app/data/` — precomputed artifacts, manifests, and cached data (not all committed)

## Spec

Build to [SPEC.md](SPEC.md).

## Local development

1. Create a Python 3.11 environment.
2. Install dependencies:
   - `pip install -r requirements.txt`
3. Copy `.env.example` to `.env` if you need overrides.
4. Start the server:
   - `uvicorn app.main:app --reload`

Health check:
- `GET /health`

## Important rules

- No fabricated data.
- No fake AI.
- No secrets in the repo.
- Heavy chemistry tooling is installed offline when needed.
- Free hosting may wipe local disk; treat cache and job store as disposable unless persisted.

## Build order

1. Repo hygiene, `/health`, and a deployable skeleton.
2. Offline Hamiltonian precompute, then VQE core and tests.
3. Quantum endpoints; frontend can start here against mocks.
4. LiH, active-space H2O, and resource estimator.
5. PubChem and RCSB lookups with validation and cache.
6. Real complexes with manifest and placeholder handling.
7. Jobs and resilience.
8. Optional AI: ESOL model first.
