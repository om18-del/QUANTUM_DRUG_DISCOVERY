# SPEC.md: Quantum Drug Discovery Explorer

Give this file to the coding agent at the start of EVERY session. If a request conflicts with it, stop and ask.

## 1. What we are building
A hosted **web app** (frontend and API deployed separately) for students and educators, with three tabs:
1. **Quantum Lab**: VQE energy calculations on small molecules, checked against exact classical results; noise and shots study; resource estimator for larger molecules.
2. **Molecule Explorer**: search a medicine (PubChem) or protein (RCSB PDB), 3D view, properties, Lipinski check.
3. **Drug-Target Studio**: curated real drug-protein complexes with interaction tables and docking validation.

It is a study and exploration tool. It does not discover drugs and is not medical advice.

## 2. Non-negotiable rules
- **No fabricated data.** Every number is computed by us or fetched from a named source. Missing data is `null` with `"status": "placeholder"`, and the UI shows "not computed".
- **No fake AI.** Any fallback prediction must return `"method": "heuristic"` (or `"model": "not trained"`) and say so.
- **No overclaiming quantum.** Live VQE only for small problems (about 12 qubits or fewer). Larger molecules get the classical result and a resource estimate.
- **Report every assumption** (library versions, endpoint shapes, field names) in your reply before moving on.
- **Never commit secrets.** Config comes from environment variables. `.gitignore` covers `.env`, caches, models and local artifacts.
- Do not add Redis, Celery or other heavy infrastructure. Use FastAPI BackgroundTasks with a SQLite job store.

## 3. Stack
- Python 3.11 (change only if an install fails, and say so). FastAPI, uvicorn, Pydantic, httpx, tenacity, diskcache.
- Quantum: `qiskit` (not qiskit-terra), `qiskit-aer`. Pin all versions. Confirm exact estimator import paths with a small probe test before relying on them. Noiseless: `StatevectorEstimator`. Noisy: Aer estimator.
- If `qiskit_algorithms` VQE does not accept V2 primitives, write the VQE loop manually (estimator + SciPy COBYLA).
- Heavy chemistry tools (PySCF, Qiskit Nature, RDKit training) run OFFLINE (WSL or Colab). The API loads their output files.
- Hosting: the API on Render, Railway or a Hugging Face Space (not Vercel, because the libraries are too large). The frontend can be on Vercel. CORS origins come from env.
- Free hosts may wipe local disk: treat cache and job store as disposable.

## 4. Repo layout
```
app/main.py, config.py
app/core/vqe.py, exact.py, resources.py
app/routers/quantum.py, molecules.py, proteins.py, complexes.py, jobs.py, ai.py
app/services/cache_manager.py, external_api.py, rdkit_helper.py
app/data/hamiltonians/*.json
app/data/curated_complexes/manifest.json
app/data/models/manifest.json
tests/  scripts/precompute_hamiltonians.py
DATA_SOURCES.md  .gitignore  requirements.txt
```

## 5. API contract (stable; frontend builds against this)
- `GET /health`
- `GET /api/quantum/h2-curve`: VQE and exact energies across bond lengths
- `POST /api/quantum/simulate` `{molecule, bond_length, shots, noise}`: returns energy, exact energy, error (mHa), convergence history
- `GET /api/quantum/hadamard-demo`
- `POST /api/quantum/resources` `{molecule}`: qubits, Hamiltonian terms, rough circuit size
- `GET /api/molecule/search?q=`, `POST /api/molecule/properties`
- `GET /api/protein/{pdb_id}`
- `GET /api/complexes/list`, `GET /api/complexes/{id}/interactions`, `GET /api/complexes/{id}/docking-validation`
- `GET /api/jobs/{id}`: state is `queued | running | done | failed`; long runs return `202` with a job id
- `POST /api/ai/resolve`, `POST /api/ai/predict-properties` (optional, built last)

Error shape for every endpoint: `{"error": {"code": str, "message": str}}`. Validate all input with Pydantic and enforce size limits (for example SMILES length).

## 6. Quantum rules
- Hamiltonians are precomputed offline per bond length into JSON. The API never calls PySCF live.
- Molecules in scope: H2 (core), LiH, active-space H2O. Larger drugs (for example aspirin): classical full-molecule result plus the resource estimator, no live VQE.
- Accuracy target: **noiseless** VQE within 1.6 mHa of exact for H2. Noisy and low-shot modes are expected to show larger error. Report the achieved error honestly and never force a pass.
- VQE energy must never be below the exact energy (variational principle). Reference check: H2 in STO-3G near 0.74 angstrom is about -1.137 Ha.

## 7. Data rules
- External calls: 10 s timeout, retries, cache-first, per-source TTLs. PubChem/RCSB: finite TTL. Lipinski by SMILES: effectively infinite. Complexes: immutable.
- Validate external responses with Pydantic before returning. Normalize SMILES. Check PDB ID format.
- Complexes: use **non-covalent** pairs for docking validation (for example HIV protease with a known inhibitor). Do not use covalent ligands such as 6LU7/N3 for redocking. Every entry has a status: `computed` or `placeholder`. Never invent RMSD or interaction values.
- Keep `DATA_SOURCES.md` current (source, licence, paper behind each method).

## 8. Optional ML (last)
One task first: ESOL solubility regression. Scaffold split (not random). Report RMSE against a mean-prediction baseline. Add an applicability warning for out-of-range molecules. Train offline; save the model with its feature pipeline and a manifest entry; pin sklearn and RDKit versions; load only our own model files. If no model is loaded, return the honest "not trained" response.

## 9. Build order (finish and verify each step before the next)
1. Repo hygiene (`.gitignore`, secrets), then `/health` and a skeleton deploy on a host.
2. Offline Hamiltonian script, then VQE core and exact solver with tests; GitHub Actions runs pytest on push.
3. Quantum endpoints (frontend starts in parallel against mocked responses).
4. LiH, active-space H2O, resource estimator.
5. PubChem and RCSB lookups with validation and cache.
6. Real complexes with manifest and placeholder handling.
7. Jobs and resilience.
8. Optional AI.

## 10. Definition of done for each step
Clean install from `requirements.txt`; server starts; new endpoints return the documented shape; tests pass locally and in CI; assumptions listed; changes committed.
