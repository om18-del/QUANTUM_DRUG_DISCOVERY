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
- `GET /api/quantum/molecules`: the registry (currently `h2`, `lih`) with `label`, `file` and `available` per entry, where `available` is whether that artifact file exists. A known-but-absent molecule is `available: false` (its curve is `503`), not a placeholder curve.
- `GET /api/quantum/h2-curve?molecule=`: VQE and exact energies across bond lengths (default `molecule=h2`, so the old path is unchanged). For a non-H2 molecule it serves the stored curve only — HF, FCI, nuclear repulsion, `n_qubits`, the `active_space` text and an exact ground energy from a NumPy diagonalisation of the stored Pauli terms, cross-checked against the stored FCI to 1e-8 Ha — with no VQE fields. An unknown `molecule` is `404 unknown_molecule`; a missing artifact is `503 hamiltonian_unavailable`; more than 8 qubits in the stored Hamiltonian is `422 too_many_qubits`.
- `POST /api/quantum/simulate` `{molecule, bond_length (or bond_length_angstrom), shots (null = statevector), noise, seed, maxiter}`: runs a real VQE (fixed Ry ansatz, SciPy COBYLA, Qiskit estimator) against the artifact's Pauli terms; returns energies, error (mHa), iterations, convergence history and seed. The ansatz is 2 qubits, so VQE is H2 only: a `molecule` that is a different registry entry is `422 vqe_not_supported_for_molecule` and names what *is* available (the curve and the exact diagonalisation).
- `GET /api/quantum/hadamard-demo`
- `POST /api/quantum/resources` `{molecule}`: qubits, Hamiltonian terms, rough circuit size
- `GET /api/molecule/search?q=`, `POST /api/molecule/properties`
- `POST /api/molecule/simulate` `{smiles, num_conformers}`: 3D conformer ensemble + drug-likeness rules. Energy is reported relative to the lowest conformer in the same run; `method` is a machine-readable code (`3d_mmff94`, `3d_uff`, `3d_unminimised`, `2d_fallback`, `unavailable`).
- `GET /api/protein/{pdb_id}`
- `GET /api/complexes/list`, `GET /api/complexes/{id}/interactions`, `GET /api/complexes/{id}/structure`, `GET /api/complexes/{id}/docking-validation`
  - `interactions` and `structure` count **unique residues**, not atom pairs (`contact_summary_basis`), and return the exact rule behind each label in `contact_rules`. `structure` also returns the literature-linked mechanism steps.
  - Every curated entry declares `target_origin` (`viral` or `human`) and `binding_site_type` (`catalytic` or `allosteric`), so a human enzyme is never presented as a virus and an allosteric inhibitor is not reported as blocking the active site.
- `GET /api/jobs/{id}`: state is `queued | running | done | failed`; long runs return `202` with a job id
- `POST /api/ai/resolve`, `POST /api/ai/predict-properties` (optional, built last)

Error shape for every endpoint: `{"error": {"code": str, "message": str}}`. Validate all input with Pydantic and enforce size limits (for example SMILES length).

### Frontend delivery

During development and in single-service deployments, FastAPI also serves the frontend:

- `GET /` and `GET /index.html` return `frontend/static/index.html`.
- `GET /app/*` serves the static assets (`app.js`, `styles.css`) with `Cache-Control: no-cache, must-revalidate`, because a browser otherwise applies heuristic freshness to a stylesheet and keeps rendering a stale file.
- `GET /docs` redirects to the OpenAPI UI at `/api/docs`.

When the frontend is deployed separately (for example on Vercel), it targets this API using the configured CORS origins and these paths are simply unused. The API routes are unchanged either way, so the contract above holds in both deployments.

## 6. Quantum rules
- Hamiltonians are precomputed offline per bond length into JSON. The API never calls PySCF live.
- Bond-length grid of the offline H2 artifact: **58 points** — every point of 0.20 to 3.00 angstrom in 0.05 steps (57 uniform points), plus **0.74** angstrom, the UI default and H2's real equilibrium separation. 0.74 is not on the 0.05 step (it would need index 10.8), so it is added explicitly. Ascending, with no duplicates.
- Bond lengths that are not on that grid: `POST /api/quantum/simulate` **snaps to the nearest stored point and discloses it**. It does not reject the request, and it **never interpolates an energy** — a number is only ever reported for a separation that was actually computed offline, because an interpolated curve between computed points would be our own construction rather than a result. The response echoes both values so a snapped answer cannot be read as an exact one: `bond_length_requested_angstrom`, `bond_length_used_angstrom`, `bond_length_is_exact` (`true` only when the two are equal), and a `note` naming the stored point that answered. An exact tie (equidistant stored points) resolves to the **shorter** bond length. `GET /api/quantum/h2-curve` is unaffected: it reports the stored points themselves.
- Molecules in scope: H2 (core), LiH, active-space H2O. Larger drugs (for example aspirin): classical full-molecule result plus the resource estimator, no live VQE.
- **Only H2 is simulated with VQE.** The single ansatz in `app/core/vqe.py` is a 2-qubit Ry+CX circuit, so a second molecule would need a second ansatz that does not exist yet. LiH (active space) serves its precomputed curve and an exact reference read/diagonalised from its own artifact; it is never reported as a VQE result. The exact reference is a dense NumPy diagonalisation and is refused above 8 qubits (`422 too_many_qubits`) rather than approximated.
- Accuracy target: **noiseless** VQE within 1.6 mHa of exact for H2. Noisy and low-shot modes are expected to show larger error. Report the achieved error honestly and never force a pass.
- VQE energy must never be below the exact energy (variational principle). Reference check: H2 in STO-3G near 0.74 angstrom is about -1.137 Ha. A run that reports below the exact energy by more than 1e-9 Ha is **refused** with a `variational_principle_violation` error instead of being served.
- **Scope of that guard: noiseless statevector runs only.** A single shot-based estimate is a statistic and can legitimately sit below the exact energy on sampling noise, so a shot run is never refused for that reason. It instead reports its uncertainty as `vqe_energy_std_ha` and keeps `within_chemical_accuracy` as `null` (chemical accuracy remains a noiseless claim). A statevector run keeps the strict 1e-9 Ha guard.
- **No winner's curse in a shot run's headline number.** The optimiser is judged by the noiseless expectation at each restart's parameters, but the reported `vqe_energy_ha` for a shot run is a **fresh re-evaluation at the final optimal parameters**, not the minimum of the noisy optimisation trace (which is biased downward by selection). Its uncertainty is the estimator's own variance metadata when the backend reports it, otherwise the mean and sample standard deviation of 10 re-evaluations at those parameters with different sampler seeds; the response names the method used (`vqe_energy_uncertainty_method`).
- **The optimiser budget is spelled out.** A run reports `n_restarts`, `maxiter_per_restart` and `total_energy_evaluations` alongside the backward-compatible `iterations` (equal to `total_energy_evaluations`), so "450 iterations" reads as 3 restarts x 150 evaluations per restart rather than an opaque total.

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
9. Phase 3 extensions (section 11) — scheduled only after the VQE core above is implemented and verified.

**The VQE core stays first.** The offline Hamiltonian artifact, the VQE loop and exact solver, the quantum endpoints, then LiH / active-space H2O and the resource estimator are the project's first priority; no Phase 3 item may be moved ahead of them (section 11, "Ordering rule").

## 10. Definition of done for each step
Clean install from `requirements.txt`; server starts; new endpoints return the documented shape; tests pass locally and in CI; assumptions listed; changes committed.

## 11. Phase 3: extensions

Phase 3 is everything that comes after the core is working. Each item below states its goal, inputs, outputs, tests, what is exact versus estimated, and its risks.

**Ordering rule.** The VQE core (B3–B7 in the project plan) stays **first** in the project order: the offline Hamiltonian artifact, the VQE loop and exact solver, the quantum endpoints, then LiH / active-space H2O and the resource estimator. Phase 3 work does not jump that queue — an item here starts only once the core is implemented and verified (section 9).

### P3-1 Resource estimator (`POST /api/quantum/resources`)
- **Goal.** For molecules too large for live VQE, report the qubit count, the Hamiltonian term count and a rough circuit size, and plot qubits against molecule size.
- **Inputs.** A molecule name or SMILES; the STO-3G orbital count for it.
- **Outputs.** The documented endpoint is **extended, not replaced**: `n_qubits`, `hamiltonian_terms`, `estimated_circuit_gates`, `estimated_circuit_depth`, `artifact`, `status`, `note`, plus the orbital count and a per-field marker saying whether the field is exact or estimated. Also the qubits-versus-size series for H2, LiH, H2O, aspirin and one larger drug.
- **Tests.** The orbital-to-qubit mapping for every molecule on the plot; the ansatz parameter formulas; qubits increasing with molecule size; the response shape of section 5; the `status`/`note` response when the input is missing.
- **Exact vs estimated.** Qubit counts are claimed **exact**, from the STO-3G orbital count. Ansatz parameter counts are exact for the stated ansatz. **Pauli term counts and circuit depth are "estimated scaling"** unless computed offline for a named small molecule — the response and the UI must say which one a number is.
- **Risks.** An orbital count that does not match the basis actually used; publishing an estimate as exact; a plot that implies accuracy for molecules we never simulate.
- **Blocker before this ships.** The orbital and qubit counts are currently **hand-derived and unverified**. A test must confirm the mapping for every molecule on the plot before the estimator ships; until that test passes, the endpoint must not claim `exact` and the figure must not be used in slides or reports (`RESULTS.md` rule).

### P3-2 Noise and shots study (simulation)
- **Goal.** Show how sampling noise and shots move the VQE energy away from the exact value.
- **Inputs.** The H2 artifact and the simulate endpoint; a shots list, noise parameters, a repeat count and an RNG seed.
- **Outputs.** Per-setting tables (mean, spread, error versus exact), a per-run convergence history, and stored JSON for the write-up.
- **Tests.** The noiseless result still meets the section 6 target; more shots reduce the spread; noisy runs are worse than noiseless; the VQE energy is never below the exact energy.
- **Exact vs estimated.** The exact energy is exact. VQE energies come from sampling, so they carry a spread — report the spread, never a single favourable run.
- **Risks.** A noise model that does not resemble real hardware; confusing device noise with our own sampling statistic; reading a seed-dependent wobble as a trend.
- **Dependency rule.** Before adding any package, check **which package actually provides a device-based noise model compatible with the pinned `qiskit` 1.4.6 and `qiskit-aer`**, and report the versions before adding it. Do not add a dependency on the assumption that one exists.

### P3-3 Redocking validation (**offline**, Colab/WSL)
- **Goal.** Independently check that a docking protocol reproduces known crystal poses.
- **Inputs.** 1VRT, 2HU4 and 1HSG with their deposited ligands; AutoDock Vina, run offline only.
- **Outputs.** Per complex: `crystal_pose_file`, `docked_pose_file`, `rmsd_angstrom`, `status` (`computed` or `failed`), and the run record — Vina version, box, exhaustiveness, seed, and the symmetry handling used.
- **Tests.** RMSD computed with the stated symmetry handling; the 2.0 angstrom success threshold applied; `failed` runs reported as failed, not omitted.
- **Exact vs estimated.** The RMSD is a measured number for the run we performed. The docked pose itself is a **prediction, not measured binding**, and must be labelled that way in the API and the UI.
- **Risks.** Using a covalent ligand (never 6LU7/N3 for redocking — section 7); symmetry handling changing the RMSD; presenting one seed as a reproducible protocol.

### P3-4 One real IBM hardware run of H2
- **Goal.** One honest hardware comparison against simulation, run offline.
- **Inputs.** An IBM Quantum token **from environment variables only, never committed**; the H2 artifact; a few bond lengths.
- **Outputs.** JSON per run with the backend name, the date and the shots, carrying the hardware error next to the simulated and exact energies.
- **Tests.** The token is absent from the repo and from the frontend; the stored JSON loads into the API; the hardware error is reported as measured and never forced to pass.
- **Exact vs estimated.** The hardware energy is measured on that backend on that date; the exact energy is exact; a single run must not be generalised into a claim about hardware accuracy.
- **Risks.** Free-tier queueing and quota; a backend that changes between submission and run; **state current free-tier limits only after checking them on IBM's site, and recheck them on the day of the run** — record the date they were checked.

### P3-5 Guided lesson mode (frontend)
- **Goal.** A short, sourced walkthrough with one question after each step.
- **Inputs.** The existing panels and endpoints; a lesson definition of steps, source links and one question per step.
- **Outputs.** A lesson panel in the frontend. No accounts, no grading, answers are not stored.
- **Tests.** Every step cites a source; each step has exactly one question; nothing is submitted or persisted.
- **Exact vs estimated.** Explanatory only, but it displays numbers from the same endpoints as the rest of the app, so the `RESULTS.md` rule applies to them too.
- **Risks.** Drifting from the sourced literature; making the lesson look like an assessment; adding accounts or tracking later without saying so.

### P3-6 Write-up discipline
- **Goal.** Every number in a slide or a report is traceable to the command that produced it.
- **Inputs.** A measured number and the command that produced it.
- **Outputs.** A row in `RESULTS.md` (number, value, command, date, file) and, where relevant, an entry in `LIMITATIONS.md`.
- **Tests.** Every number used in the write-up has a `RESULTS.md` row; every row names a reproducible command; every unverified item appears in `LIMITATIONS.md`.
- **Exact vs estimated.** Nothing here is estimated: an item is either listed with the command that produced it, or it does not appear.
- **Risks.** A number copied from memory; re-running a command that yields a different value without adding a new row; slides outrunning the two files.
