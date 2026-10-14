# API Contract

This document is the stable contract the frontend can target. The backend must not
change these shapes without updating this file and the SPEC.md build order.

All sample payloads below are sample data, not real computed results.

## Base URL

- Local development: `http://127.0.0.1:8000`
- Production: set by deployment environment.

All API routes are prefixed with `/api` unless noted otherwise.

## Common error shape

All errors use this shape:

```json
{
  "error": {
    "code": "string",
    "message": "string"
  }
}
```

Input is validated with Pydantic. Invalid input returns `422` with the same error shape where
applicable. Rate limits, timeouts, and upstream failures return error responses, not raw traces.

## Authentication

None by default. If authentication is added later, it must be documented here and in SPEC.md.

## CORS

CORS origins are configured via environment variables. The backend must not hardcode trusted
frontend origins in code.

## Frontend delivery

In this deployment FastAPI also serves the frontend, so a single process runs the whole app:

- `GET /` and `GET /index.html` return `frontend/static/index.html`.
- `GET /app/*` serves the static assets with `Cache-Control: no-cache, must-revalidate`.
  A stylesheet without an explicit directive is cached heuristically by the browser, which
  makes edits appear not to apply. The stylesheet URL is also versioned (`styles.css?v=N`).
- `GET /docs` redirects to `/api/docs` (the OpenAPI UI).

If the frontend is instead deployed separately (for example on Vercel), it calls the same API
routes and the `/`, `/index.html` and `/app/*` routes are simply unused. The API contract below
is identical in both deployments.

## Endpoints

### Health

- `GET /health`
- Returns `200` with a simple status payload.
- Used by deployment monitors and the frontend bootstrap.

Example request:
```http
GET /health HTTP/1.1
Host: 127.0.0.1:8000
```

Example response:
```json
{
  "status": "ok",
  "app": "Quantum Drug Discovery Explorer API"
}
```

### Quantum

#### GET /api/quantum/molecules

- Lists the molecules the artifact layer knows about, with whether each one's artifact file is
  actually present on disk. It reads the registry and stats two files; it computes nothing.
- `available` is `true` only when `app/data/hamiltonians/<file>` exists. A molecule with
  `available: false` is a **known** molecule whose artifact has not been generated yet — its
  curve answers `503`, never placeholder numbers.
- The registry is deliberately tiny (H2 and LiH). A molecule absent from this list is not
  "coming soon": the curve endpoint answers `404 unknown_molecule` for it.

Example request:
```http
GET /api/quantum/molecules HTTP/1.1
Host: 127.0.0.1:8000
```

Example response (LiH's artifact is not generated in this repository):
```json
{
  "molecules": [
    {
      "molecule": "h2",
      "label": "H2",
      "file": "h2_sto3g.json",
      "available": true
    },
    {
      "molecule": "lih",
      "label": "LiH (active space)",
      "file": "lih_sto3g_active.json",
      "available": false
    }
  ],
  "count": 2,
  "available_count": 1
}
```

Example response for `GET /api/quantum/h2-curve?molecule=nope` — standard error shape, `404`:
```json
{
  "error": {
    "code": "unknown_molecule",
    "message": "unknown molecule 'nope'; expected one of ['h2', 'lih']"
  }
}
```

#### GET /api/quantum/h2-curve

- Optional query parameter `molecule` (default `h2`), so the existing path keeps working
  unchanged. `h2` returns the response documented in this section; `lih` returns the
  artifact-only response documented below.
- Returns HF, exact (FCI) and VQE energies across bond lengths, built from the offline artifact
  `app/data/hamiltonians/h2_sto3g.json` (58 stored points: 57 uniform from 0.20 to 3.00 Å plus
  the UI default 0.74 Å).
- The artifact **is** in the repository. PySCF publishes no Windows wheel, so it was generated
  with `scripts/precompute_hamiltonians.ipynb` in Colab or WSL and is checked by
  `scripts/validate_hamiltonian_file.py`; this endpoint only reads it.
- `metadata` carries the documented provenance, checked against the file before anything is
  served: basis `STO-3G`, `parity` mapping with two-qubit reduction, `n_qubits: 2`, and
  `source: "precomputed offline with PySCF + qiskit-nature"`. An artifact that declares
  something else answers `500` instead of serving the mismatch.
- The parallel arrays (`bond_lengths_angstrom`, `hf_energies_ha`, `exact_energies_ha`,
  `vqe_energies_ha`, `nuclear_repulsion_energies_ha`) and the per-point `points` list answer the
  same data twice: `hf`/`exact`/`vqe` are TOTAL energies (electronic plus nuclear repulsion) as
  the artifact's `convention` field states, with the nuclear part also reported separately.
- No optimiser, shots or noise run inside this request: each point is the stored PySCF HF/FCI
  reference and the stored ground energy of the mapped qubit Hamiltonian plus nuclear
  repulsion, and `note` says so on every response. The real VQE loop lives on the simulate
  endpoint below.
- An earlier build served a hand-typed 4x4 matrix from `app/data/hamiltonians/` and labelled it
  `"status": "computed"`. Those files have been deleted. This endpoint serves either the stored
  grid or an error (503 when the file is absent, 500 when it is present but malformed); it never
  fills a gap with placeholder numbers.

Example request:
```http
GET /api/quantum/h2-curve HTTP/1.1
Host: 127.0.0.1:8000
```

Example response while the artifact is missing — standard error shape, `503`:
```json
{
  "error": {
    "code": "hamiltonian_unavailable",
    "message": "h2_sto3g.json has not been generated, so no H2 curve is available. Run scripts/precompute_hamiltonians.ipynb in Colab or WSL and save the result into app/data/hamiltonians/. This endpoint reports the curve as unavailable rather than showing placeholder numbers."
  }
}
```

Example response with the committed artifact (the parallel arrays and the `points` list
truncated to their first four of 58 entries; every number below is copied verbatim from the
live response):
```json
{
  "molecule": "H2",
  "basis": "STO-3G",
  "metadata": {
    "basis": "STO-3G",
    "mapping": "parity",
    "two_qubit_reduction": true,
    "mapping_note": "parity mapping with two-qubit reduction",
    "n_qubits": 2,
    "source": "precomputed offline with PySCF + qiskit-nature"
  },
  "artifact": "h2_sto3g.json",
  "status": "computed",
  "note": "Computed from the 58 stored points in h2_sto3g.json. Exact energies are the stored PySCF FCI reference, and VQE energies are the stored ground energy of the mapped qubit Hamiltonian plus nuclear repulsion, already minimised in the offline precompute. No optimisation, shots or noise run in this request.",
  "units": "hartree",
  "convention": "Each entry's pauli_terms describe the electronic Hamiltonian only. The smallest eigenvalue of those terms plus nuclear_repulsion_ha equals fci_energy_ha. Determine the ground energy by adding nuclear_repulsion_ha.",
  "point_count": 58,
  "bond_lengths_angstrom": [
    0.2,
    0.25,
    0.3,
    0.35
  ],
  "hf_energies_ha": [
    0.16417501206169538,
    -0.3050005337593622,
    -0.5938277585357272,
    -0.7804549020715259
  ],
  "exact_energies_ha": [
    0.15748213479836348,
    -0.3122699027007898,
    -0.6018037107656871,
    -0.7892693924044136
  ],
  "vqe_energies_ha": [
    0.15748213479836437,
    -0.312269902700788,
    -0.6018037107656862,
    -0.7892693924044119
  ],
  "nuclear_repulsion_energies_ha": [
    2.6458860546,
    2.11670884368,
    1.7639240364000002,
    1.5119348883428574
  ],
  "vqe_error_mha": [
    8.881784197001252e-13,
    1.7763568394002505e-12,
    8.881784197001252e-13,
    1.7763568394002505e-12
  ],
  "max_deviation_from_fci_ha": 2.4424906541753444e-15,
  "points": [
    {
      "bond_length_angstrom": 0.2,
      "hf_energy_ha": 0.16417501206169538,
      "fci_energy_ha": 0.15748213479836348,
      "nuclear_repulsion_ha": 2.6458860546,
      "vqe_energy_ha": 0.15748213479836437
    },
    {
      "bond_length_angstrom": 0.25,
      "hf_energy_ha": -0.3050005337593622,
      "fci_energy_ha": -0.3122699027007898,
      "nuclear_repulsion_ha": 2.11670884368,
      "vqe_energy_ha": -0.312269902700788
    },
    {
      "bond_length_angstrom": 0.3,
      "hf_energy_ha": -0.5938277585357272,
      "fci_energy_ha": -0.6018037107656871,
      "nuclear_repulsion_ha": 1.7639240364000002,
      "vqe_energy_ha": -0.6018037107656862
    },
    {
      "bond_length_angstrom": 0.35,
      "hf_energy_ha": -0.7804549020715259,
      "fci_energy_ha": -0.7892693924044136,
      "nuclear_repulsion_ha": 1.5119348883428574,
      "vqe_energy_ha": -0.7892693924044119
    }
  ]
}
```

#### GET /api/quantum/h2-curve?molecule=lih (any non-H2 registry molecule)

- **No VQE is claimed here.** For a molecule other than H2 this endpoint serves exactly two
  things: the stored artifact data, and an exact ground energy from a NumPy diagonalisation of
  the same stored Pauli terms. The response deliberately has **no** `vqe_energies_ha` and no
  `vqe_error_mha`, so nothing can be read as a variational result for a molecule we never
  optimised.
- Fields: `molecule`, `basis`, `n_qubits`, `active_space` (the free-text description written by
  the precompute script, saying exactly which orbitals were frozen and removed), `metadata`,
  `artifact`, `status`, `units`, `convention`, `point_count`, `note`, the parallel arrays
  `bond_lengths_angstrom`, `hf_energies_ha`, `fci_energies_ha`, `nuclear_repulsion_energies_ha`,
  `exact_energies_ha`, and the same rows as `points` (`bond_length_angstrom`, `hf_energy_ha`,
  `fci_energy_ha`, `nuclear_repulsion_ha`, `exact_energy_ha`).
- `exact_energies_ha` is computed here, not stored: each stored point's `pauli_terms` are
  assembled into a dense matrix and diagonalised, the nuclear repulsion is added, and the total
  is cross-checked against that point's stored `fci_energy_ha`. A disagreement beyond 1e-8 Ha
  means the file contradicts itself and answers `500` `simulation_error` instead of reporting an
  exact reference.
- **`n_qubits` is capped at 8** for that dense diagonalisation. A wider stored Hamiltonian
  answers `422 too_many_qubits`; the exact reference is refused rather than approximated, and no
  curve data is served alongside the error.
- A missing artifact answers `503 hamiltonian_unavailable` (the file name and the molecule label
  are named in the message), and an unknown `molecule` answers `404 unknown_molecule`. No example
  response is shown here because `lih_sto3g_active.json` has not been generated: this repository
  contains no LiH numbers, and none are invented for the documentation.

#### POST /api/quantum/simulate

- Body: `{ "molecule": "H2", "bond_length": number, "shots": null | number, "noise": number,
  "seed": number, "maxiter": number }`. `bond_length_angstrom` is accepted as an alias for
  `bond_length`. `seed` (0..4294967295), `maxiter` (1..10000) and `shots` (positive integer or
  null) are optional; bad values answer `422` with `invalid_seed` / `invalid_maxiter` /
  `invalid_shots`.
- **VQE is H2 only.** A `molecule` naming a registry molecule other than H2 answers `422`
  `vqe_not_supported_for_molecule` with no result fields, because this build has one fixed
  2-qubit Ry ansatz and will not pretend it fits another molecule:

```json
{
  "error": {
    "code": "vqe_not_supported_for_molecule",
    "message": "Only H2 is simulated with VQE in this build. Only the stored curve and the exact diagonalisation are available for LiH (active space)."
  }
}
```

  That molecule's data is reachable through `GET /api/quantum/h2-curve?molecule=lih` instead
  (stored curve plus an exact diagonalisation). A `molecule` that is not in the registry at all
  keeps this endpoint's existing `status: "placeholder"` answer rather than erroring.
- **A real VQE runs in this request**: a fixed `2 x [Ry on each qubit, then a CX ladder]` ansatz
  minimised by SciPy COBYLA (3 restarts, default 150 evaluations per restart) against the
  artifact's Pauli terms through a Qiskit estimator. `shots: null` (the default) selects the
  noiseless statevector estimator; a number runs the shot-based Aer estimator with that many
  shots per evaluation, and `seed` reaches it as `backend_options={"seed_simulator": ...}` as
  well as seeding the restarts, so a fixed seed reproduces a run exactly.
- The exact reference is a live NumPy diagonalisation of the same stored Pauli terms plus the
  nuclear repulsion, cross-checked against the artifact's stored `fci_energy_ha` (1e-8 Ha) before
  the run starts; `hf_energy_ha` and `fci_energy_ha` themselves are read from the artifact.
- `within_chemical_accuracy` (1.6 mHa) is reported only for noiseless runs: a shot-based
  estimate is a statistic, so it returns `null` there and the error is reported as measured.
- **The budget is spelled out:** `n_restarts`, `maxiter_per_restart` and
  `total_energy_evaluations` (equal to the backward-compatible `iterations`), so the run is read
  as "3 restarts × 150 evaluations", not an opaque total.
- **Variational guard, scoped to noiseless statevector runs:** if such a run reports an energy
  more than 1e-9 Ha below the exact ground energy, the answer is a `500` with code
  `variational_principle_violation` and **no result fields** — a physically impossible result is
  refused, not served. A shot run is **not** refused for sitting below exact on sampling noise;
  it reports that noise instead (next bullet).
- **No winner's curse in a shot run's headline number.** The optimiser picks its best restart by
  the *noiseless* expectation at each parameter set, but the reported `vqe_energy_ha` is a
  **fresh re-evaluation at the final optimal parameters**, not the minimum of the noisy
  optimisation trace. Its uncertainty is reported as `vqe_energy_std_ha` (in Ha) with
  `vqe_energy_uncertainty_method` naming the source: `estimator_variance_metadata` (the
  estimator's own variance metadata, `sqrt(variance / shots)`) or `seed_reevaluation_10x` (the
  mean and sample standard deviation of 10 re-evaluations with different sampler seeds). That
  is sampling error, not a hardware result, and the `note` says so.
- Every run carries `classical_simulation: true` and says in `note` that no quantum hardware was
  used and no quantum advantage is claimed.
- The stored grid has no point at every requested bond length, so an off-grid request snaps to
  the nearest stored point (never interpolated; an exact tie breaks toward the shorter bond
  length) and the response discloses both values: `bond_length_requested_angstrom` and
  `bond_length_used_angstrom`, with `bond_length_is_exact` telling them apart — plus the same
  three facts as `requested_bond_length`, `used_bond_length` and `snapped`. `bond_length`
  echoes the point the energies came from, because the panel badges it next to them.

Example request:
```http
POST /api/quantum/simulate HTTP/1.1
Host: 127.0.0.1:8000
Content-Type: application/json

{
  "molecule": "H2",
  "bond_length": 0.74,
  "shots": null,
  "noise": 0.0,
  "seed": 7
}
```

- `noise` is accepted for the contract and echoed back, but no noise model is applied and the
  response says so. The sampling noise that `shots` implies **is** real: the optimiser sees that
  many shots per energy evaluation.

Example response without the offline artifact (standard error shape, `503`):
```json
{
  "error": {
    "code": "hamiltonian_unavailable",
    "message": "h2_sto3g.json has not been generated, so no H2 energy is available. Run scripts/precompute_hamiltonians.ipynb in Colab or WSL and save the result into app/data/hamiltonians/. This endpoint reports the energy as unavailable rather than showing placeholder numbers."
  }
}
```

Example response for an on-grid request (0.74 Å is the stored default; `convergence_history_ha`
shows the first 4 of its 450 entries — every number is copied verbatim from the live response):
```json
{
  "molecule": "H2",
  "status": "computed",
  "units": "hartree",
  "bond_length": 0.74,
  "requested_bond_length": 0.74,
  "used_bond_length": 0.74,
  "snapped": false,
  "bond_length_requested_angstrom": 0.74,
  "bond_length_used_angstrom": 0.74,
  "bond_length_is_exact": true,
  "hf_energy_ha": -1.1167593073964255,
  "fci_energy_ha": -1.1372838344885028,
  "nuclear_repulsion_ha": 0.7151043390810812,
  "exact_energy_ha": -1.137283834488502,
  "vqe_energy_ha": -1.137283834480475,
  "error_vs_exact_ha": 8.026912468039882e-12,
  "abs_error_vs_exact_ha": 8.026912468039882e-12,
  "error_mha": 8.026912468039882e-09,
  "vqe_energy_std_ha": null,
  "vqe_energy_uncertainty_method": null,
  "within_chemical_accuracy": true,
  "chemical_accuracy_threshold_ha": 0.0016,
  "noiseless": true,
  "iterations": 450,
  "total_energy_evaluations": 450,
  "n_restarts": 3,
  "maxiter_per_restart": 150,
  "convergence_history_ha": [
    -1.0646672340673422,
    -1.0949480951674475,
    -1.035481833399952,
    -1.027491460831999
  ],
  "seed": 7,
  "seed_simulator_applied": false,
  "maxiter": 150,
  "restarts": 3,
  "ansatz": "2 x [Ry on each of 2 qubit(s), then a CX ladder]",
  "optimizer": "COBYLA",
  "shots": null,
  "noise": 0.0,
  "estimator_backend": "statevector",
  "estimator_class_name": "StatevectorEstimator",
  "classical_simulation": true,
  "note": "Real VQE run at 0.74 A: the requested 0.74 A is on the stored grid in h2_sto3g.json. 450 COBYLA energy evaluations over 3 restarts of a fixed 2 x [Ry on each of 2 qubit(s), then a CX ladder] ansatz, statevector estimator, seed 7 (noiseless statevector run). Classical simulation only: no quantum hardware was used, no quantum advantage is claimed. The exact reference diagonalises the same stored Pauli terms with NumPy; hf_energy_ha and fci_energy_ha are read from h2_sto3g.json. No noise model is applied."
}
```

Example response for an off-grid request (0.741 Å, snapped and disclosed; `convergence_history_ha`
truncated the same way):
```json
{
  "molecule": "H2",
  "status": "computed",
  "units": "hartree",
  "bond_length": 0.74,
  "requested_bond_length": 0.741,
  "used_bond_length": 0.74,
  "snapped": true,
  "bond_length_requested_angstrom": 0.741,
  "bond_length_used_angstrom": 0.74,
  "bond_length_is_exact": false,
  "hf_energy_ha": -1.1167593073964255,
  "fci_energy_ha": -1.1372838344885028,
  "nuclear_repulsion_ha": 0.7151043390810812,
  "exact_energy_ha": -1.137283834488502,
  "vqe_energy_ha": -1.137283834480475,
  "error_vs_exact_ha": 8.026912468039882e-12,
  "abs_error_vs_exact_ha": 8.026912468039882e-12,
  "error_mha": 8.026912468039882e-09,
  "vqe_energy_std_ha": null,
  "vqe_energy_uncertainty_method": null,
  "within_chemical_accuracy": true,
  "chemical_accuracy_threshold_ha": 0.0016,
  "noiseless": true,
  "iterations": 450,
  "total_energy_evaluations": 450,
  "n_restarts": 3,
  "maxiter_per_restart": 150,
  "convergence_history_ha": [
    -1.0646672340673422,
    -1.0949480951674475,
    -1.035481833399952,
    -1.027491460831999
  ],
  "seed": 7,
  "seed_simulator_applied": false,
  "maxiter": 150,
  "restarts": 3,
  "ansatz": "2 x [Ry on each of 2 qubit(s), then a CX ladder]",
  "optimizer": "COBYLA",
  "shots": null,
  "noise": 0.0,
  "estimator_backend": "statevector",
  "estimator_class_name": "StatevectorEstimator",
  "classical_simulation": true,
  "note": "0.741 A is not on the stored grid in h2_sto3g.json; answered with the nearest stored point at 0.74 A instead, never interpolated (an exact tie breaks toward the shorter bond length). 450 COBYLA energy evaluations over 3 restarts of a fixed 2 x [Ry on each of 2 qubit(s), then a CX ladder] ansatz, statevector estimator, seed 7 (noiseless statevector run). Classical simulation only: no quantum hardware was used, no quantum advantage is claimed. The exact reference diagonalises the same stored Pauli terms with NumPy; hf_energy_ha and fci_energy_ha are read from h2_sto3g.json. No noise model is applied."
}
```

Example response when the variational guard trips (captured through the real router with a
deliberately broken optimiser stub — the guard test in `tests/test_endpoints.py` produces it the
same way; a genuine run has never triggered it):
```json
{
  "error": {
    "code": "variational_principle_violation",
    "message": "the run at 0.74 A reported VQE energy -1.187283834489 Ha, below the exact -1.137283834489 Ha by 5.000e-02 Ha -- more than the 1e-09 Ha rounding allowance the variational principle allows. The result is refused instead of reported."
  }
}
```

#### GET /api/quantum/hadamard-demo

- Returns a single-qubit superposition demo.

Example request:
```http
GET /api/quantum/hadamard-demo HTTP/1.1
Host: 127.0.0.1:8000
```

Example response:
```json
{
  "circuit": "H|0>",
  "initial_state": [1.0, 0.0],
  "final_state": [0.7071, 0.7071],
  "probabilities": {
    "0": 0.5,
    "1": 0.5
  }
}
```

#### GET /api/quantum/resources?formula=&molecule=&active_orbitals=

- **What would this molecule need?** Counted from its formula for the STO-3G minimal basis. No
  Hamiltonian is built, nothing is diagonalised and no artifact is read, so this endpoint answers
  whether or not the Hamiltonian files are on disk — unlike the curve and simulate endpoints.
- Query parameters: `formula` (a plain molecular formula, e.g. `H2O` or `C9H8O4`), `molecule` (a
  registry key — `h2` or `lih`, read as the formula `H2` or `LiH`), and the optional
  `active_orbitals` (integer k, `1..n_spatial_orbitals`, which also estimates the k-orbital
  active-space problem). If both `formula` and `molecule` are given, `formula` wins and
  `requested.note` says so.
- **With no parameters it answers a usage description** instead of guessing a molecule: `status:
  "usage"`, the request shape, the supported element list and the per-atom STO-3G function
  counts.
- The counts, and the rule behind each one:
  - `n_spatial_orbitals`: one contracted basis function per occupied atomic shell — H and He 1,
    Li–Ne 5, Na–Ar 9 — summed over the formula.
  - `n_electrons`: the sum of the neutral atomic numbers. **No charge input is taken**, so this is
    the neutral count; spin and charge are not modelled (`notes.spin_and_charge`) and an odd count
    (open shell) additionally reports `notes.open_shell`.
  - `n_qubits_full = 2 * n_spatial_orbitals` (one qubit per spin orbital).
  - `n_qubits_after_parity_reduction = n_qubits_full - 2` for the parity mapping with two-qubit
    reduction (`parity_reduction_removed_qubits` is that fixed 2).
  - `statevector_memory_bytes = 2 ** n_qubits_after_parity_reduction * 16` (complex128) with
    `statevector_memory_human` as KB/MB/GB/TB, or `"more than 10^15 bytes"` above 1.126e15 bytes.
  - `hilbert_dimension` is returned **as a string**: 2**144 would lose precision as a JSON number
    in a browser, and the integer memory figure beside it has the same caveat for very large
    molecules.
- `feasibility` judges the reduced qubit count against `interactive_cap_qubits: 8` (this app's own
  limit — a **project choice**) and `laptop_statevector_ceiling_qubits: 25` (rough), giving a
  `verdict` of `simulable_in_app`, `precompute_only` or `beyond_classical_simulation` with a
  one-sentence `reason`. These are project limits and rough estimates, not measured hardware or
  algorithm bounds. `simulable_in_app` means the size fits this app's budget — **not** that a
  simulation exists for that molecule, because H2 is the only molecule with a VQE ansatz.
- `tags` maps each figure to `exact` (computed by a closed-form rule in code) or `estimated` (a
  judgement or a rough limit). The feasibility verdict, the laptop ceiling and everything in the
  active-space block are `estimated`.
- **No Pauli-term count and no circuit size is reported**, because a formula cannot give either.
  The only statement about term growth is the asymptotic "roughly N^4 in the number of orbitals",
  an order-of-magnitude note in `notes.pauli_term_scaling`.
- Errors, standard shape, `422`: `invalid_formula` (not element-symbol-and-count text),
  `unsupported_element` (names the supported symbols), `invalid_active_orbitals` (not an integer
  in `1..n_spatial_orbitals`) and `missing_formula` (`active_orbitals` with no subject).
- `POST /api/quantum/resources` remains as a body-shaped alias: `{ "formula": ... }`,
  `{ "molecule": ... }` (including the older `{"molecule": "C9H8O4"}`, read as a formula) and
  `{ "active_orbitals": k }` all answer through the same code path, and an empty body returns the
  same usage description as the parameterless GET.

Example request:
```http
GET /api/quantum/resources?formula=H2O HTTP/1.1
Host: 127.0.0.1:8000
```

Example response (live capture; H2O, 7 STO-3G spatial orbitals):
```json
{
  "status": "computed",
  "subject": {
    "input": "H2O",
    "formula": "H2O",
    "atom_counts": {
      "H": 2,
      "O": 1
    }
  },
  "basis": "STO-3G",
  "method": "computed_by_formula",
  "method_note": "Computed by formula for the STO-3G basis (no PySCF run, no simulation). The classical simulation limits below are project limits and rough estimates.",
  "n_spatial_orbitals": 7,
  "basis_functions_per_element": {
    "H": 2,
    "O": 5
  },
  "n_electrons": 10,
  "n_qubits_full": 14,
  "n_qubits_after_parity_reduction": 12,
  "parity_reduction_removed_qubits": 2,
  "statevector_memory_bytes": 65536,
  "statevector_memory_human": "64.0 KB",
  "hilbert_dimension": "4096",
  "statevector_basis": "n_qubits_after_parity_reduction",
  "bytes_per_amplitude": 16,
  "feasibility": {
    "based_on_qubits": 12,
    "interactive_cap_qubits": 8,
    "laptop_statevector_ceiling_qubits": 25,
    "verdict": "precompute_only",
    "reason": "12 qubits after the parity two-qubit reduction is above this app's 8-qubit interactive limit but inside a rough 25-qubit laptop statevector ceiling, so the numbers would have to be precomputed offline and served as an artifact.",
    "note": "These limits are project choices and rough estimates, not measured hardware or algorithm bounds."
  },
  "active_space": null,
  "tags": {
    "n_spatial_orbitals": "exact",
    "n_electrons": "exact",
    "n_qubits_full": "exact",
    "n_qubits_after_parity_reduction": "exact",
    "parity_reduction_removed_qubits": "exact",
    "statevector_memory_bytes": "exact",
    "statevector_memory_human": "exact",
    "hilbert_dimension": "exact",
    "feasibility.interactive_cap_qubits": "exact",
    "feasibility.laptop_statevector_ceiling_qubits": "estimated",
    "feasibility.verdict": "estimated",
    "notes.pauli_term_scaling": "estimated"
  },
  "notes": {
    "spin_and_charge": "Spin and charge are not modelled: the electron count is the sum of the neutral atomic numbers in the formula, and no charge input is taken.",
    "method": "Computed by formula for the STO-3G basis, not by running PySCF: no Hamiltonian was built and no diagonalisation was performed for this molecule.",
    "pauli_term_scaling": "Pauli term growth is roughly N^4 in the number of orbitals (estimated, order-of-magnitude). No specific term count is computed here, because a term count depends on the mapped Hamiltonian, which this estimator does not build.",
    "limits": "The classical simulation limits in 'feasibility' are project limits and rough estimates; they are not measured hardware or algorithm bounds.",
    "interactive_limit": "Above 8 qubits this app refuses a dense diagonalisation (too_many_qubits), and H2 is the only molecule with a VQE ansatz, so 'simulable_in_app' means the size fits the app's budget rather than that a simulation exists for this molecule.",
    "verification": "The STO-3G orbital counts are counted from the minimal-basis rule, not read from a basis file, and have not been cross-checked against PySCF on this host."
  },
  "verified_against_pyscf": false,
  "classical_simulation": true,
  "requested": {
    "formula": "H2O",
    "molecule": null,
    "resolved_formula": "H2O",
    "resolved_from": "formula"
  }
}
```

The five documented examples, as this endpoint returns them (live capture):

| `?formula=` | spatial orbitals | electrons | full qubits | reduced qubits | statevector memory | verdict |
| --- | --- | --- | --- | --- | --- | --- |
| `H2` | 2 | 2 | 4 | 2 | 64 bytes | `simulable_in_app` |
| `LiH` | 6 | 4 | 12 | 10 | 16.0 KB | `precompute_only` |
| `H2O` | 7 | 10 | 14 | 12 | 64.0 KB | `precompute_only` |
| `CH4` | 9 | 10 | 18 | 16 | 1.0 MB | `precompute_only` |
| `C9H8O4` (aspirin) | 73 | 94 | 146 | 144 | more than 10^15 bytes | `beyond_classical_simulation` |

Example request for the optional active-space block:
```http
GET /api/quantum/resources?formula=H2O&active_orbitals=4 HTTP/1.1
Host: 127.0.0.1:8000
```

Example response (the `active_space` block only; the full-molecule fields above are unchanged, and
the block carries six extra `estimated` tags keyed `active_space.*`):
```json
{
  "active_space": {
    "active_orbitals": 4,
    "n_spatial_orbitals": 7,
    "n_qubits": 8,
    "n_qubits_after_parity_reduction": 6,
    "feasibility": {
      "based_on_qubits": 6,
      "interactive_cap_qubits": 8,
      "laptop_statevector_ceiling_qubits": 25,
      "verdict": "simulable_in_app",
      "reason": "6 qubits after the parity two-qubit reduction is at or below this app's 8-qubit interactive limit, so a dense statevector calculation of this size fits the app's budget.",
      "note": "These limits are project choices and rough estimates, not measured hardware or algorithm bounds."
    },
    "note": "An active space is an approximation chosen by the user: every figure in this block describes the k-orbital problem, not the full molecule, and the frozen orbitals it leaves out change the energy.",
    "statevector_memory_bytes": 1024,
    "statevector_memory_human": "1.0 KB",
    "hilbert_dimension": "64",
    "statevector_basis": "n_qubits_after_parity_reduction",
    "bytes_per_amplitude": 16
  }
}
```

Example request with no parameters (a usage answer, not a guessed molecule):
```http
GET /api/quantum/resources HTTP/1.1
Host: 127.0.0.1:8000
```

Example response:
```json
{
  "status": "usage",
  "endpoint": "GET /api/quantum/resources",
  "message": "Pass ?formula= (for example H2O or C9H8O4) or ?molecule=h2|lih to estimate what that molecule would need. No molecule was assumed for this request.",
  "query_parameters": {
    "formula": "a plain molecular formula, e.g. H2O, CH4, C9H8O4",
    "molecule": "a known molecule key: h2, lih",
    "active_orbitals": "optional integer k: also estimate the k-orbital active-space problem"
  },
  "example": {
    "request": "?formula=H2O",
    "formula": "H2O"
  },
  "basis": "STO-3G",
  "supported_elements": [
    "H",
    "He",
    "Li",
    "Be",
    "B",
    "C",
    "N",
    "O",
    "F",
    "Ne",
    "Na",
    "Mg",
    "Al",
    "Si",
    "P",
    "S",
    "Cl",
    "Ar"
  ],
  "basis_functions_per_atom": {
    "H": 1,
    "He": 1,
    "Li": 5,
    "Be": 5,
    "B": 5,
    "C": 5,
    "N": 5,
    "O": 5,
    "F": 5,
    "Ne": 5,
    "Na": 9,
    "Mg": 9,
    "Al": 9,
    "Si": 9,
    "P": 9,
    "S": 9,
    "Cl": 9,
    "Ar": 9
  },
  "method": "computed_by_formula",
  "note": "Computed by formula for the STO-3G basis (no PySCF run, no simulation). The classical simulation limits this estimator reports are project limits and rough estimates."
}
```

Example response for `GET /api/quantum/resources?formula=Fe2O3` — standard error shape, `422`:
```json
{
  "error": {
    "code": "unsupported_element",
    "message": "No STO-3G basis-function count for element(s): Fe. Supported elements are H, He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, P, S, Cl, Ar."
  }
}
```

### Molecule

#### GET /api/molecule/search?q=

- Searches PubChem.
- Cached with a finite TTL.
- The deposited 3D conformer is requested first (PubChem's `?record_type=3d`). When PubChem
  has no 3D record for the CID (HTTP 404) or answers with a flat 2D record, the default record
  is fetched instead; the request never fails just because no 3D conformer exists.
- `sdf_dimension` reports what was actually returned: `"3D"` or `"2D"`. It is detected from
  the SDF itself (the `3D`/`2D` marker in the molfile header, confirmed against the z spread of
  the atom block), never assumed. A record whose header says 3D but whose coordinates are flat
  is reported as `"2D"`, so a flat record is never presented as a 3D structure.

Example request:
```http
GET /api/molecule/search?q=aspirin HTTP/1.1
Host: 127.0.0.1:8000
```

Example response:
```json
{
  "query": "aspirin",
  "pubchem_cid": 2244,
  "smiles": "CC(=O)Oc1ccccc1C(=O)O",
  "molecular_weight": 180.16,
  "formula": "C9H8O4",
  "sdf_dimension": "3D",
  "sdf_loaded": true,
  "note": "3D conformer from PubChem",
  "status": "computed"
}
```

#### POST /api/molecule/properties

- Body: SMILES or equivalent validated input.
- Returns RDKit-based descriptors and Lipinski check.

Example request:
```http
POST /api/molecule/properties HTTP/1.1
Host: 127.0.0.1:8000
Content-Type: application/json

{
  "smiles": "CC(=O)Oc1ccccc1C(=O)O"
}
```

Example response:
```json
{
  "smiles": "CC(=O)Oc1ccccc1C(=O)O",
  "molecular_weight": 180.16,
  "logp": 1.2,
  "h_bond_donors": 1,
  "h_bond_acceptors": 4,
  "rotatable_bonds": 3,
  "formula": "C9H8O4",
  "lipinski_rule_of_five": {
    "pass": true,
    "violations": 0
  },
  "status": "computed"
}
```

#### POST /api/molecule/simulate

- Body: `{ "smiles": "...", "num_conformers": number }`
- Returns a 3D conformer ensemble, molecular descriptors, and drug-likeness findings.
- `method` is a machine-readable code, not prose: `3d_mmff94`, `3d_uff`,
  `3d_unminimised`, `2d_fallback` or `unavailable`. The prose description is in `method_detail`.
- **Energies are relative.** `relative_energy_kcal_per_mol` is measured within one ensemble
  (the lowest minimum is `0.0`). `lowest_absolute_energy_kcal_per_mol` is returned for
  transparency with `absolute_energies_comparable: false`: an absolute force-field energy
  changes with the conformer pool, the embedding seed and the RDKit build, so it must not be
  compared between runs.
- **Populations are estimates** (`populations_are_estimates: true`). They treat the
  force-field energy as a free energy, so they ignore entropy and solvent.
- Every `findings` entry carries the published threshold and its `source`.

Example request:
```http
POST /api/molecule/simulate HTTP/1.1
Host: 127.0.0.1:8000
Content-Type: application/json

{
  "smiles": "CC(=O)Oc1ccccc1C(=O)O",
  "num_conformers": 8
}
```

Example response (truncated):
```json
{
  "smiles": "CC(=O)Oc1ccccc1C(=O)O",
  "status": "computed",
  "method": "3d_mmff94",
  "method_detail": "ETKDGv3 distance-geometry embedding of 8 conformers + MMFF94 energy minimisation. Energies are reported relative to the lowest minimum in this run; populations are Boltzmann estimates at 298.15 K.",
  "energy_basis": "relative_to_lowest_conformer_in_this_run",
  "ensemble": {
    "available": true,
    "force_field": "MMFF94",
    "conformer_count": 8,
    "energy_span_kcal_per_mol": 1.5,
    "lowest_absolute_energy_kcal_per_mol": 18.91,
    "absolute_energies_comparable": false,
    "populations_are_estimates": true
  },
  "findings": [
    {
      "name": "Lipinski rule of five",
      "verdict": "pass",
      "source": "Lipinski et al., Adv. Drug Deliv. Rev. 1997 (rule of five)"
    }
  ]
}
```

Example "looks like a PDB id" response — the endpoint refuses to guess:
```json
{
  "smiles": "1CRN",
  "status": "placeholder",
  "note": "'1CRN' looks like a PDB structure identifier rather than a SMILES string. Load it on the Protein Viewer tab, or open the Drug-Target Studio to see the drug bound to its target."
}
```

### Protein

#### GET /api/protein/{pdb_id}

- Fetches PDB metadata and structure reference from RCSB.
- PDB ID format is validated before lookup.

Example request:
```http
GET /api/protein/6LU7 HTTP/1.1
Host: 127.0.0.1:8000
```

Example response:
```json
{
  "pdb_id": "6LU7",
  "title": "SARS-CoV-2 main protease",
  "resolution": 2.16,
  "organism": "SARS-CoV-2",
  "chains": ["A"],
  "ligands": [],
  "structure_file_reference": "6LU7.pdb",
  "status": "computed"
}
```

Example placeholder response:
```json
{
  "pdb_id": "1ABC",
  "status": "placeholder"
}
```

### Complexes

#### GET /api/complexes/list

- Returns curated complexes from the manifest.

Example request:
```http
GET /api/complexes/list HTTP/1.1
Host: 127.0.0.1:8000
```

Example response:
```json
{
  "complexes": [        {
      "id": "hiv-protease-indinavir",
      "protein_pdb": "1HSG",
      "ligand_name": "Indinavir",
      "ligand_code": "MK1",
      "target_origin": "viral",
      "binding_site_type": "catalytic",
      "status": "computed"
    }
  ]
}
```

`target_origin` is `viral` or `human` and exists so a human drug target is never presented as a
virus. `binding_site_type` is `catalytic` or `allosteric`: an allosteric inhibitor is expected
to miss the catalytic residues, and the mechanism reports that rather than calling it a failure.

#### GET /api/complexes/{id}/structure

- Downloads the deposited coordinates, keeps protein plus the selected ligand, and returns
  the measured contacts and the resulting mechanism of action.
- `contact_summary` counts **unique residues**, not atom pairs (`contact_summary_basis`).
  `atom_pair_counts` is returned separately for transparency.
- `contact_rules` gives the exact rule behind every label (donor/acceptor plus an angle check
  for a hydrogen bond; charged groups for a salt bridge; C-C only for hydrophobic; anything
  else within the cutoff is a plain `close_contact`).
- `contact_labels` states, per label, whether the number is an inference or a direct
  measurement: `hydrogen_bond` and `salt_bridge` are `"geometric heuristic"` (the deposited
  file has no hydrogens and no protonation states, so identity is inferred from element and
  atom name), while `hydrophobic` and `close_contact` are `"measured geometry"`.
- **`contact_summary.hydrogen_bond` and `contact_summary.salt_bridge` are UPPER BOUNDS, not
  exact counts.** Every ligand N/O atom is treated as *both* a possible donor and a possible
  acceptor, because the deposited file carries no hydrogens and no protonation states. The set
  of genuine donors and acceptors is a subset of that permissive set, so the rule can only
  over-count; a residue reported here may not form that interaction in reality. A caller that
  needs a definite count must not read these as measurements.
- Salt bridges are evaluated before hydrogen bonds and the first match wins, so a residue can
  be reported under both types. The per-type lists **overlap** and must not be summed;
  `atom_pair_counts` is the separate, non-overlapping atom-pair tally.
- `near_miss_hydrogen_bonds` lists donor/acceptor pairs whose geometry is acceptable but whose
  distance exceeds the cutoff. This is what explains a hydrogen-bond count of 0: it is reported
  as a near miss rather than silently dropped.
- Each mechanism step is tagged `computed` or `literature`. The two are never conflated, and
  no docking score or binding free energy is claimed.
- Each residue's `atom_pairs` list is capped at 12 entries for response size. Every residue
  carries `atom_pairs_truncated` (`true` when pairs were dropped), and the top level carries
  `atom_pairs_truncated` plus `atom_pairs_truncated_note`. `atom_pair_counts` is **never**
  truncated and always holds the complete tally, so a capped list is never mistaken for the
  full set of contacts.

Example request:
```http
GET /api/complexes/hiv-protease-indinavir/structure HTTP/1.1
Host: 127.0.0.1:8000
```

Example response (truncated):
```json
{
  "complex_id": "hiv-protease-indinavir",
  "protein_pdb": "1HSG",
  "ligand_code": "MK1",
  "target_origin": "viral",
  "binding_site_type": "catalytic",
  "status": "computed",
  "binding_site_residue_count": 30,
  "contact_summary_basis": "unique_residues",
  "contact_summary": {
    "salt_bridge": 1,
    "hydrogen_bond": 5,
    "hydrophobic": 20,
    "close_contact": 28
  },
  "contact_rules": {
    "hydrogen_bond": "Donor/acceptor heavy-atom pair (N or O on both sides) with D...A <= 3.5 A AND the antecedent-D...A angle >= 90 deg."
  },
  "contact_labels": {
    "hydrogen_bond": "geometric heuristic",
    "salt_bridge": "geometric heuristic",
    "hydrophobic": "measured geometry",
    "close_contact": "measured geometry"
  },
  "atom_pairs_truncated": false,
  "atom_pairs_truncated_note": "Per-residue atom_pairs are capped at 12 per residue for response size; atom_pair_counts is the complete tally and is not truncated.",
  "near_miss_hydrogen_bonds": [
    {
      "protein_residue": "ASP B29",
      "protein_atom": "OD1",
      "ligand_atom": "O4",
      "donor_side": "ligand",
      "distance_angstrom": 3.76,
      "angle_deg": 143.9,
      "cutoff_angstrom": 3.5,
      "reason": "Donor/acceptor geometry is acceptable, but the distance is over the 3.5 A hydrogen-bond cutoff, so this is not counted as a hydrogen bond."
    }
  ],
  "mechanism": {
    "status": "computed",
    "catalytic_residues_occluded": 2,
    "catalytic_residues_curated": 2
  }
}
```

Example placeholder response when a curated entry has no deposited coordinates available:
```json
{
  "complex_id": "hiv-protease-indinavir",
  "status": "placeholder",
  "note": "Could not download the deposited coordinates for 1HSG (RCSB PDB unreachable or the entry has no PDB-format file)."
}
```

#### GET /api/complexes/{id}/interactions

- Returns interaction data if computed.
- Each of `hydrogen_bonds`, `hydrophobic_contacts`, `salt_bridges` and `close_contacts` holds
  **unique residues**, not atom pairs; each residue carries its own `atom_pairs`. The counts
  therefore match `contact_summary` in `GET /api/complexes/{id}/structure`.
- Each residue's `atom_pairs` is capped at 12 (see `GET /api/complexes/{id}/structure`). The
  per-residue `atom_pairs_truncated` flag plus the top-level `atom_pairs_truncated` and
  `atom_pairs_truncated_note` state whether the cap applied; `atom_pair_counts` is complete.
- A residue that makes more than one kind of contact appears in more than one list, so the
  lists must not be summed as if they were disjoint.
- `contact_labels` marks `hydrogen_bond` and `salt_bridge` as `"geometric heuristic"` and
  `hydrophobic` / `close_contact` as `"measured geometry"`, so a caller can tell an inference
  from a measurement without reading the prose in `method`.
- `hydrogen_bonds` and `salt_bridges` are **upper bounds**: every ligand N/O counts as both a
  possible donor and a possible acceptor, so the true set is a subset of what is listed. The
  lists overlap (a residue can hold a salt bridge and a hydrogen bond), so their lengths do not
  add up to `binding_site_residue_count`.
- `near_miss_hydrogen_bonds` carries the donor/acceptor pairs rejected only by the distance
  cutoff, which is how a zero hydrogen-bond count is explained rather than left bare.
- When nothing was computed, the interaction lists are `null` (**not** `[]`) with
  `"status": "placeholder"` and a `note` saying why.

Example request:
```http
GET /api/complexes/hiv-protease-indinavir/interactions HTTP/1.1
Host: 127.0.0.1:8000
```

Example response when there is no curated entry to work from (note the `null`s):
```json
{
  "complex_id": "1VRT",
  "hydrogen_bonds": null,
  "hydrophobic_contacts": null,
  "salt_bridges": null,
  "interacting_residues": [],
  "status": "placeholder",
  "note": "No curated entry for complex '1VRT'"
}
```

Example computed response (truncated to one residue and one near miss):
```json
{
  "complex_id": "hiv-protease-indinavir",
  "ligand_code": "MK1",
  "status": "computed",
  "contact_summary_basis": "unique_residues",
  "contact_summary": {
    "salt_bridge": 1,
    "hydrogen_bond": 5,
    "hydrophobic": 20,
    "close_contact": 28
  },
  "contact_labels": {
    "hydrogen_bond": "geometric heuristic",
    "salt_bridge": "geometric heuristic",
    "hydrophobic": "measured geometry",
    "close_contact": "measured geometry"
  },
  "hydrogen_bonds": [
    {
      "residue": "ASP",
      "chain": "B",
      "residue_number": "25",
      "residue_label": "ASP B25",
      "types": ["salt_bridge", "hydrogen_bond", "close_contact"],
      "contact_count": 11,
      "min_distance_angstrom": 2.63,
      "atom_pairs": [
        {"protein_atom": "OD1", "ligand_atom": "O2", "distance_angstrom": 2.63, "type": "hydrogen_bond"},
        {"protein_atom": "OD2", "ligand_atom": "N3", "distance_angstrom": 4.0, "type": "salt_bridge"}
      ]
    }
  ],
  "hydrophobic_contacts": [],
  "salt_bridges": [],
  "close_contacts": [],
  "near_miss_hydrogen_bonds": [
    {
      "protein_residue": "ASP B29",
      "protein_atom": "OD1",
      "ligand_atom": "O4",
      "donor_side": "ligand",
      "distance_angstrom": 3.76,
      "angle_deg": 143.9,
      "cutoff_angstrom": 3.5,
      "reason": "Donor/acceptor geometry is acceptable, but the distance is over the 3.5 A hydrogen-bond cutoff, so this is not counted as a hydrogen bond."
    }
  ]
}
```

#### GET /api/complexes/{id}/docking-validation

- **No docking is performed in this build**, for any complex, so this endpoint never returns an
  RMSD. It always reports `"status": "placeholder"` with `rmsd_angstrom: null` and a `note`
  explaining that the measured crystal contacts are reported instead of a docking score.
- There is deliberately no `"status": "computed"` example here: an RMSD of 0.0 or any other
  number would be invented, and the endpoint has no code path that produces a real one.

Example request:
```http
GET /api/complexes/hiv-protease-indinavir/docking-validation HTTP/1.1
Host: 127.0.0.1:8000
```

Example response (this is the only shape it returns):
```json
{
  "complex_id": "hiv-protease-indinavir",
  "crystal_pose_file": null,
  "docked_pose_file": null,
  "rmsd_angstrom": null,
  "status": "placeholder",
  "note": "No docking run has been performed for this complex, so there is no RMSD to report. The Drug-Target Studio instead reports measured contacts from the deposited crystal coordinates, which require no docking."
}
```

### Jobs

#### GET /api/jobs/{id}

- Returns state: `queued | running | done | failed` (SPEC section 5).
- Long-running work uses BackgroundTasks and the SQLite job store.
- **No HTTP route starts a job.** Creating a job is internal-only (`enqueue_job` in
  `app/services/jobs_service.py`); this build exposes job *status* only. A state outside the
  documented enum is rejected by the internal enqueue path rather than stored, so an unknown
  state can never appear in this response.

Example request:
```http
GET /api/jobs/abc-123 HTTP/1.1
Host: 127.0.0.1:8000
```

Example response:
```json
{
  "job_id": "abc-123",
  "state": "done",
  "created_at": "2026-10-08T12:00:00Z",
  "updated_at": "2026-10-08T12:01:00Z",
  "result": {
    "molecule": "H2",
    "error_mha": 1.0
  },
  "error": null
}
```

Example failed response:
```json
{
  "job_id": "abc-123",
  "state": "failed",
  "created_at": "2026-10-08T12:00:00Z",
  "updated_at": "2026-10-08T12:00:30Z",
  "result": null,
  "error": {
    "code": "simulation_error",
    "message": "Simulation failed"
  }
}
```

### AI

#### POST /api/ai/resolve

- Maps text to verified identifiers when possible.
- Returns "no match" honestly when it cannot resolve.

Example request:
```http
POST /api/ai/resolve HTTP/1.1
Host: 127.0.0.1:8000
Content-Type: application/json

{
  "text": "covid protease"
}
```

Example response:
```json
{
  "query": "covid protease",
  "resolved": {
    "pdb_id": "6LU7"
  },
  "method": "alias_lookup",
  "status": "resolved"
}
```

Example unresolved response:
```json
{
  "query": "blahblah",
  "resolved": null,
  "method": "alias_lookup",
  "status": "no_match"
}
```

#### POST /api/ai/predict-properties

- Optional endpoint, built last.
- If no trained model is loaded, returns `method: "heuristic"` or `model: "not trained"`.

Example request:
```http
POST /api/ai/predict-properties HTTP/1.1
Host: 127.0.0.1:8000
Content-Type: application/json

{
  "smiles": "CC(=O)Oc1ccccc1C(=O)O"
}
```

Example "not trained" response:
```json
{
  "smiles": "CC(=O)Oc1ccccc1C(=O)O",
  "model": "not trained",
  "method": "heuristic",
  "predictions": null,
  "status": "not_trained"
}
```

Example trained response:
```json
{
  "smiles": "CC(=O)Oc1ccccc1C(=O)O",
  "model": "esol-rf-v1",
  "method": "trained_model",
  "predictions": {
    "water_solubility_logS": -1.5
  },
  "applicability_warning": null,
  "status": "computed"
}
```

## Input constraints

- SMILES length is bounded.
- Query length is bounded.
- Molecule inputs are normalized before use.
- PDB IDs are validated before lookup.

## Data rules

- No fabricated numbers.
- No fake AI.
- Complexes use `status: "placeholder"` for missing data.
- Docking-validation demos should use non-covalent complexes where possible.
- External API results are cached with per-source TTLs.
- Heavy chemistry is computed offline and loaded from files.

## Versioning

If the contract changes, update:
- this file,
- SPEC.md,
- and the relevant router tests.
