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

#### GET /api/quantum/h2-curve

- Returns VQE and exact energies across bond lengths, built from the offline artifact
  `app/data/hamiltonians/h2_sto3g.json`.
- **That artifact is not in the repository yet.** PySCF publishes no Windows wheel, so it is
  generated with `scripts/precompute_hamiltonians.ipynb` in Colab or WSL.
- An earlier build served a hand-typed 4x4 matrix from `app/data/hamiltonians/` and labelled it
  `"status": "computed"`. Those files have been deleted, and the endpoint now reports the curve
  as unavailable rather than presenting placeholder numbers as computed physics.

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

Example response once the artifact exists, before the curve is wired to its format:
```json
{
  "molecule": "H2",
  "basis": "sto-3g",
  "artifact": "h2_sto3g.json",
  "bond_lengths_angstrom": [],
  "exact_energies_ha": [],
  "vqe_energies_ha": [],
  "vqe_error_mha": [],
  "status": "placeholder",
  "note": "h2_sto3g.json is present but this endpoint is not yet wired to its format, so the H2 curve cannot be drawn in this build."
}
```

#### POST /api/quantum/simulate

- Body: `{ "molecule": "...", "bond_length": number, "shots": number, "noise": number }`
- Returns energy, exact energy, error in mHa, and convergence history.
- Long runs may return `202` with a job id.

Example request:
```http
POST /api/quantum/simulate HTTP/1.1
Host: 127.0.0.1:8000
Content-Type: application/json

{
  "molecule": "H2",
  "bond_length": 0.74,
  "shots": 1024,
  "noise": 0.0
}
```

- No shot or noise simulation is performed in this build. `shots` and `noise` are accepted for the
  contract and echoed back, but the response says so rather than implying a noisy run happened.

Example response without the offline artifact (standard error shape, `503`):
```json
{
  "error": {
    "code": "hamiltonian_unavailable",
    "message": "h2_sto3g.json has not been generated, so no H2 curve is available. Run scripts/precompute_hamiltonians.ipynb in Colab or WSL and save the result into app/data/hamiltonians/. This endpoint reports the curve as unavailable rather than showing placeholder numbers."
  }
}
```

Example response once the artifact exists but the curve is not yet wired to it:
```json
{
  "molecule": "H2",
  "bond_length": 0.74,
  "exact_energy_ha": null,
  "vqe_energy_ha": null,
  "error_mha": null,
  "convergence_history_ha": null,
  "shots": 1024,
  "noise": 0.0,
  "estimator_class_name": null,
  "status": "placeholder",
  "note": "h2_sto3g.json is present but this endpoint is not yet wired to its format, so the H2 curve cannot be drawn in this build."
}
```

Example queued response:
```json
{
  "job_id": "abc-123",
  "status": "queued"
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

#### POST /api/quantum/resources

- Body: `{ "molecule": "..." }`
- Returns qubits, Hamiltonian terms, and a rough circuit size estimate.

Example request:
```http
POST /api/quantum/resources HTTP/1.1
Host: 127.0.0.1:8000
Content-Type: application/json

{
  "molecule": "H2"
}
```

The H2 numbers were previously derived from the deleted placeholder matrix, so they are withheld
until the real artifact is wired in.

Example response:
```json
{
  "molecule": "H2",
  "n_qubits": null,
  "hamiltonian_terms": null,
  "estimated_circuit_gates": null,
  "estimated_circuit_depth": null,
  "artifact": "h2_sto3g.json",
  "status": "placeholder",
  "note": "h2_sto3g.json has not been generated, so no H2 curve is available. Run scripts/precompute_hamiltonians.ipynb in Colab or WSL and save the result into app/data/hamiltonians/. This endpoint reports the curve as unavailable rather than showing placeholder numbers."
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
