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

- Returns VQE and exact energies across bond lengths.
- Hamiltonians are loaded from precomputed JSON, not computed live.

Example request:
```http
GET /api/quantum/h2-curve HTTP/1.1
Host: 127.0.0.1:8000
```

Example response:
```json
{
  "molecule": "H2",
  "basis": "sto-3g",
  "bond_lengths_angstrom": [0.2, 0.74, 2.5],
  "exact_energies_ha": [-1.0, -1.137, -0.8],
  "vqe_energies_ha": [-1.0, -1.136, -0.8],
  "vqe_error_mha": [0.0, 1.0, 5.0],
  "status": "computed"
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

Example response:
```json
{
  "molecule": "H2",
  "bond_length": 0.74,
  "exact_energy_ha": -1.137,
  "vqe_energy_ha": -1.136,
  "error_mha": 1.0,
  "convergence_history_ha": [-1.0, -1.1, -1.13, -1.136],
  "status": "computed"
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

Example response:
```json
{
  "molecule": "H2",
  "qubits": 2,
  "hamiltonian_terms": 4,
  "estimated_circuit_gates": 12,
  "estimated_circuit_depth": 6,
  "note": "Estimate only, not a full compilation"
}
```

### Molecule

#### GET /api/molecule/search?q=

- Searches PubChem.
- Cached with a finite TTL.

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
  "complexes": [
    {
      "id": "hiv-protease-indinavir",
      "protein_pdb": "1HSG",
      "ligand_name": "Indinavir",
      "status": "placeholder"
    }
  ]
}
```

#### GET /api/complexes/{id}/interactions

- Returns interaction data if computed.
- Missing data is `null` with `"status": "placeholder"`.

Example request:
```http
GET /api/complexes/hiv-protease-indinavir/interactions HTTP/1.1
Host: 127.0.0.1:8000
```

Example response:
```json
{
  "complex_id": "hiv-protease-indinavir",
  "hydrogen_bonds": [],
  "hydrophobic_contacts": [],
  "salt_bridges": [],
  "status": "placeholder"
}
```

Example computed response:
```json
{
  "complex_id": "hiv-protease-indinavir",
  "hydrogen_bonds": [
    {
      "donor_atom": "A:NE2",
      "acceptor_atom": "LIG:O2",
      "distance_angstrom": 2.9
    }
  ],
  "hydrophobic_contacts": [],
  "salt_bridges": [],
  "status": "computed"
}
```

#### GET /api/complexes/{id}/docking-validation

- Returns docking validation data if computed.
- Never returns invented RMSD or interaction values.

Example request:
```http
GET /api/complexes/hiv-protease-indinavir/docking-validation HTTP/1.1
Host: 127.0.0.1:8000
```

Example response:
```json
{
  "complex_id": "hiv-protease-indinavir",
  "crystal_pose_file": "1HSG_ligand.mol2",
  "docked_pose_file": "1HSG_ligand_docked.mol2",
  "rmsd_angstrom": 1.2,
  "status": "computed"
}
```

Example placeholder response:
```json
{
  "complex_id": "hiv-protease-indinavir",
  "status": "placeholder"
}
```

### Jobs

#### GET /api/jobs/{id}

- Returns state: `queued | running | done | failed`.
- Long-running work uses BackgroundTasks and the SQLite job store.

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
