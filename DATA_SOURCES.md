# Data Sources and Provenance

This file tracks where each piece of data or method comes from, what license applies,
and what paper or resource backs it. It is meant to stay current as the project grows.

## External APIs

- **PubChem**
  - Use: molecule search and compound metadata.
  - License: consult PubChem's current terms before commercial or large-scale use.
  - Notes: results are cached with a finite TTL. Input is normalized before lookup.

- **RCSB PDB**
  - Use: protein structure fetch and metadata.
  - License: consult RCSB's current usage policy.
  - Notes: PDB ID format is validated before download. Structures are cached with a finite TTL.

## Quantum

- **Qiskit and Qiskit Aer**
  - Use: Hamiltonian simulation, statevector estimation, noisy estimation, and VQE loop.
  - Notes: exact estimator import paths are confirmed against the pinned version before use.
  - Hamiltonians for the demo molecules are precomputed offline into JSON and loaded by the API.

- **Scientific reference point**
  - H2 in STO-3G near 0.74 angstrom is about -1.137 Ha. This is used as a sanity check for the exact solver and VQE loop.

## Molecules and properties

- **RDKit**
  - Use: SMILES normalization and validation, molecular descriptors, Lipinski rule-of-five computation.
  - Notes: descriptors are computed by us from the molecule, not fetched from an external "truth" service.

## Drug-target complexes

- Docking validation and interaction tables are served from curated artifact files.
- Every curated entry has an explicit status: `computed` or `placeholder`.
- Placeholder entries use `null` values and must not present invented RMSD or interaction data.
- Docking-validation demos should use non-covalent complexes where possible.

## AI and property prediction

- Optional module. Built last.
- If present, trained models are produced by us and stored with their preprocessing pipeline and manifest.
- The API loads only model files it produced itself.
- If no model is loaded, the endpoint returns an honest "not trained" or "heuristic" response.

## Licensing and attribution

- Keep this file updated whenever a new dataset, API, or method is added.
- Record the license, source URL, and the paper/resource behind each method.
- Do not paste external terms here from memory if they may have changed; verify current terms before relying on them.
