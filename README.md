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

## Generating the H2 Hamiltonian artifact

The Quantum Lab's H2 energy curve and its resource estimate are both built from one offline
artifact, `app/data/hamiltonians/h2_sto3g.json`. That file is **not** committed, so until you
generate it the app reports the curve and the estimate as *not available* rather than showing
placeholder numbers.

PySCF publishes no Windows wheel, so the precompute has to run where PySCF works — WSL or
Google Colab:

1. Open `scripts/precompute_hamiltonians.ipynb` in Colab (or run
   `python scripts/precompute_hamiltonians.py` inside WSL).
2. Install the pinned set the notebook lists (`pyscf`, `qiskit-nature` and their pins).
3. Run the notebook end to end. It writes `h2_sto3g.json`.
4. Save that file into `app/data/hamiltonians/` in this repository.
5. Restart the server and open the Quantum Lab tab — “Refresh curve” will draw the curve.

The viewer panels load `3Dmol-min.js` from `https://3dmol.csb.pitt.edu`, so the 3D viewers need
that host to be reachable; the API itself never makes a network call at request time for the
quantum panels.

## Important rules

- No fabricated data.
- No fake AI.
- No secrets in the repo.
- Heavy chemistry tooling is installed offline when needed.
- Free hosting may wipe local disk; treat cache and job store as disposable unless persisted.

## Reading the contact numbers

The complex analysis (`/api/complexes/{id}/structure` and `/interactions`) classifies every
protein/ligand atom pair within 4.5 Å from the deposited crystal coordinates. Two of the four
labels are measurements and two are inferences, and the API says which is which via
`contact_labels`:

- `hydrophobic` and `close_contact` are **measured geometry** — direct heavy-atom distances.
- `hydrogen_bond` and `salt_bridge` are a **geometric heuristic**. The deposited files contain
  no hydrogens and no protonation states, so hydrogen bonding is inferred from element and atom
  name, and a ligand O/N is treated as a possible anion/cation.

**Hydrogen-bond and salt-bridge counts are upper bounds, not exact counts.** Every ligand N/O
atom is treated as *both* a possible donor and a possible acceptor; the set of genuine donors and
acceptors is a subset of that, so the rule can only over-count. A pair that fails the 3.5 Å
cutoff while having acceptable geometry is reported separately in `near_miss_hydrogen_bonds`
rather than being silently dropped — that is what explains a count of 0.

Salt bridges are evaluated before hydrogen bonds and the first match wins, so a residue can be
listed under both types. The per-type lists are therefore **overlapping**, not disjoint: the
counts must not be summed, and `contact_summary` counts unique residues rather than atom pairs
(`contact_summary_basis`).

## Build order

1. Repo hygiene, `/health`, and a deployable skeleton.
2. Offline Hamiltonian precompute, then VQE core and tests.
3. Quantum endpoints; frontend can start here against mocks.
4. LiH, active-space H2O, and resource estimator.
5. PubChem and RCSB lookups with validation and cache.
6. Real complexes with manifest and placeholder handling.
7. Jobs and resilience.
8. Optional AI: ESOL model first.
