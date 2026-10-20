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

## Deployment

One container runs the API and the frontend together. The `Dockerfile` and `.dockerignore` in
the repository root build it, and [docs/DEPLOYMENT.md](docs/DEPLOYMENT.md) explains why a
container (Render or a Hugging Face Space) is the supported shape and why a serverless platform
such as Vercel is not. PySCF is not in the image: the Hamiltonians are offline artifacts.

## Generating the Hamiltonian artifacts

The Quantum Lab's H2 energy curve and its single-point simulation are both built from one
offline artifact, `app/data/hamiltonians/h2_sto3g.json`, which **is** committed (58 points,
validated by `scripts/validate_hamiltonian_file.py`). The curve endpoint only reads it. The
simulate endpoint runs a real VQE — a small Ry ansatz minimised with SciPy COBYLA through a
Qiskit estimator — against the artifact's stored Pauli terms and compares it with an exact NumPy
diagonalisation of the same operator. It is a classical simulation (no quantum hardware, no
noise model), a fixed seed makes it reproducible, and `python scripts/measure_h2_vqe.py`
reproduces the measurements listed in RESULTS.md. The resource estimator panel answers a different
question — *what would this molecule need?* — and answers it by counting from the molecular
formula in the STO-3G minimal basis (`GET /api/quantum/resources?formula=`, with `POST` kept as a
body-shaped alias). It reads no artifact and runs no simulation, it tags every figure `exact` or
`estimated`, and it reports no Pauli-term or circuit-size number at all, because a formula cannot
give one.

PySCF publishes no Windows wheel, so regenerating the artifact has to happen where PySCF works
— WSL or Google Colab:

1. Open `scripts/precompute_hamiltonians.ipynb` in Colab (or run
   `python scripts/precompute_hamiltonians.py` inside WSL).
2. Install the pinned set the notebook lists (`pyscf`, `qiskit-nature` and their pins).
3. Run the notebook end to end. It writes `h2_sto3g.json`.
4. Save that file into `app/data/hamiltonians/` in this repository (overwriting the committed
   one) and run `python scripts/validate_hamiltonian_file.py h2` on it.

LiH uses the same script with a different molecule key:
`python scripts/precompute_hamiltonians.py --molecule lih` writes
`lih_sto3g_active.json` (1.0–3.0 Å in 0.1 steps plus 1.595 Å), freezing the Li 1s core and
removing the two highest virtual spatial orbitals. **It is not in this repository and no LiH
numbers exist anywhere in this build.** Once you generate it and save it as
`app/data/hamiltonians/lih_sto3g_active.json`, `GET /api/quantum/molecules` reports it as
available and `GET /api/quantum/h2-curve?molecule=lih` serves its stored curve plus an exact
diagonalisation reference. VQE still runs for H2 only.
5. Restart the server and open the Quantum Lab tab — “Refresh curve” will redraw the curve.

The viewer panels load `3Dmol-min.js` from `https://3dmol.csb.pitt.edu`, so the 3D viewers need
that host to be reachable; the API itself never makes a network call at request time for the
quantum panels.

## Docking validation (offline redocking, Colab)

`scripts/redock_validation.py` redocks known drug-protein complexes with AutoDock Vina and
measures how close every docked pose is to the deposited crystal pose. Vina and meeko publish no
Windows wheel, so it runs where they do — Google Colab (or Linux):

1. Open `scripts/redock_validation.ipynb` in Colab.
2. Run the install cell: `pip install "requests>=2.31" "numpy>=1.26" "rdkit" "vina==1.2.7" "meeko==0.8.0"`.
3. Run the notebook end to end. It writes `redock_validation.json` for the three default targets
   (`1VRT`, `2HU4`, `1HSG`) with **one fixed setting for all of them**
   (`--exhaustiveness 16 --seed 7`, `num_modes 9`, a 22 A box on the crystal ligand centroid).
4. Save that file as `app/data/docking/redock_validation.json` in this repository and check it with
   `python scripts/validate_docking_file.py` (pass another path as its argument if you saved it
   elsewhere). The notebook's last cell downloads the file and prints the same command.

The ligand is identified from the RCSB metadata — the single non-polymer entity with at least 10
heavy atoms, so 1VRT's magnesium ion is excluded rather than guessed away — and a target whose
ligand cannot be identified with certainty is recorded as failed instead of docked. Bond orders
come from the RCSB chemical component definition, and the RMSD is heavy-atom, symmetry-aware and
computed **without re-alignment**.

What this does **not** show: the box is centred on the true crystal site, so redocking is an
**easy** test, and a Vina score is a **prediction, not measured binding** — not a Kd, Ki or IC50.
Results cover the three complexes only. Nothing is retuned per target and no run is repeated
until it passes: a target that fails is written out as `status: "failed"` with its reason.
**The result file is not in this repository — the pipeline has not been run yet** (see
LIMITATIONS.md), and Part B is what will produce it and later serve it to the app.

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
