# Deployment

How to run this project on a host. Docs only: no credentials, tokens or secrets
belong in this file, in the repository, or in the image.

## Recommendation in one line

**One Docker container built from the `Dockerfile` in the repository root**,
running `uvicorn app.main:app`. That single process serves the JSON API *and* the
frontend in `frontend/static/` (`GET /`, `/index.html`, `/app/*` in
`app/main.py`). Run it on a container host such as **Render** or a **Hugging Face
Space**, both of which have free tiers that fit a study tool like this one.

```bash
docker build -t qdde .
docker run --rm -p 8000:8000 --env-file .env qdde
# then: http://127.0.0.1:8000/health   and   http://127.0.0.1:8000/
```

`Dockerfile` and `.dockerignore` live in the repository root. The image is built
from `python:3.11-slim`, installs `requirements.txt`, copies `app/`, `frontend/`
and the top-level Markdown files, runs as a non-root user, exposes the app's port
and defines a `HEALTHCHECK` against `/health`.

## Why not Vercel serverless

**Vercel's serverless functions are not a suitable backend for this app.** Two
independent reasons, either of which is sufficient.

1. **Dependency and request limits.** The backend imports RDKit (descriptors and
   3D conformers) and Qiskit plus Qiskit Aer (the VQE path, including the
   shot-based Aer estimator). Those packages alone, measured in this checkout's
   virtual environment with `du`-style recursive sizing of
   `.venv/Lib/site-packages`, are:

   | package | installed size |
   |---|---|
   | `scipy` | ~91 MiB |
   | `qiskit-aer` | ~40 MiB |
   | `qiskit` | ~26 MiB |
   | `rdkit` | ~25 MiB |
   | `numpy` | ~22 MiB |
   | everything else in `requirements.txt` | ~220 MiB |
   | **whole site-packages tree** | **~423 MiB** |

   Older project notes said "RDKit and Qiskit exceed its 50 MB limit". The
   *conclusion* still holds, but that number is out of date: Vercel's current
   documented standard is a **250 MB uncompressed** function bundle
   (~50 MB compressed, or **500 MB uncompressed for Python**), with larger
   bundles only on Fluid Compute — see the
   [Vercel knowledge-base article](https://vercel.com/kb/guide/troubleshooting-function-250mb-limit),
   checked 2026-10-13. So the accurate statement is: the measured tree is within
   ~15 % of the Python limit with nothing else in it, and it is a fixed cost this
   app cannot trim without giving up RDKit descriptors, 3D conformers, or the real
   VQE loop. Whatever the exact ceiling, we are near it for no good reason, since
   the same code fits comfortably in a container.

2. **Shape of the application.** This is a long-running ASGI server that also
   serves static files, keeps a SQLite job store and a disk cache. Serverless
   functions have an ephemeral, effectively read-only filesystem, so the SQLite
   store and the disk cache would silently stop persisting, and every request
   would pay a cold start that re-imports Qiskit. The API also expects a single
   process holding one SQLite file — exactly what a container gives and a
   function platform does not.

Splitting the frontend onto Vercel while the API runs in a container is possible
(the CORS origins are configuration; see `CORS_ORIGINS` below), but it is not
needed: `app/main.py` already serves the frontend, so one container is simpler and
has one origin to reason about. **This build is documented and tested as the
single-container shape.**

## What is inside the image, and what is not

- **PySCF and `qiskit-nature` are NOT installed in the image.** They publish no
  usable Windows wheel and are only needed to *generate* Hamiltonians, which is
  done offline (`scripts/precompute_hamiltonians.ipynb`, Colab or WSL). The image
  contains the results of that offline step as committed JSON under
  `app/data/hamiltonians/`, and the running API only reads those files. There is
  no live PySCF call anywhere in a request path.
- **`.dockerignore` keeps secrets and local state out of the image** (`.env`,
  `.venv`, caches, logs, SQLite files, `scripts/`, `tests/`, `docs/`). `app/` and
  `frontend/` are deliberately included: they are the application.
- The disk cache and the job store live under `app/data/`, which the image makes
  writable by the non-root user. On a host with an ephemeral filesystem they are
  **disposable by design** (SPEC.md section 3): losing them costs a cache rebuild
  and drops queued jobs, never correctness.

## Configuration

All configuration comes from environment variables (`app/config.py`); none of it
is a secret for this project.

| variable | default | notes |
|---|---|---|
| `PORT` | `8000` | Read by the container's `CMD`. Render sets it; Hugging Face Spaces expects `7860`. |
| `HOST` | `127.0.0.1` | The app's own default. The container overrides it with `--host 0.0.0.0`, which is required outside a container. |
| `CORS_ORIGINS` | empty | Comma-separated origins. Empty is correct for the single-container shape (same origin). Set it only when the frontend is deployed separately. |
| `CACHE_DIR` | `app/data/cache` | Disposable on an ephemeral filesystem. |
| `JOB_STORE_PATH` | `app/data/jobs.sqlite` | Disposable, as above. |
| `H2_HAMILTONIAN_DIR` | `app/data/hamiltonians` | Where the offline artifacts are read from. |

`.env.example` lists the template; copy it to `.env` locally and never commit it.
No API keys are required by this build.

## Health check and smoke test

The image defines:

```bash
python -c "import os, sys, urllib.request; port = os.environ.get('PORT', '8000'); \
  sys.exit(0 if urllib.request.urlopen('http://127.0.0.1:' + port + '/health', timeout=4).status == 200 else 1)"
```

After deploying, check these, in order:

1. `GET /health` → `200 {"status": "ok", ...}`.
2. `GET /` → `200 text/html` and `GET /app/styles.css` → `200 text/css` with
   `Cache-Control: no-cache`. The API and the frontend ship together, so a
   healthy API with a broken frontend means the `frontend/` copy is at fault.
3. `GET /api/quantum/molecules` → both registry entries, with `available: true`
   for H2 (its artifact is committed) and `available: false` for LiH unless the
   user generated `lih_sto3g_active.json` offline and committed it.
4. `GET /api/quantum/h2-curve` → `200` with 58 stored points and
   `status: "computed"`. A `503 hamiltonian_unavailable` here means the image was
   built without `app/data/hamiltonians/`.
5. `POST /api/quantum/simulate {"bond_length": 0.74, "seed": 7}` → `200` with a
   VQE energy and `within_chemical_accuracy: true`. This is the heaviest request
   in the app (a real COBYLA optimisation through a Qiskit estimator); if a free
   tier is too small for it, that shows up here first.

## Free-tier notes

- Free container tiers sleep after inactivity, so the first request after a sleep
  pays a cold start (importing Qiskit and Aer). That is a latency cost, not a
  correctness one.
- Pick a plan with enough memory for the VQE request: the shot-based path builds
  an Aer estimator per restart and the statevector path builds a 4x4 observable,
  which is small, but the interpreter plus Qiskit is not.
- Nothing here needs a GPU, a quantum backend, or a paid service, and no quantum
  advantage is claimed anywhere in this project (see `LIMITATIONS.md`).

## Not solved here

- **CI has never run on a runner.** The workflow file exists; "tests pass in CI"
  is unverified and is listed in `LIMITATIONS.md`.
- **No CI/CD to the host is configured.** Deploys are manual: build the image,
  push or point the host at the repository.
- **The Hermitian matrix, the artifacts and the results are unchanged by
  deploying.** Deploying this app does not make any stored number more accurate.

## Future work (not implemented)

None of the following exists in this build. They are not part of the results, no
endpoint serves them, and no number anywhere in the app is produced by them.

- **Fragment-based scaling.** Splitting a larger molecule into fragments, solving
  each small fragment (where the existing 2-qubit-style ansatz may still apply)
  and recombining the energies. This is the standard route past "H2 only", but it
  needs a fragmentation scheme, a consistent energy-combination rule, and a
  validation study against full-molecule references before it may report anything.
- **ML property prediction.** A learned model (for example solubility) would be a
  **prediction**, not a measurement: it would have to be labelled as such in the
  API and the UI, ship with the training data split and an RMSE against a
  mean-prediction baseline, and answer honestly when no model is loaded. Nothing
  of this is trained or loaded today, and SPEC.md section 8 keeps it last.
- **De novo generation.** Generating new molecules is not attempted. This project
  does not discover drugs; it is a study tool.

Adding any of these changes what the app can claim, so each one needs its own
tests, its own `LIMITATIONS.md` entry and its own `RESULTS.md` rows — the rules in
SPEC.md section 11 and `RESULTS.md` apply as much to a new feature as to a new
number.
