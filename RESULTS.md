# RESULTS.md

Every measured number that appears in a slide, a report or the README is listed here with the command that produced it.

## Rule

- **No number appears in slides or reports unless it is listed in the table below.**

## Numbers

| Number | Value | Command that produced it | Date | File |
| --- | --- | --- | --- | --- |
| H2 VQE vs exact, noiseless run at 0.74 Å (seed 7) | VQE -1.137283834 Ha, exact -1.137283834 Ha, \|VQE - exact\| 8.0e-12 Ha, 3 restarts x 450 evaluations | `python scripts/measure_h2_vqe.py --bond-lengths 0.74 --seed 7` | 2026-10-13 | scripts/measure_h2_vqe.py |
| H2 VQE vs exact, 8192-shot run at 0.74 Å (seed 7) | VQE -1.137005090 Ha ± 2.461 mHa, exact -1.137283834 Ha, \|VQE - exact\| 0.279 mHa, 3 x 208 evaluations; uncertainty from `estimator_variance_metadata` (sqrt(variance/shots)) | `python scripts/measure_h2_vqe.py --bond-lengths 0.74 --seed 7 --shots 8192` | 2026-10-13 | scripts/measure_h2_vqe.py |
| H2 VQE vs exact at 0.5 / 0.74 / 1.5 / 2.5 Å, noiseless, seed 7 | \|VQE - exact\| = 2.7e-10 / 8.0e-12 / 6.0e-13 / 5.4e-4 Ha — all below the 1.6 mHa target | `python scripts/measure_h2_vqe.py` | 2026-10-13 | scripts/measure_h2_vqe.py |
| Worst H2 VQE error over all 58 stored points, noiseless, seed 7 | 1.016e-3 Ha at 2.65 Å (target 1.6e-3 Ha); every point inside chemical accuracy, no variational-principle refusal | `python scripts/measure_h2_vqe.py --all` | 2026-10-13 | scripts/measure_h2_vqe.py |
| STO-3G resource estimates for five molecules (spatial orbitals / electrons / full qubits / reduced qubits / statevector memory / verdict) | H2 2/2/4/2, 64 bytes, `simulable_in_app`; LiH 6/4/12/10, 16.0 KB, `precompute_only`; H2O 7/10/14/12, 64.0 KB, `precompute_only`; CH4 9/10/18/16, 1.0 MB, `precompute_only`; C9H8O4 (aspirin) 73/94/146/144, more than 10^15 bytes, `beyond_classical_simulation` | `python scripts/estimate_resources_examples.py` | 2026-10-14 | scripts/estimate_resources_examples.py |
| STO-3G orbital counts cross-checked with PySCF `mol.nao` (the five formulas above) | H2 2, LiH 6, H2O 7, CH4 9, C9H8O4 73 spatial orbitals — all match the formula-derived counts; electrons 2/4/10/10/94. Run in a Colab session, **not reproducible in this repository** (no PySCF wheel here) and no PySCF version was recorded | PySCF `mol.nao` with the STO-3G basis, reported for this project | not recorded | not in this repository; recorded in `PYSCF_VERIFIED_COUNTS`, `app/core/resources.py` |

## First real VQE measurement (2026-10-13)

The simulate endpoint's own code path (`app/services/quantum_service.py`), noiseless statevector
estimator, fixed seed 7, 3 COBYLA restarts × 150 evaluations (450 evaluations in total); exact =
NumPy diagonalisation of the artifact's Pauli terms plus nuclear repulsion:

| bond length | HF | VQE | exact | abs error (Ha) | restarts x evaluations | seed |
| --- | --- | --- | --- | --- | --- | --- |
| 0.50 A | -1.042996275 | -1.055159794 | -1.055159794 | 2.746e-10 | 3 x 450 | 7 |
| 0.74 A | -1.116759307 | -1.137283834 | -1.137283834 | 8.027e-12 | 3 x 450 | 7 |
| 1.50 A | -0.910873555 | -0.998149353 | -0.998149353 | 6.013e-13 | 3 x 450 | 7 |
| 2.50 A | -0.702943600 | -0.935511640 | -0.936054920 | 5.433e-04 | 3 x 450 | 7 |

Reproduce with `python scripts/measure_h2_vqe.py` (or `--all` for the full 58-point grid). The
0.74 Å run lands on H2's equilibrium reference of about -1.137 Ha within chemical accuracy, as
SPEC.md section 6 requires.

## First shot-based VQE measurement (2026-10-13)

The same code path with the shot-based Aer estimator. The headline energy is a fresh
re-evaluation at the final optimal parameters (not the minimum of the noisy trace) and its
uncertainty is the estimator's own variance metadata, reported as `sqrt(variance / shots)` —
sampling error, not a hardware result, and `within_chemical_accuracy` stays `null` by design:

| bond length | HF | VQE | std (mHa) | exact | abs error (mHa) | restarts x evaluations | seed |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0.74 A | -1.116759307 | -1.137005090 | 2.461 | -1.137283834 | 0.279 | 3 x 208 | 7 |

Reproduce with `python scripts/measure_h2_vqe.py --bond-lengths 0.74 --seed 7 --shots 8192`.

## First STO-3G resource estimates (2026-10-14)

These five rows are **counted from the molecular formula**, not simulated and not fetched from a
database: spatial orbitals from the STO-3G minimal-basis rule (H/He 1 function, Li–Ne 5,
Na–Ar 9), electrons as the sum of the neutral atomic numbers, one qubit per spin orbital, minus
the fixed two-qubit saving of the parity mapping with two-qubit reduction, and a complex128
statevector of `2 ** n_qubits * 16` bytes. Every figure is tagged `exact` or `estimated` in the
endpoint's `tags` dict. Reproduce with `python scripts/estimate_resources_examples.py` (it calls
`GET /api/quantum/resources?formula=` in-process, so no server is needed).

| formula | spatial orbitals | electrons | full qubits | reduced qubits | statevector memory | verdict |
| --- | --- | --- | --- | --- | --- | --- |
| `H2` | 2 | 2 | 4 | 2 | 64 bytes | `simulable_in_app` |
| `LiH` | 6 | 4 | 12 | 10 | 16.0 KB | `precompute_only` |
| `H2O` | 7 | 10 | 14 | 12 | 64.0 KB | `precompute_only` |
| `CH4` | 9 | 10 | 18 | 16 | 1.0 MB | `precompute_only` |
| `C9H8O4` (aspirin) | 73 | 94 | 146 | 144 | more than 10^15 bytes | `beyond_classical_simulation` |

The verdicts come from `interactive_cap_qubits: 8` (this app's own limit, a project choice) and
`laptop_statevector_ceiling_qubits: 25` (a rough estimate); they are project limits and rough
estimates, not measured hardware or algorithm bounds. `simulable_in_app` means the size fits the
app's budget, **not** that a simulation was run — none of these five molecules is simulated by
this build, and H2 is the only molecule with a VQE ansatz. No Pauli-term count or circuit size is
reported for any of them, because a formula cannot give one.

**Orbital counts for these five molecules were cross-checked with PySCF `mol.nao` in Colab; other
formulas are counted by the same rule but not machine-confirmed.** The five payloads therefore
carry `verified_against_pyscf: true`, every other formula `false` (see LIMITATIONS.md). No PySCF
version is recorded, because none was captured, and the check is not re-run in this repository.

## Redocking validation: results pending the Colab run (no numbers yet)

`scripts/redock_validation.py` redocks 1VRT, 2HU4 and 1HSG offline and would write its per-pose
RMSDs and Vina scores to `app/data/docking/redock_validation.json`, which is produced by the
Colab run of `scripts/redock_validation.ipynb` (`python scripts/redock_validation.py OUT.json`,
by default `--exhaustiveness 16 --seed 7`). **That run has not happened and the file is not in
this repository, so there is no redocking number to list above** — no RMSD, no Vina score,
nothing in the API and nothing in the UI. The pipeline itself has been written and tested, but a
pipeline is not a result.

When it runs, every number that reaches a slide or a report gets a row in the table above with
the command that produced it. Two things must be said alongside any such row: a redocking pass
means the search rediscovered a pose **that was already inside a 22 A box centred on the true
crystal site**, which is an easy test and says nothing about virtual screening or affinity
prediction; and a Vina score is a **prediction, not measured binding**, never a Kd, Ki or IC50.

**Recording rule (applies to the first run and every rerun).** "The first complete run of the
redocking pipeline at the fixed settings is the record, whatever it shows, including failures.
Any rerun needs a dated log entry stating why; both results stay on file. Parameters are not
tuned per target to obtain a pass." `tests/test_result_file_guardrails.py` validates the committed
file as soon as it exists, so a result cannot arrive here without the schema and honesty checks
having run.
