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
