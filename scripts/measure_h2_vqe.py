"""Run the simulate endpoint's VQE at named bond lengths and print a table.

Usage::

    python scripts/measure_h2_vqe.py                  # 0.5, 0.74, 1.5, 2.5 A
    python scripts/measure_h2_vqe.py --all            # every stored point
    python scripts/measure_h2_vqe.py --bond-lengths 0.74 1.1 --seed 11
    python scripts/measure_h2_vqe.py --bond-lengths 0.74 --shots 8192

Every printed number comes from ``compute_h2_simulate`` -- the same code path
the API serves -- which reads the offline artifact and runs the VQE loop; this
script computes and invents nothing itself. Exits non-zero if any point misses
SPEC.md section 6's 1.6 mHa noiseless target, or if a run raises (including the
variational-principle guard).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from app.services.quantum_service import compute_h2_curve, compute_h2_simulate  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description="Measure VQE against exact for stored H2 points")
    parser.add_argument(
        "--bond-lengths",
        type=float,
        nargs="+",
        default=[0.5, 0.74, 1.5, 2.5],
        help="bond lengths to run (each snapped to the stored grid by the service)",
    )
    parser.add_argument("--seed", type=int, default=7, help="RNG seed for the VQE restarts")
    parser.add_argument("--all", action="store_true", help="run every stored bond length")
    parser.add_argument(
        "--shots",
        type=int,
        default=None,
        help="run mode: omit for the noiseless statevector run, or a shot count (e.g. 8192)",
    )
    args = parser.parse_args()

    lengths = compute_h2_curve()["bond_lengths_angstrom"] if args.all else args.bond_lengths

    if args.shots is not None:
        header = "| bond length (A) | HF | VQE | std (mHa) | exact | abs error (mHa) | restarts x evaluations | seed |"
    else:
        header = "| bond length (A) | HF | VQE | exact | abs error (Ha) | restarts x evaluations | seed |"
    print(header)
    print("| " + " | ".join(["---"] * (header.count("|") - 1)) + " |")

    failures: list[float] = []
    for bond in lengths:
        result = compute_h2_simulate(bond, shots=args.shots, seed=args.seed)
        error = abs(result["error_vs_exact_ha"])
        budget = f"{result['n_restarts']} x {result['total_energy_evaluations']}"
        if args.shots is not None:
            std_mha = (
                result["vqe_energy_std_ha"] * 1000.0
                if result["vqe_energy_std_ha"] is not None
                else float("nan")
            )
            print(
                f"| {result['used_bond_length']:.2f} A"
                f" | {result['hf_energy_ha']:.9f}"
                f" | {result['vqe_energy_ha']:.9f}"
                f" | {std_mha:.3f}"
                f" | {result['exact_energy_ha']:.9f}"
                f" | {error * 1000.0:.3f}"
                f" | {budget}"
                f" | {result['seed']} |"
            )
        else:
            print(
                f"| {result['used_bond_length']:.2f} A"
                f" | {result['hf_energy_ha']:.9f}"
                f" | {result['vqe_energy_ha']:.9f}"
                f" | {result['exact_energy_ha']:.9f}"
                f" | {error:.3e}"
                f" | {budget}"
                f" | {result['seed']} |"
            )
        if error > result["chemical_accuracy_threshold_ha"] and args.shots is None:
            failures.append(result["used_bond_length"])

    if failures:
        print(
            f"exceeded the 1.6 mHa target at: {', '.join(f'{b:g} A' for b in failures)}",
            file=sys.stderr,
        )
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
