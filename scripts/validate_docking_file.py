#!/usr/bin/env python
"""Standalone validator for the offline redocking-validation result file.

Usage:
    python scripts/validate_docking_file.py [path]

Defaults to ``app/data/docking/redock_validation.json`` -- the file the Colab run
writes, which is **not** in this repository (Part A of the redocking task does not
run Vina, so no result exists yet). Exits 0 only when every check passes, and
non-zero with a precise list of problems otherwise. It never writes, repairs or
fabricates the file: a missing file is reported as missing.

Standard library only, and it deliberately does not import ``app`` or the
redocking script: this is an independent check of a file that was produced
elsewhere (a Colab session), so it re-derives every relation it can instead of
trusting the writer.

What it checks -- the shapes and the honesty rules:

* ``schema_version``, ``generated_at``, ``tool_versions`` (python, vina, meeko,
  rdkit, numpy), ``pass_threshold_angstrom`` (2.0), ``pass_rule``, ``disclaimers``
  and ``targets`` are present, and the schema version is one this validator knows.
* ``1VRT``, ``2HU4`` and ``1HSG`` are all present and no other target is: the
  disclaimers claim three complexes, so a file with more of them contradicts
  itself.
* Every ``ok`` target has 1..9 modes, ranks 1..n, scores in ascending order,
  finite scores and RMSDs, RMSD >= 0, ``top1_rmsd`` equal to mode rank 1's RMSD,
  ``best_mode_rmsd`` equal to the smallest mode RMSD, and
  ``passed == (top1_rmsd < pass_threshold_angstrom)``.
* Every ``failed`` target has ``passed: false``, a non-empty ``failure_reason``
  and **no** claimed RMSD (a failed target must not carry a pass metric).
* One fixed setting for every target: every ``ok`` target carries the same
  ``vina_params`` (with ``num_modes`` 9) and the same 22 A box. Per-target
  parameter tuning is the thing this rule exists to catch, so it is a failure,
  not a warning.
* The four required disclaimer statements are present: scores are predictions and
  not measured binding; redocking is an easy test; the results cover three
  complexes only; and no parameter was retuned to obtain a pass.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from datetime import datetime
from pathlib import Path
from typing import Any, Sequence

SCHEMA_VERSION = "1.0"
EXPECTED_PDB_IDS: tuple[str, ...] = ("1VRT", "2HU4", "1HSG")
EXPECTED_PASS_THRESHOLD = 2.0
EXPECTED_PASS_RULE = "top-ranked pose RMSD < 2.0 A"
EXPECTED_BOX_SIZE: tuple[float, float, float] = (22.0, 22.0, 22.0)
MAX_MODES = 9
NUM_MODES = 9
# The values are compared as stored; 1e-9 is far below the 3-decimal precision the
# writer records and far above float noise.
TOLERANCE = 1e-9

DEFAULT_RESULT = (
    Path(__file__).resolve().parents[1] / "app" / "data" / "docking" / "redock_validation.json"
)

TOP_LEVEL_FIELDS = (
    "schema_version",
    "generated_at",
    "tool_versions",
    "pass_threshold_angstrom",
    "pass_rule",
    "disclaimers",
    "targets",
)

TOOL_VERSION_FIELDS = ("python", "vina", "meeko", "rdkit", "numpy")

TARGET_FIELDS = (
    "pdb_id",
    "ligand_code",
    "receptor_chains_used",
    "box_center",
    "box_size",
    "vina_params",
    "modes",
    "top1_rmsd",
    "best_mode_rmsd",
    "passed",
    "status",
    "failure_reason",
)

MODE_FIELDS = ("rank", "score_kcal_per_mol", "rmsd_to_crystal_angstrom")

VINA_PARAM_FIELDS = ("exhaustiveness", "seed", "num_modes")

# Each entry is (what the statement has to say, lowercase fragments of which at
# least one must appear). A statement that says none of them is a problem.
REQUIRED_DISCLAIMER_STATEMENTS: tuple[tuple[str, tuple[str, ...]], ...] = (
    (
        "scores are predictions, not measured binding",
        ("not measured binding", "not measured affinity", "not measured", "not an experimental"),
    ),
    (
        "redocking is an easy test",
        ("easy test", "easy redocking", "easy test:", "an easy"),
    ),
    (
        "the results cover three complexes only",
        ("complexes only", "3 complexes only", "three complexes only"),
    ),
    (
        "no parameter was retuned to obtain a pass",
        ("retuned", "re-tuned", "tuned per target", "no tuning"),
    ),
)


class Problems:
    """Collects every problem instead of stopping at the first one."""

    def __init__(self) -> None:
        self.items: list[str] = []

    def add(self, message: str) -> None:
        self.items.append(message)

    def require(self, condition: bool, message: str) -> bool:
        if not condition:
            self.add(message)
        return bool(condition)


def is_finite_number(value: Any) -> bool:
    """True for a real, finite number -- and False for bools and strings."""
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return False
    return math.isfinite(float(value))


def is_sequence_of_size(value: Any, size: int) -> bool:
    return isinstance(value, (list, tuple)) and len(value) == size


def check_top_level(payload: Any, problems: Problems) -> dict[str, Any] | None:
    if not isinstance(payload, dict):
        problems.add(f"the file's top level must be a JSON object, found {type(payload).__name__}")
        return None

    for field in TOP_LEVEL_FIELDS:
        problems.require(field in payload, f"missing top-level field {field!r}")

    version = payload.get("schema_version")
    problems.require(
        version == SCHEMA_VERSION,
        f"schema_version must be {SCHEMA_VERSION!r}, found {version!r}",
    )

    generated_at = payload.get("generated_at")
    if problems.require(
        isinstance(generated_at, str) and generated_at.strip() != "",
        f"generated_at must be a non-empty string, found {generated_at!r}",
    ):
        try:
            datetime.fromisoformat(generated_at.replace("Z", "+00:00"))
        except ValueError:
            problems.add(
                f"generated_at {generated_at!r} is not an ISO-8601 timestamp; it must "
                "record when the run happened"
            )

    versions = payload.get("tool_versions")
    if problems.require(
        isinstance(versions, dict) and bool(versions),
        f"tool_versions must be a non-empty object, found {versions!r}",
    ):
        for field in TOOL_VERSION_FIELDS:
            value = versions.get(field)
            problems.require(
                isinstance(value, str) and value.strip() != "",
                f"tool_versions is missing a usable {field!r} (found {value!r})",
            )

    threshold = payload.get("pass_threshold_angstrom")
    if problems.require(
        is_finite_number(threshold),
        f"pass_threshold_angstrom must be a finite number, found {threshold!r}",
    ):
        problems.require(
            abs(float(threshold) - EXPECTED_PASS_THRESHOLD) <= TOLERANCE,
            f"pass_threshold_angstrom must be {EXPECTED_PASS_THRESHOLD}, found {threshold!r}",
        )

    rule = payload.get("pass_rule")
    problems.require(
        rule == EXPECTED_PASS_RULE,
        f"pass_rule must be {EXPECTED_PASS_RULE!r}, found {rule!r}",
    )

    check_disclaimers(payload.get("disclaimers"), problems)

    targets = payload.get("targets")
    if not problems.require(
        isinstance(targets, list) and bool(targets),
        f"targets must be a non-empty list, found {type(targets).__name__}",
    ):
        return None
    return payload


def check_disclaimers(disclaimers: Any, problems: Problems) -> None:
    if not problems.require(
        isinstance(disclaimers, list) and bool(disclaimers),
        f"disclaimers must be a non-empty list, found {type(disclaimers).__name__}",
    ):
        return
    for index, disclaimer in enumerate(disclaimers):
        problems.require(
            isinstance(disclaimer, str) and disclaimer.strip() != "",
            f"disclaimers[{index}] must be a non-empty string, found {disclaimer!r}",
        )
    text = " ".join(str(item).lower() for item in disclaimers)
    for statement, fragments in REQUIRED_DISCLAIMER_STATEMENTS:
        problems.require(
            any(fragment in text for fragment in fragments),
            f"the disclaimers do not state that {statement} (looked for any of "
            f"{list(fragments)}); each of these statements is required",
        )


def check_expected_targets(targets: Sequence[Any], problems: Problems) -> None:
    pdb_ids = [target.get("pdb_id") for target in targets if isinstance(target, dict)]
    missing = [pdb_id for pdb_id in EXPECTED_PDB_IDS if pdb_id not in pdb_ids]
    problems.require(
        not missing,
        f"the expected PDB ID(s) {missing} are missing; the run must cover "
        f"{', '.join(EXPECTED_PDB_IDS)}",
    )
    unexpected = [pdb_id for pdb_id in pdb_ids if pdb_id not in EXPECTED_PDB_IDS]
    problems.require(
        not unexpected,
        f"unexpected PDB ID(s) {unexpected} are present; the disclaimers say the "
        f"results cover {len(EXPECTED_PDB_IDS)} complexes only "
        f"({', '.join(EXPECTED_PDB_IDS)}), so the file contradicts itself",
    )


def check_mode(target_id: str, index: int, mode: Any, problems: Problems, *, ok: bool) -> None:
    where = f"{target_id} modes[{index}]"
    if not isinstance(mode, dict):
        problems.add(f"{where} must be an object, found {type(mode).__name__}")
        return
    for field in MODE_FIELDS:
        problems.require(field in mode, f"{where} is missing {field!r}")

    rank = mode.get("rank")
    problems.require(
        isinstance(rank, int) and not isinstance(rank, bool) and rank == index + 1,
        f"{where} must have rank {index + 1}, found {rank!r}",
    )

    score = mode.get("score_kcal_per_mol")
    problems.require(
        is_finite_number(score),
        f"{where}.score_kcal_per_mol must be a finite number, found {score!r}",
    )

    rmsd = mode.get("rmsd_to_crystal_angstrom")
    if not ok:
        problems.require(
            rmsd is None,
            f"{where}.rmsd_to_crystal_angstrom must be null in a failed target (found "
            f"{rmsd!r}); a failed target must not report an RMSD",
        )
        return
    if rmsd is None:
        problems.add(
            f"{where}.rmsd_to_crystal_angstrom is null in an ok target; every docked "
            "mode of a completed target must carry its measured RMSD"
        )
        return
    if problems.require(
        is_finite_number(rmsd),
        f"{where}.rmsd_to_crystal_angstrom must be a finite number, found {rmsd!r}",
    ):
        problems.require(
            float(rmsd) >= 0.0,
            f"{where}.rmsd_to_crystal_angstrom must be >= 0, found {rmsd!r}",
        )


def check_ok_target(
    target: dict[str, Any], problems: Problems, *, threshold: float = EXPECTED_PASS_THRESHOLD
) -> None:
    target_id = str(target.get("pdb_id"))

    ligand_code = target.get("ligand_code")
    problems.require(
        isinstance(ligand_code, str) and ligand_code.strip() != "",
        f"{target_id}: ligand_code must be a non-empty string, found {ligand_code!r}",
    )

    chains = target.get("receptor_chains_used")
    if problems.require(
        isinstance(chains, list) and bool(chains),
        f"{target_id}: receptor_chains_used must be a non-empty list, found {chains!r}",
    ):
        for chain in chains:
            problems.require(
                isinstance(chain, str) and chain.strip() != "",
                f"{target_id}: receptor_chains_used entries must be non-empty strings, "
                f"found {chain!r}",
            )

    center = target.get("box_center")
    if problems.require(
        is_sequence_of_size(center, 3),
        f"{target_id}: box_center must be three numbers, found {center!r}",
    ):
        for axis, value in zip("xyz", center):
            problems.require(
                is_finite_number(value),
                f"{target_id}: box_center.{axis} must be a finite number, found {value!r}",
            )

    check_box_size(target_id, target.get("box_size"), problems)
    check_vina_params(target_id, target.get("vina_params"), problems)

    modes = target.get("modes")
    if problems.require(
        isinstance(modes, list),
        f"{target_id}: modes must be a list, found {type(modes).__name__}",
    ):
        problems.require(
            1 <= len(modes) <= MAX_MODES,
            f"{target_id}: modes must hold 1..{MAX_MODES} docked poses, found {len(modes)}",
        )
        for index, mode in enumerate(modes):
            check_mode(target_id, index, mode, problems, ok=True)
        scores = [
            mode["score_kcal_per_mol"]
            for mode in modes
            if isinstance(mode, dict) and is_finite_number(mode.get("score_kcal_per_mol"))
        ]
        problems.require(
            all(
                earlier <= later + TOLERANCE
                for earlier, later in zip(scores, scores[1:])
            ),
            f"{target_id}: modes must be sorted by ascending score, found "
            f"{[mode.get('score_kcal_per_mol') for mode in modes if isinstance(mode, dict)]}",
        )

        rmsds = [
            float(mode["rmsd_to_crystal_angstrom"])
            for mode in modes
            if isinstance(mode, dict) and is_finite_number(mode.get("rmsd_to_crystal_angstrom"))
        ]
    else:
        rmsds = []

    top1 = target.get("top1_rmsd")
    if problems.require(
        is_finite_number(top1),
        f"{target_id}: top1_rmsd must be a finite number for an ok target, found {top1!r}",
    ):
        problems.require(
            float(top1) >= 0.0,
            f"{target_id}: top1_rmsd must be >= 0, found {top1!r}",
        )
        if isinstance(modes, list) and modes and isinstance(modes[0], dict):
            first = modes[0].get("rmsd_to_crystal_angstrom")
            if is_finite_number(first):
                problems.require(
                    abs(float(top1) - float(first)) <= TOLERANCE,
                    f"{target_id}: top1_rmsd {top1!r} is not the rank-1 mode's RMSD "
                    f"{first!r}",
                )
        expected_passed = float(top1) < threshold
        problems.require(
            target.get("passed") is expected_passed,
            f"{target_id}: passed must be {expected_passed} for "
            f"top1_rmsd {top1!r} (pass rule {EXPECTED_PASS_RULE}), found "
            f"{target.get('passed')!r}",
        )

    best = target.get("best_mode_rmsd")
    if problems.require(
        is_finite_number(best),
        f"{target_id}: best_mode_rmsd must be a finite number for an ok target, "
        f"found {best!r}",
    ):
        problems.require(
            float(best) >= 0.0,
            f"{target_id}: best_mode_rmsd must be >= 0, found {best!r}",
        )
        if rmsds:
            problems.require(
                abs(float(best) - min(rmsds)) <= TOLERANCE,
                f"{target_id}: best_mode_rmsd {best!r} is not the smallest mode RMSD "
                f"{min(rmsds)!r}",
            )

    problems.require(
        target.get("failure_reason") is None,
        f"{target_id}: an ok target must have failure_reason null, found "
        f"{target.get('failure_reason')!r}",
    )


def check_failed_target(target: dict[str, Any], problems: Problems) -> None:
    target_id = str(target.get("pdb_id"))

    problems.require(
        target.get("passed") is False,
        f"{target_id}: a failed target must have passed false, found "
        f"{target.get('passed')!r}",
    )
    reason = target.get("failure_reason")
    problems.require(
        isinstance(reason, str) and reason.strip() != "",
        f"{target_id}: a failed target must record why it failed in failure_reason, "
        f"found {reason!r}",
    )
    for field in ("top1_rmsd", "best_mode_rmsd"):
        problems.require(
            target.get(field) is None,
            f"{target_id}: a failed target must not report {field} (found "
            f"{target.get(field)!r}); a failure is a failure, not a partial result",
        )

    modes = target.get("modes")
    if modes is None:
        return
    if problems.require(
        isinstance(modes, list),
        f"{target_id}: modes must be a list, found {type(modes).__name__}",
    ):
        problems.require(
            len(modes) <= MAX_MODES,
            f"{target_id}: a failed target must not report more than {MAX_MODES} modes, "
            f"found {len(modes)}",
        )
        for index, mode in enumerate(modes):
            check_mode(target_id, index, mode, problems, ok=False)
        scores = [
            mode["score_kcal_per_mol"]
            for mode in modes
            if isinstance(mode, dict) and is_finite_number(mode.get("score_kcal_per_mol"))
        ]
        problems.require(
            all(
                earlier <= later + TOLERANCE
                for earlier, later in zip(scores, scores[1:])
            ),
            f"{target_id}: modes must be sorted by ascending score",
        )


def check_box_size(target_id: str, box_size: Any, problems: Problems) -> None:
    if not problems.require(
        is_sequence_of_size(box_size, 3),
        f"{target_id}: box_size must be three numbers, found {box_size!r}",
    ):
        return
    for axis, value, expected in zip("xyz", box_size, EXPECTED_BOX_SIZE):
        if problems.require(
            is_finite_number(value),
            f"{target_id}: box_size.{axis} must be a finite number, found {value!r}",
        ):
            problems.require(
                abs(float(value) - expected) <= TOLERANCE,
                f"{target_id}: box_size.{axis} must be {expected}, found {value!r}; every "
                "target uses the same documented (easy, box-on-the-true-site) setup",
            )


def check_vina_params(target_id: str, params: Any, problems: Problems) -> None:
    if not problems.require(
        isinstance(params, dict) and bool(params),
        f"{target_id}: vina_params must be a non-empty object, found {params!r}",
    ):
        return
    for field in VINA_PARAM_FIELDS:
        problems.require(field in params, f"{target_id}: vina_params is missing {field!r}")
    exhaustiveness = params.get("exhaustiveness")
    problems.require(
        isinstance(exhaustiveness, int)
        and not isinstance(exhaustiveness, bool)
        and exhaustiveness >= 1,
        f"{target_id}: vina_params.exhaustiveness must be a positive integer, found "
        f"{exhaustiveness!r}",
    )
    seed = params.get("seed")
    problems.require(
        isinstance(seed, int) and not isinstance(seed, bool),
        f"{target_id}: vina_params.seed must be an integer, found {seed!r}",
    )
    num_modes = params.get("num_modes")
    problems.require(
        num_modes == NUM_MODES,
        f"{target_id}: vina_params.num_modes must be {NUM_MODES}, found {num_modes!r}",
    )


def check_single_setting(targets: Sequence[dict[str, Any]], problems: Problems) -> None:
    """Every ok target must carry the same parameters: one setting for all of them."""
    ok_targets = [t for t in targets if t.get("status") == "ok"]
    if len(ok_targets) < 2:
        return
    first = ok_targets[0]
    reference = {
        "vina_params": first.get("vina_params"),
        "box_size": first.get("box_size"),
    }
    for target in ok_targets[1:]:
        for field, value in reference.items():
            problems.require(
                target.get(field) == value,
                f"{target['pdb_id']}: {field} {target.get(field)!r} differs from "
                f"{first['pdb_id']}'s {value!r}; one fixed setting is used for every "
                "target and nothing is retuned per target",
            )


def validate(payload: Any) -> list[str]:
    """Every problem found in ``payload`` (empty list means it is valid)."""
    problems = Problems()
    top = check_top_level(payload, problems)
    if top is None:
        return problems.items

    targets = top["targets"]
    check_expected_targets(targets, problems)

    # The task's rule is passed == (top1_rmsd < pass_threshold) with the file's
    # own threshold; that threshold is separately required to be 2.0 above.
    stored_threshold = top.get("pass_threshold_angstrom")
    threshold = (
        float(stored_threshold) if is_finite_number(stored_threshold) else EXPECTED_PASS_THRESHOLD
    )

    for index, target in enumerate(targets):
        if not isinstance(target, dict):
            problems.add(f"targets[{index}] must be an object, found {type(target).__name__}")
            continue
        target_id = f"targets[{index}]"
        for field in TARGET_FIELDS:
            if field not in target:
                problems.add(f"{target_id} (pdb_id {target.get('pdb_id')!r}) is missing {field!r}")
        status = target.get("status")
        if status not in ("ok", "failed"):
            problems.add(
                f"{target_id} (pdb_id {target.get('pdb_id')!r}): status must be 'ok' or "
                f"'failed', found {status!r}"
            )
            continue
        if status == "ok":
            check_ok_target(target, problems, threshold=threshold)
        else:
            check_failed_target(target, problems)

    check_single_setting([t for t in targets if isinstance(t, dict)], problems)
    return problems.items


def load(path: Path) -> Any:
    """Read the JSON file, or exit with the reason it cannot be read."""
    if not path.exists():
        print(f"FAIL: {path} does not exist.", file=sys.stderr)
        print(
            "The redocking result file is produced by the Colab run of "
            "scripts/redock_validation.py; it is not created by this validator.",
            file=sys.stderr,
        )
        raise SystemExit(1)
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        print(f"FAIL: {path} is not valid JSON ({exc}).", file=sys.stderr)
        raise SystemExit(1) from exc


def summarise(payload: Any) -> str:
    lines = [
        f"{'pdb_id':<8} {'ligand':<8} {'status':<7} {'passed':<7} {'top1_rmsd':>10} "
        f"{'best_rmsd':>10} {'modes':>5}",
    ]
    for target in payload["targets"]:
        top1 = target.get("top1_rmsd")
        best = target.get("best_mode_rmsd")
        lines.append(
            f"{str(target.get('pdb_id')):<8} {str(target.get('ligand_code') or '-'):<8} "
            f"{str(target.get('status')):<7} {str(bool(target.get('passed'))):<7} "
            f"{'-' if top1 is None else f'{float(top1):.3f}':>10} "
            f"{'-' if best is None else f'{float(best):.3f}':>10} "
            f"{len(target.get('modes') or []):>5}"
        )
    return "\n".join(lines)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Validate the offline redocking-validation result file: schema and "
            "honesty rules, one fixed Vina setting for every target, and failures "
            "recorded as failures."
        )
    )
    parser.add_argument(
        "path",
        nargs="?",
        default=str(DEFAULT_RESULT),
        help=f"result file to validate (default: {DEFAULT_RESULT})",
    )
    args = parser.parse_args(argv)

    path = Path(args.path)
    payload = load(path)
    problems = validate(payload)

    if problems:
        print(f"FAIL: {len(problems)} problem(s) in {path}", file=sys.stderr)
        for index, problem in enumerate(problems, start=1):
            print(f"  {index:2d}. {problem}", file=sys.stderr)
        return 1

    print(f"OK: every check passed for {path}")
    print(f"pass rule: {payload['pass_rule']} (threshold {payload['pass_threshold_angstrom']} A)")
    print(f"generated: {payload['generated_at']}  tools: {payload['tool_versions']}")
    print(summarise(payload))
    print(
        "Reminder: these are docking predictions for an easy redocking test (box on "
        "the true crystal site), not measured binding and not a screening result."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
