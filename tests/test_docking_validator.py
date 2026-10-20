"""Tests for ``scripts/validate_docking_file.py``.

The real artifact (``app/data/docking/redock_validation.json``) is produced by an
offline Colab run of ``scripts/redock_validation.py`` and is **not** in this
repository: Part A of the redocking task writes the pipeline, not the result. So
these tests exercise the validator against hand-written payloads in a temporary
directory, and almost all of what matters is what the validator **rejects**.

The validator is run as a subprocess rather than imported, because its contract
is a process exit code plus a message on stderr -- importing it would skip exactly
the part the caller depends on (``tests/test_redock_validation.py`` imports it for
the payload-shape cross-check).

Every payload here is a fixture. None of these numbers is a docking result, and
nothing in this file writes to the repository.
"""

from __future__ import annotations

import copy
import json
import subprocess
import sys
from pathlib import Path
from typing import Any, Callable

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
VALIDATOR = REPO_ROOT / "scripts" / "validate_docking_file.py"

DISCLAIMERS = [
    "Scores are docking predictions, not measured binding: they are not a Kd, Ki or IC50.",
    "Redocking is an easy test: the box is on the true crystal site, so a low RMSD "
    "says the search can rediscover a pose that is already inside the box.",
    "These results cover 3 complexes only (1VRT, 2HU4, 1HSG).",
    "No parameter was retuned per target and no run was repeated until it passed.",
]


def _mode(rank: int, score: float, rmsd: float) -> dict[str, Any]:
    return {"rank": rank, "score_kcal_per_mol": score, "rmsd_to_crystal_angstrom": rmsd}


def _ok_target(pdb_id: str, ligand: str, modes: list[dict[str, Any]]) -> dict[str, Any]:
    rmsds = [mode["rmsd_to_crystal_angstrom"] for mode in modes]
    return {
        "pdb_id": pdb_id,
        "ligand_code": ligand,
        "ligand_identification": {"selection_rule": "metadata rule"},
        "ligand_definition": {"smiles": "C"},
        "ligand_instance": {"chain": "A"},
        "ligand_preparation": {"heavy_atom_count": 20},
        "receptor_chains_used": ["A"],
        "receptor_preparation": {"atom_count": 2962},
        "box_center": [0.383, 81.705, 109.195],
        "box_size": [22.0, 22.0, 22.0],
        "vina_params": {
            "exhaustiveness": 16,
            "seed": 7,
            "num_modes": 9,
            "energy_range_kcal_per_mol": 3.0,
            "seed_used": 7,
        },
        "modes": modes,
        "top1_rmsd": rmsds[0],
        "best_mode_rmsd": min(rmsds),
        "passed": rmsds[0] < 2.0,
        "status": "ok",
        "failure_reason": None,
        "failure_type": None,
        "inputs": {"receptor_pdbqt_sha256": "0" * 64},
        "notes": [],
        "crystal_pose_self_check_rmsd_angstrom": 0.0,
    }


def _failed_target(pdb_id: str, reason: str = "no non-polymer entity passes the ligand rule") -> dict[str, Any]:
    return {
        "pdb_id": pdb_id,
        "ligand_code": None,
        "ligand_identification": None,
        "ligand_definition": None,
        "ligand_instance": None,
        "ligand_preparation": None,
        "receptor_chains_used": [],
        "receptor_preparation": None,
        "box_center": None,
        "box_size": [22.0, 22.0, 22.0],
        "vina_params": {
            "exhaustiveness": 16,
            "seed": 7,
            "num_modes": 9,
            "energy_range_kcal_per_mol": 3.0,
        },
        "modes": [],
        "top1_rmsd": None,
        "best_mode_rmsd": None,
        "passed": False,
        "status": "failed",
        "failure_reason": reason,
        "failure_type": "TargetFailure",
        "inputs": None,
        "notes": [],
        "crystal_pose_self_check_rmsd_angstrom": None,
    }


def _payload() -> dict[str, Any]:
    """A good file: two docked targets (one passing, one not) and one failure."""
    return {
        "schema_version": "1.0",
        "generated_at": "2026-10-16T09:30:00Z",
        "generated_by": "scripts/redock_validation.py",
        "tool_versions": {
            "python": "3.11.12",
            "vina": "1.2.7",
            "meeko": "0.8.0",
            "rdkit": "2026.03.6",
            "numpy": "2.0.2",
        },
        "pass_threshold_angstrom": 2.0,
        "pass_rule": "top-ranked pose RMSD < 2.0 A",
        "box_size_angstrom": [22.0, 22.0, 22.0],
        "setup": {"kind": "redocking of the deposited crystal pose"},
        "vina_params": {"exhaustiveness": 16, "seed": 7, "num_modes": 9},
        "disclaimers": list(DISCLAIMERS),
        "targets": [
            _ok_target(
                "1VRT",
                "NVP",
                [
                    _mode(1, -9.514, 1.25),
                    _mode(2, -9.101, 1.31),
                    _mode(3, -8.812, 4.02),
                ],
            ),
            _ok_target(
                "2HU4",
                "G39",
                [
                    _mode(1, -7.204, 2.75),
                    _mode(2, -7.010, 2.60),
                ],
            ),
            _failed_target("1HSG"),
        ],
    }


def _write(tmp_path: Path, payload: Any, name: str = "redock_validation.json") -> Path:
    path = tmp_path / name
    if isinstance(payload, str):
        path.write_text(payload, encoding="utf-8")
    else:
        path.write_text(json.dumps(payload), encoding="utf-8")
    return path


def _run(path: Path) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(VALIDATOR), str(path)],
        cwd=str(REPO_ROOT),
        capture_output=True,
        text=True,
    )


def test_a_good_file_passes(tmp_path: Path) -> None:
    """Positive control: without this, the failure cases below prove nothing."""
    result = _run(_write(tmp_path, _payload()))

    assert result.returncode == 0, result.stderr
    assert "OK: every check passed" in result.stdout
    # The summary and the reminder must be printed, not merely asserted.
    assert "1VRT" in result.stdout and "NVP" in result.stdout
    assert "top-ranked pose RMSD < 2.0 A" in result.stdout
    assert "not measured binding" in result.stdout


def test_a_good_file_with_all_three_docked_as_ok_passes(tmp_path: Path) -> None:
    payload = _payload()
    payload["targets"][2] = _ok_target("1HSG", "MK1", [_mode(1, -11.5, 0.83)])
    result = _run(_write(tmp_path, payload))

    assert result.returncode == 0, result.stderr


def test_a_failed_target_is_allowed_and_reported(tmp_path: Path) -> None:
    payload = _payload()
    payload["targets"][0] = _failed_target("1VRT", "the ligand code could not be identified")
    result = _run(_write(tmp_path, payload))

    assert result.returncode == 0, result.stderr
    assert "failed" in result.stdout


def test_missing_file_fails_and_is_not_created(tmp_path: Path) -> None:
    missing = tmp_path / "not_written_yet.json"

    result = _run(missing)

    assert result.returncode != 0
    assert "does not exist" in result.stderr
    assert not missing.exists(), "the validator must not create or invent the result file"


def test_malformed_json_fails(tmp_path: Path) -> None:
    result = _run(_write(tmp_path, "{not json at all"))

    assert result.returncode != 0
    assert "not valid JSON" in result.stderr


def test_usage_error_for_extra_arguments(tmp_path: Path) -> None:
    result = subprocess.run(
        [sys.executable, str(VALIDATOR), "a.json", "b.json"],
        cwd=str(REPO_ROOT),
        capture_output=True,
        text=True,
    )

    assert result.returncode == 2
    assert "usage:" in result.stderr


# --------------------------------------------------------------------------- #
# one broken rule at a time
# --------------------------------------------------------------------------- #


def _wrong_schema(payload: dict[str, Any]) -> dict[str, Any]:
    payload["schema_version"] = "2.0"
    return payload


def _drop_top_level_field(payload: dict[str, Any]) -> dict[str, Any]:
    del payload["disclaimers"]
    return payload


def _drop_a_tool_version(payload: dict[str, Any]) -> dict[str, Any]:
    del payload["tool_versions"]["numpy"]
    return payload


def _wrong_threshold(payload: dict[str, Any]) -> dict[str, Any]:
    payload["pass_threshold_angstrom"] = 5.0
    return payload


def _wrong_pass_rule(payload: dict[str, Any]) -> dict[str, Any]:
    payload["pass_rule"] = "best pose RMSD < 2.0 A"
    return payload


def _not_a_timestamp(payload: dict[str, Any]) -> dict[str, Any]:
    payload["generated_at"] = "yesterday"
    return payload


def _empty_disclaimers(payload: dict[str, Any]) -> dict[str, Any]:
    payload["disclaimers"] = []
    return payload


def _drop_the_easy_test_statement(payload: dict[str, Any]) -> dict[str, Any]:
    payload["disclaimers"] = [d for d in payload["disclaimers"] if "easy test" not in d]
    return payload


def _drop_the_three_complexes_statement(payload: dict[str, Any]) -> dict[str, Any]:
    payload["disclaimers"] = [d for d in payload["disclaimers"] if "complexes only" not in d]
    return payload


def _drop_the_not_measured_binding_statement(payload: dict[str, Any]) -> dict[str, Any]:
    payload["disclaimers"] = [d for d in payload["disclaimers"] if "not measured binding" not in d]
    return payload


def _drop_the_no_tuning_statement(payload: dict[str, Any]) -> dict[str, Any]:
    payload["disclaimers"] = [d for d in payload["disclaimers"] if "retuned" not in d]
    return payload


def _drop_a_target(payload: dict[str, Any]) -> dict[str, Any]:
    payload["targets"] = [t for t in payload["targets"] if t["pdb_id"] != "2HU4"]
    return payload


def _add_an_unexpected_target(payload: dict[str, Any]) -> dict[str, Any]:
    extra = copy.deepcopy(payload["targets"][0])
    extra["pdb_id"] = "6LU7"
    payload["targets"].append(extra)
    return payload


def _no_modes(payload: dict[str, Any]) -> dict[str, Any]:
    target = payload["targets"][0]
    target["modes"] = []
    target["top1_rmsd"] = None
    target["best_mode_rmsd"] = None
    target["passed"] = False
    return payload


def _ten_modes(payload: dict[str, Any]) -> dict[str, Any]:
    target = payload["targets"][0]
    target["modes"] = [_mode(rank, -12.0 + rank, 1.0) for rank in range(1, 11)]
    return payload


def _wrong_rank(payload: dict[str, Any]) -> dict[str, Any]:
    payload["targets"][0]["modes"][1]["rank"] = 7
    return payload


def _unsorted_scores(payload: dict[str, Any]) -> dict[str, Any]:
    modes = payload["targets"][0]["modes"]
    modes[0]["score_kcal_per_mol"], modes[1]["score_kcal_per_mol"] = (
        modes[1]["score_kcal_per_mol"],
        modes[0]["score_kcal_per_mol"],
    )
    return payload


def _negative_rmsd(payload: dict[str, Any]) -> dict[str, Any]:
    payload["targets"][0]["modes"][0]["rmsd_to_crystal_angstrom"] = -0.5
    return payload


def _non_finite_score(payload: dict[str, Any]) -> dict[str, Any]:
    payload["targets"][0]["modes"][0]["score_kcal_per_mol"] = float("nan")
    return payload


def _non_finite_rmsd(payload: dict[str, Any]) -> dict[str, Any]:
    payload["targets"][0]["modes"][0]["rmsd_to_crystal_angstrom"] = float("inf")
    return payload


def _null_rmsd_in_an_ok_target(payload: dict[str, Any]) -> dict[str, Any]:
    payload["targets"][0]["modes"][0]["rmsd_to_crystal_angstrom"] = None
    return payload


def _top1_disagrees_with_the_rank_one_mode(payload: dict[str, Any]) -> dict[str, Any]:
    payload["targets"][0]["top1_rmsd"] = 1.20
    return payload


def _best_is_not_the_minimum(payload: dict[str, Any]) -> dict[str, Any]:
    # The modes have RMSDs 1.25 / 1.31 / 4.02, so 1.30 is not the smallest.
    payload["targets"][0]["best_mode_rmsd"] = 1.30
    return payload


def _passed_contradicts_the_rmsd(payload: dict[str, Any]) -> dict[str, Any]:
    payload["targets"][1]["passed"] = True
    return payload


def _failed_but_passing(payload: dict[str, Any]) -> dict[str, Any]:
    payload["targets"][2]["passed"] = True
    return payload


def _failed_without_a_reason(payload: dict[str, Any]) -> dict[str, Any]:
    payload["targets"][2]["failure_reason"] = ""
    return payload


def _failed_with_a_pass_metric(payload: dict[str, Any]) -> dict[str, Any]:
    payload["targets"][2]["top1_rmsd"] = 1.1
    return payload


def _failed_with_a_measured_mode(payload: dict[str, Any]) -> dict[str, Any]:
    target = payload["targets"][2]
    target["modes"] = [_mode(1, -8.0, 1.1)]
    return payload


def _per_target_tuning(payload: dict[str, Any]) -> dict[str, Any]:
    payload["targets"][1]["vina_params"]["seed"] = 99
    return payload


def _wrong_num_modes(payload: dict[str, Any]) -> dict[str, Any]:
    payload["targets"][0]["vina_params"]["num_modes"] = 11
    return payload


def _different_box(payload: dict[str, Any]) -> dict[str, Any]:
    payload["targets"][0]["box_size"] = [30.0, 30.0, 30.0]
    return payload


def _box_not_three_numbers(payload: dict[str, Any]) -> dict[str, Any]:
    payload["targets"][0]["box_size"] = [22.0, 22.0]
    return payload


def _missing_target_field(payload: dict[str, Any]) -> dict[str, Any]:
    del payload["targets"][0]["box_center"]
    return payload


def _missing_mode_field(payload: dict[str, Any]) -> dict[str, Any]:
    del payload["targets"][0]["modes"][0]["score_kcal_per_mol"]
    return payload


def _bad_status(payload: dict[str, Any]) -> dict[str, Any]:
    payload["targets"][0]["status"] = "computed"
    return payload


def _no_ligand_code(payload: dict[str, Any]) -> dict[str, Any]:
    payload["targets"][0]["ligand_code"] = None
    return payload


def _no_receptor_chains(payload: dict[str, Any]) -> dict[str, Any]:
    payload["targets"][0]["receptor_chains_used"] = []
    return payload


def _not_a_box_center(payload: dict[str, Any]) -> dict[str, Any]:
    payload["targets"][0]["box_center"] = ["0.383", "81.705", "109.195"]
    return payload


INVALID_PAYLOADS: list[tuple[Callable[[dict[str, Any]], dict[str, Any]], str]] = [
    (_wrong_schema, "schema_version must be"),
    (_drop_top_level_field, "missing top-level field 'disclaimers'"),
    (_drop_a_tool_version, "tool_versions is missing a usable 'numpy'"),
    (_wrong_threshold, "pass_threshold_angstrom must be 2.0"),
    (_wrong_pass_rule, "pass_rule must be"),
    (_not_a_timestamp, "not an ISO-8601 timestamp"),
    (_empty_disclaimers, "disclaimers must be a non-empty list"),
    (_drop_the_easy_test_statement, "easy test"),
    (_drop_the_three_complexes_statement, "three complexes only"),
    (_drop_the_not_measured_binding_statement, "not measured binding"),
    (_drop_the_no_tuning_statement, "no parameter was retuned"),
    (_drop_a_target, "expected PDB ID(s) ['2HU4'] are missing"),
    (_add_an_unexpected_target, "unexpected PDB ID(s) ['6LU7']"),
    (_no_modes, "modes must hold 1..9 docked poses"),
    (_ten_modes, "modes must hold 1..9 docked poses"),
    (_wrong_rank, "must have rank 2"),
    (_unsorted_scores, "sorted by ascending score"),
    (_negative_rmsd, "must be >= 0"),
    (_non_finite_score, "score_kcal_per_mol must be a finite number"),
    (_non_finite_rmsd, "rmsd_to_crystal_angstrom must be a finite number"),
    (_null_rmsd_in_an_ok_target, "every docked mode of a completed target must carry"),
    (_top1_disagrees_with_the_rank_one_mode, "is not the rank-1 mode's RMSD"),
    (_best_is_not_the_minimum, "is not the smallest mode RMSD"),
    (_passed_contradicts_the_rmsd, "passed must be False"),
    (_failed_but_passing, "a failed target must have passed false"),
    (_failed_without_a_reason, "must record why it failed"),
    (_failed_with_a_pass_metric, "must not report top1_rmsd"),
    (_failed_with_a_measured_mode, "a failed target must not report an RMSD"),
    (_per_target_tuning, "one fixed setting is used for every target"),
    (_wrong_num_modes, "num_modes must be 9"),
    (_different_box, "box_size.x must be 22.0"),
    (_box_not_three_numbers, "box_size must be three numbers"),
    (_missing_target_field, "is missing 'box_center'"),
    (_missing_mode_field, "modes[0] is missing 'score_kcal_per_mol'"),
    (_bad_status, "status must be 'ok' or 'failed'"),
    (_no_ligand_code, "ligand_code must be a non-empty string"),
    (_no_receptor_chains, "receptor_chains_used must be a non-empty list"),
    (_not_a_box_center, "box_center.x must be a finite number"),
]


@pytest.mark.parametrize(
    ("mutate", "expected"),
    INVALID_PAYLOADS,
    ids=[mutate.__name__.lstrip("_") for mutate, _ in INVALID_PAYLOADS],
)
def test_each_broken_rule_fails_with_a_clear_message(
    tmp_path: Path, mutate: Callable[[dict[str, Any]], dict[str, Any]], expected: str
) -> None:
    payload = mutate(_payload())

    result = _run(_write(tmp_path, payload))

    assert result.returncode != 0, f"expected a non-zero exit for {expected!r}"
    assert "FAIL:" in result.stderr
    assert expected in result.stderr, f"expected {expected!r} in:\n{result.stderr}"


def test_every_problem_is_reported_not_just_the_first(tmp_path: Path) -> None:
    payload = _payload()
    payload["pass_threshold_angstrom"] = 5.0
    payload["targets"][0]["top1_rmsd"] = 4.0

    result = _run(_write(tmp_path, payload))

    assert result.returncode != 0
    assert "problem(s) in" in result.stderr
    assert "pass_threshold_angstrom must be 2.0" in result.stderr
    assert "top1_rmsd" in result.stderr
