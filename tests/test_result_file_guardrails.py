"""Guardrail: a committed result file must pass its own validator.

The offline pipelines write JSON into ``app/data/``, and the validators are the
only thing standing between those files and the app: nothing else re-derives
their numbers. These tests make that check automatic. They do not create,
repair or fabricate a file -- a missing file is skipped with the reason, and a
**skip is not a pass**: it means there is nothing in the repository to check yet.

* ``app/data/docking/redock_validation.json`` is produced by the Colab run of
  ``scripts/redock_validation.py`` and is not in this repository, so its test
  skips here. The moment the file is committed, the test starts running.
* ``app/data/hamiltonians/*.json`` holds the offline PySCF artifacts. Every file
  that exists is validated; the committed H2 artifact is checked today.

The validators are imported from ``scripts/`` rather than run as subprocesses so
a failure can report every problem it found, not just an exit code.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
DOCKING_FILE = REPO_ROOT / "app" / "data" / "docking" / "redock_validation.json"
DOCKING_VALIDATOR = REPO_ROOT / "scripts" / "validate_docking_file.py"
HAMILTONIAN_DIR = REPO_ROOT / "app" / "data" / "hamiltonians"
HAMILTONIAN_VALIDATOR = REPO_ROOT / "scripts" / "validate_hamiltonian_file.py"


def _load_module(name: str, path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None, f"cannot load {path}"
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _hamiltonian_files() -> list[Path]:
    return sorted(HAMILTONIAN_DIR.glob("*.json"))


def test_docking_result_file_passes_its_validator_if_present() -> None:
    if not DOCKING_FILE.is_file():
        pytest.skip(
            f"{DOCKING_FILE.relative_to(REPO_ROOT)} does not exist: the redocking pipeline has "
            "not been run here (it needs Vina and meeko, which publish no Windows wheel), so "
            "there is no docking result to validate. Skip is not a pass."
        )

    validator = _load_module("validate_docking_file", DOCKING_VALIDATOR)
    try:
        payload: Any = json.loads(DOCKING_FILE.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:  # pragma: no cover - a broken file must fail
        pytest.fail(f"{DOCKING_FILE.relative_to(REPO_ROOT)} is not valid JSON: {exc}")

    problems = validator.validate(payload)
    assert not problems, (
        f"{DOCKING_FILE.relative_to(REPO_ROOT)} is committed but fails "
        f"{DOCKING_VALIDATOR.name}:\n"
        + "\n".join(f"  {index}. {problem}" for index, problem in enumerate(problems, start=1))
    )


@pytest.mark.parametrize(
    "artifact", _hamiltonian_files(), ids=lambda path: path.name or "none"
)
def test_hamiltonian_artifacts_pass_their_validator(artifact: Path, capsys: pytest.CaptureFixture[str]) -> None:
    if not artifact.is_file():
        pytest.skip(f"{artifact} does not exist; nothing to validate. Skip is not a pass.")

    validator = _load_module("validate_hamiltonian_file", HAMILTONIAN_VALIDATOR)
    exit_code = validator.validate(artifact)
    captured = capsys.readouterr()
    assert exit_code == 0, (
        f"{artifact.relative_to(REPO_ROOT)} is committed but fails "
        f"{HAMILTONIAN_VALIDATOR.name}:\n{captured.err or captured.out}"
    )
    assert "OK: every check passed." in captured.out


def test_hamiltonian_guardrail_is_not_vacuous(tmp_path: Path) -> None:
    """The guard above only means something if the validator can reject a file.

    A copy of the committed H2 artifact with one field deliberately broken must
    fail, so a green guard is evidence rather than a test that always passes.
    """
    existing = _hamiltonian_files()
    if not existing:
        pytest.skip("no committed Hamiltonian artifact to copy; nothing to corrupt. Skip is not a pass.")

    validator = _load_module("validate_hamiltonian_file", HAMILTONIAN_VALIDATOR)
    payload = json.loads(existing[0].read_text(encoding="utf-8"))
    assert validator.validate(existing[0]) == 0, "the original artifact must pass first"

    payload["n_qubits"] = 99
    broken = tmp_path / "broken_artifact.json"
    broken.write_text(json.dumps(payload), encoding="utf-8")

    assert validator.validate(broken) == 1, "a file with a wrong qubit count must be rejected"


def test_docking_result_guardrail_is_not_vacuous() -> None:
    """The docking validator is exercised even though no result file exists yet.

    A minimal payload that looks like a result -- three targets, plausible modes
    -- but drops every disclaimer must be rejected. This is the same rule the
    guard above enforces on a real file.
    """
    validator = _load_module("validate_docking_file", DOCKING_VALIDATOR)
    empty = {
        "schema_version": validator.SCHEMA_VERSION,
        "generated_at": "2026-10-18T00:00:00+00:00",
        "tool_versions": {
            "python": "3.11",
            "vina": "1.2.7",
            "meeko": "0.8.0",
            "rdkit": "2023.09",
            "numpy": "1.26",
        },
        "pass_threshold_angstrom": 2.0,
        "pass_rule": validator.EXPECTED_PASS_RULE,
        "disclaimers": [],
        "targets": [],
    }

    assert validator.validate(empty), "a payload with no disclaimers and no targets must be rejected"
