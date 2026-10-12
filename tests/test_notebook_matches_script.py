"""Keep the Colab notebook's embedded script in sync with the real script.

``scripts/precompute_hamiltonians.ipynb`` contains a ``%%writefile`` cell whose body
is meant to be a byte-for-byte copy of ``scripts/precompute_hamiltonians.py``. The
two drifted apart once already -- the cell lost its ``mol = gto.M(...)`` block while
the script kept it -- which would have made the notebook write a broken script.

The notebook also has to be internally consistent about *where* the script lives. A
fresh Colab session has no ``scripts/`` directory, so a ``%%writefile`` target of
``scripts/precompute_hamiltonians.py`` fails with FileNotFoundError while the run
cell executes ``precompute_hamiltonians.py`` from the working directory. That
specific mismatch shipped once, so it is asserted here too.

The checks are deliberately not vacuous: along with the real files they are exercised
against a faithful copy (must pass) and against mutated copies in a temporary
directory (must fail), so a broken comparison cannot pass by comparing nothing.
"""

from __future__ import annotations

import json
import re
from itertools import zip_longest
from pathlib import Path
from typing import Callable

import pytest

ROOT = Path(__file__).resolve().parents[1]
NOTEBOOK_PATH = ROOT / "scripts" / "precompute_hamiltonians.ipynb"
SCRIPT_PATH = ROOT / "scripts" / "precompute_hamiltonians.py"

EXPECTED_CELL_COUNT = 11
WRITEFILE_MAGIC = "%%writefile"
RUN_COMMAND_RE = re.compile(r"^\s*!\s*python\s+(?P<script>\S+)")


def _normalise_newlines(text: str) -> str:
    """Collapse CRLF/CR to LF so only real content differences remain."""
    return text.replace("\r\n", "\n").replace("\r", "\n")


def _line_or_missing(lines: list[str], index: int) -> str:
    return lines[index] if 0 <= index < len(lines) else "<no such line>"


def _writefile_cell_index(notebook: dict) -> int:
    """Index of the notebook's single ``%%writefile`` cell."""
    matches = [
        i
        for i, cell in enumerate(notebook["cells"])
        if "".join(cell.get("source", [])).startswith(WRITEFILE_MAGIC)
    ]
    assert len(matches) == 1, (
        f"expected exactly one {WRITEFILE_MAGIC!r} cell in the notebook, found {len(matches)}"
    )
    return matches[0]


def _writefile_cell_body(notebook: dict, index: int) -> str:
    """The cell body: the ``%%writefile`` line stripped, line endings normalised."""
    source = "".join(notebook["cells"][index]["source"])
    assert source.startswith(WRITEFILE_MAGIC), source[:60]
    return _normalise_newlines(source.split("\n", 1)[1])


def _writefile_target(notebook: dict) -> str:
    """The path the ``%%writefile`` cell writes to, as written in the magic line."""
    index = _writefile_cell_index(notebook)
    magic_line = "".join(notebook["cells"][index]["source"]).split("\n", 1)[0]
    parts = magic_line.split()
    assert len(parts) == 2, f"{WRITEFILE_MAGIC} cell must be '{WRITEFILE_MAGIC} <path>', got {magic_line!r}"
    return parts[1]


def _run_cell_script(notebook: dict) -> str:
    """The script reference the ``!python ...`` run cell executes, *as written*."""
    runs: list[str] = []
    for cell in notebook["cells"]:
        if cell.get("cell_type") != "code":
            continue
        for line in "".join(cell.get("source", [])).split("\n"):
            match = RUN_COMMAND_RE.match(line)
            if match:
                runs.append(match.group("script"))
    assert len(runs) == 1, (
        f"expected exactly one '!python <script> ...' run cell, found {runs}"
    )
    return runs[0]


def _script_ref(reference: str) -> str:
    """A script reference as written, with a leading './' removed.

    Only the redundant './' is normalised. A directory component is *not* stripped:
    ``scripts/precompute_hamiltonians.py`` and ``precompute_hamiltonians.py`` name
    different files, and conflating them is what let the Colab breakage ship.
    """
    reference = reference.strip()
    while reference.startswith("./"):
        reference = reference[2:]
    return reference


def _verify_cell_index(notebook: dict) -> int:
    """Index of the notebook's "Verify and download" cell, found by its download call."""
    matches = [
        i
        for i, cell in enumerate(notebook["cells"])
        if cell.get("cell_type") == "code" and "files.download" in "".join(cell.get("source", []))
    ]
    assert len(matches) == 1, (
        f"expected exactly one code cell calling files.download, found {len(matches)}"
    )
    return matches[0]


def assert_notebook_writefile_cell_matches_script(
    notebook_path: Path = NOTEBOOK_PATH, script_path: Path = SCRIPT_PATH
) -> None:
    """Raise ``AssertionError`` naming both files when they have drifted apart."""
    notebook = json.loads(notebook_path.read_text(encoding="utf-8"))
    index = _writefile_cell_index(notebook)
    body = _writefile_cell_body(notebook, index)
    script = _normalise_newlines(script_path.read_text(encoding="utf-8"))

    if body == script:
        return

    body_lines = body.splitlines()
    script_lines = script.splitlines()
    first = next(
        (
            n
            for n, (a, b) in enumerate(zip_longest(body_lines, script_lines), start=1)
            if a != b
        ),
        1,
    )
    raise AssertionError(
        f"{notebook_path} (cell {index}, its {WRITEFILE_MAGIC} cell) and {script_path} have "
        f"drifted: {len(body_lines)} body line(s) vs {len(script_lines)} script line(s), first "
        f"difference at line {first}.\n"
        f"Update {notebook_path} so the {WRITEFILE_MAGIC} cell matches {script_path} exactly "
        f"(the script is the source of truth), or update {script_path} if the cell is the one "
        f"that is correct. Keep both, and use the same line endings.\n"
        f"  notebook: {_line_or_missing(body_lines, first - 1)!r}\n"
        f"  script:   {_line_or_missing(script_lines, first - 1)!r}"
    )


def assert_notebook_runs_the_file_it_writes(
    notebook_path: Path = NOTEBOOK_PATH,
) -> None:
    """Raise ``AssertionError`` when the ``%%writefile`` target is not the file that runs.

    The two references are compared as written (only a redundant leading './' is
    normalised). ``scripts/precompute_hamiltonians.py`` and ``precompute_hamiltonians.py``
    are treated as different, because that mismatch is exactly what broke a fresh
    Colab session: the writefile target went into a directory that does not exist
    while the run cell executed the working-directory copy.
    """
    notebook = json.loads(notebook_path.read_text(encoding="utf-8"))
    target_raw = _writefile_target(notebook)
    run_raw = _run_cell_script(notebook)

    if _script_ref(target_raw) == _script_ref(run_raw):
        return

    directory_hint = ""
    if "/" in target_raw or "\\" in target_raw:
        directory_hint = (
            f" The {WRITEFILE_MAGIC} target {target_raw!r} also carries a directory component; "
            "a fresh Colab session has no such directory, so name the script in the working "
            "directory instead and keep the repository path only where the script is described."
        )

    raise AssertionError(
        f"{notebook_path} is internally inconsistent: the {WRITEFILE_MAGIC} cell writes "
        f"{target_raw!r} but the run cell executes {run_raw!r}.\n"
        f"Update {notebook_path} so the {WRITEFILE_MAGIC} target and the '!python ...' run cell "
        f"name the same working-directory-relative file."
        + directory_hint
    )


MAGIC_PREFIXES = ("!", "%")


def _strip_magic_lines(source: str) -> str:
    """Remove IPython ``!``/``%`` lines so what remains is plain Python."""
    return "\n".join(
        line for line in source.split("\n") if not line.lstrip().startswith(MAGIC_PREFIXES)
    )


def assert_every_code_cell_parses(notebook_path: Path = NOTEBOOK_PATH) -> None:
    """Raise ``AssertionError`` on the first code cell that is not valid Python.

    ``!``/``%`` lines are excluded first, so anything reported here is real Python
    source -- a stray leading character from a bad copy/paste or a leaked diff
    marker (a ``-`` before ``import json, pathlib`` in the verify cell shipped
    once) shows up as a ``SyntaxError``.
    """
    notebook = json.loads(notebook_path.read_text(encoding="utf-8"))
    for index, cell in enumerate(notebook["cells"]):
        if cell.get("cell_type") != "code":
            continue
        source = _strip_magic_lines("".join(cell.get("source", [])))
        try:
            compile(source, f"{notebook_path.name}[cell {index}]", "exec")
        except SyntaxError as exc:
            lines = source.split("\n")
            offender = (
                lines[exc.lineno - 1]
                if exc.lineno is not None and exc.lineno <= len(lines)
                else "<no such line>"
            )
            raise AssertionError(
                f"{notebook_path} code cell {index} does not parse "
                f"({type(exc).__name__}): {exc.msg} at line "
                f"{exc.lineno} ({offender!r}). Lines starting with '!' or '%' are excluded "
                f"before compiling, so this is broken Python source -- check for a stray "
                f"leading character (a leaked diff marker) or a bad paste."
            ) from exc


def test_notebook_is_valid_json_with_eleven_cells() -> None:
    try:
        notebook = json.loads(NOTEBOOK_PATH.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        pytest.fail(f"{NOTEBOOK_PATH} is not valid JSON: {exc}")

    assert isinstance(notebook, dict)
    cells = notebook["cells"]
    assert len(cells) == EXPECTED_CELL_COUNT, (
        f"{NOTEBOOK_PATH} should contain {EXPECTED_CELL_COUNT} cells, found {len(cells)}"
    )


def test_writefile_cell_body_matches_the_script() -> None:
    assert_notebook_writefile_cell_matches_script()


def test_writefile_target_matches_the_run_cell() -> None:
    assert_notebook_runs_the_file_it_writes()


def test_every_code_cell_parses() -> None:
    assert_every_code_cell_parses()


def _copy_pair(tmp_path: Path) -> tuple[Path, Path]:
    """A faithful copy of both files, so mutations never touch the real ones."""
    notebook_path = tmp_path / NOTEBOOK_PATH.name
    script_path = tmp_path / SCRIPT_PATH.name
    notebook_path.write_bytes(NOTEBOOK_PATH.read_bytes())
    script_path.write_bytes(SCRIPT_PATH.read_bytes())
    return notebook_path, script_path


def test_a_faithful_copy_passes(tmp_path: Path) -> None:
    """Control: both checks pass on an unmodified copy."""
    notebook_path, script_path = _copy_pair(tmp_path)
    assert_notebook_writefile_cell_matches_script(notebook_path, script_path)
    assert_notebook_runs_the_file_it_writes(notebook_path)


def _drop_final_body_line(source: list[str]) -> list[str]:
    return source[:-1]


def _append_marker_to_first_body_line(source: list[str]) -> list[str]:
    mutated = list(source)
    mutated[1] = mutated[1].rstrip("\n") + "  # drifted\n"
    return mutated


def _add_extra_trailing_newline(source: list[str]) -> list[str]:
    mutated = list(source)
    mutated[-1] = mutated[-1] + "\n"
    return mutated


MUTATIONS = [
    pytest.param(_drop_final_body_line, id="missing-final-line"),
    pytest.param(_append_marker_to_first_body_line, id="altered-line"),
    pytest.param(_add_extra_trailing_newline, id="extra-trailing-newline"),
]


def _write_mutated_notebook(
    tmp_path: Path, mutate: Callable[[list[str]], list[str]]
) -> Path:
    """Copy the notebook into *tmp_path* with the writefile cell body mutated."""
    notebook_path, _ = _copy_pair(tmp_path)
    notebook = json.loads(NOTEBOOK_PATH.read_text(encoding="utf-8"))
    index = _writefile_cell_index(notebook)
    notebook["cells"][index]["source"] = mutate(list(notebook["cells"][index]["source"]))
    notebook_path.write_text(
        json.dumps(notebook, indent=2) + "\n", encoding="utf-8", newline="\n"
    )
    return notebook_path


@pytest.mark.parametrize("mutate", MUTATIONS)
def test_a_drifted_copy_fails(
    tmp_path: Path, mutate: Callable[[list[str]], list[str]]
) -> None:
    """The body check must fail on a drifted copy, or it proves nothing."""
    notebook_path = _write_mutated_notebook(tmp_path, mutate)
    script_path = tmp_path / SCRIPT_PATH.name

    with pytest.raises(AssertionError) as excinfo:
        assert_notebook_writefile_cell_matches_script(notebook_path, script_path)

    message = str(excinfo.value)
    # The message must name both files so it is clear which one to update.
    assert str(notebook_path) in message
    assert str(script_path) in message


def test_writefile_target_mismatch_is_detected(tmp_path: Path) -> None:
    """The exact regression that broke Colab must be caught by the path check.

    Mutating the ``%%writefile`` target to the directory-qualified name leaves the
    cell *body* untouched (the body check still passes) but breaks the notebook's
    internal consistency, which is what the path check exists to catch.
    """
    notebook_path, script_path = _copy_pair(tmp_path)
    notebook = json.loads(NOTEBOOK_PATH.read_text(encoding="utf-8"))
    index = _writefile_cell_index(notebook)
    notebook["cells"][index]["source"][0] = "%%writefile scripts/precompute_hamiltonians.py\n"
    notebook_path.write_text(
        json.dumps(notebook, indent=2) + "\n", encoding="utf-8", newline="\n"
    )

    # The body is unchanged, so the body check alone would not notice.
    assert_notebook_writefile_cell_matches_script(notebook_path, script_path)

    with pytest.raises(AssertionError) as excinfo:
        assert_notebook_runs_the_file_it_writes(notebook_path)

    message = str(excinfo.value)
    assert "scripts/precompute_hamiltonians.py" in message
    assert "precompute_hamiltonians.py" in message
    assert str(notebook_path) in message


def test_a_stray_leading_dash_is_detected(tmp_path: Path) -> None:
    """A leaked diff marker in the verify cell must fail the parse check.

    The stray '-' that produced ``SyntaxError: invalid syntax`` on ``-import json,
    pathlib`` is reproduced exactly: one character prepended to the cell's first
    line, in a temporary copy.
    """
    notebook_path, _ = _copy_pair(tmp_path)
    notebook = json.loads(NOTEBOOK_PATH.read_text(encoding="utf-8"))
    index = _verify_cell_index(notebook)
    notebook["cells"][index]["source"][0] = "-" + notebook["cells"][index]["source"][0]
    first_line = notebook["cells"][index]["source"][0].rstrip("\n")
    notebook_path.write_text(
        json.dumps(notebook, indent=2) + "\n", encoding="utf-8", newline="\n"
    )

    with pytest.raises(AssertionError) as excinfo:
        assert_every_code_cell_parses(notebook_path)

    message = str(excinfo.value)
    assert f"cell {index}" in message
    assert first_line in message
    assert "SyntaxError".lower() in message.lower()
