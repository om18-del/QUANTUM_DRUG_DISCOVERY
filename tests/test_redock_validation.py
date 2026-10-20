"""Tests for the offline redocking pipeline in ``scripts/redock_validation.py``.

Nothing here needs the network, AutoDock Vina or meeko: the script is expected to
import on a host where none of them exist (that is ``test_import_safety``), and
everything else is either a pure helper fed hand-written PDB/PDBQT text or an
RDKit-gated check of the preparation and RMSD path, which skip where RDKit is
blocked (``tests/availability.py``).

The RMSD checks are the point of the file. They use a **synthetic** pose built
from the crystal coordinates of a hand-written ligand block, which is the only
way to prove three separate claims without running a docking program:

* the atom correspondence between a pose and the ligand molecule is verified
  rather than assumed (a permuted pose must fail, not return a number),
* the RMSD does **not** re-align (a 2.5 A translation must come back as 2.5 A),
* the RMSD **is** symmetry-aware (a benzene ring written with its atoms rotated
  by one position must come back as 0, while the atom-by-atom distance is 1.39 A).

No number in this file is a docking result and none of it is written to any
result file.
"""

from __future__ import annotations

import importlib.util
import json
import math
import subprocess
import sys
import textwrap
from pathlib import Path
from typing import Any

import pytest

from tests.availability import require_rdkit_chem

ROOT = Path(__file__).resolve().parents[1]
REDOCK_SCRIPT = ROOT / "scripts" / "redock_validation.py"
VALIDATOR_SCRIPT = ROOT / "scripts" / "validate_docking_file.py"


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None, f"cannot load {path}"
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


redock = _load("redock_validation", REDOCK_SCRIPT)
validator = _load("validate_docking_file", VALIDATOR_SCRIPT)


# --------------------------------------------------------------------------- #
# import safety: the script must import with none of its heavy dependencies
# --------------------------------------------------------------------------- #

BLOCKED_MODULES = ("vina", "meeko", "rdkit", "requests", "numpy")

IMPORT_SAFETY_PROGRAM = textwrap.dedent(
    """
    import importlib.util
    import sys

    BLOCKED = %r


    class Blocker:
        \"\"\"Refuse to import the heavy dependencies, as a Windows host without
        Vina/meeko/RDKit would.\"\"\"

        def find_spec(self, fullname, path=None, target=None):
            if fullname.split(".")[0] in BLOCKED:
                raise ImportError(f"{fullname} is blocked in this test")
            return None


    sys.meta_path.insert(0, Blocker())
    for name in BLOCKED:
        sys.modules.pop(name, None)

    spec = importlib.util.spec_from_file_location("redock_validation", sys.argv[1])
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    leaked = sorted(name for name in sys.modules if name.split(".")[0] in BLOCKED)
    print("LEAKED:", leaked)
    print("TARGETS:", ",".join(module.EXPECTED_TARGETS))
    print("THRESHOLD:", module.PASS_THRESHOLD_ANGSTROM)
    print("BOX:", module.BOX_SIZE_ANGSTROM)
    print("DISCLAIMERS:", len(module.DISCLAIMERS))

    try:
        module._require("vina")
    except SystemExit as exc:
        message = str(exc)
        print("REQUIRE_SAYS_PIP:", module.PIP_INSTALL_LINE in message)
        print("REQUIRE_SAYS_WINDOWS:", "Windows" in message)
    """
) % (BLOCKED_MODULES,)


def test_import_safety_without_vina_meeko_rdkit_requests_numpy() -> None:
    """The module imports, and imports none of its heavy dependencies, without them.

    This is the whole point of keeping every heavy import inside a function: the
    Windows CI runner has no Vina, no meeko and (sometimes) no usable RDKit, and
    the test suite has to stay green there.
    """
    result = subprocess.run(
        [sys.executable, "-c", IMPORT_SAFETY_PROGRAM, str(REDOCK_SCRIPT)],
        capture_output=True,
        text=True,
        cwd=str(ROOT),
    )

    assert result.returncode == 0, f"stderr:\n{result.stderr}"
    assert "LEAKED: []" in result.stdout, result.stdout
    assert "TARGETS: 1VRT,2HU4,1HSG" in result.stdout
    assert f"THRESHOLD: {redock.PASS_THRESHOLD_ANGSTROM}" in result.stdout
    assert "DISCLAIMERS: 5" in result.stdout
    # The error path must tell the user what to install and why, not just fail.
    assert "REQUIRE_SAYS_PIP: True" in result.stdout
    assert "REQUIRE_SAYS_WINDOWS: True" in result.stdout


def test_no_heavy_import_at_module_level_in_the_source() -> None:
    """A readable guard: no top-level ``import requests/rdkit/vina/meeko/numpy``."""
    offenders = []
    for line in REDOCK_SCRIPT.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if not line[:1].strip():  # indented -> inside a function
            continue
        for module in BLOCKED_MODULES:
            if stripped.startswith(f"import {module}") or stripped.startswith(
                f"from {module}"
            ):
                offenders.append(stripped)
    assert not offenders, f"module-level heavy import(s): {offenders}"


# --------------------------------------------------------------------------- #
# hand-written fixtures
# --------------------------------------------------------------------------- #

# Nicotinic acid (C6H5NO2, 9 heavy atoms), coordinates embedded from a single
# deterministic RDKit embedding. Written as HETATM in chain L, residue 1.
NICOTINIC_BLOCK = textwrap.dedent(
    """\
    HETATM    1  O1  LIG L   1       2.767   0.703   0.457  1.00  0.00           O
    HETATM    2  C1  LIG L   1       1.865  -0.138  -0.180  1.00  0.00           C
    HETATM    3  O2  LIG L   1       2.337  -1.053  -0.899  1.00  0.00           O
    HETATM    4  C2  LIG L   1       0.421   0.016  -0.037  1.00  0.00           C
    HETATM    5  C3  LIG L   1      -0.386  -0.878  -0.717  1.00  0.00           C
    HETATM    6  C4  LIG L   1      -1.756  -0.700  -0.555  1.00  0.00           C
    HETATM    7  C5  LIG L   1      -2.252   0.316   0.242  1.00  0.00           C
    HETATM    8  N1  LIG L   1      -1.426   1.177   0.897  1.00  0.00           N
    HETATM    9  C6  LIG L   1      -0.084   1.026   0.757  1.00  0.00           C
    """
)
NICOTINIC_SMILES = "OC(=O)c1cccnc1"

# Benzene on an exact hexagon of radius 1.39 A, so a ring permutation is exactly
# a 60 degree rotation of the molecule.
BENZENE_RADIUS = 1.39
BENZENE_ATOMS = [
    (
        BENZENE_RADIUS * math.cos(math.radians(60 * index)),
        BENZENE_RADIUS * math.sin(math.radians(60 * index)),
        0.0,
    )
    for index in range(6)
]
BENZENE_BLOCK = "".join(
    f"HETATM{index + 1:5d}  C   BEN B   1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C\n"
    for index, (x, y, z) in enumerate(BENZENE_ATOMS)
)
BENZENE_SMILES = "c1ccccc1"

PROTEIN_AND_LIGAND = textwrap.dedent(
    """\
    ATOM      1  N   GLY A   1       0.000   0.000   0.000  1.00  0.00           N
    ATOM      2  CA  GLY A   1       1.000   0.000   0.000  1.00  0.00           C
    ATOM      3  N   GLY B   1      40.000  40.000  40.000  1.00  0.00           N
    ATOM      4  C   GLY C   1      10.000  10.000  10.000  1.00  0.00           C
    HETATM    5  C1  LIG A 500       0.500   0.500   0.500  1.00  0.00           C
    HETATM    6  O1  LIG A 500       1.500   0.500   0.500  1.00  0.00           O
    HETATM    7  O   HOH A 600       9.000   9.000   9.000  1.00  0.00           O
    HETATM    8  MG  MG  A 601      20.000  20.000  20.000  1.00  0.00          MG
    """
)

# The same residue in two chains with alternate locations: the second instance
# has an 'A' and a 'B' conformer plus one atom with no altLoc at all.
INSTANCES_BLOCK = textwrap.dedent(
    """\
    ATOM      1  N   GLY A   1       0.000   0.000   0.000  1.00  0.00           N
    HETATM    2  C1  LIG A 500       0.500   0.500   0.500  1.00  0.00           C
    HETATM    3  C2  LIG A 500       1.500   0.500   0.500  1.00  0.00           C
    HETATM    4  C1 ALIG B 500       5.500   5.500   5.500  1.00  0.00           C
    HETATM    5  C2 ALIG B 500       6.500   5.500   5.500  1.00  0.00           C
    HETATM    6  C3 BLIG B 500       7.500   5.500   5.500  1.00  0.00           C
    HETATM    7  C4  LIG B 500       8.500   5.500   5.500  1.00  0.00           C
    """
)


def pdbqt_line(
    serial: int,
    name: str,
    x: float,
    y: float,
    z: float,
    atom_type: str,
    *,
    residue: str = "LIG",
    chain: str = "A",
    resseq: int = 901,
) -> str:
    """One PDBQT ATOM record with the AutoDock type in the trailing column."""
    prefix = (
        f"ATOM  {serial:5d} {name:<4} {residue:>3} {chain}{resseq:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00     0.000"
    )
    return prefix.ljust(78) + f"{atom_type:>2}"


ELEMENT_TYPE = {"C": "C", "N": "N", "O": "OA", "S": "SA", "H": "HD"}


def pose_from_molecule(mol, *, shift=(0.0, 0.0, 0.0), order=None) -> str:
    """A PDBQT pose holding ``mol``'s heavy atoms, in molecule order or ``order``."""
    conformer = mol.GetConformer()
    indices = list(range(mol.GetNumAtoms())) if order is None else list(order)
    lines = []
    for serial, index in enumerate(indices, start=1):
        atom = mol.GetAtomWithIdx(index)
        position = conformer.GetAtomPosition(index)
        lines.append(
            pdbqt_line(
                serial,
                f"{atom.GetSymbol()}{index + 1}",
                position.x + shift[0],
                position.y + shift[1],
                position.z + shift[2],
                ELEMENT_TYPE[atom.GetSymbol()],
            )
        )
    return "\n".join(lines) + "\n"


def _nicotinic_ligand():
    """(heavy_mol, prep_mol, info) for the hand-written nicotinic acid block."""
    instance = redock.ligand_instances(NICOTINIC_BLOCK, "LIG")[0]
    definition = {
        "comp_id": "LIG",
        "smiles": NICOTINIC_SMILES,
        "smiles_field": "synthetic test definition",
        "atom_count_heavy": 9,
    }
    return redock.prepare_crystal_ligand("LIG", instance, definition)


# --------------------------------------------------------------------------- #
# pure helpers
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(
    ("raw", "expected"),
    [
        ("1VRT,2HU4,1HSG", ["1VRT", "2HU4", "1HSG"]),
        ("1vrt, 2hu4", ["1VRT", "2HU4"]),
        ("1VRT", ["1VRT"]),
    ],
)
def test_parse_targets(raw: str, expected: list[str]) -> None:
    assert redock.parse_targets(raw) == expected


@pytest.mark.parametrize("raw", ["", "  ", "1VR", "ABCDE", "1VRT,1VRT", "1VRT,,1HSG-"])
def test_parse_targets_rejects_bad_lists(raw: str) -> None:
    with pytest.raises(ValueError):
        redock.parse_targets(raw)


def test_ligand_centroid_is_the_mean_of_the_ligand_atoms() -> None:
    atoms = [
        {"x": 0.0, "y": 0.0, "z": 0.0},
        {"x": 2.0, "y": 4.0, "z": -6.0},
    ]
    assert redock.ligand_centroid(atoms) == [1.0, 2.0, -3.0]


def test_ligand_centroid_refuses_an_empty_ligand() -> None:
    with pytest.raises(redock.TargetFailure):
        redock.ligand_centroid([])


def test_pdb_atoms_reads_one_model_only() -> None:
    text = (
        PROTEIN_AND_LIGAND
        + "ENDMDL\n"
        + "ATOM      9  N   GLY A   9      99.000  99.000  99.000  1.00  0.00           N\n"
    )
    atoms = redock.pdb_atoms(text, record="ATOM")
    assert [atom["serial"] for atom in atoms] == ["1", "2", "3", "4"]


def test_ligand_instances_finds_every_instance_in_file_order() -> None:
    instances = redock.ligand_instances(INSTANCES_BLOCK, "LIG")
    # The second instance keeps its 'A' conformer (C1, C2) and the shared atom
    # with no altLoc (C4); the 'B' conformer (C3) is dropped.
    assert [(i["chain"], i["resseq"], len(i["atoms"])) for i in instances] == [
        ("A", "500", 2),
        ("B", "500", 3),
    ]


def test_ligand_instances_prefers_the_first_altloc_and_reports_it() -> None:
    """A residue with altLocs 'A'/'B' and one atom without any.

    The blank-altLoc atom is kept because a deposited residue routinely leaves
    shared atoms without an altLoc; the 'A' conformer is chosen for the rest and
    the choice is recorded on the instance.
    """
    instance = redock.ligand_instances(INSTANCES_BLOCK, "LIG")[1]
    assert instance["altloc"] == "A"
    assert [atom["name"] for atom in instance["atoms"]] == ["C1", "C2", "C4"]


def test_ligand_instances_records_a_blank_altloc_for_a_residue_without_one() -> None:
    instance = redock.ligand_instances(INSTANCES_BLOCK, "LIG")[0]
    assert instance["altloc"] == " "
    assert [atom["name"] for atom in instance["atoms"]] == ["C1", "C2"]


def test_ligand_instances_falls_back_to_the_only_altloc_present() -> None:
    text = (
        "HETATM    1  C1 BLIG A 500       5.500   5.500   5.500  1.00  0.00           C\n"
        "HETATM    2  C2 BLIG A 500       6.500   5.500   5.500  1.00  0.00           C\n"
    )
    instance = redock.ligand_instances(text, "LIG")[0]
    assert instance["altloc"] == "B"
    assert len(instance["atoms"]) == 2


def test_ligand_instances_counts_deposited_hydrogens() -> None:
    text = (
        "HETATM    1  C1  LIG A 500       5.500   5.500   5.500  1.00  0.00           C\n"
        "HETATM    2  H1  LIG A 500       8.500   5.500   5.500  1.00  0.00           H\n"
    )
    assert redock.ligand_instances(text, "LIG")[0]["hydrogen_atoms_in_deposit"] == 1


def test_ligand_instances_ignores_other_residues() -> None:
    assert redock.ligand_instances(PROTEIN_AND_LIGAND, "ZZZ") == []


def test_chains_near_ligand_keeps_only_the_chains_touching_the_ligand() -> None:
    instances = redock.ligand_instances(PROTEIN_AND_LIGAND, "LIG")
    chains = redock.chains_near_ligand(PROTEIN_AND_LIGAND, instances[0]["atoms"])
    # Chain B is 40 A away and chain C is 10 A away; the default cutoff is 5 A.
    assert chains == ["A"]


def test_chains_near_ligand_cutoff_is_inclusive_and_documented() -> None:
    instances = redock.ligand_instances(PROTEIN_AND_LIGAND, "LIG")
    wide = redock.chains_near_ligand(
        PROTEIN_AND_LIGAND, instances[0]["atoms"], cutoff=20.0
    )
    assert wide == ["A", "C"]


def test_chains_near_ligand_returns_nothing_when_no_chain_is_close() -> None:
    instances = redock.ligand_instances(PROTEIN_AND_LIGAND, "LIG")
    assert (
        redock.chains_near_ligand(PROTEIN_AND_LIGAND, instances[0]["atoms"], cutoff=0.1)
        == []
    )


def test_receptor_pdb_text_keeps_only_the_selected_protein_chains() -> None:
    text = redock.receptor_pdb_text(PROTEIN_AND_LIGAND, ["A"])
    atoms = redock.pdb_atoms(text, record="ATOM")
    assert [atom["serial"] for atom in atoms] == ["1", "2"]
    # Every HETATM goes: waters, ions and the ligand itself.
    assert redock.pdb_atoms(text, record="HETATM") == []
    assert text.rstrip().endswith("END")
    assert "LIG" not in text and "HOH" not in text


def test_parse_pdbqt_atoms_reads_the_trailing_autodock_type() -> None:
    text = (
        pdbqt_line(1, "C1", 1.0, 2.0, 3.0, "A")
        + "\n"
        + pdbqt_line(2, "O1", 4.0, 5.0, 6.0, "OA")
        + "\n"
    )
    atoms = redock.parse_pdbqt_atoms(text)
    assert [(a["serial"], a["type"]) for a in atoms] == [("1", "A"), ("2", "OA")]
    assert (atoms[0]["x"], atoms[0]["y"], atoms[0]["z"]) == (1.0, 2.0, 3.0)


def test_parse_pdbqt_atoms_reads_a_two_letter_type() -> None:
    text = pdbqt_line(1, "CL1", 1.0, 2.0, 3.0, "Cl") + "\n"
    assert redock.parse_pdbqt_atoms(text)[0]["type"] == "Cl"


def test_parse_pdbqt_atoms_reports_an_unreadable_type_as_empty() -> None:
    # The type column is cut off, so nothing alphabetic is left to read.
    line = pdbqt_line(1, "C1", 1.0, 2.0, 3.0, "C")[:77].rstrip()
    atoms = redock.parse_pdbqt_atoms(line + "\n")
    assert atoms[0]["type"] == ""
    assert redock.element_of_pdbqt_type(atoms[0]["type"]) is None


def test_parse_pdbqt_atoms_falls_back_to_the_last_alphabetic_field() -> None:
    line = pdbqt_line(1, "C1", 1.0, 2.0, 3.0, "C")[:77].rstrip() + " A"
    assert redock.parse_pdbqt_atoms(line + "\n")[0]["type"] == "A"


def test_parse_pdbqt_models_splits_vina_pose_models() -> None:
    text = (
        f"MODEL 1\n{pdbqt_line(1, 'C1', 1.0, 2.0, 3.0, 'C')}\nENDMDL\n"
        f"MODEL 2\n{pdbqt_line(1, 'C1', 4.0, 5.0, 6.0, 'C')}\n"
        f"{pdbqt_line(2, 'O1', 7.0, 8.0, 9.0, 'OA')}\nENDMDL\n"
    )
    models = redock.parse_pdbqt_models(text)
    assert [len(model) for model in models] == [1, 2]
    assert models[1][0]["x"] == 4.0


def test_parse_pdbqt_models_treats_a_file_without_models_as_one_pose() -> None:
    models = redock.parse_pdbqt_models(pdbqt_line(1, "C1", 1.0, 2.0, 3.0, "C") + "\n")
    assert len(models) == 1 and len(models[0]) == 1


def test_count_pdbqt_types_reports_unreadable_types() -> None:
    atoms = [{"type": "OA"}, {"type": "OA"}, {"type": ""}]
    assert redock.count_pdbqt_types(atoms) == {"OA": 2, "unreadable": 1}


@pytest.mark.parametrize(
    ("atom_type", "element"),
    [("C", "C"), ("A", "C"), ("N", "N"), ("NA", "N"), ("OA", "O"), ("SA", "S"), ("HD", "H")],
)
def test_element_of_pdbqt_type_maps_the_documented_types(atom_type: str, element: str) -> None:
    assert redock.element_of_pdbqt_type(atom_type) == element


def test_element_of_pdbqt_type_is_none_for_an_unknown_type() -> None:
    assert redock.element_of_pdbqt_type("Q") is None
    assert redock.element_of_pdbqt_type("") is None


# --------------------------------------------------------------------------- #
# target records: a failure is a failure, and the payload says so
# --------------------------------------------------------------------------- #


def test_empty_target_record_has_every_key_and_no_result() -> None:
    record = redock.empty_target_record("1VRT", exhaustiveness=16, seed=7)
    assert set(validator.TARGET_FIELDS) <= set(record)
    assert record["status"] == "failed"
    assert record["passed"] is False
    assert record["modes"] == []
    assert record["top1_rmsd"] is None
    assert record["best_mode_rmsd"] is None
    assert record["failure_reason"] is None
    assert record["box_size"] == list(redock.BOX_SIZE_ANGSTROM)


def test_run_target_records_a_target_failure_instead_of_raising(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    def explode(pdb_id: str):
        raise redock.TargetFailure(f"metadata for {pdb_id} is unavailable")

    monkeypatch.setattr(redock, "entry_metadata", explode)
    record = redock.run_target("1VRT", exhaustiveness=16, seed=7, work_dir=tmp_path)

    assert record["status"] == "failed"
    assert record["passed"] is False
    assert record["failure_type"] == "TargetFailure"
    assert record["failure_reason"] == "metadata for 1VRT is unavailable"
    assert record["modes"] == [] and record["top1_rmsd"] is None


def test_run_target_records_an_unexpected_exception_as_a_failure(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """A bug is still written out as a failed target, named as itself."""

    def explode(pdb_id: str):
        raise RuntimeError("boom")

    monkeypatch.setattr(redock, "entry_metadata", explode)
    record = redock.run_target("1HSG", exhaustiveness=16, seed=7, work_dir=tmp_path)

    assert record["status"] == "failed"
    assert record["failure_type"] == "RuntimeError"
    assert record["failure_reason"] == "unexpected RuntimeError: boom"


def test_format_summary_table_shows_failures_and_truncates_long_reasons() -> None:
    ok = redock.empty_target_record("1VRT", exhaustiveness=16, seed=7)
    ok.update(
        {
            "ligand_code": "NVP",
            "receptor_chains_used": ["A"],
            "modes": [
                {"rank": 1, "score_kcal_per_mol": -9.5, "rmsd_to_crystal_angstrom": 1.25}
            ],
            "top1_rmsd": 1.25,
            "best_mode_rmsd": 1.25,
            "passed": True,
            "status": "ok",
        }
    )
    failed = redock.empty_target_record("2HU4", exhaustiveness=16, seed=7)
    failed["failure_reason"] = "x" * 200

    table = redock.format_summary_table([ok, failed])

    assert "1VRT" in table and "NVP" in table and "1.250" in table
    assert "2HU4" in table and "failed" in table
    assert "..." in table
    assert max(len(line) for line in table.splitlines()) < 200


def _synthetic_payload() -> dict[str, Any]:
    """A payload shaped like the writer's, with hand-written, consistent numbers."""
    records = []
    for pdb_id, ligand, top1, scores in (
        ("1VRT", "NVP", 1.25, [-9.5, -9.1, -8.8]),
        ("2HU4", "G39", 2.75, [-7.2, -7.0, -6.9]),
    ):
        record = redock.empty_target_record(pdb_id, exhaustiveness=16, seed=7)
        record.update(
            {
                "ligand_code": ligand,
                "receptor_chains_used": ["A"],
                "box_center": [1.0, 2.0, 3.0],
                "modes": [
                    {
                        "rank": rank,
                        "score_kcal_per_mol": score,
                        "rmsd_to_crystal_angstrom": round(top1 + 0.1 * (rank - 1), 3),
                    }
                    for rank, score in enumerate(scores, start=1)
                ],
            }
        )
        record["top1_rmsd"] = record["modes"][0]["rmsd_to_crystal_angstrom"]
        record["best_mode_rmsd"] = min(
            mode["rmsd_to_crystal_angstrom"] for mode in record["modes"]
        )
        record["passed"] = record["top1_rmsd"] < redock.PASS_THRESHOLD_ANGSTROM
        record["status"] = "ok"
        records.append(record)

    failed = redock.empty_target_record("1HSG", exhaustiveness=16, seed=7)
    failed["failure_reason"] = "no non-polymer entity passes the ligand rule"
    failed["failure_type"] = "TargetFailure"
    records.append(failed)

    return redock.build_payload(
        records,
        versions={"python": "3.11.9", "numpy": "1.26.4", "rdkit": "2026.03.6", "vina": "1.2.7", "meeko": "0.8.0"},
        exhaustiveness=16,
        seed=7,
    )


def test_build_payload_shape_and_provenance() -> None:
    records = [
        redock.empty_target_record(pdb_id, exhaustiveness=16, seed=7)
        for pdb_id in redock.EXPECTED_TARGETS
    ]
    payload = redock.build_payload(
        records,
        versions={"python": "3.11.9", "numpy": "1", "rdkit": "1", "vina": "1", "meeko": "1"},
        exhaustiveness=16,
        seed=7,
    )

    for field in validator.TOP_LEVEL_FIELDS:
        assert field in payload, field
    assert payload["schema_version"] == redock.SCHEMA_VERSION
    assert payload["pass_threshold_angstrom"] == 2.0
    assert payload["pass_rule"] == "top-ranked pose RMSD < 2.0 A"
    assert payload["box_size_angstrom"] == [22.0, 22.0, 22.0]
    assert payload["vina_params"]["one_setting_for_all_targets"] is True
    assert payload["vina_params"]["num_modes"] == 9
    assert payload["generated_at"].endswith("Z")
    # The records are carried through, not recomputed.
    assert payload["targets"] == records
    # The easy-setup note must be stated in the file itself, not only in the docs.
    assert "easy" in payload["setup"]["note"].lower()
    assert "true" in payload["setup"]["box_center"]


def test_build_payload_disclaimers_state_every_required_thing() -> None:
    payload = _synthetic_payload()
    text = " ".join(payload["disclaimers"]).lower()
    assert "not measured binding" in text
    assert "easy test" in text
    assert "complexes only" in text
    assert "retuned" in text
    assert "prediction" in text


def test_writer_and_validator_agree_on_the_contract() -> None:
    """The two scripts must not drift apart on the shared rules."""
    assert redock.SCHEMA_VERSION == validator.SCHEMA_VERSION
    assert redock.EXPECTED_TARGETS == validator.EXPECTED_PDB_IDS
    assert redock.PASS_THRESHOLD_ANGSTROM == validator.EXPECTED_PASS_THRESHOLD
    assert redock.PASS_RULE == validator.EXPECTED_PASS_RULE
    assert redock.BOX_SIZE_ANGSTROM == validator.EXPECTED_BOX_SIZE
    assert redock.NUM_MODES == validator.NUM_MODES
    assert redock.NUM_MODES == validator.MAX_MODES


def test_a_writer_payload_passes_the_validator() -> None:
    """The closest thing to an end-to-end check without Vina: writer shape vs rules."""
    payload = _synthetic_payload()
    problems = validator.validate(payload)
    assert problems == [], problems


def test_the_writer_payload_check_is_not_vacuous() -> None:
    payload = _synthetic_payload()
    payload["targets"][0]["top1_rmsd"] = 4.0  # would be reported as a pass otherwise
    problems = validator.validate(payload)
    assert problems, "the validator must notice a broken payload"
    assert any("rank-1 mode's RMSD" in problem for problem in problems), problems


def test_the_writer_payload_round_trips_through_json() -> None:
    payload = _synthetic_payload()
    reloaded = json.loads(json.dumps(payload))
    assert validator.validate(reloaded) == []
    assert reloaded == payload


# --------------------------------------------------------------------------- #
# RDKit-gated: ligand preparation and the RMSD path
# --------------------------------------------------------------------------- #


def test_prepare_crystal_ligand_takes_bond_orders_from_the_definition() -> None:
    require_rdkit_chem()
    from rdkit import Chem

    heavy, prep, info = _nicotinic_ligand()

    assert info["heavy_atom_count"] == 9
    assert info["deposited_matches_definition_smiles"] is True
    # C6H5NO2: five hydrogens added on top of the nine heavy atoms.
    assert info["atoms_with_hydrogens"] == 14
    assert prep.GetNumAtoms() == 14
    assert heavy.GetNumAtoms() == 9
    # The assigned molecule is the definition, not a geometry guess.
    assert Chem.MolToSmiles(heavy, isomericSmiles=False) == Chem.MolToSmiles(
        Chem.MolFromSmiles(NICOTINIC_SMILES), isomericSmiles=False
    )


def test_prepare_crystal_ligand_refuses_a_definition_that_does_not_cover_the_ligand() -> None:
    require_rdkit_chem()
    instance = redock.ligand_instances(NICOTINIC_BLOCK, "LIG")[0]
    wrong = {"comp_id": "LIG", "smiles": "c1ccccc1", "atom_count_heavy": 6}

    with pytest.raises(redock.TargetFailure) as excinfo:
        redock.prepare_crystal_ligand("LIG", instance, wrong)

    message = str(excinfo.value)
    assert "LIG" in message
    assert "heavy" in message or "definition" in message


def test_pose_rmsd_is_zero_at_the_crystal_coordinates() -> None:
    require_rdkit_chem()
    heavy, prep, _ = _nicotinic_ligand()

    pose = redock.parse_pdbqt_atoms(pose_from_molecule(heavy))
    rmsd = redock.pose_rmsd_angstrom(prep, heavy, pose, where="test")

    assert rmsd < 1e-9


def test_pose_rmsd_does_not_re_align_the_pose() -> None:
    """A pure translation must come back as its own size.

    This is what separates ``rdMolAlign.CalcRMS`` (in place) from
    ``GetBestRMS``: had the pose been fitted onto the crystal pose first, a 2.5 A
    shift would have measured as 0.
    """
    require_rdkit_chem()
    heavy, prep, _ = _nicotinic_ligand()

    pose = redock.parse_pdbqt_atoms(pose_from_molecule(heavy, shift=(2.5, 0.0, 0.0)))
    rmsd = redock.pose_rmsd_angstrom(prep, heavy, pose, where="test")

    assert rmsd == pytest.approx(2.5, abs=1e-6)


def test_pose_rmsd_is_symmetry_aware() -> None:
    """A ring written with its atoms rotated by one position is the same pose.

    Every atom of this benzene is a carbon, so nothing but symmetry handling can
    tell the two apart: the atom-by-atom distance is the 1.39 A ring radius run
    as a chord, while the symmetry-aware answer is zero.
    """
    require_rdkit_chem()
    instance = redock.ligand_instances(BENZENE_BLOCK, "BEN")[0]
    definition = {"comp_id": "BEN", "smiles": BENZENE_SMILES, "atom_count_heavy": 6}
    heavy, prep, info = redock.prepare_crystal_ligand("BEN", instance, definition)
    assert info["heavy_atom_count"] == 6

    rotated = [index + 1 if index < 5 else 0 for index in range(6)]
    pose = redock.parse_pdbqt_atoms(pose_from_molecule(heavy, order=rotated))
    rmsd = redock.pose_rmsd_angstrom(prep, heavy, pose, where="test")

    naive = math.sqrt(
        sum(
            (position[axis] - BENZENE_ATOMS[(index + 1) % 6][axis]) ** 2
            for index, position in enumerate(BENZENE_ATOMS)
            for axis in range(3)
        )
        / 6
    )
    assert naive == pytest.approx(1.39, abs=0.01)
    assert rmsd < 1e-6


def test_pose_mapping_rejects_a_pose_whose_atoms_do_not_match_the_molecule() -> None:
    """A permuted pose must fail loudly instead of returning a wrong number.

    The oxygen and the nitrogen swap places in the file, so the pose's atom 3 is
    typed as nitrogen where the molecule's atom 3 is oxygen.
    """
    require_rdkit_chem()
    heavy, prep, _ = _nicotinic_ligand()

    order = list(range(9))
    order[2], order[7] = order[7], order[2]
    pose = redock.parse_pdbqt_atoms(pose_from_molecule(heavy, order=order))

    with pytest.raises(redock.TargetFailure) as excinfo:
        redock.pose_rmsd_angstrom(prep, heavy, pose, where="docked mode 1")

    message = str(excinfo.value)
    assert "docked mode 1" in message
    assert "atom correspondence is not trustworthy" in message


def test_pose_mapping_rejects_an_unknown_autodock_type() -> None:
    require_rdkit_chem()
    heavy, prep, _ = _nicotinic_ligand()

    pose = redock.parse_pdbqt_atoms(pose_from_molecule(heavy))
    pose[4]["type"] = "Q"
    with pytest.raises(redock.TargetFailure) as excinfo:
        redock.pose_rmsd_angstrom(prep, heavy, pose, where="test")

    assert "type-to-element table" in str(excinfo.value)


def test_pose_mapping_rejects_a_pose_with_the_wrong_atom_count() -> None:
    require_rdkit_chem()
    heavy, prep, _ = _nicotinic_ligand()

    pose = redock.parse_pdbqt_atoms(pose_from_molecule(heavy))[:-1]
    with pytest.raises(redock.TargetFailure) as excinfo:
        redock.pose_rmsd_angstrom(prep, heavy, pose, where="test")

    assert "heavy atoms" in str(excinfo.value)


def test_crystal_pose_self_check_is_zero_for_the_prepared_ligand() -> None:
    """The check the script runs before docking: the ligand PDBQT is the crystal pose."""
    require_rdkit_chem()
    heavy, prep, _ = _nicotinic_ligand()

    pose_text = pose_from_molecule(heavy)
    rmsd = redock.crystal_pose_self_check(prep, heavy, pose_text)

    assert rmsd < 1e-9


def test_crystal_pose_self_check_catches_a_broken_correspondence() -> None:
    require_rdkit_chem()
    heavy, prep, _ = _nicotinic_ligand()

    order = list(range(9))
    order[0], order[1] = order[1], order[0]
    with pytest.raises(redock.TargetFailure):
        redock.crystal_pose_self_check(prep, heavy, pose_from_molecule(heavy, order=order))
