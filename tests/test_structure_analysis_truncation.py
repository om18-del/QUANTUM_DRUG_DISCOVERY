"""Tests for the per-residue atom_pairs cap (atom_pairs_truncated)."""

from __future__ import annotations

import math

from fastapi.testclient import TestClient

from app.main import app
from app.services.structure_analysis import ATOM_PAIRS_PER_RESIDUE_CAP, analyse_complex

client = TestClient(app)


def _atom(
    record: str,
    serial: int,
    name: str,
    res: str,
    chain: str,
    resseq: int,
    x: float,
    y: float,
    z: float,
    element: str,
) -> str:
    """Format one ATOM/HETATM line with the columns the PDB spec requires."""
    return (
        f"{record:<6}{serial:>5} {name:<4} {res:>3} {chain}{resseq:>4}    "
        f"{x:>8.3f}{y:>8.3f}{z:>8.3f}{1.0:>6.2f}{20.0:>6.2f}{'':>10}{element:>2}"
    )


def _crowded_complex(n_ligand_atoms: int) -> str:
    """One protein residue whose single carbon contacts N ligand carbons.

    Every ligand carbon sits 3.0 A from the protein carbon, so each is a
    hydrophobic contact and the residue accumulates one atom pair per ligand atom.
    """
    lines = [_atom("ATOM", 1, "CB", "ASP", "A", 25, 0.0, 0.0, 0.0, "C")]
    for i in range(n_ligand_atoms):
        angle = 2.0 * math.pi * i / n_ligand_atoms
        lines.append(
            _atom(
                "HETATM",
                100 + i,
                f"C{i + 1}",
                "MK1",
                "B",
                400,
                3.0 * math.cos(angle),
                3.0 * math.sin(angle),
                0.0,
                "C",
            )
        )
    lines.append("END")
    return "\n".join(lines)


def test_atom_pairs_are_truncated_and_flagged() -> None:
    n = ATOM_PAIRS_PER_RESIDUE_CAP + 1
    analysis = analyse_complex(_crowded_complex(n), "MK1")

    assert analysis["status"] == "computed"
    residue = analysis["interacting_residues"][0]
    assert residue["residue_label"] == "ASP A25"
    assert len(residue["atom_pairs"]) == ATOM_PAIRS_PER_RESIDUE_CAP
    assert residue["atom_pairs_truncated"] is True

    # The tally stays complete even though the list is capped.
    assert analysis["atom_pair_counts"]["hydrophobic"] == n
    assert analysis["atom_pairs_truncated"] is True
    assert analysis["note"] and "capped" in analysis["note"]


def test_atom_pairs_below_the_cap_are_not_flagged() -> None:
    n = ATOM_PAIRS_PER_RESIDUE_CAP - 3
    analysis = analyse_complex(_crowded_complex(n), "MK1")

    residue = analysis["interacting_residues"][0]
    assert len(residue["atom_pairs"]) == n
    assert residue["atom_pairs_truncated"] is False
    assert analysis["atom_pairs_truncated"] is False
    assert analysis["note"] is None
    # The cap is documented even when it did not bite.
    assert analysis["atom_pairs_truncated_note"]


def test_structure_endpoint_exposes_truncation(monkeypatch) -> None:
    monkeypatch.setattr(
        "app.services.complexes_service.http_get_bytes_cached",
        lambda *_args, **_kwargs: _crowded_complex(ATOM_PAIRS_PER_RESIDUE_CAP + 1).encode("utf-8"),
    )

    response = client.get("/api/complexes/hiv-protease-indinavir/structure")

    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "computed"
    assert data["atom_pairs_truncated"] is True
    assert data["note"] and "capped" in data["note"]
