from __future__ import annotations

import json
from types import SimpleNamespace

import pytest

from fastapi.testclient import TestClient

from app.main import app
from app.services.molecule_simulation import METHOD_CODES
from app.services.quantum_service import H2_HAMILTONIAN_ARTIFACT
from tests.availability import (
    require_live_network,
    require_rdkit_descriptors,
)

client = TestClient(app)


def _assert_error_shape(response, expected_code: str, expected_status: int = 422) -> dict:
    assert response.status_code == expected_status, response.text
    data = response.json()
    assert "error" in data
    assert data["error"]["code"] == expected_code
    assert "message" in data["error"]
    return data


def test_health() -> None:
    response = client.get("/health")
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "ok"
    assert data["app"] == "Quantum Drug Discovery Explorer API"


def _point_quantum_service_at(tmp_path, monkeypatch) -> None:
    """Point the quantum service at a directory we control.

    Both branches below are then deterministic, whether or not a real
    h2_sto3g.json happens to exist in the working tree.
    """
    monkeypatch.setattr(
        "app.services.quantum_service.get_settings",
        lambda: SimpleNamespace(h2_hamiltonian_dir=tmp_path),
    )


def test_h2_curve_reports_unavailable_when_artifact_is_missing(tmp_path, monkeypatch) -> None:
    """No fabricated matrix may be served, so the curve is honestly unavailable."""
    _point_quantum_service_at(tmp_path, monkeypatch)

    response = client.get("/api/quantum/h2-curve")

    data = _assert_error_shape(response, "hamiltonian_unavailable", 503)
    message = data["error"]["message"]
    assert H2_HAMILTONIAN_ARTIFACT in message
    assert "not been generated" in message
    assert "Colab" in message
    # This endpoint really does draw a curve, so the word belongs here -- which
    # is exactly why the shared sentence used to leak it into the simulate panel.
    assert "no H2 curve is available" in message
    assert "curve" in message.lower()


def test_h2_curve_is_placeholder_when_artifact_exists_but_is_not_wired(
    tmp_path, monkeypatch
) -> None:
    (tmp_path / H2_HAMILTONIAN_ARTIFACT).write_text("{}", encoding="utf-8")
    _point_quantum_service_at(tmp_path, monkeypatch)

    response = client.get("/api/quantum/h2-curve")

    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "placeholder"
    assert data["bond_lengths_angstrom"] == []
    assert data["exact_energies_ha"] == []
    assert data["note"]


def test_hadamard_demo() -> None:
    response = client.get("/api/quantum/hadamard-demo")
    assert response.status_code == 200
    data = response.json()
    assert data["circuit"] == "H|0>"
    assert data["probabilities"]["0"] == pytest.approx(0.5, abs=1e-9)
    assert data["probabilities"]["1"] == pytest.approx(0.5, abs=1e-9)


def test_simulate_missing_body() -> None:
    response = client.post("/api/quantum/simulate", content="not json")
    _assert_error_shape(response, "invalid_body", 400)


def test_simulate_invalid_field() -> None:
    response = client.post(
        "/api/quantum/simulate",
        content=json.dumps({"bond_length": "bad"}),
    )
    _assert_error_shape(response, "invalid_bond_length", 422)


def test_simulate_h2_reports_unavailable_without_the_artifact(tmp_path, monkeypatch) -> None:
    """Without the offline artifact there is no honest energy to report."""
    _point_quantum_service_at(tmp_path, monkeypatch)

    response = client.post(
        "/api/quantum/simulate",
        content=json.dumps({"bond_length": 0.74, "shots": 1024, "noise": 0.0}),
    )

    data = _assert_error_shape(response, "hamiltonian_unavailable", 503)
    message = data["error"]["message"]
    assert H2_HAMILTONIAN_ARTIFACT in message
    # The simulate panel reports one energy point. It must name that energy, not
    # the curve that a different endpoint draws.
    assert "no H2 energy is available" in message
    assert "curve" not in message.lower(), message


def test_simulate_h2_is_placeholder_when_artifact_exists_but_is_not_wired(
    tmp_path, monkeypatch
) -> None:
    (tmp_path / H2_HAMILTONIAN_ARTIFACT).write_text("{}", encoding="utf-8")
    _point_quantum_service_at(tmp_path, monkeypatch)

    response = client.post(
        "/api/quantum/simulate",
        content=json.dumps({"bond_length": 0.74, "shots": 1024, "noise": 0.0}),
    )

    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "placeholder"
    assert data["exact_energy_ha"] is None
    assert data["vqe_energy_ha"] is None
    assert data["shots"] == 1024
    assert data["noise"] == 0.0
    assert data["convergence_history_ha"] is None
    # Same rule on the other branch: a present-but-unwired artifact still must
    # not describe this endpoint in terms of the curve.
    assert data["note"]
    assert "curve" not in data["note"].lower(), data["note"]


def test_simulate_unsupported_molecule_is_placeholder() -> None:
    response = client.post(
        "/api/quantum/simulate",
        content=json.dumps({"bond_length": 0.74, "shots": 1, "noise": 0.0, "molecule": "caffeine"}),
    )
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "placeholder"
    assert data["vqe_energy_ha"] is None


def test_resources_h2_is_placeholder_without_the_artifact(tmp_path, monkeypatch) -> None:
    """These numbers were previously derived from the fabricated matrix."""
    _point_quantum_service_at(tmp_path, monkeypatch)

    response = client.post(
        "/api/quantum/resources",
        content=json.dumps({"molecule": "H2"}),
    )

    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "placeholder"
    assert data["n_qubits"] is None
    assert data["hamiltonian_terms"] is None
    assert data["note"]
    # It is a resource estimate, not a curve, so the panel must not be told to
    # look at one.
    assert "no H2 resource estimate is available" in data["note"]
    assert "curve" not in data["note"].lower(), data["note"]


def test_resources_unsupported_is_placeholder() -> None:
    response = client.post(
        "/api/quantum/resources",
        content=json.dumps({"molecule": "caffeine"}),
    )
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "placeholder"
    assert "n_qubits" not in data


def test_molecule_search_live_pubchem() -> None:
    require_live_network()


def test_molecule_search_unknown_name_returns_compound_not_found() -> None:
    """A name PubChem does not have is a 404 answer, not a retry and not a 200.

    Replaces ``test_molecule_search_returns_placeholder_when_pubchem_fails``:
    A1 separates "not found" from "not working", so the honest response is an
    error envelope with one clean, URL-free message.
    """
    query = "xyzzynonexistentsearchterm"
    response = client.get(f"/api/molecule/search?q={query}")

    assert response.status_code == 404
    data = response.json()
    assert set(data) == {"error"}, f"expected a single error envelope, got {set(data)}"
    assert data["error"]["code"] == "compound_not_found"
    assert data["error"]["message"] == (
        f"No compound found for '{query}'. Check the spelling."
    )
    message = data["error"]["message"].lower()
    for leak in ("http://", "https://", "pubchem", "traceback"):
        assert leak not in message


def test_molecule_properties_live_rdkit() -> None:
    """The endpoint must return real RDKit descriptors, not a placeholder."""
    require_rdkit_descriptors()
    response = client.post(
        "/api/molecule/properties",
        content=json.dumps({"smiles": "CC(=O)Oc1ccccc1C(=O)O"}),
    )
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "computed", data.get("note")
    assert data["formula"] == "C9H8O4"
    assert data["h_bond_donors"] == 1
    assert data["h_bond_acceptors"] == 4
    assert data["rotatable_bonds"] == 3
    assert data["lipinski_rule_of_five"]["violations"] == 0


def test_molecule_properties_returns_placeholder_when_rdkit_fails(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def _rdkit_unavailable(*_args: object, **_kwargs: object) -> None:
        raise RuntimeError("RDKit not available in this build")

    monkeypatch.setattr(
        "app.services.rdkit_helper.rdkit_properties_from_smiles", _rdkit_unavailable
    )
    response = client.post(
        "/api/molecule/properties",
        content=json.dumps({"smiles": "CC(=O)Oc1ccccc1C(=O)O"}),
    )
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "placeholder"
    assert data["smiles"] is not None
    assert data["note"] is not None
    assert "RDKit" in data["note"]


def test_molecule_properties_invalid_smiles_is_placeholder() -> None:
    response = client.post(
        "/api/molecule/properties",
        content=json.dumps({"smiles": "not a valid smiles"}),
    )
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "placeholder"
    assert data["note"]


def _pdb_atom(
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


# A physically plausible fragment: the backbone amide N of ASP A25 donates a
# hydrogen bond to a ligand oxygen placed in the direction the N-H bond points
# (the N -> CA direction is +x, so the acceptor sits on -x), and the ASP
# carboxylate oxygen is 3.11 A from a ligand nitrogen for the salt bridge. The
# geometry is what makes the classification unambiguous, and it is independent
# of any network round-trip or the machine's RDKit build.
PDB_FIXTURE = "\n".join([
    _pdb_atom("ATOM", 1, "N", "ASP", "A", 25, 0.000, 0.000, 0.000, "N"),
    _pdb_atom("ATOM", 2, "CA", "ASP", "A", 25, 1.460, 0.000, 0.000, "C"),
    _pdb_atom("ATOM", 3, "OD1", "ASP", "A", 25, 1.200, 1.400, 0.000, "O"),
    _pdb_atom("HETATM", 9, "O1", "MK1", "B", 400, -2.800, 0.200, 0.000, "O"),
    _pdb_atom("HETATM", 10, "N2", "MK1", "B", 400, 1.500, 4.500, 0.000, "N"),
    _pdb_atom("HETATM", 11, "C1", "MK1", "B", 400, 6.000, 6.000, 6.000, "C"),
    "END",
])


def test_complexes_list_reports_real_ligand_codes() -> None:
    response = client.get("/api/complexes/list")
    assert response.status_code == 200
    data = response.json()
    assert "complexes" in data
    assert data["count"] > 0
    codes = {item["id"]: item["ligand_code"] for item in data["complexes"]}
    # The manifest must carry the HETATM codes actually deposited in each entry.
    assert codes["hiv-protease-indinavir"] == "MK1"
    assert codes["ace-lisinopril"] == "LPR"
    # Codes verified against the deposited files for the viral additions.
    assert codes["hiv-rt-nevirapine"] == "NVP"
    assert codes["influenza-neuraminidase-oseltamivir"] == "G39"


def test_far_donor_acceptor_pair_is_a_near_miss_not_a_hydrogen_bond() -> None:
    """Ideal donor/acceptor geometry beyond 3.5 A must not be called a hydrogen bond.

    This is the nevirapine situation in 1VRT: N14...O Lys101 is 4.05 A at 164 deg.
    """
    from app.services.structure_analysis import HBOND_MAX, analyse_complex

    pdb = "\n".join(
        [
            _pdb_atom("ATOM", 1, "N", "ASP", "A", 25, 0.000, 0.000, 0.000, "N"),
            _pdb_atom("ATOM", 2, "CA", "ASP", "A", 25, 1.460, 0.000, 0.000, "C"),
            _pdb_atom("HETATM", 9, "O1", "MK1", "B", 400, -4.050, 0.100, 0.000, "O"),
            "END",
        ]
    )
    analysis = analyse_complex(pdb, "MK1")

    assert analysis["status"] == "computed"
    assert analysis["contact_summary"]["hydrogen_bond"] == 0

    misses = analysis["near_miss_hydrogen_bonds"]
    assert len(misses) == 1
    miss = misses[0]
    assert miss["protein_residue"] == "ASP A25"
    assert miss["ligand_atom"] == "O1"
    assert miss["distance_angstrom"] > HBOND_MAX
    assert miss["distance_angstrom"] == pytest.approx(4.05, abs=0.05)
    assert miss["angle_deg"] is not None and miss["angle_deg"] > 150.0
    assert "cutoff" in miss["reason"]


def test_contact_labels_mark_hydrogen_bonds_as_a_heuristic() -> None:
    from app.services.structure_analysis import CONTACT_LABELS

    assert CONTACT_LABELS["hydrogen_bond"] == "geometric heuristic"
    assert CONTACT_LABELS["salt_bridge"] == "geometric heuristic"
    assert CONTACT_LABELS["hydrophobic"] == "measured geometry"


def test_complexes_label_human_versus_viral_targets() -> None:
    """ACE must be presented as a human enzyme, never as a virus."""
    response = client.get("/api/complexes/list")
    data = response.json()
    by_id = {item["id"]: item for item in data["complexes"]}

    assert by_id["ace-lisinopril"]["target_origin"] == "human"
    assert "not a virus" in by_id["ace-lisinopril"]["display_name"].lower()

    for viral_id in (
        "hiv-protease-indinavir",
        "hiv-rt-nevirapine",
        "influenza-neuraminidase-oseltamivir",
    ):
        assert by_id[viral_id]["target_origin"] == "viral"

    assert by_id["hiv-rt-nevirapine"]["binding_site_type"] == "allosteric"
    assert by_id["hiv-protease-indinavir"]["binding_site_type"] == "catalytic"


def test_allosteric_mechanism_does_not_claim_active_site_blockade() -> None:
    """An allosteric inhibitor is expected to miss the catalytic residues."""
    from app.services.structure_analysis import build_mechanism

    analysis = {
        "status": "computed",
        "ligand_code": "NVP",
        "ligand_atom_count": 20,
        "binding_site_residue_count": 1,
        "contact_summary": {
            "hydrogen_bond": 0,
            "salt_bridge": 0,
            "hydrophobic": 5,
            "close_contact": 9,
        },
        "interacting_residues": [
            {
                "residue": "TYR",
                "chain": "A",
                "residue_number": "181",
                "residue_label": "TYR A181",
                "types": ["hydrophobic"],
                "contact_count": 3,
                "min_distance_angstrom": 3.04,
                "atom_pairs": [],
            }
        ],
        "binding_site_residues": ["TYR A181"],
    }
    meta = {
        "binding_site_type": "allosteric",
        "target_origin": "viral",
        "binding_pocket_residues": [{"label": "Tyr181", "match": "TYR A181"}],
        "catalytic_residues": [
            {"label": "Asp110 (catalytic)", "match": "ASP A110", "role": "catalytic"}
        ],
        "target_class": "HIV-1 reverse transcriptase",
        "drug_class": "NNRTI",
        "mechanism_summary": "summary",
        "consequence": "consequence",
    }

    mechanism = build_mechanism(analysis, meta)
    blockade = next(s for s in mechanism["steps"] if s["step"] == "blockade")

    assert blockade["value"] == "allosteric"
    assert "allosteric" in blockade["title"].lower()
    assert "not the active site" in blockade["title"].lower()
    assert "Tyr181" in blockade["detail"]
    # The catalytic residue is curated but must be reported as NOT in contact.
    assert mechanism["catalytic_residues_occluded"] == 0
    assert mechanism["catalytic_residues_curated"] == 1
    assert mechanism["binding_site_type"] == "allosteric"
    assert mechanism["catalytic_residues"][0]["in_contact"] is False


def test_catalytic_mechanism_still_reports_occlusion() -> None:
    """The catalytic branch must be unchanged by the allosteric branch."""
    from app.services.structure_analysis import build_mechanism

    analysis = {
        "status": "computed",
        "ligand_code": "MK1",
        "ligand_atom_count": 45,
        "binding_site_residue_count": 1,
        "contact_summary": {
            "hydrogen_bond": 1,
            "salt_bridge": 0,
            "hydrophobic": 3,
            "close_contact": 4,
        },
        "interacting_residues": [
            {
                "residue": "ASP",
                "chain": "A",
                "residue_number": "25",
                "residue_label": "ASP A25",
                "types": ["hydrogen_bond"],
                "contact_count": 1,
                "min_distance_angstrom": 2.77,
                "atom_pairs": [],
            }
        ],
        "binding_site_residues": ["ASP A25"],
    }
    meta = {
        "binding_site_type": "catalytic",
        "catalytic_residues": [
            {"label": "Asp25", "match": "ASP A25", "role": "catalytic aspartate"}
        ],
        "mechanism_summary": "summary",
        "consequence": "consequence",
    }

    mechanism = build_mechanism(analysis, meta)
    blockade = next(s for s in mechanism["steps"] if s["step"] == "blockade")
    assert blockade["value"] == "occluded"
    assert mechanism["catalytic_residues_occluded"] == 1


def test_complex_interactions_computed_from_coordinates(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        "app.services.complexes_service.http_get_bytes_cached",
        lambda *_args, **_kwargs: PDB_FIXTURE.encode("utf-8"),
    )
    response = client.get("/api/complexes/hiv-protease-indinavir/interactions")
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "computed", data.get("note")
    assert data["ligand_code"] == "MK1"
    assert data["contact_summary"]["hydrogen_bond"] >= 1
    assert data["contact_summary"]["salt_bridge"] >= 1
    assert "ASP A25" in data["binding_site_residues"]


def test_complex_interactions_placeholder_when_pdb_unavailable(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def _unreachable(*_args: object, **_kwargs: object) -> None:
        raise RuntimeError("RCSB unreachable")

    monkeypatch.setattr(
        "app.services.complexes_service.http_get_bytes_cached", _unreachable
    )
    response = client.get("/api/complexes/hiv-protease-indinavir/interactions")
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "placeholder"
    assert data["complex_id"] == "hiv-protease-indinavir"
    assert data["note"] is not None


def test_complex_structure_exposes_mechanism_and_coordinates(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        "app.services.complexes_service.http_get_bytes_cached",
        lambda *_args, **_kwargs: PDB_FIXTURE.encode("utf-8"),
    )
    response = client.get("/api/complexes/hiv-protease-indinavir/structure")
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "computed", data.get("note")
    # The viewer payload is the trimmed protein + ligand, not the raw deposit.
    assert "MK1" in data["structure_file_bytes"]
    assert data["structure_file_reference"] == "1HSG.pdb"
    mechanism = data["mechanism"]
    assert mechanism["status"] == "computed"
    steps = {step["step"] for step in mechanism["steps"]}
    assert {"target", "occupy", "anchor", "blockade", "consequence"} <= steps
    asp25 = next(c for c in mechanism["catalytic_residues"] if c["match"] == "ASP A25")
    assert asp25["in_contact"] is True
    assert asp25["min_distance_angstrom"] is not None


def test_complex_structure_unknown_id_is_placeholder() -> None:
    response = client.get("/api/complexes/not-a-real-complex/structure")
    assert response.status_code == 200
    assert response.json()["status"] == "placeholder"


def test_complex_docking_validation_placeholder() -> None:
    response = client.get("/api/complexes/hiv-protease-indinavir/docking-validation")
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "placeholder"
    assert data["complex_id"] == "hiv-protease-indinavir"
    assert data["rmsd_angstrom"] is None


def test_invalid_pdb_id() -> None:
    response = client.get("/api/protein/ABC")
    _assert_error_shape(response, "invalid_pdb_id", 400)


def test_protein_live_rcsb() -> None:
    require_live_network()


def test_protein_returns_placeholder_when_rcsb_fails(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def _rcsb_unreachable(*_args: object, **_kwargs: object) -> None:
        raise RuntimeError("RCSB unreachable")

    monkeypatch.setattr("app.services.rscb_service.pdb_entry_json", _rcsb_unreachable)
    response = client.get("/api/protein/6LU7")
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "placeholder"
    assert data["pdb_id"] == "6LU7"
    assert data["note"] is not None
    assert "RCSB" in data["note"]


def test_ai_resolve_no_match() -> None:
    response = client.post(
        "/api/ai/resolve",
        content=json.dumps({"query": "blahblah"}),
    )
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "no_match"
    assert data["resolved"] is None


def test_ai_resolve_scaffold_example() -> None:
    response = client.post(
        "/api/ai/resolve",
        content=json.dumps({"query": "covid protease"}),
    )
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "resolved"


def test_ai_predict_properties_not_trained() -> None:
    response = client.post(
        "/api/ai/predict-properties",
        content=json.dumps({"smiles": "CC(=O)Oc1ccccc1C(=O)O"}),
    )
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "not_trained"
    assert data["model"] == "not trained"
    assert data["method"] == "heuristic"
    assert data["predictions"] is None


def test_molecule_simulate_requires_smiles() -> None:
    response = client.post("/api/molecule/simulate", content=json.dumps({}))
    _assert_error_shape(response, "invalid_smiles", 422)


def test_molecule_simulate_unparseable_smiles_is_placeholder() -> None:
    response = client.post(
        "/api/molecule/simulate",
        content=json.dumps({"smiles": "not a valid smiles"}),
    )
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "placeholder"
    assert data["note"] is not None


def test_molecule_simulate_shape_with_mocked_chemistry(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Pin the simulation payload contract without depending on this machine's
    RDKit build, which varies with the OS application-control policy."""
    properties = {
        "molecular_weight": 180.16,
        "logp": 1.31,
        "h_bond_donors": 1,
        "h_bond_acceptors": 4,
        "rotatable_bonds": 3,
        "formula": "C9H8O4",
        "exact_mass": 180.04,
        "tpsa": 63.6,
        "heavy_atom_count": 13,
        "ring_count": 1,
        "aromatic_ring_count": 1,
        "fraction_csp3": 0.111,
        "veber": {"pass": True, "rotatable_bonds": 3, "tpsa": 63.6},
        "lipinski_rule_of_five": {"pass": True, "violations": 0},
    }
    conformers = [
        {
            "index": 0,
            "relative_energy_kcal_per_mol": 0.0,
            "minimised": True,
            "population": 0.7,
            "molblock": None,
        },
        {
            "index": 1,
            "relative_energy_kcal_per_mol": 1.0,
            "minimised": True,
            "population": 0.3,
            "molblock": None,
        },
    ]
    monkeypatch.setattr(
        "app.services.rdkit_helper.normalize_smiles_safe", lambda s: s
    )
    monkeypatch.setattr(
        "app.services.rdkit_helper.rdkit_properties_from_smiles", lambda s: properties
    )
    monkeypatch.setattr(
        "app.services.rdkit_helper.smiles_to_conformer_ensemble",
        lambda s, num_conformers=8, max_iters=400, temperature_k=298.15: {
            "available": True,
            "method": "mocked ensemble",
            "force_field": "MMFF94",
            "temperature_k": 298.15,
            "conformer_count": 2,
            "energy_span_kcal_per_mol": 1.0,
            "lowest_absolute_energy_kcal_per_mol": 10.0,
            "energy_reference": "relative_to_lowest_conformer_in_this_run",
            "absolute_energies_comparable": False,
            "absolute_energy_note": "absolute energies are not comparable between runs",
            "populations_are_estimates": True,
            "population_note": "populations are Boltzmann estimates",
            "conformers": conformers,
            "lowest_energy_index": 0,
            "lowest_energy_relative_kcal_per_mol": 0.0,
            "lowest_energy_molblock": "MOLBLOCK",
            "note": None,
        },
    )
    monkeypatch.setattr(
        "app.services.rdkit_helper.smiles_to_3d_structure",
        lambda s, num_conformers=8, max_iters=400: {
            "available": True,
            "is_3d": True,
            "molblock": "MOLBLOCK",
            "method": "mocked 3D",
            "force_field": "MMFF94",
            "relative_energy_kcal_per_mol": 0.0,
            "atom_count": 21,
            "conformer_count": 2,
            "note": None,
        },
    )
    monkeypatch.setattr(
        "app.services.rdkit_helper.smiles_to_2d_structure",
        lambda s: {"available": True, "molblock": "TWOD", "method": "mocked 2D", "atom_count": 13, "note": None},
    )

    response = client.post(
        "/api/molecule/simulate",
        content=json.dumps({"smiles": "CC(=O)Oc1ccccc1C(=O)O", "num_conformers": 2}),
    )
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "computed"
    assert data["structure"]["is_3d"] is True
    assert data["structure"]["molblock"] == "MOLBLOCK"
    assert data["ensemble"]["conformer_count"] == 2
    assert data["depiction"]["molblock"] == "TWOD"
    # Populations must be reported, and the sourced rules must be evaluated.
    assert [c["population"] for c in data["ensemble"]["conformers"]] == [0.7, 0.3]
    # A5: geometry method is a machine-readable code, energies are relative, and
    # populations are labelled as estimates.
    assert data["method"] in METHOD_CODES
    assert data["method_detail"]
    assert data["energy_basis"] == "relative_to_lowest_conformer_in_this_run"
    assert data["ensemble"]["absolute_energies_comparable"] is False
    assert data["ensemble"]["populations_are_estimates"] is True
    assert data["ensemble"]["population_note"]
    assert [c["relative_energy_kcal_per_mol"] for c in data["ensemble"]["conformers"]] == [0.0, 1.0]
    assert data["findings"], "drug-likeness findings should be present"
    assert all(f["source"] for f in data["findings"])
    steps = {step["step"] for step in data["simulation_steps"]}
    assert {"parse", "descriptors", "embed", "minimise", "population", "interpret"} <= steps
