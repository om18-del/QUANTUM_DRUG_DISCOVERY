"""Precompute real H2 / LiH STO-3G qubit Hamiltonians, offline.

This script must be run where PySCF works (Linux/WSL/Colab). PySCF publishes no
Windows wheel and ``qiskit-nature`` declares PySCF with ``sys_platform != "win32"``,
so there is deliberately no Windows path here.

Every coefficient in the output file is produced by the computation below. Nothing
is typed in by hand. The script refuses to write a file unless, for every bond
length, the smallest eigenvalue of the stored qubit operator plus the stored
nuclear repulsion reproduces the PySCF FCI energy to 1e-6 Ha.

Run:
    python scripts/precompute_hamiltonians.py                    # H2  -> h2_sto3g.json
    python scripts/precompute_hamiltonians.py --molecule lih     # LiH -> lih_sto3g_active.json
    python scripts/precompute_hamiltonians.py --molecule lih out.json
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import sys
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any

import numpy as np

# Pins. These are stated in one place so a mismatch is loud rather than silent.
#   qiskit-nature 0.8.0 requires qiskit>=1.4,<3.0  -> satisfied by qiskit 1.4.6
#   qiskit-nature 0.8.0 requires numpy>=2          -> Colab ships numpy 2.x
PINNED_VERSIONS = {
    "qiskit": "1.4.6",
    "qiskit-nature": "0.8.0",
    "pyscf": "2.14.0",
}

BASIS = "sto-3g"
MAPPING = "parity"
TWO_QUBIT_REDUCTION = True

# H2 in STO-3G: 4 spin orbitals, 2 electrons -> parity mapping with the
# two-qubit reduction leaves 2 qubits. This is known for H2 and is asserted below.
# The LiH active space is smaller than the full basis, so its qubit count is
# recorded from the mapped operator, never assumed.
H2_EXPECTED_QUBITS = 2

# The H2 grid: every point of 0.20 .. 3.00 angstrom inclusive in 0.05 steps (57
# uniform points), plus 0.74 angstrom -- the UI default and H2's real equilibrium
# bond length, which the 0.05 step misses (0.74 would need index 10.8). 58 points
# in total, sorted ascending with no duplicates. Built from integers so
# floating-point drift cannot silently change the grid, and the set union makes a
# duplicate impossible even if 0.74 ever lands exactly on the step.
H2_UNIFORM_BOND_LENGTHS: list[float] = [round(0.20 + 0.05 * i, 2) for i in range(57)]
H2_BOND_LENGTHS: list[float] = sorted({*H2_UNIFORM_BOND_LENGTHS, 0.74})
assert len(H2_BOND_LENGTHS) == 58, f"expected 58 bond lengths, got {len(H2_BOND_LENGTHS)}"

# The LiH grid: 1.0 .. 3.0 angstrom in 0.1 steps (21 uniform points), plus 1.595
# angstrom near the equilibrium separation, which the 0.1 step misses.
LIH_UNIFORM_BOND_LENGTHS: list[float] = [round(1.0 + 0.1 * i, 1) for i in range(21)]
LIH_BOND_LENGTHS: list[float] = sorted({*LIH_UNIFORM_BOND_LENGTHS, 1.595})

TOLERANCE_HA = 1e-6

# One entry per supported molecule. The grid is declared here and written into the
# file so the validator can read it back instead of assuming a point count.
MOLECULES: dict[str, dict[str, Any]] = {
    "h2": {
        "label": "H2",
        "file": "h2_sto3g.json",
        "template": "H 0.0 0.0 0.0; H 0.0 0.0 {r}",
        "bond_lengths": H2_BOND_LENGTHS,
        "freeze_core": False,
        "remove_highest_virtuals": 0,
        "expected_qubits": H2_EXPECTED_QUBITS,
    },
    "lih": {
        "label": "LiH",
        "file": "lih_sto3g_active.json",
        "template": "Li 0.0 0.0 0.0; H 0.0 0.0 {r}",
        "bond_lengths": LIH_BOND_LENGTHS,
        # Freeze the Li 1s core and drop the two highest virtual spatial orbitals.
        # What that leaves is described from the transformer's own numbers and
        # stored in the file as `active_space`.
        "freeze_core": True,
        "remove_highest_virtuals": 2,
        "expected_qubits": None,
    },
}

DEFAULT_MOLECULE = "h2"
DEFAULT_OUTPUT_DIR = Path(__file__).resolve().parents[1] / "app" / "data" / "hamiltonians"


def resolve_output_path(molecule_key: str, explicit: str | None) -> Path:
    """Where to write the JSON.

    Precedence: an explicit argument, then ``QDDE_HAMILTONIAN_OUT``, then the repo
    location. The override exists because in Colab the script does not sit inside
    the repository, so the repo-relative default would be wrong.
    """
    if explicit:
        return Path(explicit).resolve()
    override = os.environ.get("QDDE_HAMILTONIAN_OUT")
    if override:
        return Path(override).resolve()
    return DEFAULT_OUTPUT_DIR / MOLECULES[molecule_key]["file"]


def _require_dependencies() -> dict[str, str]:
    """Import the stack and report the versions actually present.

    A missing dependency or a version mismatch is raised here rather than worked
    around, so the output file can always be attributed to exact versions.
    """
    try:
        from pyscf import fci, gto, scf  # noqa: F401
        from qiskit.quantum_info import SparsePauliOp  # noqa: F401
        from qiskit_nature.second_q.drivers import PySCFDriver  # noqa: F401
        from qiskit_nature.second_q.mappers import ParityMapper  # noqa: F401
    except Exception as exc:  # noqa: BLE001 - report exactly what failed
        raise SystemExit(
            "Cannot import the quantum-chemistry stack.\n"
            f"  platform : {platform.platform()}\n"
            f"  python   : {sys.version.split()[0]}\n"
            f"  error    : {type(exc).__name__}: {exc}\n"
            "PySCF has no Windows wheel; run this in WSL or Colab."
        ) from exc

    found: dict[str, str] = {}
    for name in PINNED_VERSIONS:
        try:
            found[name] = version(name)
        except PackageNotFoundError:
            found[name] = "MISSING"

    problems = [
        f"{name}: expected {expected}, found {found[name]}"
        for name, expected in PINNED_VERSIONS.items()
        if found[name] != expected
    ]
    if problems:
        raise SystemExit(
            "Dependency versions do not match the pins.\n  "
            + "\n  ".join(problems)
            + "\nInstall the pinned set (see the Colab notebook) before running."
        )
    return found


def _pauli_terms(operator: Any) -> list[dict[str, Any]]:
    """Convert a SparsePauliOp into JSON-friendly Pauli terms."""
    labels = operator.paulis.to_labels()
    terms: list[dict[str, Any]] = []
    for label, coeff in zip(labels, operator.coeffs):
        value = complex(coeff)
        if abs(value.imag) > 1e-10:
            raise ValueError(f"non-real coefficient for {label}: {value}")
        terms.append({"pauli": str(label), "coefficient": float(value.real)})
    # Largest contribution first; every term is kept.
    terms.sort(key=lambda t: (-abs(t["coefficient"]), t["pauli"]))
    return terms


def _ground_energy(operator: Any) -> float:
    """Smallest eigenvalue of the (Hermitian) qubit operator."""
    matrix = np.asarray(operator.to_matrix())
    return float(np.linalg.eigvalsh(matrix)[0].real)


def _fermionic_operator(problem: Any) -> Any:
    """Fetch the electronic fermionic operator across qiskit-nature API shapes."""
    ops = problem.second_q_ops()
    if isinstance(ops, dict):
        # 0.8 returns a mapping keyed by operator name.
        if "ElectronicEnergy" in ops:
            return ops["ElectronicEnergy"]
        return next(iter(ops.values()))
    return ops[0]


def _pyscf_reference(bond_length: float, atom: str) -> tuple[Any, float, float]:
    """PySCF's own RHF and FCI reference energies for one geometry.

    PySCF is used directly for HF and FCI. Its ``gto.M`` accepts the plain
    ``"Angstrom"`` string; the qiskit-nature driver below needs the
    ``DistanceUnit`` enum instead, because 0.8.0 stores the unit as a
    ``DistanceUnit`` and reads ``self._unit.value`` (raising an AttributeError
    when it is given a bare string).
    """
    from pyscf import fci as pyscf_fci
    from pyscf import gto, scf

    mol = gto.M(atom=atom, basis=BASIS, unit="Angstrom", verbose=0)
    mf = scf.RHF(mol)
    mf.verbose = 0
    mf.kernel()
    hf_energy = float(mf.e_tot)

    cisolver = pyscf_fci.FCI(mol, mf.mo_coeff)
    cisolver.verbose = 0
    fci_energy, _ = cisolver.kernel()
    return mol, hf_energy, float(fci_energy)


def build_point(bond_length: float, molecule_key: str) -> tuple[dict[str, Any], str]:
    """Compute one bond length: HF, FCI, the mapped Hamiltonian and the active space.

    Returns the stored point and the free-text active-space description (the same
    for every bond length of a molecule, so the caller records it once at the top
    level). Nothing about the active space is assumed: the description is built
    from the frozen/removed orbital indices and from what the transformed problem
    reports.
    """
    config = MOLECULES[molecule_key]
    atom = config["template"].format(r=bond_length)

    from qiskit.quantum_info import SparsePauliOp
    from qiskit_nature.second_q.drivers import PySCFDriver
    from qiskit_nature.second_q.mappers import ParityMapper
    from qiskit_nature.units import DistanceUnit

    mol, hf_energy, fci_energy = _pyscf_reference(bond_length, atom)

    problem = PySCFDriver(atom=atom, basis=BASIS, unit=DistanceUnit.ANGSTROM).run()
    nuclear_repulsion = float(problem.nuclear_repulsion_energy)

    if abs(nuclear_repulsion - float(mol.energy_nuc())) > 1e-6:
        raise ValueError(
            f"nuclear repulsion disagrees between PySCF and Qiskit Nature at "
            f"{bond_length} A: {nuclear_repulsion} vs {mol.energy_nuc()}"
        )

    n_remove = int(config["remove_highest_virtuals"])
    if config["freeze_core"] or n_remove:
        from qiskit_nature.second_q.transformers import FreezeCoreTransformer

        # ``remove_orbitals`` uses zero-based spatial-orbital indices of the FULL
        # problem. Counting back from the AO count removes the highest virtuals
        # without hard-coding indices that depend on the basis.
        n_ao = int(mol.nao)
        remove_orbitals = [n_ao - n_remove + i for i in range(n_remove)] if n_remove else None
        problem = FreezeCoreTransformer(
            freeze_core=bool(config["freeze_core"]),
            remove_orbitals=remove_orbitals,
        ).transform(problem)

        # "What was frozen" is derived, not described by hand: the transformer
        # reports the active space it produced, so the number of core spatial
        # orbitals it dropped is the full-basis count minus the orbitals we asked
        # for minus what is left. If that bookkeeping does not add up, say so
        # instead of writing a description that may be wrong.
        active_spatial = int(problem.num_spatial_orbitals)
        frozen_count = (n_ao - n_remove) - active_spatial
        if frozen_count < 0:
            raise ValueError(
                f"the transformed problem reports {active_spatial} active spatial orbitals, but "
                f"{n_ao} basis orbitals minus {n_remove} removed leaves {n_ao - n_remove}; the "
                "active-space bookkeeping does not add up, so the active space cannot be "
                "described or stored."
            )
        frozen_text = (
            f"froze {frozen_count} core spatial orbital(s) "
            f"(zero-based indices {list(range(frozen_count))})"
            if frozen_count
            else "froze no core orbital"
        )
        removed_text = (
            f"removed the {n_remove} highest virtual spatial orbital(s) "
            f"(zero-based indices {remove_orbitals})"
            if n_remove
            else "removed no virtual orbital"
        )
        active_space = (
            f"{config['label']} STO-3G: {n_ao} spatial orbitals in the full basis; {frozen_text} "
            f"and {removed_text}; active space = {active_spatial} spatial / "
            f"{problem.num_spin_orbitals} spin orbitals, {sum(problem.num_particles)} active "
            "electrons."
        )
        fermionic_op = problem.hamiltonian.second_q_op()
        mapper = ParityMapper(num_particles=problem.num_particles)
    else:
        active_space = (
            f"{config['label']} STO-3G: full space (no core frozen, no virtual orbitals "
            f"removed); {problem.num_spatial_orbitals} spatial / {problem.num_spin_orbitals} spin "
            f"orbitals, {sum(problem.num_particles)} electrons."
        )
        fermionic_op = _fermionic_operator(problem)
        # H2: 2 electrons in the (1, 1) sector, the configuration this project's
        # VQE code is built for.
        mapper = ParityMapper(num_particles=(1, 1))

    operator = mapper.map(fermionic_op)
    n_qubits = int(operator.num_qubits)

    expected = config["expected_qubits"]
    if expected is not None and n_qubits != expected:
        raise ValueError(
            f"expected {expected} qubits from parity + two-qubit reduction at {bond_length} A, "
            f"got {n_qubits}. The mapper configuration assumed by app/core vqe code no longer holds."
        )

    energy_total = _ground_energy(operator)

    # Fix the convention so that, for every point,
    #     smallest eigenvalue of the stored operator + nuclear_repulsion == FCI.
    # Qiskit Nature may or may not fold nuclear repulsion into the constant term;
    # detect which, rather than assuming.
    if abs(energy_total - fci_energy) <= TOLERANCE_HA:
        nuclear_term_present = True
    elif abs(energy_total + nuclear_repulsion - fci_energy) <= TOLERANCE_HA:
        nuclear_term_present = False
    else:
        raise ValueError(
            f"qubit operator does not reproduce FCI at {bond_length} A: "
            f"diag={energy_total:.12f}, diag+nuc={energy_total + nuclear_repulsion:.12f}, "
            f"fci={fci_energy:.12f}"
        )

    if nuclear_term_present:
        identity = SparsePauliOp(["I" * n_qubits], coeffs=[-nuclear_repulsion])
        operator = operator + identity
        energy_total = _ground_energy(operator)

    deviation = abs(energy_total + nuclear_repulsion - fci_energy)
    if deviation > TOLERANCE_HA:
        raise ValueError(
            f"post-adjustment check failed at {bond_length} A: deviation {deviation:.3e} Ha"
        )

    point = {
        "bond_length_angstrom": bond_length,
        "nuclear_repulsion_ha": nuclear_repulsion,
        "hf_energy_ha": hf_energy,
        "fci_energy_ha": fci_energy,
        "qubit_energy_ha": energy_total,
        "qubit_energy_plus_nuclear_ha": energy_total + nuclear_repulsion,
        "deviation_from_fci_ha": deviation,
        "n_qubits": n_qubits,
        "pauli_terms": _pauli_terms(operator),
        "nuclear_repulsion_in_qubit_operator": False,
        "energy_span_two_state_ha": None,
    }
    return point, active_space


def main() -> None:
    parser = argparse.ArgumentParser(description="Precompute STO-3G qubit Hamiltonians offline.")
    parser.add_argument(
        "out_path",
        nargs="?",
        default=None,
        help="output JSON (defaults to the repository artifact path for the molecule)",
    )
    parser.add_argument("--molecule", choices=sorted(MOLECULES), default=DEFAULT_MOLECULE)
    args = parser.parse_args()

    molecule_key = args.molecule
    config = MOLECULES[molecule_key]
    bond_lengths: list[float] = config["bond_lengths"]

    versions = _require_dependencies()
    print("versions:", ", ".join(f"{k}={v}" for k, v in versions.items()))
    print(f"platform: {platform.platform()}  python {sys.version.split()[0]}")
    print(f"molecule: {config['label']} ({molecule_key})")
    print(f"bond lengths: {len(bond_lengths)} points, {bond_lengths[0]} to {bond_lengths[-1]} A")

    points: list[dict[str, Any]] = []
    worst = 0.0
    active_space: str | None = None
    for bond_length in bond_lengths:
        point, point_active_space = build_point(bond_length, molecule_key)
        if active_space is None:
            active_space = point_active_space
        elif point_active_space != active_space:
            raise SystemExit("active-space description changed between bond lengths")
        points.append(point)
        worst = max(worst, point["deviation_from_fci_ha"])
        print(
            f"  r={bond_length:6.3f} A  HF={point['hf_energy_ha']:.6f}  "
            f"FCI={point['fci_energy_ha']:.6f}  "
            f"qubit+nuc={point['qubit_energy_plus_nuclear_ha']:.6f}  "
            f"|dev|={point['deviation_from_fci_ha']:.2e}  "
            f"terms={len(point['pauli_terms'])}"
        )

    if worst > TOLERANCE_HA:
        raise SystemExit(f"refusing to write: worst deviation {worst:.3e} Ha exceeds {TOLERANCE_HA}")

    n_qubits = int(points[0]["n_qubits"])
    payload = {
        "molecule": config["label"],
        "basis": BASIS,
        "mapping": MAPPING,
        "two_qubit_reduction": TWO_QUBIT_REDUCTION,
        "n_qubits": n_qubits,
        "units": "hartree",
        "active_space": active_space,
        "generated_by": "scripts/precompute_hamiltonians.py",
        "library_versions": versions,
        "platform": platform.platform(),
        "convention": (
            "Each entry's pauli_terms describe the electronic Hamiltonian only. "
            "The smallest eigenvalue of those terms plus nuclear_repulsion_ha equals "
            "fci_energy_ha. Determine the ground energy by adding nuclear_repulsion_ha."
        ),
        "max_deviation_from_fci_ha": worst,
        "point_count": len(points),
        # Declared so the validator reads the grid from the file instead of
        # assuming a point count or a step.
        "bond_length_grid": list(bond_lengths),
        "points": points,
    }

    destination = resolve_output_path(molecule_key, args.out_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"\nwrote {destination}")
    print(f"points={len(points)}  worst deviation={worst:.3e} Ha  n_qubits={n_qubits}")
    print(f"Copy this file into app/data/hamiltonians/{config['file']} in the repository.")


if __name__ == "__main__":
    main()
