"""Precompute real H2 STO-3G qubit Hamiltonians, offline.

This script must be run where PySCF works (Linux/WSL/Colab). PySCF publishes no
Windows wheel and ``qiskit-nature`` declares PySCF with ``sys_platform != "win32"``,
so there is deliberately no Windows path here.

Every coefficient in the output file is produced by the computation below. Nothing
is typed in by hand. The script refuses to write a file unless, for every bond
length, the smallest eigenvalue of the stored qubit operator plus the stored
nuclear repulsion reproduces the PySCF FCI energy to 1e-6 Ha.

Run:
    python scripts/precompute_hamiltonians.py

Output:
    app/data/hamiltonians/h2_sto3g.json
"""

from __future__ import annotations

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
EXPECTED_QUBITS = 2

# The grid: every point of 0.20 .. 3.00 angstrom inclusive in 0.05 steps (57
# uniform points), plus 0.74 angstrom -- the UI default and H2's real equilibrium
# bond length, which the 0.05 step misses (0.74 would need index 10.8). 58 points
# in total, sorted ascending with no duplicates. Built from integers so
# floating-point drift cannot silently change the grid, and the set union makes a
# duplicate impossible even if 0.74 ever lands exactly on the step.
UNIFORM_BOND_LENGTHS: list[float] = [round(0.20 + 0.05 * i, 2) for i in range(57)]
BOND_LENGTHS: list[float] = sorted({*UNIFORM_BOND_LENGTHS, 0.74})
assert len(BOND_LENGTHS) == 58, f"expected 58 bond lengths, got {len(BOND_LENGTHS)}"

TOLERANCE_HA = 1e-6

DEFAULT_OUTPUT_PATH = (
    Path(__file__).resolve().parents[1] / "app" / "data" / "hamiltonians" / "h2_sto3g.json"
)


def output_path() -> Path:
    """Where to write the JSON.

    Precedence: an explicit argv path, then ``QDDE_HAMILTONIAN_OUT``, then the
    repo location. The override exists because in Colab the script does not sit
    inside the repository, so the repo-relative default would be wrong.
    """
    if len(sys.argv) > 1:
        return Path(sys.argv[1]).resolve()
    override = os.environ.get("QDDE_HAMILTONIAN_OUT")
    if override:
        return Path(override).resolve()
    return DEFAULT_OUTPUT_PATH


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


def build_point(bond_length: float) -> dict[str, Any]:
    """Compute one bond length: HF, FCI, and the 2-qubit Hamiltonian."""
    from pyscf import fci as pyscf_fci
    from pyscf import gto, scf
    from qiskit.quantum_info import SparsePauliOp
    from qiskit_nature.second_q.drivers import PySCFDriver
    from qiskit_nature.units import DistanceUnit
    from qiskit_nature.second_q.mappers import ParityMapper

    # PySCF is used directly for the HF and FCI reference energies. Its ``gto.M``
    # accepts the plain ``"Angstrom"`` string; the qiskit-nature driver below needs
    # the DistanceUnit enum instead, because 0.8.0 stores the unit as a
    # DistanceUnit and reads ``self._unit.value`` (raising an AttributeError when it
    # is given a bare string).
    mol = gto.M(
        atom=f"H 0.0 0.0 0.0; H 0.0 0.0 {bond_length}",
        basis=BASIS,
        unit="Angstrom",
        verbose=0,
    )
    mf = scf.RHF(mol)
    mf.verbose = 0
    mf.kernel()
    hf_energy = float(mf.e_tot)

    cisolver = pyscf_fci.FCI(mol, mf.mo_coeff)
    cisolver.verbose = 0
    fci_energy, _ = cisolver.kernel()
    fci_energy = float(fci_energy)

    # Qiskit Nature builds the second-quantised operator from the same basis.
    driver = PySCFDriver(
        atom=f"H 0.0 0.0 0.0; H 0.0 0.0 {bond_length}",
        basis=BASIS,
        unit=DistanceUnit.ANGSTROM,
    )
    problem = driver.run()
    nuclear_repulsion = float(problem.nuclear_repulsion_energy)

    if abs(nuclear_repulsion - float(mol.energy_nuc())) > 1e-6:
        raise ValueError(
            f"nuclear repulsion disagrees between PySCF and Qiskit Nature at "
            f"{bond_length} A: {nuclear_repulsion} vs {mol.energy_nuc()}"
        )

    fermionic_op = _fermionic_operator(problem)

    # H2 in STO-3G: 4 spin orbitals, 2 electrons -> parity mapping with the
    # two-qubit reduction leaves 2 qubits.
    mapper = ParityMapper(num_particles=(1, 1))
    operator = mapper.map(fermionic_op)

    if operator.num_qubits != EXPECTED_QUBITS:
        raise ValueError(
            f"expected {EXPECTED_QUBITS} qubits from parity + two-qubit reduction "
            f"at {bond_length} A, got {operator.num_qubits}. "
            "The mapper configuration assumed by app/core vqe code no longer holds."
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
        identity = SparsePauliOp(["I" * EXPECTED_QUBITS], coeffs=[-nuclear_repulsion])
        operator = operator + identity
        energy_total = _ground_energy(operator)

    deviation = abs(energy_total + nuclear_repulsion - fci_energy)
    if deviation > TOLERANCE_HA:
        raise ValueError(
            f"post-adjustment check failed at {bond_length} A: deviation {deviation:.3e} Ha"
        )

    return {
        "bond_length_angstrom": bond_length,
        "nuclear_repulsion_ha": nuclear_repulsion,
        "hf_energy_ha": hf_energy,
        "fci_energy_ha": fci_energy,
        "qubit_energy_ha": energy_total,
        "qubit_energy_plus_nuclear_ha": energy_total + nuclear_repulsion,
        "deviation_from_fci_ha": deviation,
        "n_qubits": int(operator.num_qubits),
        "pauli_terms": _pauli_terms(operator),
        "nuclear_repulsion_in_qubit_operator": False,
        "energy_span_two_state_ha": None,
    }


def main() -> None:
    versions = _require_dependencies()
    print("versions:", ", ".join(f"{k}={v}" for k, v in versions.items()))
    print(f"platform: {platform.platform()}  python {sys.version.split()[0]}")
    print(f"bond lengths: {len(BOND_LENGTHS)} points, {BOND_LENGTHS[0]} to {BOND_LENGTHS[-1]} A")

    points: list[dict[str, Any]] = []
    worst = 0.0
    for bond_length in BOND_LENGTHS:
        point = build_point(bond_length)
        points.append(point)
        worst = max(worst, point["deviation_from_fci_ha"])
        print(
            f"  r={bond_length:4.2f} A  HF={point['hf_energy_ha']:.6f}  "
            f"FCI={point['fci_energy_ha']:.6f}  "
            f"qubit+nuc={point['qubit_energy_plus_nuclear_ha']:.6f}  "
            f"|dev|={point['deviation_from_fci_ha']:.2e}  "
            f"terms={len(point['pauli_terms'])}"
        )

    if worst > TOLERANCE_HA:
        raise SystemExit(f"refusing to write: worst deviation {worst:.3e} Ha exceeds {TOLERANCE_HA}")

    payload = {
        "molecule": "H2",
        "basis": BASIS,
        "mapping": MAPPING,
        "two_qubit_reduction": TWO_QUBIT_REDUCTION,
        "n_qubits": EXPECTED_QUBITS,
        "units": "hartree",
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
        "points": points,
    }

    destination = output_path()
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"\nwrote {destination}")
    print(f"points={len(points)}  worst deviation={worst:.3e} Ha")
    print(
        "Copy this file into app/data/hamiltonians/h2_sto3g.json in the repository."
    )


if __name__ == "__main__":
    main()
