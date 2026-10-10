# LIMITATIONS.md

What this build does **not** verify. Anything listed here is either unmeasured, unavailable, or known to be weaker than it looks; it must not be presented as a result.

## Known limitations

- **The Colab/WSL Hamiltonian script has never been executed.** `scripts/precompute_hamiltonians.py` and `scripts/precompute_hamiltonians.ipynb` have not been run, so `app/data/hamiltonians/` holds only `.gitkeep`. No H2 curve numbers exist, and the Quantum Lab reports the curve as unavailable rather than showing placeholder values.
- **RDKit tests are skipped where RDKit is blocked.** On a host where the RDKit DLL is blocked, the RDKit-dependent tests skip instead of running. A skip is not a pass, and the conformer/property paths are unverified on that host.
- **The ligand N/O donor/acceptor rule is permissive.** Hydrogen-bond contact counts come from a geometric heuristic over deposited heavy atoms (no hydrogens, no protonation states), so they are **upper bounds**, not observed hydrogen bonds. The same applies to the salt-bridge rule.
- **CI has never run on a runner.** `.github/workflows/ci.yml` exists, but no job has executed on a GitHub-hosted runner, so "tests pass in CI" is unverified.
- **Docking is not built.** There is no docking pipeline, no Vina run and no RMSD. The Drug-Target Studio reports measured contacts from deposited crystal coordinates, which need no docking.
