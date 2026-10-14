"""Print the five documented STO-3G resource estimates.

Run from the repository root (no server needed -- the endpoint is called through
the in-process TestClient):

    python scripts/estimate_resources_examples.py

The table it prints is the one quoted in RESULTS.md and docs/API_CONTRACT.md.
Every figure comes from ``GET /api/quantum/resources?formula=`` and is counted by
formula in the STO-3G minimal basis; nothing here is typed in by hand, and no
simulation is run for any of these molecules.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Running ``python scripts/this_file.py`` puts scripts/ on sys.path, not the
# repository root, so the app package has to be reachable explicitly.
REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from fastapi.testclient import TestClient  # noqa: E402  (after sys.path)

from app.main import app  # noqa: E402  (after sys.path)

EXAMPLES = ("H2", "LiH", "H2O", "CH4", "C9H8O4")

HEADER = (
    "| formula | spatial orbitals | electrons | full qubits | reduced qubits | "
    "statevector memory | verdict |"
)


def main() -> int:
    client = TestClient(app)
    print(HEADER)
    print("| --- | --- | --- | --- | --- | --- | --- |")

    for formula in EXAMPLES:
        response = client.get("/api/quantum/resources", params={"formula": formula})
        if response.status_code != 200:
            print(
                f"{formula}: HTTP {response.status_code} {response.text}",
                file=sys.stderr,
            )
            return 1
        payload = response.json()
        print(
            f"| `{formula}` | {payload['n_spatial_orbitals']} | {payload['n_electrons']} | "
            f"{payload['n_qubits_full']} | {payload['n_qubits_after_parity_reduction']} | "
            f"{payload['statevector_memory_human']} | {payload['feasibility']['verdict']} |"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
