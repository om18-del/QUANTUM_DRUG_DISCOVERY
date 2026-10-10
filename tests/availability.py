"""Availability helpers for tests that need an external dependency.

Several tests in this suite were previously skipped with an unconditional
``pytest.skip(...)`` whose message claimed a dependency was missing. That is
misleading: it hides whether the code works. The helpers here skip only when the
dependency genuinely is not usable, and they say why.
"""

from __future__ import annotations

import os

import pytest

LIVE_NETWORK_ENV = "QDDE_LIVE_NETWORK"
_TRUTHY = {"1", "true", "yes", "on"}


def live_network_enabled() -> bool:
    return os.environ.get(LIVE_NETWORK_ENV, "").strip().lower() in _TRUTHY


def require_live_network() -> None:
    """Skip unless live-network integration tests were explicitly enabled."""
    if not live_network_enabled():
        pytest.skip(
            f"live network test disabled by default; set {LIVE_NETWORK_ENV}=1 to run it"
        )


def rdkit_reason() -> str | None:
    """Return why ``rdkit.Chem`` is unusable here, or ``None`` when it works."""
    try:
        from rdkit import Chem  # noqa: F401
    except Exception as exc:  # noqa: BLE001 - report whatever the OS/RDKit raises
        return f"{type(exc).__name__}: {exc}"
    return None


def require_rdkit_chem() -> None:
    reason = rdkit_reason()
    if reason is not None:
        pytest.skip(f"rdkit.Chem is not usable in this environment ({reason})")


def require_rdkit_descriptors() -> None:
    """Descriptors need both ``rdkit.Chem`` and ``rdMolDescriptors``."""
    reason = rdkit_reason()
    if reason is not None:
        pytest.skip(f"rdkit.Chem is not usable in this environment ({reason})")
    from app.services.rdkit_helper import rdkit_capabilities

    caps = rdkit_capabilities()
    if not caps.get("descriptors"):
        pytest.skip(f"RDKit descriptors unavailable in this environment (capabilities={caps})")
