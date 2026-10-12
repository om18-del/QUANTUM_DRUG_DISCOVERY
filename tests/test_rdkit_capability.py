"""Tests for the RDKit capability cache: failures must not be cached forever."""

from __future__ import annotations

import builtins

from app.services import rdkit_helper


def _reset_capability_cache(monkeypatch) -> None:
    # A real rdkit.Chem may already be in sys.modules from another test; the import
    # blocker below is what makes ``from rdkit import Chem`` fail deterministically
    # regardless of test order.
    monkeypatch.setattr(rdkit_helper, "_CAPABILITY_CACHE", None)


def _block_rdkit_import(monkeypatch):
    """Make ``import rdkit`` raise ImportError, as a blocked install would."""
    real_import = builtins.__import__

    def blocked(name, *args, **kwargs):
        if name == "rdkit":
            raise ImportError("simulated: rdkit import blocked")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", blocked)
    return real_import


def test_failed_import_is_reported_and_not_cached(monkeypatch) -> None:
    _reset_capability_cache(monkeypatch)
    _block_rdkit_import(monkeypatch)

    first = rdkit_helper.rdkit_compat()

    assert first.get("Chem") is None
    assert "ImportError" in (first.get("note") or "")
    # The failure must not be cached, or a later attempt could never succeed.
    assert rdkit_helper._CAPABILITY_CACHE is None


def test_failed_import_then_success(monkeypatch) -> None:
    _reset_capability_cache(monkeypatch)
    real_import = _block_rdkit_import(monkeypatch)

    assert rdkit_helper.rdkit_compat().get("Chem") is None

    # The same call succeeds as soon as the import resolves again.
    monkeypatch.setattr(builtins, "__import__", real_import)

    assert rdkit_helper.rdkit_compat().get("Chem") is not None


def test_unavailable_reason_includes_the_import_error(monkeypatch) -> None:
    _reset_capability_cache(monkeypatch)
    _block_rdkit_import(monkeypatch)

    reason = rdkit_helper.rdkit_unavailable_reason()

    assert "not importable" in reason
    assert "ImportError" in reason


def test_successful_capability_is_cached(monkeypatch) -> None:
    sentinel = object()
    monkeypatch.setattr(rdkit_helper, "_CAPABILITY_CACHE", {"Chem": sentinel})

    assert rdkit_helper.rdkit_compat()["Chem"] is sentinel
