from __future__ import annotations

from pathlib import Path
from typing import Any

from app.config import get_settings


def resolve_query(query: str) -> dict[str, Any]:
    """Resolve a short text query to a verified identifier when possible.

    This build has no alias database, so most queries return ``no_match``
    honestly. Only a couple of curried scaffold examples are recognized so the
    contract shape is demonstrable without inventing real mappings.
    """
    normalized = query.strip().lower()
    if normalized == "covid protease":
        return {
            "query": query,
            "resolved": {"pdb_id": "6LU7"},
            "method": "alias_lookup",
            "status": "resolved",
            "note": "Recognized scaffold example; no live alias DB in this build",
        }
    if normalized == "hiv protease":
        return {
            "query": query,
            "resolved": {"complex_id": "hiv-protease-indinavir"},
            "method": "alias_lookup",
            "status": "resolved",
            "note": "Recognized scaffold example; no live alias DB in this build",
        }
    return {
        "query": query,
        "resolved": None,
        "method": "alias_lookup",
        "status": "no_match",
        "note": "No live alias lookup in this build",
    }


def predict_properties(smiles: str, models_manifest_path: Path) -> dict[str, Any]:
    """Property prediction stub.

    Until a trained model is loaded, every response is honest: ``model:
    "not trained"`` with ``method: "heuristic"``.
    """
    del models_manifest_path  # reserved for model loading
    settings = get_settings()
    return {
        "smiles": smiles,
        "model": "not trained",
        "method": "heuristic",
        "predictions": None,
        "applicability_warning": None,
        "status": "not_trained",
        "note": "No trained model loaded in this build",
    }
