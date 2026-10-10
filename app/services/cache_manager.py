from __future__ import annotations

from typing import Any


def get_cached(key: str) -> Any:
    """Return cached value for key, or None.

    This build has no functional cache backend. The cache manager is a stub so
    external-api wiring can call it later without special casing import paths.
    """
    del key
    return None


def set_cached(key: str, value: Any, ttl_seconds: int) -> None:
    """Store a value in the cache.

    This build has no functional cache backend.
    """
    del key, value, ttl_seconds
