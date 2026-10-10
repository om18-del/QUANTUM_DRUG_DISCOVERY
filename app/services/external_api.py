from __future__ import annotations

import json
import logging
from datetime import timedelta
from pathlib import Path
from typing import Any

import diskcache
import httpx
from tenacity import (
    retry,
    retry_if_exception_type,
    stop_after_attempt,
    wait_exponential,
)

from app.config import get_settings

logger = logging.getLogger(__name__)

PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
RSCB_PDB_BASE = "https://files.rcsb.org"


def _cache() -> diskcache.Cache:
    settings = get_settings()
    settings.cache_dir.mkdir(parents=True, exist_ok=True)
    return diskcache.Cache(str(settings.cache_dir))


def _cache_key(prefix: str, *parts: str) -> str:
    return f"{prefix}:{':'.join(parts)}"


def _cache_ttl(settings, source: str) -> int:
    if source == "pubchem":
        return settings.cache_ttl_pubchem_hours * 3600
    if source == "pdb":
        return settings.cache_ttl_pdb_hours * 3600
    if source == "smiles_descriptors":
        return settings.cache_ttl_smiles_descriptors_hours * 3600
    return 3600


class UpstreamError(RuntimeError):
    """Wrapped upstream API failure for tenacity retry eligibility."""


class NotFoundError(RuntimeError):
    """The upstream has no such record (HTTP 404/410).

    This is an answer, not a failure, so it is never retried. Its message is
    deliberately free of URLs and of library detail, because it is shown to the
    person asking for the compound.
    """


def _is_retryable_status(status: int) -> bool:
    """Only a status that could change by asking again is worth retrying.

    404/410 mean "no such record" and 400 means our request is wrong; resending
    those cannot alter the answer. A rate limit (429) or a server error (5xx) can,
    so those are retried with backoff.
    """
    return status == 429 or status >= 500


def time_sleep(seconds: float) -> None:
    """Sleep via the real ``time`` module, so tests can patch one symbol."""
    import time

    time.sleep(seconds)


def _http_get_json(
    client: httpx.Client,
    url: str,
    source: str,
    retries: int,
    backoff_seconds: float,
) -> Any:
    last: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            response = client.get(url, timeout=10.0, follow_redirects=True)
        except httpx.TimeoutException as exc:
            last = exc
            logger.warning("Upstream %s timeout, attempt %d/%d", source, attempt, retries)
            if attempt < retries:
                time_sleep(backoff_seconds * (2 ** (attempt - 1)))
            continue
        except httpx.NetworkError as exc:
            last = exc
            logger.warning("Upstream %s network error, attempt %d/%d", source, attempt, retries)
            if attempt < retries:
                time_sleep(backoff_seconds * (2 ** (attempt - 1)))
            continue

        status = response.status_code
        if status >= 400:
            if status in (404, 410):
                # "Not found" is an answer. Stop immediately: retrying a name the
                # registry does not have only burns time, and the message must not
                # carry the URL or the HTTP library's wording.
                logger.info(
                    "Upstream %s has no record for that request (HTTP %d)", source, status
                )
                raise NotFoundError(f"{source} has no record for that request")
            if status == 400:
                raise UpstreamError(f"{source} rejected the request as malformed (HTTP 400)")
            if not _is_retryable_status(status):
                raise UpstreamError(f"{source} returned HTTP {status}")

            last = UpstreamError(f"HTTP {status}")
            logger.warning(
                "Upstream %s returned HTTP %s, attempt %d/%d",
                source,
                status,
                attempt,
                retries,
            )
            if attempt < retries:
                time_sleep(backoff_seconds * (2 ** (attempt - 1)))
            continue

        try:
            return response.json()
        except ValueError as exc:
            raise UpstreamError(f"{source} returned a response that is not valid JSON") from exc

    if last is not None:
        raise UpstreamError(f"{source} failed after {retries} attempts: {last}") from last
    raise UpstreamError(f"{source} failed after {retries} attempts")


def http_get_json_cached(
    url: str,
    source: str = "pubchem",
    retries: int | None = None,
    backoff_seconds: float = 1.0,
) -> Any:
    settings = get_settings()
    if retries is None:
        retries = settings.external_max_retries
    cache = _cache()
    key = _cache_key(source, url)
    ttl = _cache_ttl(settings, source)
    try:
        value = cache.get(key)
        if value is not None:
            logger.debug("Cache hit: %s", key[:80])
            return value
    except Exception as exc:
        logger.warning("Cache read failed for %s: %s", key[:80], exc)

    value = _http_get_json(
        httpx.Client(timeout=10.0, follow_redirects=True),
        url,
        source,
        retries,
        backoff_seconds,
    )
    try:
        cache.set(key, value, expire=ttl)
    except Exception as exc:
        logger.warning("Cache write failed for %s: %s", key[:80], exc)
    return value


def http_get_bytes_cached(
    url: str,
    source: str = "pubchem",
    retries: int | None = None,
    backoff_seconds: float = 1.0,
) -> bytes:
    settings = get_settings()
    if retries is None:
        retries = settings.external_max_retries
    cache = _cache()
    key = _cache_key(source, url)
    ttl = _cache_ttl(settings, source)
    try:
        value = cache.get(key)
        if isinstance(value, bytes):
            logger.debug("Bytes cache hit: %s", key[:80])
            return value
    except Exception as exc:
        logger.warning("Bytes cache read failed for %s: %s", key[:80], exc)

    last: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            response = httpx.Client(timeout=10.0, follow_redirects=True).get(
                url, timeout=10.0, follow_redirects=True
            )
        except httpx.TimeoutException as exc:
            last = exc
            logger.warning("Upstream %s bytes timeout, attempt %d/%d", source, attempt, retries)
            if attempt < retries:
                time_sleep(backoff_seconds * (2 ** (attempt - 1)))
            continue
        except httpx.NetworkError as exc:
            last = exc
            logger.warning(
                "Upstream %s bytes network error, attempt %d/%d", source, attempt, retries
            )
            if attempt < retries:
                time_sleep(backoff_seconds * (2 ** (attempt - 1)))
            continue

        status = response.status_code
        if status >= 400:
            # Same policy as the JSON path: 404/410 are answers, 400 is our fault,
            # and only a rate limit or a server error is worth asking again.
            if status in (404, 410):
                logger.info(
                    "Upstream %s bytes: no record for that request (HTTP %d)", source, status
                )
                raise NotFoundError(f"{source} has no record for that request")
            if status == 400:
                raise UpstreamError(f"{source} rejected the request as malformed (HTTP 400)")
            if not _is_retryable_status(status):
                raise UpstreamError(f"{source} returned HTTP {status}")

            last = UpstreamError(f"HTTP {status}")
            logger.warning(
                "Upstream %s bytes returned HTTP %s, attempt %d/%d",
                source,
                status,
                attempt,
                retries,
            )
            if attempt < retries:
                time_sleep(backoff_seconds * (2 ** (attempt - 1)))
            continue

        value = response.content
        break
    else:
        if last is not None:
            raise UpstreamError(f"{source} bytes failed after {retries} attempts: {last}") from last
        raise UpstreamError(f"{source} bytes failed after {retries} attempts")

    try:
        cache.set(key, value, expire=ttl)
    except Exception as exc:
        logger.warning("Bytes cache write failed for %s: %s", key[:80], exc)
    return value


def invalidate_cached_url(url: str, source: str = "pubchem") -> None:
    cache = _cache()
    key = _cache_key(source, url)
    try:
        cache.delete(key)
    except Exception as exc:
        logger.warning("Cache invalidation failed for %s: %s", key[:80], exc)
