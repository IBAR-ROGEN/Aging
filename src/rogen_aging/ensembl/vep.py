"""Path-based JSON cache for Ensembl VEP lookups over :class:`EnsemblClient`.

Keeps the on-disk filenames used by the July annotation pipeline and the
LA-SNP VEP script so warm production caches remain valid.
"""

from __future__ import annotations

import json
import logging
from collections.abc import Callable
from pathlib import Path
from typing import Any

from rogen_aging.ensembl.client import EnsemblApiError, EnsemblClient

LOG = logging.getLogger(__name__)

JsonPayload = Any


def as_vep_list(payload: JsonPayload) -> list[dict[str, Any]] | None:
    """Normalize a VEP JSON response to a list of variant dicts.

    Args:
        payload: Raw JSON (list, dict, or ``None``).

    Returns:
        A list of VEP records, or ``None`` when ``payload`` is empty/missing.
    """
    if payload is None:
        return None
    if isinstance(payload, list):
        return [item for item in payload if isinstance(item, dict)]
    if isinstance(payload, dict):
        return [payload]
    return None


def _read_json_file(path: Path) -> JsonPayload | None:
    if not path.is_file():
        return None
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def _write_json_file(path: Path, payload: JsonPayload) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)


def _cached_fetch(
    cache_file: Path,
    *,
    cache_only: bool,
    fetch: Callable[[], JsonPayload | None],
) -> JsonPayload | None:
    cached = _read_json_file(cache_file)
    if cached is not None:
        return cached
    if cache_only:
        LOG.warning("cache miss (cache-only mode): %s", cache_file.name)
        return None
    payload = fetch()
    if payload is None:
        return None
    _write_json_file(cache_file, payload)
    return payload


def fetch_vep_id_cached(
    client: EnsemblClient,
    variant_id: str,
    cache_file: Path,
    *,
    extra_params: dict[str, Any] | None = None,
    cache_only: bool = False,
) -> list[dict[str, Any]] | None:
    """Return cached or live VEP-by-id JSON.

    Args:
        client: Configured Ensembl REST client.
        variant_id: rsID or other VEP identifier.
        cache_file: JSON cache path for this identifier.
        extra_params: Extra VEP query parameters (e.g. ``hgvs``, ``canonical``).
        cache_only: When true, never contact Ensembl on a cache miss.

    Returns:
        VEP records, or ``None`` on 404 / exhausted retries / cache-only miss.
    """

    def _fetch() -> JsonPayload | None:
        try:
            return client.get_vep_id(variant_id, extra_params=extra_params)
        except EnsemblApiError as exc:
            LOG.warning("VEP id lookup failed for %s: %s", variant_id, exc)
            return None

    return as_vep_list(_cached_fetch(cache_file, cache_only=cache_only, fetch=_fetch))


def fetch_vep_region_cached(
    client: EnsemblClient,
    region: str,
    cache_file: Path,
    *,
    extra_params: dict[str, Any] | None = None,
    cache_only: bool = False,
) -> list[dict[str, Any]] | None:
    """Return cached or live VEP-by-region JSON.

    Args:
        client: Configured Ensembl REST client.
        region: Region string ``chrom:start-end/ALT`` (alt already quoted).
        cache_file: JSON cache path for this region.
        extra_params: Extra VEP query parameters.
        cache_only: When true, never contact Ensembl on a cache miss.

    Returns:
        VEP records, or ``None`` on 404 / exhausted retries / cache-only miss.
    """

    def _fetch() -> JsonPayload | None:
        try:
            return client.get_vep_region(region, extra_params=extra_params)
        except EnsemblApiError as exc:
            LOG.warning("VEP region lookup failed for %s: %s", region, exc)
            return None

    return as_vep_list(_cached_fetch(cache_file, cache_only=cache_only, fetch=_fetch))
