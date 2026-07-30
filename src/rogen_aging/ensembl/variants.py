"""Cached Ensembl variant lookups (SQLite or JSON) over :class:`EnsemblClient`."""

from __future__ import annotations

import logging
from collections.abc import Iterable, Sequence
from pathlib import Path
from typing import Any, Literal

from rogen_aging.ensembl.cache import (
    CacheBackend,
    ResponseCache,
    cache_key_for,
    open_cache,
)
from rogen_aging.ensembl.client import EnsemblClient

LOG = logging.getLogger(__name__)

VariantEndpoint = Literal["variation", "vep"]


def fetch_variant_cached(
    variant_id: str,
    *,
    client: EnsemblClient,
    cache: ResponseCache,
    endpoint: VariantEndpoint = "variation",
    refresh: bool = False,
    vep_params: dict[str, Any] | None = None,
    phenotypes: bool = False,
) -> Any | None:
    """Fetch one variant from Ensembl, reusing a local cache entry when present.

    On a cache hit the live API is not contacted. On a miss the response is
    written back so subsequent calls for the same identifier skip the network.

    Args:
        variant_id: dbSNP-style identifier (e.g. ``rs429358``).
        client: Configured :class:`EnsemblClient` (release 116 / GRCh38).
        cache: JSON or SQLite :class:`ResponseCache`.
        endpoint: ``variation`` (default) or ``vep``.
        refresh: When true, ignore any existing cache entry and re-query.
        vep_params: Extra VEP query parameters when ``endpoint='vep'``.
        phenotypes: Include phenotypes for ``variation`` lookups.

    Returns:
        Parsed Ensembl JSON (dict for variation, list for VEP), or ``None``
        when the variant is unknown (HTTP 404). Negative results are also
        cached so repeated misses stay offline.
    """
    key_params: dict[str, Any] | None = None
    if endpoint == "vep" and vep_params:
        key_params = dict(vep_params)
    elif endpoint == "variation" and phenotypes:
        key_params = {"phenotypes": 1}

    key = cache_key_for(endpoint, variant_id, params=key_params)

    if not refresh:
        cached = cache.get(key)
        if cached is not None:
            LOG.debug("Ensembl cache hit for %s (%s)", variant_id, endpoint)
            if isinstance(cached, dict) and cached.get("__ensembl_missing__") is True:
                return None
            return cached

    LOG.debug("Ensembl cache miss for %s (%s)", variant_id, endpoint)
    if endpoint == "variation":
        payload: Any | None = client.get_variation(variant_id, phenotypes=phenotypes)
    elif endpoint == "vep":
        payload = client.get_vep_id(variant_id, extra_params=vep_params)
    else:
        raise ValueError(f"Unsupported endpoint: {endpoint!r}")

    if payload is None:
        cache.set(key, {"__ensembl_missing__": True, "id": variant_id})
    else:
        cache.set(key, payload)
    return payload


def lookup_variants_cached(
    variant_ids: Iterable[str],
    *,
    cache_path: Path | str,
    backend: CacheBackend | None = None,
    client: EnsemblClient | None = None,
    endpoint: VariantEndpoint = "variation",
    refresh: bool = False,
    vep_params: dict[str, Any] | None = None,
    phenotypes: bool = False,
    skip_empty: bool = True,
) -> dict[str, Any | None]:
    """Resolve many variant IDs via Ensembl with a shared local cache.

    This is the seamless integration point: open a SQLite or JSON cache,
    reuse already-processed variants, and only query Ensembl for misses.
    Rate limiting and exponential backoff are handled by ``client``.

    Args:
        variant_ids: Identifiers to resolve (duplicates are collapsed).
        cache_path: Directory (JSON files) or ``.sqlite``/``.db`` file path.
        backend: Force ``json`` or ``sqlite``; otherwise inferred from path.
        client: Optional shared client; one is created and closed when omitted.
        endpoint: ``variation`` or ``vep``.
        refresh: Re-query all IDs, overwriting cache entries.
        vep_params: Extra VEP query parameters when ``endpoint='vep'``.
        phenotypes: Include phenotypes for ``variation`` lookups.
        skip_empty: Drop blank / whitespace-only identifiers.

    Returns:
        Mapping from each unique input identifier to its Ensembl payload
        (or ``None`` if not found).
    """
    ordered: list[str] = []
    seen: set[str] = set()
    for raw in variant_ids:
        if raw is None:
            continue
        token = str(raw).strip()
        if skip_empty and not token:
            continue
        if token in seen:
            continue
        seen.add(token)
        ordered.append(token)

    owns_client = client is None
    ensembl = client or EnsemblClient()
    cache = open_cache(cache_path, backend=backend)
    results: dict[str, Any | None] = {}
    try:
        for variant_id in ordered:
            results[variant_id] = fetch_variant_cached(
                variant_id,
                client=ensembl,
                cache=cache,
                endpoint=endpoint,
                refresh=refresh,
                vep_params=vep_params,
                phenotypes=phenotypes,
            )
    finally:
        cache.close()
        if owns_client:
            ensembl.close()

    return results


def grch38_locus_from_variation(payload: dict[str, Any] | None) -> tuple[str, int] | None:
    """Extract primary-assembly GRCh38 ``(chromosome, position)`` from a variation payload."""
    if not payload:
        return None
    mappings = payload.get("mappings") or []
    chromosomal: list[dict[str, Any]] = []
    for mapping in mappings:
        if mapping.get("assembly_name") != "GRCh38":
            continue
        # Prefer chromosome; allow missing coord_system from some archives.
        if mapping.get("coord_system") not in (None, "chromosome"):
            continue
        seq_region = str(mapping.get("seq_region_name", ""))
        if not seq_region or seq_region.upper().endswith("_PATCH"):
            continue
        chromosomal.append(mapping)

    preferred = [m for m in chromosomal if m.get("coord_system") == "chromosome"] or chromosomal
    if not preferred:
        return None
    preferred.sort(
        key=lambda item: (
            len(str(item.get("seq_region_name", ""))),
            str(item.get("seq_region_name", "")),
        )
    )
    chosen = preferred[0]
    chrom = chosen.get("seq_region_name")
    start = chosen.get("start")
    if chrom is None or start is None:
        return None
    return str(chrom), int(start)


def lookup_grch38_loci_cached(
    variant_ids: Sequence[str],
    *,
    cache_path: Path | str,
    backend: CacheBackend | None = None,
    client: EnsemblClient | None = None,
    refresh: bool = False,
) -> dict[str, tuple[str, int] | None]:
    """Resolve rsIDs to GRCh38 loci using cached ``/variation`` responses."""
    payloads = lookup_variants_cached(
        variant_ids,
        cache_path=cache_path,
        backend=backend,
        client=client,
        endpoint="variation",
        refresh=refresh,
    )
    return {
        variant_id: grch38_locus_from_variation(
            payload if isinstance(payload, dict) else None
        )
        for variant_id, payload in payloads.items()
    }
