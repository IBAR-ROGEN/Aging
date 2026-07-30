"""Ensembl REST helpers pinned to release 116 (GRCh38).

Provides an HTTP client with polite pacing, ``Retry-After`` handling, and
exponential backoff, plus SQLite/JSON caching for variant lookups.
"""

from rogen_aging.ensembl.cache import (
    JsonFileCache,
    SqliteCache,
    cache_key_for,
    open_cache,
)
from rogen_aging.ensembl.client import (
    DEFAULT_BASE_URL,
    ENSEMBL_ASSEMBLY,
    ENSEMBL_RELEASE,
    EnsemblApiError,
    EnsemblClient,
)
from rogen_aging.ensembl.variants import (
    fetch_variant_cached,
    grch38_locus_from_variation,
    lookup_grch38_loci_cached,
    lookup_variants_cached,
)

__all__ = [
    "DEFAULT_BASE_URL",
    "ENSEMBL_ASSEMBLY",
    "ENSEMBL_RELEASE",
    "EnsemblApiError",
    "EnsemblClient",
    "JsonFileCache",
    "SqliteCache",
    "cache_key_for",
    "fetch_variant_cached",
    "grch38_locus_from_variation",
    "lookup_grch38_loci_cached",
    "lookup_variants_cached",
    "open_cache",
]
