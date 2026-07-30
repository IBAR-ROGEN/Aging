"""Ensembl REST HTTP client pinned to release 116 (GRCh38) with polite retries."""

from __future__ import annotations

import logging
import time
from typing import Any
from urllib.parse import quote, urlencode, urljoin

import requests

LOG = logging.getLogger(__name__)

# Ensembl release 116 (June 2026) archive REST — GRCh38 primary assembly.
ENSEMBL_RELEASE = 116
ENSEMBL_ASSEMBLY = "GRCh38"
DEFAULT_BASE_URL = "https://jun2026.rest.ensembl.org"
DEFAULT_SPECIES = "human"

# Ensembl recommends staying well under the global rate limit; ~3 req/s is polite.
DEFAULT_MIN_INTERVAL_SEC = 0.34
DEFAULT_TIMEOUT_SEC = 30.0
DEFAULT_MAX_RETRIES = 5
DEFAULT_BACKOFF_BASE_SEC = 0.5
DEFAULT_BACKOFF_CAP_SEC = 60.0

RETRYABLE_STATUS_CODES = frozenset({429, 500, 502, 503, 504})


class EnsemblApiError(RuntimeError):
    """Raised when an Ensembl REST request fails after retries."""

    def __init__(
        self,
        message: str,
        *,
        status_code: int | None = None,
        url: str | None = None,
    ) -> None:
        super().__init__(message)
        self.status_code = status_code
        self.url = url


def _parse_retry_after(header_value: str | None, fallback: float) -> float:
    """Parse ``Retry-After`` as seconds; fall back when missing or non-numeric."""
    if header_value is None:
        return fallback
    try:
        return max(float(header_value), 0.0)
    except ValueError:
        # HTTP-date form is rare for Ensembl; treat as unusable.
        return fallback


def _backoff_seconds(attempt: int, *, base: float, cap: float) -> float:
    """Exponential backoff for attempt ``1..N`` (attempt 1 → ``base``)."""
    return min(cap, base * (2 ** max(attempt - 1, 0)))


class EnsemblClient:
    """Thin Ensembl REST wrapper with pacing, ``Retry-After``, and backoff.

    Targets Ensembl **release 116** on **GRCh38** via the June 2026 archive
    REST host by default. Override ``base_url`` only when you intentionally need
    a different archive or the live ``rest.ensembl.org`` mirror.
    """

    def __init__(
        self,
        *,
        base_url: str = DEFAULT_BASE_URL,
        release: int = ENSEMBL_RELEASE,
        assembly: str = ENSEMBL_ASSEMBLY,
        species: str = DEFAULT_SPECIES,
        min_interval_sec: float = DEFAULT_MIN_INTERVAL_SEC,
        timeout_sec: float = DEFAULT_TIMEOUT_SEC,
        max_retries: int = DEFAULT_MAX_RETRIES,
        backoff_base_sec: float = DEFAULT_BACKOFF_BASE_SEC,
        backoff_cap_sec: float = DEFAULT_BACKOFF_CAP_SEC,
        session: requests.Session | None = None,
        user_agent: str = (
            "rogen-aging-ensembl/0.1 "
            "(ROGEN; academic research; Ensembl REST 116/GRCh38)"
        ),
    ) -> None:
        if release != ENSEMBL_RELEASE:
            LOG.warning(
                "Client configured for Ensembl release %s; module defaults target %s.",
                release,
                ENSEMBL_RELEASE,
            )
        self.base_url = base_url.rstrip("/") + "/"
        self.release = release
        self.assembly = assembly
        self.species = species
        self.min_interval_sec = min_interval_sec
        self.timeout_sec = timeout_sec
        self.max_retries = max_retries
        self.backoff_base_sec = backoff_base_sec
        self.backoff_cap_sec = backoff_cap_sec
        self._owns_session = session is None
        self._session = session or requests.Session()
        self._session.headers.setdefault("Content-Type", "application/json")
        self._session.headers.setdefault("Accept", "application/json")
        self._session.headers.setdefault("User-Agent", user_agent)
        self._last_request_end: float | None = None

    def close(self) -> None:
        """Close the owned HTTP session, if any."""
        if self._owns_session:
            self._session.close()

    def __enter__(self) -> EnsemblClient:
        return self

    def __exit__(self, *_exc: object) -> None:
        self.close()

    def _pace(self) -> None:
        if self._last_request_end is None:
            return
        elapsed = time.monotonic() - self._last_request_end
        wait = self.min_interval_sec - elapsed
        if wait > 0:
            time.sleep(wait)

    def build_url(self, path: str, params: dict[str, Any] | None = None) -> str:
        """Join ``path`` onto the base URL and optionally append query params."""
        clean = path.lstrip("/")
        url = urljoin(self.base_url, clean)
        if params:
            query = urlencode({k: v for k, v in params.items() if v is not None}, doseq=True)
            if query:
                sep = "&" if "?" in url else "?"
                url = f"{url}{sep}{query}"
        return url

    def request_json(
        self,
        method: str,
        path: str,
        *,
        params: dict[str, Any] | None = None,
        json_body: Any | None = None,
        allow_404: bool = False,
    ) -> Any | None:
        """Perform an HTTP request and return parsed JSON.

        Respects the minimum inter-request interval, retries transient failures
        with exponential backoff, and honours ``Retry-After`` on ``429`` /
        overload responses.

        Args:
            method: HTTP method (``GET``, ``POST``, …).
            path: Path relative to the Ensembl REST base URL.
            params: Optional query string parameters.
            json_body: Optional JSON body for ``POST``/``PUT``.
            allow_404: When true, return ``None`` for HTTP 404 instead of raising.

        Returns:
            Parsed JSON payload, or ``None`` when ``allow_404`` and the resource
            is missing.

        Raises:
            EnsemblApiError: When retries are exhausted or a non-retryable error
                occurs.
        """
        url = self.build_url(path, params)
        attempt = 0
        while True:
            attempt += 1
            self._pace()
            try:
                response = self._session.request(
                    method.upper(),
                    url,
                    json=json_body,
                    timeout=self.timeout_sec,
                )
            except requests.RequestException as exc:
                if attempt > self.max_retries:
                    raise EnsemblApiError(
                        f"Ensembl request failed after {self.max_retries} attempts: {exc}",
                        url=url,
                    ) from exc
                sleep_s = _backoff_seconds(
                    attempt,
                    base=self.backoff_base_sec,
                    cap=self.backoff_cap_sec,
                )
                LOG.warning(
                    "Ensembl network error (%s); retry %s/%s in %.2fs",
                    exc,
                    attempt,
                    self.max_retries,
                    sleep_s,
                )
                time.sleep(sleep_s)
                continue
            finally:
                self._last_request_end = time.monotonic()

            if response.status_code == 404 and allow_404:
                return None

            if response.status_code in RETRYABLE_STATUS_CODES:
                if attempt > self.max_retries:
                    raise EnsemblApiError(
                        f"Ensembl HTTP {response.status_code} after {self.max_retries} attempts",
                        status_code=response.status_code,
                        url=url,
                    )
                fallback = _backoff_seconds(
                    attempt,
                    base=self.backoff_base_sec,
                    cap=self.backoff_cap_sec,
                )
                sleep_s = _parse_retry_after(response.headers.get("Retry-After"), fallback)
                # Prefer the server's Retry-After, but never sleep less than backoff floor.
                sleep_s = max(sleep_s, min(fallback, self.backoff_cap_sec))
                remaining = response.headers.get("X-RateLimit-Remaining")
                LOG.warning(
                    "Ensembl HTTP %s (Retry-After=%s, remaining=%s); sleeping %.2fs",
                    response.status_code,
                    response.headers.get("Retry-After"),
                    remaining,
                    sleep_s,
                )
                time.sleep(sleep_s)
                continue

            if not response.ok:
                raise EnsemblApiError(
                    f"Ensembl HTTP {response.status_code}: {response.text[:300]}",
                    status_code=response.status_code,
                    url=url,
                )

            if not response.content:
                return None
            return response.json()

    def get_json(
        self,
        path: str,
        *,
        params: dict[str, Any] | None = None,
        allow_404: bool = False,
    ) -> Any | None:
        """GET ``path`` and return parsed JSON."""
        return self.request_json("GET", path, params=params, allow_404=allow_404)

    def post_json(
        self,
        path: str,
        *,
        json_body: Any,
        params: dict[str, Any] | None = None,
        allow_404: bool = False,
    ) -> Any | None:
        """POST JSON to ``path`` and return the parsed response."""
        return self.request_json(
            "POST",
            path,
            params=params,
            json_body=json_body,
            allow_404=allow_404,
        )

    def get_variation(self, variant_id: str, *, phenotypes: bool = False) -> dict[str, Any] | None:
        """Fetch a variation record by rsID / variant name (GRCh38 mappings).

        Endpoint: ``GET /variation/{species}/{id}``.
        """
        params: dict[str, Any] = {"content-type": "application/json"}
        if phenotypes:
            params["phenotypes"] = 1
        path = f"variation/{quote(self.species, safe='')}/{quote(variant_id, safe='')}"
        payload = self.get_json(path, params=params, allow_404=True)
        if payload is None:
            return None
        if not isinstance(payload, dict):
            raise EnsemblApiError(
                f"Unexpected variation payload type for {variant_id!r}: {type(payload).__name__}",
                url=path,
            )
        return payload

    def get_vep_id(
        self,
        variant_id: str,
        *,
        extra_params: dict[str, Any] | None = None,
    ) -> list[dict[str, Any]] | None:
        """Run VEP for a single variant identifier.

        Endpoint: ``GET /vep/{species}/id/{id}``.
        """
        params: dict[str, Any] = {"content-type": "application/json"}
        if extra_params:
            params.update(extra_params)
        path = f"vep/{quote(self.species, safe='')}/id/{quote(variant_id, safe='')}"
        payload = self.get_json(path, params=params, allow_404=True)
        if payload is None:
            return None
        if isinstance(payload, dict):
            return [payload]
        if not isinstance(payload, list):
            raise EnsemblApiError(
                f"Unexpected VEP payload type for {variant_id!r}: {type(payload).__name__}",
                url=path,
            )
        return payload

    def info_data(self) -> dict[str, Any]:
        """Return ``/info/data`` (release numbers served by this REST host)."""
        payload = self.get_json("info/data", params={"content-type": "application/json"})
        if not isinstance(payload, dict):
            raise EnsemblApiError("Unexpected /info/data payload")
        return payload
