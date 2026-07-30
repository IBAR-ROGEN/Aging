"""Local SQLite and JSON-file caches for Ensembl REST responses."""

from __future__ import annotations

import hashlib
import json
import sqlite3
import time
from pathlib import Path
from typing import Any, Literal, Protocol

CacheBackend = Literal["json", "sqlite"]


class ResponseCache(Protocol):
    """Minimal read/write interface for cached Ensembl JSON payloads."""

    def get(self, key: str) -> Any | None:
        """Return the cached payload for ``key``, or ``None`` on a miss."""

    def set(self, key: str, value: Any) -> None:
        """Store ``value`` under ``key``."""

    def delete(self, key: str) -> None:
        """Remove ``key`` if present."""

    def clear(self) -> None:
        """Drop all cached entries."""

    def close(self) -> None:
        """Release any held resources."""


def cache_key_for(endpoint: str, identifier: str, *, params: dict[str, Any] | None = None) -> str:
    """Build a stable cache key for an Ensembl request.

    Args:
        endpoint: Logical endpoint name (e.g. ``variation``, ``vep``).
        identifier: Variant or feature identifier (e.g. ``rs429358``).
        params: Optional query parameters that affect the response.

    Returns:
        A filesystem- and SQLite-safe cache key string.
    """
    base = f"{endpoint}:{identifier.strip().lower()}"
    if not params:
        return base
    encoded = json.dumps(params, sort_keys=True, separators=(",", ":"), default=str)
    digest = hashlib.sha1(encoded.encode("utf-8")).hexdigest()[:12]
    return f"{base}:{digest}"


class JsonFileCache:
    """One JSON file per cache key under a directory."""

    def __init__(self, cache_dir: Path | str) -> None:
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    def _path_for(self, key: str) -> Path:
        digest = hashlib.sha1(key.encode("utf-8")).hexdigest()
        return self.cache_dir / f"{digest}.json"

    def get(self, key: str) -> Any | None:
        path = self._path_for(key)
        if not path.is_file():
            return None
        with path.open(encoding="utf-8") as handle:
            envelope = json.load(handle)
        if not isinstance(envelope, dict) or "payload" not in envelope:
            return envelope
        return envelope["payload"]

    def set(self, key: str, value: Any) -> None:
        path = self._path_for(key)
        envelope = {
            "key": key,
            "cached_at": time.time(),
            "payload": value,
        }
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", encoding="utf-8") as handle:
            json.dump(envelope, handle, indent=2, sort_keys=True)
            handle.write("\n")

    def delete(self, key: str) -> None:
        path = self._path_for(key)
        path.unlink(missing_ok=True)

    def clear(self) -> None:
        for path in self.cache_dir.glob("*.json"):
            path.unlink(missing_ok=True)

    def close(self) -> None:
        return None


class SqliteCache:
    """SQLite-backed response cache keyed by request identity."""

    def __init__(self, db_path: Path | str) -> None:
        self.db_path = Path(db_path)
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        self._conn = sqlite3.connect(self.db_path)
        self._conn.execute(
            """
            CREATE TABLE IF NOT EXISTS ensembl_cache (
                cache_key TEXT PRIMARY KEY,
                payload TEXT NOT NULL,
                cached_at REAL NOT NULL
            )
            """
        )
        self._conn.commit()

    def get(self, key: str) -> Any | None:
        row = self._conn.execute(
            "SELECT payload FROM ensembl_cache WHERE cache_key = ?",
            (key,),
        ).fetchone()
        if row is None:
            return None
        return json.loads(row[0])

    def set(self, key: str, value: Any) -> None:
        payload = json.dumps(value, sort_keys=True)
        self._conn.execute(
            """
            INSERT INTO ensembl_cache (cache_key, payload, cached_at)
            VALUES (?, ?, ?)
            ON CONFLICT(cache_key) DO UPDATE SET
                payload = excluded.payload,
                cached_at = excluded.cached_at
            """,
            (key, payload, time.time()),
        )
        self._conn.commit()

    def delete(self, key: str) -> None:
        self._conn.execute("DELETE FROM ensembl_cache WHERE cache_key = ?", (key,))
        self._conn.commit()

    def clear(self) -> None:
        self._conn.execute("DELETE FROM ensembl_cache")
        self._conn.commit()

    def close(self) -> None:
        self._conn.close()


def open_cache(
    path: Path | str,
    *,
    backend: CacheBackend | None = None,
) -> ResponseCache:
    """Open a JSON-directory or SQLite cache, inferring the backend from ``path``.

    Args:
        path: Directory for JSON files, or a ``.sqlite``/``.db`` file path.
        backend: Force ``json`` or ``sqlite``. When omitted, directories and
            non-``.sqlite``/``.db`` paths use JSON; otherwise SQLite.

    Returns:
        A :class:`ResponseCache` implementation.
    """
    resolved = Path(path)
    chosen = backend
    if chosen is None:
        suffix = resolved.suffix.lower()
        if suffix in {".sqlite", ".db", ".sqlite3"}:
            chosen = "sqlite"
        else:
            chosen = "json"

    if chosen == "sqlite":
        return SqliteCache(resolved)
    return JsonFileCache(resolved)
