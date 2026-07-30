"""Unit tests for the Ensembl REST wrapper and local caches."""

from __future__ import annotations

from typing import Any
from unittest.mock import MagicMock

import pytest
import requests

from rogen_aging.ensembl import (
    ENSEMBL_ASSEMBLY,
    ENSEMBL_RELEASE,
    EnsemblApiError,
    EnsemblClient,
    JsonFileCache,
    SqliteCache,
    fetch_variant_cached,
    grch38_locus_from_variation,
    lookup_variants_cached,
    open_cache,
)
from rogen_aging.ensembl.client import _backoff_seconds, _parse_retry_after


def test_release_and_assembly_defaults() -> None:
    assert ENSEMBL_RELEASE == 116
    assert ENSEMBL_ASSEMBLY == "GRCh38"
    client = EnsemblClient()
    assert "jun2026" in client.base_url
    client.close()


def test_parse_retry_after_and_backoff() -> None:
    assert _parse_retry_after("2.5", 1.0) == 2.5
    assert _parse_retry_after(None, 1.0) == 1.0
    assert _parse_retry_after("soon", 1.0) == 1.0
    assert _backoff_seconds(1, base=0.5, cap=60.0) == 0.5
    assert _backoff_seconds(3, base=0.5, cap=60.0) == 2.0
    assert _backoff_seconds(10, base=0.5, cap=3.0) == 3.0


def test_json_and_sqlite_cache_roundtrip(tmp_path: Any) -> None:
    json_cache = JsonFileCache(tmp_path / "json_cache")
    json_cache.set("variation:rs1", {"id": "rs1"})
    assert json_cache.get("variation:rs1") == {"id": "rs1"}
    assert json_cache.get("missing") is None
    json_cache.delete("variation:rs1")
    assert json_cache.get("variation:rs1") is None
    json_cache.close()

    sqlite_cache = SqliteCache(tmp_path / "ensembl.sqlite")
    sqlite_cache.set("variation:rs2", {"id": "rs2", "mappings": []})
    assert sqlite_cache.get("variation:rs2")["id"] == "rs2"
    sqlite_cache.clear()
    assert sqlite_cache.get("variation:rs2") is None
    sqlite_cache.close()

    inferred = open_cache(tmp_path / "auto.sqlite")
    inferred.set("k", {"ok": True})
    assert inferred.get("k") == {"ok": True}
    inferred.close()


def test_client_respects_retry_after(monkeypatch: pytest.MonkeyPatch) -> None:
    sleeps: list[float] = []
    monkeypatch.setattr("rogen_aging.ensembl.client.time.sleep", sleeps.append)

    responses = [
        MagicMock(
            status_code=429,
            headers={"Retry-After": "1.5", "X-RateLimit-Remaining": "0"},
            ok=False,
            text="rate limited",
            content=b"",
        ),
        MagicMock(
            status_code=200,
            headers={},
            ok=True,
            content=b'{"releases":[116]}',
            json=MagicMock(return_value={"releases": [116]}),
        ),
    ]
    session = MagicMock()
    session.headers = {}
    session.request.side_effect = responses

    client = EnsemblClient(session=session, min_interval_sec=0.0, max_retries=3)
    payload = client.info_data()
    assert payload == {"releases": [116]}
    assert sleeps and sleeps[0] >= 1.5
    client.close()


def test_client_exponential_backoff_on_network_error(monkeypatch: pytest.MonkeyPatch) -> None:
    sleeps: list[float] = []
    monkeypatch.setattr("rogen_aging.ensembl.client.time.sleep", sleeps.append)

    session = MagicMock()
    session.headers = {}
    session.request.side_effect = [
        requests.ConnectionError("boom"),
        MagicMock(
            status_code=200,
            headers={},
            ok=True,
            content=b'{"ok":true}',
            json=MagicMock(return_value={"ok": True}),
        ),
    ]

    client = EnsemblClient(
        session=session,
        min_interval_sec=0.0,
        max_retries=3,
        backoff_base_sec=0.25,
    )
    assert client.get_json("info/ping") == {"ok": True}
    assert sleeps == [0.25]
    client.close()


def test_client_raises_after_exhausted_retries(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr("rogen_aging.ensembl.client.time.sleep", lambda _s: None)
    session = MagicMock()
    session.headers = {}
    session.request.return_value = MagicMock(
        status_code=503,
        headers={},
        ok=False,
        text="unavailable",
        content=b"",
    )
    client = EnsemblClient(session=session, min_interval_sec=0.0, max_retries=2)
    with pytest.raises(EnsemblApiError) as exc_info:
        client.get_json("variation/human/rs1")
    assert exc_info.value.status_code == 503
    client.close()


def test_fetch_variant_cached_skips_duplicate_network(tmp_path: Any) -> None:
    session = MagicMock()
    session.headers = {}
    session.request.return_value = MagicMock(
        status_code=200,
        headers={},
        ok=True,
        content=b'{"name":"rs429358","mappings":[]}',
        json=MagicMock(return_value={"name": "rs429358", "mappings": []}),
    )
    client = EnsemblClient(session=session, min_interval_sec=0.0)
    cache = SqliteCache(tmp_path / "variants.sqlite")

    first = fetch_variant_cached("rs429358", client=client, cache=cache)
    second = fetch_variant_cached("rs429358", client=client, cache=cache)

    assert first == second == {"name": "rs429358", "mappings": []}
    assert session.request.call_count == 1
    cache.close()
    client.close()


def test_lookup_variants_cached_json_backend(tmp_path: Any, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr("rogen_aging.ensembl.client.time.sleep", lambda _s: None)

    payloads = {
        "rs1": {"name": "rs1"},
        "rs2": None,
    }

    class FakeClient(EnsemblClient):
        def get_variation(self, variant_id: str, *, phenotypes: bool = False) -> dict[str, Any] | None:
            return payloads[variant_id]

    results = lookup_variants_cached(
        ["rs1", "rs2", "rs1", ""],
        cache_path=tmp_path / "cache",
        backend="json",
        client=FakeClient(min_interval_sec=0.0),
    )
    assert set(results) == {"rs1", "rs2"}
    assert results["rs1"] == {"name": "rs1"}
    assert results["rs2"] is None

    # Second pass should be fully offline.
    results2 = lookup_variants_cached(
        ["rs1", "rs2"],
        cache_path=tmp_path / "cache",
        backend="json",
        client=FakeClient(min_interval_sec=0.0),
    )
    assert results2 == results


def test_grch38_locus_from_variation() -> None:
    payload = {
        "mappings": [
            {
                "assembly_name": "GRCh37",
                "coord_system": "chromosome",
                "seq_region_name": "19",
                "start": 1,
            },
            {
                "assembly_name": "GRCh38",
                "coord_system": "chromosome",
                "seq_region_name": "19",
                "start": 44908684,
            },
            {
                "assembly_name": "GRCh38",
                "coord_system": "chromosome",
                "seq_region_name": "19_PATCH",
                "start": 99,
            },
        ]
    }
    assert grch38_locus_from_variation(payload) == ("19", 44908684)
    assert grch38_locus_from_variation(None) is None
