#!/usr/bin/env python3
"""Map locked QTD dataset IDs to QTS study folders via eQTL Catalogue HTTPS listings.

Reads ``analysis/eqtl/eqtl_datasets_used.tsv`` (does not modify it). Lists
``https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/`` then each QTS
folder by parsing HTML directory indexes, at most one request every 2 seconds.
Builds QTD→QTS only from those listings (no ID-number inference). Writes
``analysis/eqtl/eqtl_datasets_used_with_qts.tsv`` plus provenance.

Offline modes (no network):
  --from-saved-html DIR   read local HTML files from DIR
  --dry-run               use the synthetic fixture under tests/fixtures/

Source: https://www.ebi.ac.uk/eqtl/Data_access/

Example:
    uv run python scripts/eqtl/map_qts_to_qtd.py --dry-run
    uv run python scripts/eqtl/map_qts_to_qtd.py   # live; run in a terminal that can reach EBI
"""

from __future__ import annotations

import hashlib
import json
import re
import sys
import time
import urllib.error
import urllib.request
from datetime import UTC, datetime
from html.parser import HTMLParser
from importlib import metadata
from pathlib import Path
from typing import Any

import polars as pl
import typer

from rogen_aging.config import find_repo_root

REPO_ROOT = find_repo_root()

DEFAULT_DATASETS = REPO_ROOT / "analysis" / "eqtl" / "eqtl_datasets_used.tsv"
DEFAULT_OUT = REPO_ROOT / "analysis" / "eqtl" / "eqtl_datasets_used_with_qts.tsv"
DEFAULT_PROVENANCE = REPO_ROOT / "analysis" / "eqtl" / "eqtl_datasets_used_with_qts_provenance.json"
DEFAULT_SYNTHETIC_HTML_DIR = REPO_ROOT / "tests" / "fixtures" / "eqtl_catalogue_synthetic" / "html"

HTTPS_SUMSTATS_URL = "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/"
ACCESS_DOCS_URL = "https://www.ebi.ac.uk/eqtl/Data_access/"
MIN_REQUEST_INTERVAL_S = 2.0
EXPECTED_N = 14
FILE_ALL = ".all.tsv.gz"
FILE_CC = ".cc.tsv.gz"

_HREF_RE = re.compile(r"""href=["']([^"']+)["']""", re.IGNORECASE)


class _HrefParser(HTMLParser):
    """Collect href attribute values from an HTML directory listing."""

    def __init__(self) -> None:
        super().__init__()
        self.hrefs: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag.lower() != "a":
            return
        for key, value in attrs:
            if key.lower() == "href" and value:
                self.hrefs.append(value)


def _stop(message: str) -> None:
    typer.echo(message, err=True)
    raise typer.Exit(code=1)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _git_commit() -> str:
    import subprocess

    completed = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=REPO_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    return completed.stdout.strip()


def _package_version(name: str) -> str:
    return metadata.version(name)


def _repo_relative(path: Path) -> str:
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(REPO_ROOT))
    except ValueError:
        return str(resolved)


def parse_directory_hrefs(html_text: str) -> list[str]:
    """Return cleaned basenames from hrefs in an Apache-style directory index."""
    parser = _HrefParser()
    parser.feed(html_text)
    hrefs = list(parser.hrefs)
    if not hrefs:
        hrefs = _HREF_RE.findall(html_text)
    names: list[str] = []
    for href in hrefs:
        item = href.strip()
        if not item or item.startswith("?") or item.startswith("#"):
            continue
        if item in (".", "..", "./", "../"):
            continue
        base = item.rstrip("/").split("/")[-1]
        if not base or base in (".", ".."):
            continue
        names.append(base)
    return names


def load_locked_datasets(path: Path) -> pl.DataFrame:
    assert path.is_file(), f"MISSING INPUT: {path}"
    frame = pl.read_csv(path, separator="\t", infer_schema_length=0)
    assert "dataset_id" in frame.columns, f"missing dataset_id column: {path}"
    assert (
        frame.height == EXPECTED_N
    ), f"DATASET COUNT DISAGREES: computed={frame.height} ground_truth={EXPECTED_N}"
    ids = [str(value) for value in frame.get_column("dataset_id").to_list()]
    assert len(set(ids)) == EXPECTED_N, f"dataset_id not unique: {ids}"
    for dataset_id in ids:
        assert dataset_id.startswith("QTD"), f"expected QTD id, got {dataset_id}"
    return frame


class ListingSource:
    """Abstract directory listing (live HTTPS or local HTML)."""

    def list_names(self, url_or_key: str) -> list[str]:
        raise NotImplementedError


class RateLimitedHttpsSource(ListingSource):
    """Fetch HTTPS directory indexes with a 2-second minimum interval."""

    def __init__(self, *, timeout: float = 120.0) -> None:
        self.timeout = timeout
        self._last_request_monotonic = 0.0
        self.fetched_urls: list[str] = []

    def _throttle(self) -> None:
        now = time.monotonic()
        wait = MIN_REQUEST_INTERVAL_S - (now - self._last_request_monotonic)
        if wait > 0:
            time.sleep(wait)
        self._last_request_monotonic = time.monotonic()

    def fetch_text(self, url: str) -> str:
        self._throttle()
        self.fetched_urls.append(url)
        request = urllib.request.Request(
            url,
            headers={"User-Agent": "rogen-aging-eqtl-mapper/1.0"},
        )
        try:
            with urllib.request.urlopen(request, timeout=self.timeout) as response:
                raw = response.read()
        except urllib.error.URLError as exc:
            _stop(
                "HTTPS LISTING FAILED. Stopping. Do not try another source.\n"
                f"  url={url}\n"
                f"  error={type(exc).__name__}: {exc}"
            )
        except TimeoutError as exc:
            _stop(
                "HTTPS LISTING FAILED. Stopping. Do not try another source.\n"
                f"  url={url}\n"
                f"  error={type(exc).__name__}: {exc}"
            )
        return raw.decode("utf-8", errors="replace")

    def list_names(self, url_or_key: str) -> list[str]:
        return parse_directory_hrefs(self.fetch_text(url_or_key))


class SavedHtmlSource(ListingSource):
    """Read local HTML directory listings (synthetic or saved from EBI)."""

    def __init__(self, html_dir: Path) -> None:
        assert html_dir.is_dir(), f"MISSING INPUT: {html_dir}"
        self.html_dir = html_dir
        self.reads: list[str] = []

    def _read(self, filename: str) -> str:
        path = self.html_dir / filename
        assert path.is_file(), f"MISSING HTML FIXTURE: {path}"
        self.reads.append(str(path))
        return path.read_text(encoding="utf-8")

    def list_sumstats(self) -> list[str]:
        return parse_directory_hrefs(self._read("sumstats_index.html"))

    def list_qts(self, qts_id: str) -> list[str]:
        return parse_directory_hrefs(self._read(f"{qts_id}_index.html"))

    def list_qtd_files(self, qtd_id: str) -> list[str]:
        return parse_directory_hrefs(self._read(f"{qtd_id}_files.html"))

    def list_names(self, url_or_key: str) -> list[str]:
        # unused in saved mode; kept for interface parity
        raise AssertionError(f"SavedHtmlSource.list_names not used for key={url_or_key}")


def choose_file_type_from_names(qtd_id: str, names: list[str]) -> str:
    basenames = set(names)
    all_name = f"{qtd_id}{FILE_ALL}"
    cc_name = f"{qtd_id}{FILE_CC}"
    if all_name in basenames:
        return FILE_ALL
    if cc_name in basenames:
        return FILE_CC
    assert False, f"Neither {all_name} nor {cc_name} in listing; names={sorted(basenames)}"


def map_datasets_https(
    locked: pl.DataFrame,
    source: RateLimitedHttpsSource,
) -> tuple[pl.DataFrame, dict[str, Any]]:
    wanted = [str(value) for value in locked.get_column("dataset_id").to_list()]
    wanted_set = set(wanted)

    top_names = source.list_names(HTTPS_SUMSTATS_URL)
    qts_ids = sorted(name for name in top_names if name.startswith("QTS"))
    assert qts_ids, f"NO QTS FOLDERS under {HTTPS_SUMSTATS_URL}"
    typer.echo(f"Listed {len(qts_ids)} QTS folders under {HTTPS_SUMSTATS_URL}")

    found: dict[str, list[str]] = {qtd: [] for qtd in wanted}
    listing_log: list[dict[str, Any]] = []
    for qts_id in qts_ids:
        url = f"{HTTPS_SUMSTATS_URL}{qts_id}/"
        qtd_ids = sorted(name for name in source.list_names(url) if name.startswith("QTD"))
        listing_log.append({"qts_id": qts_id, "n_qtd": len(qtd_ids)})
        typer.echo(f"  {qts_id}: {len(qtd_ids)} QTD folders")
        for qtd_id in qtd_ids:
            if qtd_id in wanted_set:
                found[qtd_id].append(qts_id)

    missing = [qtd for qtd in wanted if not found[qtd]]
    if missing:
        _stop(
            "QTD NOT FOUND IN HTTPS LISTINGS. Stopping.\n"
            + "\n".join(f"  - {qtd}" for qtd in missing)
        )
    ambiguous = {qtd: qts_list for qtd, qts_list in found.items() if len(qts_list) != 1}
    assert not ambiguous, f"QTD mapped to multiple QTS: {ambiguous}"

    rows: list[dict[str, str]] = []
    for qtd_id in wanted:
        qts_id = found[qtd_id][0]
        file_url_dir = f"{HTTPS_SUMSTATS_URL}{qts_id}/{qtd_id}/"
        names = source.list_names(file_url_dir)
        file_type = choose_file_type_from_names(qtd_id, names)
        rows.append(
            {
                "dataset_id": qtd_id,
                "qts_id": qts_id,
                "file_type": file_type,
                "file_url": f"{file_url_dir}{qtd_id}{file_type}",
            }
        )
    return _finalize_mapped(locked, rows, qts_ids, listing_log)


def map_datasets_saved_html(
    locked: pl.DataFrame,
    source: SavedHtmlSource,
) -> tuple[pl.DataFrame, dict[str, Any]]:
    wanted = [str(value) for value in locked.get_column("dataset_id").to_list()]
    wanted_set = set(wanted)

    top_names = source.list_sumstats()
    qts_ids = sorted(name for name in top_names if name.startswith("QTS"))
    assert qts_ids, f"NO QTS FOLDERS in {source.html_dir}/sumstats_index.html"

    found: dict[str, list[str]] = {qtd: [] for qtd in wanted}
    listing_log: list[dict[str, Any]] = []
    for qts_id in qts_ids:
        qtd_ids = sorted(name for name in source.list_qts(qts_id) if name.startswith("QTD"))
        listing_log.append({"qts_id": qts_id, "n_qtd": len(qtd_ids)})
        typer.echo(f"  {qts_id}: {len(qtd_ids)} QTD folders")
        for qtd_id in qtd_ids:
            if qtd_id in wanted_set:
                found[qtd_id].append(qts_id)

    missing = [qtd for qtd in wanted if not found[qtd]]
    if missing:
        _stop(
            "QTD NOT FOUND IN SAVED HTML LISTINGS. Stopping.\n"
            + "\n".join(f"  - {qtd}" for qtd in missing)
        )
    ambiguous = {qtd: qts_list for qtd, qts_list in found.items() if len(qts_list) != 1}
    assert not ambiguous, f"QTD mapped to multiple QTS: {ambiguous}"

    rows: list[dict[str, str]] = []
    for qtd_id in wanted:
        qts_id = found[qtd_id][0]
        names = source.list_qtd_files(qtd_id)
        file_type = choose_file_type_from_names(qtd_id, names)
        rows.append(
            {
                "dataset_id": qtd_id,
                "qts_id": qts_id,
                "file_type": file_type,
                "file_url": (f"{HTTPS_SUMSTATS_URL}{qts_id}/{qtd_id}/{qtd_id}{file_type}"),
            }
        )
    return _finalize_mapped(locked, rows, qts_ids, listing_log)


def _finalize_mapped(
    locked: pl.DataFrame,
    rows: list[dict[str, str]],
    qts_ids: list[str],
    listing_log: list[dict[str, Any]],
) -> tuple[pl.DataFrame, dict[str, Any]]:
    mapped = pl.DataFrame(rows)
    assert (
        mapped.height == EXPECTED_N
    ), f"MAPPED COUNT DISAGREES: computed={mapped.height} ground_truth={EXPECTED_N}"
    out = locked.join(mapped, on="dataset_id", how="left")
    assert out.filter(pl.col("qts_id").is_null()).height == 0
    assert out.height == EXPECTED_N
    meta = {
        "n_qts_listed": len(qts_ids),
        "qts_ids_listed": qts_ids,
        "per_qts_qtd_counts": listing_log,
        "n_mapped": int(out.height),
    }
    return out, meta


def write_outputs(
    *,
    locked: pl.DataFrame,
    table: pl.DataFrame,
    listing_meta: dict[str, Any],
    datasets: Path,
    out: Path,
    provenance_out: Path,
    mode: str,
    started: str,
    extra_inputs: dict[str, Any] | None = None,
) -> None:
    original_cols = list(locked.columns)
    write_cols = original_cols + ["qts_id", "file_type"]
    assert set(write_cols).issubset(set(table.columns))
    out.parent.mkdir(parents=True, exist_ok=True)
    table.select(write_cols).write_csv(out, separator="\t")

    provenance: dict[str, Any] = {
        "generated_at_utc": started,
        "access_date_utc": started[:10],
        "access_docs_url": ACCESS_DOCS_URL,
        "https_sumstats_url": HTTPS_SUMSTATS_URL,
        "mode": mode,
        "rate_limit_seconds": MIN_REQUEST_INTERVAL_S,
        "catalogue_release_or_version": "not provided by source",
        "git_commit": _git_commit(),
        "tool_versions": {
            "python": sys.version.split()[0],
            "polars": _package_version("polars"),
            "typer": _package_version("typer"),
        },
        "inputs": {
            "eqtl_datasets_used": {
                "path": _repo_relative(datasets),
                "sha256": _sha256(datasets),
                "n_rows": locked.height,
                "modified": False,
            }
        },
        "outputs": {
            "table": {
                "path": _repo_relative(out),
                "sha256": _sha256(out),
                "n_rows": table.height,
                "columns": write_cols,
            }
        },
        "listing": listing_meta,
        "mapped_rows": table.select(write_cols + ["file_url"]).to_dicts(),
    }
    if extra_inputs:
        provenance["inputs"].update(extra_inputs)
    provenance_out.parent.mkdir(parents=True, exist_ok=True)
    provenance_out.write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    typer.echo(f"Wrote {out}")
    typer.echo(f"Wrote {provenance_out}")
    typer.echo(f"OK: mapped {table.height}/{EXPECTED_N} QTD→QTS ({mode})")


def main(
    datasets: Path = typer.Option(DEFAULT_DATASETS, "--datasets"),
    out: Path = typer.Option(DEFAULT_OUT, "--out"),
    provenance_out: Path = typer.Option(DEFAULT_PROVENANCE, "--provenance-out"),
    from_saved_html: Path | None = typer.Option(
        None,
        "--from-saved-html",
        help="Directory of local HTML indexes (no network).",
    ),
    dry_run: bool = typer.Option(
        False,
        "--dry-run",
        help="Use synthetic HTML fixtures under tests/fixtures/ (no network).",
    ),
) -> None:
    """Map locked QTDs to QTS IDs from HTTPS HTML listings or local fixtures."""
    started = datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")
    assert (
        datasets.resolve() != out.resolve()
    ), "REFUSING TO OVERWRITE locked datasets TSV; --out must differ"
    if dry_run and from_saved_html is not None:
        _stop("Pass only one of --dry-run or --from-saved-html.")

    locked = load_locked_datasets(datasets)
    typer.echo(f"Locked datasets: {datasets} n={locked.height}")

    if dry_run:
        html_dir = DEFAULT_SYNTHETIC_HTML_DIR
        assert html_dir.is_dir(), f"MISSING SYNTHETIC FIXTURE DIR: {html_dir}"
        source = SavedHtmlSource(html_dir)
        table, listing_meta = map_datasets_saved_html(locked, source)
        write_outputs(
            locked=locked,
            table=table,
            listing_meta=listing_meta,
            datasets=datasets,
            out=out,
            provenance_out=provenance_out,
            mode="dry-run-synthetic",
            started=started,
            extra_inputs={
                "synthetic_html_dir": {
                    "path": _repo_relative(html_dir),
                    "label": "synthetic",
                    "files_read": source.reads,
                }
            },
        )
        return

    if from_saved_html is not None:
        source = SavedHtmlSource(from_saved_html)
        table, listing_meta = map_datasets_saved_html(locked, source)
        write_outputs(
            locked=locked,
            table=table,
            listing_meta=listing_meta,
            datasets=datasets,
            out=out,
            provenance_out=provenance_out,
            mode="from-saved-html",
            started=started,
            extra_inputs={
                "saved_html_dir": {
                    "path": str(from_saved_html),
                    "files_read": source.reads,
                }
            },
        )
        return

    # Live HTTPS mode — caller must run this in a terminal that can reach EBI.
    live = RateLimitedHttpsSource()
    table, listing_meta = map_datasets_https(locked, live)
    write_outputs(
        locked=locked,
        table=table,
        listing_meta=listing_meta,
        datasets=datasets,
        out=out,
        provenance_out=provenance_out,
        mode="live-https",
        started=started,
        extra_inputs={"fetched_urls": live.fetched_urls},
    )


if __name__ == "__main__":
    typer.run(main)
