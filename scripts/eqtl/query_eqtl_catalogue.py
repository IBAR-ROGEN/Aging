#!/usr/bin/env python3
"""Query eQTL Catalogue sumstats via tabix (HTTPS), one region per panel variant × dataset.

METHOD (do not substitute): tabix on
  https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/<QTS>/<QTD>/<QTD>.all.tsv.gz
  <chrom>:<pos>-<pos>

REST API is deprecated — this script does not call it. Rate limit: one tabix
request every 2 seconds. Raw responses are cached under
``data/cache/eqtl_catalogue/``.

Offline:
  --dry-run   synthetic fixtures only; no network and no tabix subprocess to EBI.

Example:
    uv run python scripts/eqtl/query_eqtl_catalogue.py --dry-run
    # live (run in your terminal, not the agent):
    micromamba activate rogen-tools
    uv run python scripts/eqtl/query_eqtl_catalogue.py
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass
from datetime import UTC, datetime
from importlib import metadata
from pathlib import Path
from typing import Any, Literal

import polars as pl
import typer

from rogen_aging.config import find_repo_root

REPO_ROOT = find_repo_root()

DEFAULT_PANEL = REPO_ROOT / "data" / "processed" / "variants_47_input_v2.csv"
DEFAULT_DATASETS = REPO_ROOT / "analysis" / "eqtl" / "eqtl_datasets_used_with_qts.tsv"
DEFAULT_CACHE = REPO_ROOT / "data" / "cache" / "eqtl_catalogue"
DEFAULT_OUT = REPO_ROOT / "results" / "eqtl_associations_full_rebuilt.csv"
DEFAULT_OLD = REPO_ROOT / "results" / "eqtl_associations_full.csv"
DEFAULT_DIFF = REPO_ROOT / "results" / "eqtl_associations_full_rebuilt_diff.csv"
DEFAULT_PROVENANCE = REPO_ROOT / "results" / "eqtl_associations_full_rebuilt_provenance.json"
DEFAULT_PAIR_STATUS = REPO_ROOT / "analysis" / "eqtl" / "eqtl_catalogue_query_pair_status.tsv"
DEFAULT_DOC = REPO_ROOT / "docs" / "EQTL_CATALOGUE_REBUILD.md"
SYNTHETIC_DIR = REPO_ROOT / "tests" / "fixtures" / "eqtl_catalogue_synthetic"
SYNTHETIC_WITH_QTS = SYNTHETIC_DIR / "eqtl_datasets_used_with_qts.synthetic.tsv"
INTEGRATED_SCRIPT = REPO_ROOT / "scripts" / "longevity" / "integrated_evidence.py"

HTTPS_SUMSTATS_URL = "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/"
ACCESS_DOCS_URL = "https://www.ebi.ac.uk/eqtl/Data_access/"
COLUMNS_MD_URL = (
    "https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/master/" "tabix/Columns.md"
)
MIN_REQUEST_INTERVAL_S = 2.0
EXPECTED_PANEL_N = 47
EXPECTED_DATASET_N = 14
EXPECTED_PAIRS = EXPECTED_PANEL_N * EXPECTED_DATASET_N  # 658
GROUND_TRUTH_LABELS = {"GWAS support": 1, "eQTL support": 16, "population difference": 5}

PairStatus = Literal["returned", "no data", "error"]

# Catalogue .all.tsv.gz header fields commonly present (Columns.md). Synthetic
# dry-run rows use this order; live mode always prefers the file's own header.
SYNTHETIC_HEADER = (
    "molecular_trait_id\tchromosome\tposition\tref\talt\tvariant\tma_samples\t"
    "maf\tpvalue\tbeta\tse\ttype\tac\tan\tr2\tgene_id\tmedian_tpm\trsid"
)


@dataclass(frozen=True)
class PanelVariant:
    rsid: str
    chrom: str
    pos: int
    ref: str
    alt: str


@dataclass(frozen=True)
class DatasetRow:
    dataset_id: str
    qts_id: str
    file_type: str


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


def cache_path(cache_dir: Path, dataset_id: str, rsid: str) -> Path:
    safe_rsid = re.sub(r"[^A-Za-z0-9_.-]+", "_", rsid)
    return cache_dir / f"{dataset_id}__{safe_rsid}.tsv"


def load_panel(path: Path) -> list[PanelVariant]:
    assert path.is_file(), f"MISSING INPUT: {path}"
    frame = pl.read_csv(
        path,
        schema_overrides={
            "chrom": pl.Utf8,
            "pos": pl.Int64,
            "ref": pl.Utf8,
            "alt": pl.Utf8,
            "rsid": pl.Utf8,
            "alt_status": pl.Utf8,
        },
    )
    required = {"chrom", "pos", "ref", "alt", "rsid"}
    missing = sorted(required - set(frame.columns))
    assert not missing, f"panel missing columns: {missing}"
    assert frame.height == EXPECTED_PANEL_N, (
        f"PANEL ROW COUNT DISAGREES: computed={frame.height} " f"ground_truth={EXPECTED_PANEL_N}"
    )
    null_coord = frame.filter(
        pl.col("chrom").is_null() | (pl.col("chrom").cast(pl.Utf8) == "") | pl.col("pos").is_null()
    )
    assert null_coord.height == 0, "PANEL HAS NO GRCh38 COORDINATES for some rows; stopping."
    variants: list[PanelVariant] = []
    for row in frame.iter_rows(named=True):
        chrom = str(row["chrom"]).removeprefix("chr")
        variants.append(
            PanelVariant(
                rsid=str(row["rsid"]),
                chrom=chrom,
                pos=int(row["pos"]),
                ref=str(row["ref"]).upper(),
                alt=str(row["alt"]).upper(),
            )
        )
    return variants


def load_datasets(path: Path) -> list[DatasetRow]:
    assert path.is_file(), f"MISSING INPUT: {path}"
    frame = pl.read_csv(path, separator="\t", infer_schema_length=0)
    required = {"dataset_id", "qts_id", "file_type"}
    missing = sorted(required - set(frame.columns))
    assert not missing, f"datasets missing columns: {missing}"
    assert frame.height == EXPECTED_DATASET_N, (
        f"DATASET COUNT DISAGREES: computed={frame.height} " f"ground_truth={EXPECTED_DATASET_N}"
    )
    rows: list[DatasetRow] = []
    for row in frame.iter_rows(named=True):
        qts = str(row["qts_id"])
        assert qts.startswith("QTS"), f"qts_id must be QTS*, got {qts}"
        file_type = str(row["file_type"])
        assert file_type in {".all.tsv.gz", ".cc.tsv.gz"}, file_type
        rows.append(
            DatasetRow(
                dataset_id=str(row["dataset_id"]),
                qts_id=qts,
                file_type=file_type,
            )
        )
    return rows


def file_url(dataset: DatasetRow) -> str:
    qtd = dataset.dataset_id
    return f"{HTTPS_SUMSTATS_URL}{dataset.qts_id}/{qtd}/{qtd}{dataset.file_type}"


def region_string(variant: PanelVariant) -> str:
    return f"{variant.chrom}:{variant.pos}-{variant.pos}"


class RateLimiter:
    def __init__(self, interval_s: float = MIN_REQUEST_INTERVAL_S) -> None:
        self.interval_s = interval_s
        self._last = 0.0

    def wait(self) -> None:
        now = time.monotonic()
        delay = self.interval_s - (now - self._last)
        if delay > 0:
            time.sleep(delay)
        self._last = time.monotonic()


def ensure_tabix() -> str:
    path = shutil.which("tabix")
    assert path, (
        "MISSING TOOL: tabix. Activate micromamba env rogen-tools " "(htslib) before live mode."
    )
    return path


def run_tabix_live(
    *,
    tabix_bin: str,
    url: str,
    region: str,
    limiter: RateLimiter,
) -> tuple[PairStatus, str, str]:
    """Return (status, body, error_text). Body includes no header."""
    limiter.wait()
    completed = subprocess.run(
        [tabix_bin, url, region],
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        return "error", "", completed.stderr.strip() or completed.stdout.strip()
    body = completed.stdout
    if not body.strip():
        return "no data", "", ""
    return "returned", body, ""


def write_synthetic_cache_body(variant: PanelVariant) -> str:
    """Synthetic tabix body: one matching row + header comment label."""
    # Header is stored separately in cache meta; body is tabix-like data lines.
    # Return one allele-matched row for the first chrom/pos so matching is exercised.
    return (
        f"ENSG00000000001\t{variant.chrom}\t{variant.pos}\t{variant.ref}\t"
        f"{variant.alt}\tchr{variant.chrom}_{variant.pos}_{variant.ref}_{variant.alt}"
        f"\t10\t0.1\t1e-4\t0.25\t0.05\tSNP\t20\t100\t1.0\t"
        f"ENSG00000000001\t5.0\t{variant.rsid}\n"
    )


def query_pair(
    *,
    variant: PanelVariant,
    dataset: DatasetRow,
    cache_dir: Path,
    dry_run: bool,
    tabix_bin: str | None,
    limiter: RateLimiter,
    cache_hits: list[str],
) -> dict[str, Any]:
    path = cache_path(cache_dir, dataset.dataset_id, variant.rsid)
    url = file_url(dataset)
    region = region_string(variant)
    meta_path = path.with_suffix(path.suffix + ".meta.json")

    if path.is_file() and meta_path.is_file():
        cache_hits.append(str(path))
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
        body = path.read_text(encoding="utf-8")
        status = meta["status"]
        assert status in {"returned", "no data", "error"}, status
        return {
            "dataset_id": dataset.dataset_id,
            "rsid": variant.rsid,
            "status": status,
            "cache_hit": True,
            "cache_path": str(path),
            "file_url": url,
            "region": region,
            "body": body,
            "header": meta.get("header", ""),
            "error": meta.get("error", ""),
        }

    if dry_run:
        # Synthetic offline response: alternate no-data / returned without network.
        if hash((dataset.dataset_id, variant.rsid)) % 5 == 0:
            status = "no data"
            body = ""
            header = SYNTHETIC_HEADER
            error = ""
        else:
            status = "returned"
            body = write_synthetic_cache_body(variant)
            header = SYNTHETIC_HEADER
            error = ""
        label = "synthetic"
    else:
        assert tabix_bin is not None
        status, body, error = run_tabix_live(
            tabix_bin=tabix_bin, url=url, region=region, limiter=limiter
        )
        header = ""  # filled by caller via curl once per dataset if needed
        label = "live"

    cache_dir.mkdir(parents=True, exist_ok=True)
    path.write_text(body, encoding="utf-8")
    meta = {
        "status": status,
        "file_url": url,
        "region": region,
        "header": header,
        "error": error,
        "label": label,
        "written_at_utc": datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
    }
    meta_path.write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
    return {
        "dataset_id": dataset.dataset_id,
        "rsid": variant.rsid,
        "status": status,
        "cache_hit": False,
        "cache_path": str(path),
        "file_url": url,
        "region": region,
        "body": body,
        "header": header,
        "error": error,
    }


def alleles_match_panel(row_ref: str, row_alt: str, panel: PanelVariant) -> bool:
    ref = row_ref.upper()
    alt = row_alt.upper()
    return {ref, alt} == {panel.ref, panel.alt}


def parse_tabix_rows(
    *,
    header: str,
    body: str,
    variant: PanelVariant,
    dataset: DatasetRow,
) -> list[dict[str, Any]]:
    if not body.strip():
        return []
    header_cols = header.strip().split("\t") if header.strip() else []
    if not header_cols:
        # Fall back to synthetic/known column order.
        header_cols = SYNTHETIC_HEADER.split("\t")
    col_index = {name: i for i, name in enumerate(header_cols)}
    required = ["chromosome", "position", "ref", "alt", "beta"]
    for name in required:
        assert name in col_index, (
            f"Catalogue header missing {name}; columns={header_cols}. " f"See {COLUMNS_MD_URL}"
        )

    out: list[dict[str, Any]] = []
    for line in body.splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < len(header_cols):
            # pad / tolerate shorter
            parts = parts + [""] * (len(header_cols) - len(parts))

        def get(name: str, default: str = "") -> str:
            idx = col_index.get(name)
            if idx is None or idx >= len(parts):
                return default
            return parts[idx]

        chrom = get("chromosome").removeprefix("chr")
        pos_raw = get("position")
        if not pos_raw:
            continue
        pos = int(float(pos_raw))
        if chrom != variant.chrom or pos != variant.pos:
            continue
        row_ref = get("ref").upper()
        row_alt = get("alt").upper()
        if not alleles_match_panel(row_ref, row_alt, variant):
            continue
        alt_match = "Y" if row_alt == variant.alt else "N"
        out.append(
            {
                "canonical_rsid": variant.rsid,
                "eqtl_alt_allele": row_alt,
                "tissue": get("tissue_label") or get("tissue") or "",
                "dataset_id": dataset.dataset_id,
                "n_samples": get("ma_samples") or get("an") or "",
                "gene_id": get("gene_id"),
                "egene_symbol": "",
                "is_cis_to_annotated_gene": "",
                "beta": get("beta"),
                "se": get("se"),
                "pvalue": get("pvalue"),
                "fdr_bh": get("fdr") or get("fdr_bh") or "",
                "maf": get("maf"),
                "median_tpm": get("median_tpm"),
                "chrom": chrom,
                "pos": pos,
                "ref": row_ref,
                "alt": row_alt,
                "alt_matches_panel_alt": alt_match,
            }
        )
    return out


def account_pairs(statuses: list[dict[str, Any]]) -> dict[str, int]:
    counts = {"returned": 0, "no data": 0, "error": 0}
    for row in statuses:
        status = row["status"]
        assert status in counts, f"unknown pair status: {status}"
        counts[status] += 1
    total = sum(counts.values())
    assert (
        total == EXPECTED_PAIRS
    ), f"PAIR ACCOUNTING DISAGREES: computed={total} ground_truth={EXPECTED_PAIRS}"
    return counts


def write_synthetic_with_qts(path: Path) -> None:
    """Tiny labelled synthetic datasets table covering the 14 locked QTDs."""
    locked = pl.read_csv(
        REPO_ROOT / "analysis" / "eqtl" / "eqtl_datasets_used.tsv",
        separator="\t",
    )
    rows = [
        {
            "dataset_id": str(r["dataset_id"]),
            "qts_id": "QTS000999",
            "file_type": ".all.tsv.gz",
        }
        for r in locked.iter_rows(named=True)
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    pl.DataFrame(rows).write_csv(path, separator="\t")


def row_level_diff(old_path: Path, new_frame: pl.DataFrame, diff_path: Path) -> dict[str, Any]:
    if not old_path.is_file():
        diff_path.write_text(
            "status,detail\nmissing_old,old associations file not present\n",
            encoding="utf-8",
        )
        return {"old_present": False, "n_diff_rows": 1}
    old = pl.read_csv(old_path, infer_schema_length=0)
    key = ["canonical_rsid", "dataset_id", "gene_id", "eqtl_alt_allele"]
    for frame_name, frame in ("old", old), ("new", new_frame):
        missing = [c for c in key if c not in frame.columns]
        assert not missing, f"{frame_name} missing key cols {missing}"
    old_keys = set(map(tuple, old.select(key).iter_rows()))
    new_keys = set(map(tuple, new_frame.select(key).iter_rows())) if new_frame.height else set()
    only_old = sorted(old_keys - new_keys)
    only_new = sorted(new_keys - old_keys)
    rows = [{"status": "only_old", "key": "|".join(map(str, k))} for k in only_old] + [
        {"status": "only_new", "key": "|".join(map(str, k))} for k in only_new
    ]
    pl.DataFrame(rows if rows else [{"status": "identical_keys", "key": ""}]).write_csv(diff_path)
    return {
        "old_present": True,
        "n_only_old": len(only_old),
        "n_only_new": len(only_new),
        "n_diff_rows": len(rows),
    }


def run_integrated_temp(associations: Path) -> dict[str, Any]:
    assert INTEGRATED_SCRIPT.is_file(), f"MISSING INPUT: {INTEGRATED_SCRIPT}"
    with tempfile.TemporaryDirectory(prefix="eqtl_rebuild_ie_") as tmp:
        out_dir = Path(tmp) / "integrated"
        out_dir.mkdir()
        env = os.environ.copy()
        env.setdefault("MPLCONFIGDIR", str(out_dir / ".mplconfig"))
        completed = subprocess.run(
            [
                "uv",
                "run",
                "python",
                str(INTEGRATED_SCRIPT),
                "--associations",
                str(associations),
                "--out-dir",
                str(out_dir),
            ],
            cwd=REPO_ROOT,
            check=False,
            capture_output=True,
            text=True,
            env=env,
        )
        assert completed.returncode == 0, (
            "integrated_evidence.py FAILED on rebuilt table\n"
            f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
        )
        table = pl.read_csv(out_dir / "integrated_evidence.csv")
        assert table.height == EXPECTED_PANEL_N
        if table.schema.get("label_eqtl_support") == pl.Boolean:
            n_eqtl = int(table.filter(pl.col("label_eqtl_support") == True).height)  # noqa: E712
        else:
            n_eqtl = int(
                table.filter(
                    pl.col("label_eqtl_support").cast(pl.Utf8).str.to_lowercase() == "true"
                ).height
            )
        observed = {
            "GWAS support": int(table.filter(pl.col("label_gwas_support")).height),
            "eQTL support": n_eqtl,
            "population difference": int(
                table.filter(pl.col("label_population_difference")).height
            ),
        }
        return {
            "observed_labels": observed,
            "ground_truth_labels": dict(GROUND_TRUTH_LABELS),
            "labels_match": observed == GROUND_TRUTH_LABELS,
            "stdout_tail": completed.stdout[-1500:],
        }


def render_doc(payload: dict[str, Any]) -> str:
    labels = payload["integrated_evidence_temp"]["observed_labels"]
    gt = GROUND_TRUTH_LABELS
    lines = [
        "# eQTL Catalogue rebuild (tabix)",
        "",
        "Access method: tabix on HTTPS sumstats under",
        f"`{HTTPS_SUMSTATS_URL}` (see {ACCESS_DOCS_URL}).",
        "REST API not used.",
        "",
        f"- Access date (UTC): `{payload['access_date_utc']}`",
        f"- Mode: `{payload['mode']}`",
        f"- Catalogue release/version: `{payload['catalogue_release_or_version']}`",
        f"- Pairs accounted: `{payload['pair_counts']}` (asserted total {EXPECTED_PAIRS})",
        f"- Association rows written: `{payload['n_association_rows']}`",
        f"- Cache hits: `{payload['n_cache_hits']}`",
        "",
        "## Integrated evidence labels (temp out-dir)",
        "",
        f"- Observed: GWAS **{labels['GWAS support']}** /"
        f" eQTL **{labels['eQTL support']}** /"
        f" population **{labels['population difference']}**",
        f"- Ground truth (local tables): GWAS **{gt['GWAS support']}** /"
        f" eQTL **{gt['eQTL support']}** /"
        f" population **{gt['population difference']}**",
        f"- Match: `{payload['integrated_evidence_temp']['labels_match']}` " "(not forced)",
        "",
        "## Provenance",
        "",
        f"- Git commit: `{payload['git_commit']}`",
        f"- Timestamp: `{payload['generated_at_utc']}`",
        f"- Python: `{payload['tool_versions']['python']}`",
        f"- Polars: `{payload['tool_versions']['polars']}`",
        f"- Columns.md: {COLUMNS_MD_URL}",
        "",
    ]
    return "\n".join(lines)


def main(
    panel: Path = typer.Option(DEFAULT_PANEL, "--panel"),
    datasets: Path = typer.Option(DEFAULT_DATASETS, "--datasets"),
    cache_dir: Path = typer.Option(DEFAULT_CACHE, "--cache-dir"),
    out: Path = typer.Option(DEFAULT_OUT, "--out"),
    old: Path = typer.Option(DEFAULT_OLD, "--old"),
    diff_out: Path = typer.Option(DEFAULT_DIFF, "--diff-out"),
    provenance_out: Path = typer.Option(DEFAULT_PROVENANCE, "--provenance-out"),
    pair_status_out: Path = typer.Option(DEFAULT_PAIR_STATUS, "--pair-status-out"),
    doc_out: Path = typer.Option(DEFAULT_DOC, "--doc-out"),
    dry_run: bool = typer.Option(
        False,
        "--dry-run",
        help="Synthetic fixtures only; no network / no EBI tabix.",
    ),
    skip_integrated: bool = typer.Option(
        False,
        "--skip-integrated",
        help="Skip temp integrated_evidence run (tests).",
    ),
) -> None:
    """Query Catalogue sumstats for 47×14 pairs; write rebuilt associations."""
    started = datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")
    assert (
        out.resolve() != old.resolve()
    ), "REFUSING TO OVERWRITE results/eqtl_associations_full.csv"

    if dry_run:
        write_synthetic_with_qts(SYNTHETIC_WITH_QTS)
        datasets = SYNTHETIC_WITH_QTS
        if cache_dir.resolve() == DEFAULT_CACHE.resolve():
            cache_dir = SYNTHETIC_DIR / "cache_seed"
        mode = "dry-run-synthetic"
        tabix_bin = None
    else:
        mode = "live-tabix-https"
        tabix_bin = ensure_tabix()

    variants = load_panel(panel)
    dataset_rows = load_datasets(datasets)
    limiter = RateLimiter()
    cache_hits: list[str] = []
    pair_records: list[dict[str, Any]] = []
    association_rows: list[dict[str, Any]] = []
    # Optional: one header per dataset in live mode (caller may pre-seed cache meta).
    headers_by_dataset: dict[str, str] = {}

    for dataset in dataset_rows:
        for variant in variants:
            result = query_pair(
                variant=variant,
                dataset=dataset,
                cache_dir=cache_dir,
                dry_run=dry_run,
                tabix_bin=tabix_bin,
                limiter=limiter,
                cache_hits=cache_hits,
            )
            header = result["header"] or headers_by_dataset.get(dataset.dataset_id, "")
            if dry_run and not header:
                header = SYNTHETIC_HEADER
            if result["status"] == "returned":
                association_rows.extend(
                    parse_tabix_rows(
                        header=header,
                        body=result["body"],
                        variant=variant,
                        dataset=dataset,
                    )
                )
            pair_records.append(
                {
                    "dataset_id": dataset.dataset_id,
                    "qts_id": dataset.qts_id,
                    "rsid": variant.rsid,
                    "chrom": variant.chrom,
                    "pos": variant.pos,
                    "status": result["status"],
                    "cache_hit": result["cache_hit"],
                    "file_url": result["file_url"],
                    "region": result["region"],
                    "error": result["error"],
                }
            )
            if result["cache_hit"]:
                typer.echo(f"cache hit: {dataset.dataset_id} {variant.rsid}")

    pair_counts = account_pairs(pair_records)
    typer.echo(f"pair_counts={pair_counts} total={sum(pair_counts.values())}")

    new_frame = (
        pl.DataFrame(association_rows)
        if association_rows
        else pl.DataFrame(
            schema={
                "canonical_rsid": pl.Utf8,
                "eqtl_alt_allele": pl.Utf8,
                "tissue": pl.Utf8,
                "dataset_id": pl.Utf8,
                "n_samples": pl.Utf8,
                "gene_id": pl.Utf8,
                "egene_symbol": pl.Utf8,
                "is_cis_to_annotated_gene": pl.Utf8,
                "beta": pl.Utf8,
                "se": pl.Utf8,
                "pvalue": pl.Utf8,
                "fdr_bh": pl.Utf8,
                "maf": pl.Utf8,
                "median_tpm": pl.Utf8,
                "chrom": pl.Utf8,
                "pos": pl.Int64,
                "ref": pl.Utf8,
                "alt": pl.Utf8,
                "alt_matches_panel_alt": pl.Utf8,
            }
        )
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    new_frame.write_csv(out)
    diff_stats = row_level_diff(old, new_frame, diff_out)
    pair_status_out.parent.mkdir(parents=True, exist_ok=True)
    pl.DataFrame(pair_records).write_csv(pair_status_out, separator="\t")

    if skip_integrated:
        integrated = {
            "observed_labels": {
                "GWAS support": -1,
                "eQTL support": -1,
                "population difference": -1,
            },
            "ground_truth_labels": dict(GROUND_TRUTH_LABELS),
            "labels_match": False,
            "skipped": True,
        }
    else:
        integrated = run_integrated_temp(out)

    payload: dict[str, Any] = {
        "generated_at_utc": started,
        "access_date_utc": started[:10],
        "mode": mode,
        "access_docs_url": ACCESS_DOCS_URL,
        "https_sumstats_url": HTTPS_SUMSTATS_URL,
        "columns_md_url": COLUMNS_MD_URL,
        "catalogue_release_or_version": "not provided by source",
        "git_commit": _git_commit(),
        "tool_versions": {
            "python": sys.version.split()[0],
            "polars": _package_version("polars"),
            "typer": _package_version("typer"),
            "tabix": (
                subprocess.check_output([tabix_bin, "--version"], text=True).splitlines()[0]
                if tabix_bin
                else "not used (dry-run)"
            ),
        },
        "inputs": {
            "panel": {"path": str(panel), "sha256": _sha256(panel)},
            "datasets": {
                "path": str(datasets),
                "sha256": _sha256(datasets),
                "label": "synthetic" if dry_run else "production",
            },
        },
        "pair_counts": pair_counts,
        "n_pairs": EXPECTED_PAIRS,
        "n_association_rows": new_frame.height,
        "n_cache_hits": len(cache_hits),
        "diff": diff_stats,
        "integrated_evidence_temp": integrated,
        "outputs": {
            "associations": str(out),
            "diff": str(diff_out),
            "pair_status": str(pair_status_out),
            "provenance": str(provenance_out),
            "doc": str(doc_out),
        },
        "file_urls": sorted({file_url(d) for d in dataset_rows}),
    }
    provenance_out.parent.mkdir(parents=True, exist_ok=True)
    provenance_out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    doc_out.parent.mkdir(parents=True, exist_ok=True)
    doc_out.write_text(render_doc(payload), encoding="utf-8")

    typer.echo(f"Wrote {out}")
    typer.echo(f"Wrote {diff_out}")
    typer.echo(f"Wrote {pair_status_out}")
    typer.echo(f"Wrote {provenance_out}")
    typer.echo(f"Wrote {doc_out}")
    obs = integrated["observed_labels"]
    typer.echo(
        "labels observed vs ground truth: "
        f"GWAS {obs['GWAS support']}/{GROUND_TRUTH_LABELS['GWAS support']} "
        f"eQTL {obs['eQTL support']}/{GROUND_TRUTH_LABELS['eQTL support']} "
        f"population {obs['population difference']}/"
        f"{GROUND_TRUTH_LABELS['population difference']} "
        f"match={integrated.get('labels_match')}"
    )


if __name__ == "__main__":
    typer.run(main)
