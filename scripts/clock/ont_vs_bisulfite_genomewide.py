#!/usr/bin/env python3
"""Genome-wide ONT versus bisulfite agreement for the public HG002 / GM24385 sample.

Reads coverage thresholds from the committed chr20 dry-run metrics JSON and
reuses them exactly (no new tuning). Compares methylation fractions and clock
outputs from ONT bedMethyl versus bisulfite Bismark ``.cov`` on the same sample.

This is an input-agreement check. It does not estimate accuracy against
chronological age. A single sample cannot support an age-accuracy claim.

Example:
    uv run python scripts/clock/ont_vs_bisulfite_genomewide.py \\
      --ont-bedmethyl /path/to/genomewide.bedMethyl \\
      --bisulfite-cov data/dryrun/ont_vs_bisulfite/CpG.gz.bismark.cov.gz
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import subprocess
import sys
from datetime import UTC, datetime
from importlib import metadata
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
import pandas as pd
import polars as pl
import typer
from scipy.stats import pearsonr

from rogen_aging.clock.benchmark import (
    OUR_CLOCK_LABEL,
    missing_features_before,
    predict_pyaging,
)
from rogen_aging.clock.evaluate import build_feature_matrix, load_model
from rogen_aging.config import find_repo_root

REPO_ROOT = find_repo_root()
REQUIRED_PYAGING = "0.5.7"
MIN_PEARSON_N = 30
PYAGING_CLOCKS: tuple[str, ...] = ("horvath2013", "hannum")
TRANSFER_CLOCKS: tuple[str, ...] = (*PYAGING_CLOCKS, OUR_CLOCK_LABEL)

DEFAULT_CHR20_METRICS = REPO_ROOT / "analysis" / "dryrun" / "hg002_ont_methylation_metrics.json"
DEFAULT_BISULFITE = REPO_ROOT / "data" / "dryrun" / "ont_vs_bisulfite" / "CpG.gz.bismark.cov.gz"
DEFAULT_ZHOU = REPO_ROOT / "data" / "dryrun" / "reference" / "HM450.hg38.manifest.tsv.gz"
DEFAULT_MODEL = REPO_ROOT / "models" / "gse40279_hannum450k_elasticnet.joblib"
DEFAULT_PYAGING_CACHE = REPO_ROOT / "pyaging_data"
DEFAULT_OUT_DIR = REPO_ROOT / "analysis" / "clock"
FORBIDDEN_BAM = "gm24385_mod_2021.09/extra_analysis/all.bam"
CHR20_ONLY_REFUSAL = (
    "STOP: --ont-bedmethyl looks region-restricted (chr20 only). "
    "This script requires genome-wide ONT bedMethyl for the same public sample "
    "as the chr20 dry run. Do not substitute the chr20 window file."
)

app = typer.Typer(add_completion=False, no_args_is_help=True)


def _load_chr20_helpers() -> ModuleType:
    path = REPO_ROOT / "scripts" / "clock" / "ont_clock_site_coverage.py"
    assert path.is_file(), f"STOP: missing helper script: {path}"
    name = "ont_clock_site_coverage"
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None, f"STOP: cannot load {path}"
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _git_commit() -> str:
    result = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=REPO_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip()


def _repo_rel(path: Path) -> str:
    try:
        return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()
    except ValueError:
        return path.as_posix()


def load_min_depth_from_chr20_metrics(path: Path) -> int:
    """Return the chr20 dry-run ``min_depth``; refuse if levels disagree."""
    assert path.is_file(), f"STOP: missing chr20 metrics JSON: {path}"
    payload = json.loads(path.read_text(encoding="utf-8"))
    levels = payload.get("levels")
    assert isinstance(levels, list) and levels, f"STOP: no levels in {path}"
    depths: list[int] = []
    for level in levels:
        assert isinstance(level, dict), f"STOP: bad level entry in {path}"
        assert "min_depth" in level, f"STOP: level missing min_depth in {path}: {level}"
        depths.append(int(level["min_depth"]))
    assert (
        len(set(depths)) == 1
    ), f"STOP: chr20 metrics min_depth values disagree: {depths}. No new tuning."
    return depths[0]


def parse_modkit_combined_bed(path: Path) -> pl.DataFrame:
    """Parse a modkit ``--cpg --combine-strands`` bedMethyl into count columns.

    Expected columns (BED9+): chrom, start, end, code (``m``/``h``), …, N_valid,
    percent, N_mod, …. Rows for ``m`` and ``h`` at the same cytosine are summed.
    """
    assert path.is_file(), f"STOP: missing ONT bedMethyl: {path}"
    raw = pl.read_csv(
        path,
        separator="\t",
        has_header=False,
        infer_schema_length=1000,
        ignore_errors=False,
    )
    assert raw.width >= 12, f"STOP: modkit bedMethyl needs >= 12 columns, got {raw.width}: {path}"
    named = raw.select(
        [
            pl.col("column_1").cast(pl.Utf8).alias("chrom"),
            pl.col("column_2").cast(pl.Int64).alias("start"),
            pl.col("column_4").cast(pl.Utf8).alias("code"),
            pl.col("column_10").cast(pl.Int64).alias("n_valid"),
            pl.col("column_12").cast(pl.Int64).alias("n_mod"),
        ]
    ).filter(pl.col("code").is_in(["m", "h"]))
    assert named.height > 0, f"STOP: no m/h rows in ONT bedMethyl: {path}"
    pivoted = (
        named.group_by(["chrom", "start"])
        .agg(
            pl.col("n_valid").max().alias("n_valid"),
            pl.col("n_mod").filter(pl.col("code") == "m").sum().alias("n_5mc"),
            pl.col("n_mod").filter(pl.col("code") == "h").sum().alias("n_5hmc"),
        )
        .with_columns(
            [
                pl.col("n_5mc").fill_null(0).cast(pl.Int64),
                pl.col("n_5hmc").fill_null(0).cast(pl.Int64),
            ]
        )
    )
    return pivoted


def assert_genomewide_ont(pileup: pl.DataFrame) -> None:
    """Refuse a chr20-only ONT table. Genome-wide input is required."""
    chroms = set(pileup["chrom"].unique().to_list())
    if chroms <= {"chr20", "20"}:
        raise AssertionError(CHR20_ONLY_REFUSAL)
    assert (
        "chr1" in chroms or "1" in chroms
    ), f"STOP: ONT bedMethyl has no chr1; chroms={sorted(chroms)[:20]}"


def load_bismark_cov(path: Path, helpers: ModuleType) -> pl.DataFrame:
    """Load genome-wide Bismark ``.cov`` and merge strands."""
    assert path.is_file(), f"STOP: missing bisulfite cov: {path}"
    assert FORBIDDEN_BAM not in path.as_posix(), f"STOP: refused path: {path}"
    raw = pl.read_csv(
        path,
        separator="\t",
        has_header=False,
        columns=[1, 2, 5, 6],
        new_columns=["chrom", "start", "n_meth", "n_unmeth"],
        schema_overrides={
            "chrom": pl.Utf8,
            "start": pl.Int64,
            "n_meth": pl.Int64,
            "n_unmeth": pl.Int64,
        },
    )
    chroms = set(raw["chrom"].unique().to_list())
    assert (
        "chr1" in chroms or "1" in chroms
    ), f"STOP: bisulfite file has no chr1; refusing. chroms sample={sorted(chroms)[:20]}"
    if "chr1" not in chroms and "1" in chroms:
        raise AssertionError(
            "STOP: bisulfite contigs are bare integers (1), not chr1; refusing to join"
        )
    return helpers.merge_bismark_cpg_strands(raw)


def load_zhou_hm450(path: Path, helpers: ModuleType) -> pl.DataFrame:
    """Zhou HM450 hg38 probes as ``probe_id``, ``chrom``, 0-based ``start``."""
    assert path.is_file(), f"STOP: missing Zhou HM450 manifest: {path}"
    manifest = helpers.load_zhou_manifest(path)
    return (
        manifest.filter(pl.col("CpG_chrm").is_not_null() & pl.col("CpG_beg").is_not_null())
        .select(
            [
                pl.col("probeID").alias("probe_id"),
                pl.col("CpG_chrm").alias("chrom"),
                pl.col("CpG_beg").cast(pl.Int64).alias("start"),
            ]
        )
        .unique(subset=["probe_id"])
    )


def clock_probe_ids(
    *,
    pyaging_cache: Path,
    model_path: Path,
) -> dict[str, list[str]]:
    """Return Illumina probe IDs for each transfer clock."""
    from pyaging.predict import _pred as pyaging_pred

    load_clock = getattr(pyaging_pred, "load_clock")
    pyaging_cache.mkdir(parents=True, exist_ok=True)
    out: dict[str, list[str]] = {}
    for name in PYAGING_CLOCKS:
        model = load_clock(name, dir=str(pyaging_cache))
        feats = [
            str(feat)
            for feat in (getattr(model, "features", None) or [])
            if str(feat).startswith("cg") or str(feat).startswith("ch.")
        ]
        assert feats, f"STOP: pyaging clock {name!r} has no CpG features"
        out[name] = sorted(feats)
    retrain = load_model(model_path)
    names = getattr(retrain, "feature_names_in_", None)
    assert names is not None, f"STOP: model missing feature_names_in_: {model_path}"
    out[OUR_CLOCK_LABEL] = [
        str(name) for name in names if str(name).startswith("cg") or str(name).startswith("ch.")
    ]
    assert out[OUR_CLOCK_LABEL], f"STOP: no cg* features in {model_path}"
    return out


def pearson_report(ont: list[float], bisulfite: list[float]) -> dict[str, Any]:
    """Pearson r with explicit n; ``not estimated`` when n < 30."""
    assert len(ont) == len(bisulfite)
    n = len(ont)
    if n < MIN_PEARSON_N:
        return {
            "n": n,
            "pearson_r": "not estimated",
            "mean_ont_minus_bisulfite": (
                float(np.mean(np.asarray(ont) - np.asarray(bisulfite))) if n else None
            ),
        }
    corr = float(pearsonr(ont, bisulfite).statistic)
    delta = float(np.mean(np.asarray(ont, dtype=float) - np.asarray(bisulfite, dtype=float)))
    return {"n": n, "pearson_r": corr, "mean_ont_minus_bisulfite": delta}


def compare_level(
    *,
    level: str,
    sites: pl.DataFrame,
    ont_frac: pl.DataFrame,
    bs_frac: pl.DataFrame,
    min_depth: int,
) -> dict[str, Any]:
    """Join site table to ONT and bisulfite fractions; report Pearson."""
    keys = ["chrom", "start"]
    ont_cov = ont_frac.filter(pl.col("n_valid") >= min_depth)
    joined = (
        sites.select(keys)
        .unique()
        .join(ont_cov.select([*keys, "ont", "n_valid"]), on=keys, how="inner")
        .join(bs_frac.select([*keys, "bisulfite"]), on=keys, how="inner")
        .filter(pl.col("ont").is_not_nan() & pl.col("bisulfite").is_not_nan())
    )
    report = pearson_report(joined["ont"].to_list(), joined["bisulfite"].to_list())
    return {
        "level": level,
        "min_depth": min_depth,
        "n_sites_in_level": sites.select(keys).unique().height,
        "n_ont_depth_ge_min": ont_cov.join(
            sites.select(keys).unique(), on=keys, how="inner"
        ).height,
        **report,
    }


def betas_from_joined(
    joined: pl.DataFrame,
    probe_ids: list[str],
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Build a one-row beta frame for ``probe_ids`` from probe_id/fraction table."""
    by_probe = {
        str(row["probe_id"]): float(row["fraction"])
        for row in joined.select(["probe_id", "fraction"]).iter_rows(named=True)
    }
    present = [pid for pid in probe_ids if pid in by_probe]
    missing = [pid for pid in probe_ids if pid not in by_probe]
    row: dict[str, float] = {pid: by_probe[pid] for pid in present}
    frame = pd.DataFrame([row])
    # Placeholder age only so pyaging/rogen loaders that require the column work.
    # Never used as an accuracy target in this script.
    frame["chronological_age"] = np.nan
    return frame, {
        "n_clock_sites": len(probe_ids),
        "n_present": len(present),
        "n_missing": len(missing),
        "fraction_missing": (len(missing) / len(probe_ids)) if probe_ids else None,
        "missing_probe_ids": missing,
    }


def score_transfer(
    *,
    ont_betas: pd.DataFrame,
    bs_betas: pd.DataFrame,
    model_path: Path,
    pyaging_cache: Path,
    probe_sets: dict[str, list[str]],
) -> dict[str, Any]:
    """Clock ages from ONT betas versus bisulfite betas. Agreement only."""
    pyaging_version = metadata.version("pyaging")
    assert (
        pyaging_version == REQUIRED_PYAGING
    ), f"STOP: pyaging {pyaging_version} != required {REQUIRED_PYAGING}"

    ont_missing_before = missing_features_before(ont_betas, PYAGING_CLOCKS, cache_dir=pyaging_cache)
    bs_missing_before = missing_features_before(bs_betas, PYAGING_CLOCKS, cache_dir=pyaging_cache)
    ont_preds, ont_uns = predict_pyaging(ont_betas, PYAGING_CLOCKS, cache_dir=pyaging_cache)
    bs_preds, bs_uns = predict_pyaging(bs_betas, PYAGING_CLOCKS, cache_dir=pyaging_cache)

    model = load_model(model_path)
    ont_x, ont_align = build_feature_matrix(ont_betas, model)
    bs_x, bs_align = build_feature_matrix(bs_betas, model)
    ont_retrain = float(np.asarray(model.predict(ont_x), dtype=float)[0])
    bs_retrain = float(np.asarray(model.predict(bs_x), dtype=float)[0])

    clocks: dict[str, Any] = {}
    for name in PYAGING_CLOCKS:
        ont_age = float(ont_preds[name].iloc[0])
        bs_age = float(bs_preds[name].iloc[0])
        ont_miss = list(ont_uns[name].get("missing_features") or [])
        bs_miss = list(bs_uns[name].get("missing_features") or [])
        n_features = len(probe_sets[name])
        assert n_features > 0, f"STOP: empty feature list for {name}"
        clocks[name] = {
            "ont_age": ont_age,
            "bisulfite_age": bs_age,
            "ont_minus_bisulfite_age": ont_age - bs_age,
            "ont_missing_before_imputation": {
                "n": len(ont_missing_before[name]),
                "fraction": len(ont_missing_before[name]) / n_features,
                "probe_ids": ont_missing_before[name],
            },
            "bisulfite_missing_before_imputation": {
                "n": len(bs_missing_before[name]),
                "fraction": len(bs_missing_before[name]) / n_features,
                "probe_ids": bs_missing_before[name],
            },
            "ont_missing_features_after_pyaging": ont_miss,
            "bisulfite_missing_features_after_pyaging": bs_miss,
            "ont_percent_na": ont_uns[name].get("percent_na"),
            "bisulfite_percent_na": bs_uns[name].get("percent_na"),
            "imputation": "pyaging df_to_adata imputer_strategy=knn",
            "n_features": n_features,
        }

    clocks[OUR_CLOCK_LABEL] = {
        "ont_age": ont_retrain,
        "bisulfite_age": bs_retrain,
        "ont_minus_bisulfite_age": ont_retrain - bs_retrain,
        "ont_missing_imputed": {
            "n": int(ont_align["n_missing_imputed"]),
            "fraction": (
                float(ont_align["n_missing_imputed"]) / float(ont_align["n_clock_cpgs"])
                if ont_align["n_clock_cpgs"]
                else None
            ),
            "probe_ids": list(ont_align["missing_cpg_ids"]),
        },
        "bisulfite_missing_imputed": {
            "n": int(bs_align["n_missing_imputed"]),
            "fraction": (
                float(bs_align["n_missing_imputed"]) / float(bs_align["n_clock_cpgs"])
                if bs_align["n_clock_cpgs"]
                else None
            ),
            "probe_ids": list(bs_align["missing_cpg_ids"]),
        },
        "imputation": "rogen ElasticNet pipeline training imputer",
        "n_clock_cpgs": int(ont_align["n_clock_cpgs"]),
    }
    return {
        "pyaging_version": pyaging_version,
        "note": (
            "Clock ages compare ONT input vs bisulfite input on one sample. "
            "This is agreement between inputs, not accuracy against chronological age."
        ),
        "clocks": clocks,
    }


@app.command()
def main(
    ont_bedmethyl: Path = typer.Option(
        ...,
        "--ont-bedmethyl",
        exists=False,
        help="Genome-wide ONT modkit bedMethyl (required; no chr20 substitute).",
    ),
    bisulfite_cov: Path = typer.Option(
        DEFAULT_BISULFITE,
        "--bisulfite-cov",
        help="Genome-wide Bismark CpG cov for the same public sample.",
    ),
    chr20_metrics: Path = typer.Option(
        DEFAULT_CHR20_METRICS,
        "--chr20-metrics",
        help="Committed chr20 metrics JSON; min_depth is read from here.",
    ),
    zhou_manifest: Path = typer.Option(DEFAULT_ZHOU, "--zhou-manifest"),
    model: Path = typer.Option(DEFAULT_MODEL, "--model"),
    pyaging_cache: Path = typer.Option(DEFAULT_PYAGING_CACHE, "--pyaging-cache"),
    out_dir: Path = typer.Option(DEFAULT_OUT_DIR, "--out-dir"),
) -> None:
    """Run genome-wide ONT vs bisulfite CpG and clock-transfer checks."""
    started = datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")
    helpers = _load_chr20_helpers()

    assert ont_bedmethyl.is_file(), (
        f"STOP: missing genome-wide ONT bedMethyl: {ont_bedmethyl}. "
        "Only the chr20 window file exists in this checkout "
        "(data/dryrun/ont_vs_bisulfite/modkit_cpg_chr20.bed). "
        "Do not pass that file here."
    )
    assert bisulfite_cov.is_file(), f"STOP: missing bisulfite cov: {bisulfite_cov}"
    assert model.is_file(), f"STOP: missing retrained Hannum model: {model}"
    assert FORBIDDEN_BAM not in ont_bedmethyl.as_posix()
    assert FORBIDDEN_BAM not in bisulfite_cov.as_posix()

    min_depth = load_min_depth_from_chr20_metrics(chr20_metrics)
    typer.echo(f"min_depth={min_depth} (from {_repo_rel(chr20_metrics)}; no new tuning)")

    pyaging_version = metadata.version("pyaging")
    assert (
        pyaging_version == REQUIRED_PYAGING
    ), f"STOP: pyaging {pyaging_version} != required {REQUIRED_PYAGING}"

    ont_counts = parse_modkit_combined_bed(ont_bedmethyl)
    assert_genomewide_ont(ont_counts)
    ont_frac = helpers.ont_modification_fraction(ont_counts).select(
        [
            "chrom",
            "start",
            "n_valid",
            pl.col("fraction").alias("ont"),
        ]
    )

    bs_merged = load_bismark_cov(bisulfite_cov, helpers)
    bs_frac = bs_merged.select(["chrom", "start", pl.col("fraction").alias("bisulfite")])

    hm450 = load_zhou_hm450(zhou_manifest, helpers)
    probe_sets = clock_probe_ids(pyaging_cache=pyaging_cache, model_path=model)
    all_clock_probes = sorted({pid for ids in probe_sets.values() for pid in ids})

    # All CpGs: intersection of ONT (any depth) and bisulfite positions, then
    # Pearson on sites that also pass ONT min_depth (same filter spirit as chr20).
    all_sites = (
        ont_frac.select(["chrom", "start"])
        .join(bs_frac.select(["chrom", "start"]), on=["chrom", "start"], how="inner")
        .unique()
    )
    levels: list[dict[str, Any]] = [
        compare_level(
            level="all_cpgs",
            sites=all_sites,
            ont_frac=ont_frac,
            bs_frac=bs_frac,
            min_depth=min_depth,
        ),
        compare_level(
            level="hm450",
            sites=hm450,
            ont_frac=ont_frac,
            bs_frac=bs_frac,
            min_depth=min_depth,
        ),
    ]

    hm450_keyed = hm450.filter(pl.col("probe_id").is_in(all_clock_probes))
    for clock_name, probes in probe_sets.items():
        clock_sites = hm450.filter(pl.col("probe_id").is_in(probes))
        levels.append(
            compare_level(
                level=f"clock_{clock_name}",
                sites=clock_sites,
                ont_frac=ont_frac,
                bs_frac=bs_frac,
                min_depth=min_depth,
            )
        )

    # Beta matrices for transfer: HM450 probes covered at min_depth on each side.
    ont_probe = (
        hm450_keyed.join(
            ont_frac.filter(pl.col("n_valid") >= min_depth),
            on=["chrom", "start"],
            how="inner",
        )
        .select(["probe_id", pl.col("ont").alias("fraction")])
        .unique(subset=["probe_id"])
    )
    bs_probe = (
        hm450_keyed.join(bs_frac, on=["chrom", "start"], how="inner")
        .select(["probe_id", pl.col("bisulfite").alias("fraction")])
        .unique(subset=["probe_id"])
    )
    ont_frame, ont_cov_report = betas_from_joined(ont_probe, all_clock_probes)
    bs_frame, bs_cov_report = betas_from_joined(bs_probe, all_clock_probes)

    transfer = score_transfer(
        ont_betas=ont_frame,
        bs_betas=bs_frame,
        model_path=model,
        pyaging_cache=pyaging_cache,
        probe_sets=probe_sets,
    )
    # Attach coverage-before-clock reports (sites absent before imputation).
    transfer["coverage_before_clock"] = {
        "ont": ont_cov_report,
        "bisulfite": bs_cov_report,
        "min_depth": min_depth,
    }

    out_dir.mkdir(parents=True, exist_ok=True)
    metrics_path = out_dir / "ont_vs_bisulfite_genomewide_metrics.json"
    provenance_path = out_dir / "ont_vs_bisulfite_genomewide_provenance.json"

    metrics: dict[str, Any] = {
        "sample": "HG002 / NA24385 / GM24385 (public ONT open data + GM24385 bisulfite)",
        "scope": "genome-wide",
        "min_depth": min_depth,
        "min_depth_source": _repo_rel(chr20_metrics),
        "min_pearson_n": MIN_PEARSON_N,
        "levels": levels,
        "transfer_check": transfer,
        "prohibited": ("Do not report age accuracy from this single-sample transfer check."),
    }
    provenance: dict[str, Any] = {
        "generated_at_utc": started,
        "git_commit": _git_commit(),
        "python": sys.version.split()[0],
        "pyaging_version": pyaging_version,
        "polars_version": pl.__version__,
        "numpy_version": np.__version__,
        "required_pyaging": REQUIRED_PYAGING,
        "inputs": {
            "ont_bedmethyl": {
                "path": _repo_rel(ont_bedmethyl),
                "sha256": _sha256(ont_bedmethyl),
            },
            "bisulfite_cov": {
                "path": _repo_rel(bisulfite_cov),
                "sha256": _sha256(bisulfite_cov),
            },
            "chr20_metrics": {
                "path": _repo_rel(chr20_metrics),
                "sha256": _sha256(chr20_metrics),
            },
            "zhou_manifest": {
                "path": _repo_rel(zhou_manifest),
                "sha256": _sha256(zhou_manifest),
            },
            "model": {
                "path": _repo_rel(model),
                "sha256": _sha256(model),
            },
        },
        "outputs": {
            "metrics": _repo_rel(metrics_path),
            "provenance": _repo_rel(provenance_path),
        },
        "thresholds": {
            "min_depth": min_depth,
            "source": "chr20 dry-run metrics levels[*].min_depth (no new tuning)",
            "min_pearson_n": MIN_PEARSON_N,
        },
    }

    metrics_path.write_text(json.dumps(metrics, indent=2) + "\n", encoding="utf-8")
    provenance_path.write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    typer.echo(f"Wrote {_repo_rel(metrics_path)}")
    typer.echo(f"Wrote {_repo_rel(provenance_path)}")
    for level in levels:
        typer.echo(
            f"{level['level']}: n={level['n']} pearson_r={level['pearson_r']!r} "
            f"mean_ont_minus_bisulfite={level.get('mean_ont_minus_bisulfite')}"
        )
    for name in TRANSFER_CLOCKS:
        row = transfer["clocks"][name]
        typer.echo(
            f"transfer {name}: ont={row['ont_age']:.4f} "
            f"bisulfite={row['bisulfite_age']:.4f} "
            f"delta={row['ont_minus_bisulfite_age']:.4f}"
        )


if __name__ == "__main__":
    app()
