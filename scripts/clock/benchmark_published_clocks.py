#!/usr/bin/env python3
"""Benchmark published pyaging clocks against the retrained Hannum model.

Cohorts: GSE87571 (all samples) and the GSE40279 held-out test split only.
Never scores GSE40279 training samples with the retrained model.

Example:
    uv run python scripts/clock/benchmark_published_clocks.py
"""

from __future__ import annotations

import json
from datetime import UTC, datetime
from importlib import metadata
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer

from rogen_aging.clock.benchmark import (
    OUR_CLOCK_LABEL,
    TARGET_CLOCKS,
    age_acceleration,
    bonferroni_adjust,
    clock_feature_union,
    list_pyaging_clocks,
    load_gse40279_heldout_cohort,
    load_gse87571_cohort,
    metrics_for_clock,
    plot_accel_by_sex,
    plot_predicted_vs_age,
    predict_pyaging,
    predict_rogen_hannum,
    run_release_check,
    sex_acceleration_model,
    training_overlap_notes,
    try_include_grimage2,
)
from rogen_aging.config import find_repo_root

REPO_ROOT = find_repo_root()
DEFAULT_OUT = REPO_ROOT / "analysis" / "clocks_benchmark"
DEFAULT_MODEL = REPO_ROOT / "models" / "gse40279_hannum450k_elasticnet.joblib"
DEFAULT_GSE40279 = REPO_ROOT / "data" / "methylation" / "GSE40279_processed.parquet"
DEFAULT_GSE40279_MATRIX = REPO_ROOT / "data" / "geo" / "GSE40279_series_matrix.txt.gz"
DEFAULT_GSE87571_MATRIX = REPO_ROOT / "data" / "geo" / "GSE87571_series_matrix.txt.gz"
DEFAULT_GEO_CACHE = REPO_ROOT / "data" / "geo"
DEFAULT_PYAGING_CACHE = REPO_ROOT / "pyaging_data"


def _write_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False)


def main(
    out_dir: Path = typer.Option(DEFAULT_OUT, "--out-dir", help="Aggregate output directory."),
    model: Path = typer.Option(DEFAULT_MODEL, "--model", help="Retrained Hannum joblib."),
    gse40279_parquet: Path = typer.Option(DEFAULT_GSE40279, "--gse40279-parquet"),
    gse40279_matrix: Path = typer.Option(DEFAULT_GSE40279_MATRIX, "--gse40279-matrix"),
    gse87571_matrix: Path = typer.Option(DEFAULT_GSE87571_MATRIX, "--gse87571-matrix"),
    geo_cache_dir: Path = typer.Option(DEFAULT_GEO_CACHE, "--geo-cache-dir"),
    pyaging_cache: Path = typer.Option(DEFAULT_PYAGING_CACHE, "--pyaging-cache"),
    skip_grimage2: bool = typer.Option(False, "--skip-grimage2", help="Force-skip GrimAge2."),
) -> None:
    """Run the published-clock vs retrained-Hannum benchmark."""
    started = datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")
    out_dir.mkdir(parents=True, exist_ok=True)
    pyaging_version = metadata.version("pyaging")
    all_clock_names = list_pyaging_clocks()
    typer.echo(f"pyaging {pyaging_version}; {len(all_clock_names)} clocks listed")

    clocks = list(TARGET_CLOCKS)
    grimage_skip_reason: str | None = None
    if skip_grimage2:
        grimage_skip_reason = "Skipped by --skip-grimage2."
    else:
        clocks, grimage_skip_reason = try_include_grimage2(clocks)
        # Probe: attempt will happen later; if it fails we remove it

    typer.echo(f"Loading feature union for: {clocks}")
    restrict_set = clock_feature_union(clocks, cache_dir=pyaging_cache)
    # Keep retrained Hannum probes so scoring is not all-imputed.
    from rogen_aging.clock.evaluate import load_model

    our_model = load_model(model)
    our_feats = getattr(our_model, "feature_names_in_", None)
    if our_feats is None and hasattr(our_model, "named_steps"):
        for step in reversed(list(our_model.named_steps.values())):
            our_feats = getattr(step, "feature_names_in_", None)
            if our_feats is not None:
                break
    for feat in list(our_feats) if our_feats is not None else []:
        text = str(feat)
        if text.startswith("cg") or text.startswith("ch."):
            restrict_set.add(text)
    restrict = sorted(restrict_set)
    typer.echo(
        f"Restricting to {len(restrict)} Illumina probes "
        f"(pyaging clocks + retrained Hannum features)"
    )

    typer.echo("Loading GSE87571…")
    gse87571 = load_gse87571_cohort(
        geo_cache_dir=geo_cache_dir,
        series_matrix=gse87571_matrix if gse87571_matrix.is_file() else None,
        restrict_cpgs=restrict,
    )
    typer.echo(f"GSE87571: {gse87571.frame.shape[0]} samples × {gse87571.frame.shape[1]} columns")

    typer.echo("Loading GSE40279 held-out test split…")
    gse40279 = load_gse40279_heldout_cohort(
        parquet_path=gse40279_parquet,
        series_matrix=gse40279_matrix,
        restrict_cpgs=restrict,
    )
    typer.echo(
        f"GSE40279_heldout: {gse40279.frame.shape[0]} samples × {gse40279.frame.shape[1]} columns"
    )

    cohorts = [gse87571, gse40279]
    metric_rows: list[dict[str, Any]] = []
    sex_rows: list[dict[str, Any]] = []
    per_sample_parts: list[pd.DataFrame] = []
    uns_by_cohort: dict[str, Any] = {}
    clocks_used = list(clocks)
    clock_metadata: dict[str, Any] = {}

    for cohort in cohorts:
        typer.echo(f"Predicting pyaging clocks on {cohort.name}…")
        try:
            preds, uns = predict_pyaging(
                cohort.frame,
                clocks_used,
                cache_dir=pyaging_cache,
            )
        except Exception as exc:
            if "grimage2" in clocks_used:
                typer.echo(f"GrimAge2 failed on {cohort.name}: {exc}")
                grimage_skip_reason = (
                    f"GrimAge2 predict_age failed with available inputs "
                    f"(sex/age + CpGs): {type(exc).__name__}: {exc}"
                )
                clocks_used = [c for c in clocks_used if c != "grimage2"]
                preds, uns = predict_pyaging(
                    cohort.frame,
                    clocks_used,
                    cache_dir=pyaging_cache,
                )
            else:
                raise
        uns_by_cohort[cohort.name] = uns
        for name, payload in uns.items():
            if name not in clock_metadata and isinstance(payload.get("metadata"), dict):
                clock_metadata[name] = payload["metadata"]

        typer.echo(f"Predicting retrained Hannum on {cohort.name}…")
        our = predict_rogen_hannum(cohort.frame, model)
        combined = preds.copy()
        combined[OUR_CLOCK_LABEL] = our.to_numpy()
        age = cohort.frame["chronological_age"].to_numpy(dtype=float)
        female = cohort.frame["female"].to_numpy(dtype=float)
        sex = cohort.frame["sex"].astype(str).to_numpy()
        ethnicity = (
            cohort.frame["ethnicity"].astype(object).to_numpy()
            if "ethnicity" in cohort.frame.columns
            else np.array([None] * len(cohort.frame))
        )

        sample_table = pd.DataFrame(
            {
                "sample_id": cohort.frame.index.astype(str),
                "cohort": cohort.name,
                "chronological_age": age,
                "sex": sex,
                "female": female,
                "ethnicity": ethnicity,
            },
            index=cohort.frame.index,
        )
        all_clocks = [*clocks_used, OUR_CLOCK_LABEL]
        for clock in all_clocks:
            values = combined[clock].to_numpy(dtype=float)
            sample_table[clock] = values
            accel = age_acceleration(values, age)
            sample_table[f"{clock}_accel"] = accel
            metric_rows.append(
                metrics_for_clock(
                    clock=clock,
                    cohort=cohort.name,
                    chronological_age=age,
                    predicted=values,
                )
            )

        # Sex models (Bonferroni across clocks within cohort)
        sex_p: list[float] = []
        sex_partial: list[dict[str, Any]] = []
        for clock in all_clocks:
            accel = sample_table[f"{clock}_accel"].to_numpy(dtype=float)
            fit = sex_acceleration_model(accel=accel, female=female, age=age)
            sex_p.append(float(fit["sex_p"]))
            sex_partial.append({"clock": clock, "cohort": cohort.name, **fit})
        adjusted = bonferroni_adjust(sex_p, n_tests=len(all_clocks))
        for row, p_adj in zip(sex_partial, adjusted, strict=True):
            row["sex_p_bonferroni"] = p_adj
            row["n_tests_bonferroni"] = len(all_clocks)
            sex_rows.append(row)

        # Ethnicity descriptive (GSE40279 only)
        if cohort.name == "GSE40279_heldout":
            eth = sample_table.copy()
            eth_rows = []
            for clock in all_clocks:
                grouped = (
                    eth.groupby("ethnicity", dropna=False)[f"{clock}_accel"]
                    .agg(["count", "mean", "std"])
                    .reset_index()
                )
                grouped.insert(0, "clock", clock)
                grouped.insert(1, "cohort", cohort.name)
                eth_rows.append(grouped)
            ethnicity_table = pd.concat(eth_rows, ignore_index=True)
            _write_csv(ethnicity_table, out_dir / "ethnicity_accel_descriptive.csv")

        per_sample_parts.append(sample_table.reset_index(drop=True))

    per_sample = pd.concat(per_sample_parts, ignore_index=True)
    metrics = pd.DataFrame(metric_rows)
    sex_table = pd.DataFrame(sex_rows)
    all_clocks_final = [*clocks_used, OUR_CLOCK_LABEL]

    # Acceleration correlation matrices per cohort
    corr_frames: list[pd.DataFrame] = []
    for cohort_name in per_sample["cohort"].unique():
        sub = per_sample.loc[per_sample["cohort"] == cohort_name]
        accel_cols = [f"{c}_accel" for c in all_clocks_final]
        corr = sub[accel_cols].corr()
        corr.insert(0, "cohort", cohort_name)
        corr.insert(1, "clock", [c.replace("_accel", "") for c in corr.index])
        corr_frames.append(corr.reset_index(drop=True))
    corr_table = pd.concat(corr_frames, ignore_index=True)

    metrics_path = out_dir / "metrics_by_clock_cohort.csv"
    sex_path = out_dir / "sex_acceleration_models.csv"
    corr_path = out_dir / "acceleration_correlation.csv"
    sample_path = out_dir / "per_sample_predictions.csv"
    accel_rows: list[dict[str, Any]] = []
    for cohort_name, sub in per_sample.groupby("cohort"):
        for clock in all_clocks_final:
            col = f"{clock}_accel"
            vals = sub[col].astype(float)
            accel_rows.append(
                {
                    "clock": clock,
                    "cohort": cohort_name,
                    "n": int(len(vals)),
                    "accel_mean": float(vals.mean()),
                    "accel_sd": float(vals.std(ddof=1)),
                    "accel_mean_male": float(sub.loc[sub["sex"] == "M", col].astype(float).mean()),
                    "accel_mean_female": float(
                        sub.loc[sub["sex"] == "F", col].astype(float).mean()
                    ),
                }
            )
    accel_table = pd.DataFrame(accel_rows)
    accel_path = out_dir / "acceleration_by_clock_cohort.csv"
    _write_csv(metrics, metrics_path)
    _write_csv(sex_table, sex_path)
    _write_csv(corr_table, corr_path)
    _write_csv(accel_table, accel_path)
    # per-sample may be large but still aggregate-ish; keep for transparency
    _write_csv(per_sample, sample_path)

    for cohort_name in ("GSE87571", "GSE40279_heldout"):
        plot_predicted_vs_age(
            per_sample,
            all_clocks_final,
            cohort_name,
            out_dir / f"predicted_vs_age_{cohort_name}.png",
        )
        plot_accel_by_sex(
            per_sample,
            all_clocks_final,
            cohort_name,
            out_dir / f"acceleration_by_sex_{cohort_name}.png",
        )

    report_dir = out_dir / "release_checks"
    release_checks = [
        run_release_check(metrics_path, report_dir / "metrics_by_clock_cohort.json"),
        run_release_check(sex_path, report_dir / "sex_acceleration_models.json"),
        run_release_check(corr_path, report_dir / "acceleration_correlation.json"),
        run_release_check(accel_path, report_dir / "acceleration_by_clock_cohort.json"),
    ]

    overlap = training_overlap_notes(clock_metadata)
    provenance = {
        "generated_at_utc": started,
        "finished_at_utc": datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "pyaging_version": pyaging_version,
        "pyaging_clocks_available": all_clock_names,
        "pyaging_clocks_available_count": len(all_clock_names),
        "target_clocks_requested": list(TARGET_CLOCKS),
        "phenoage_note": (
            "pyaging clock name 'phenoage' is clinical Phenotypic Age (blood biomarkers). "
            "This methylation benchmark uses 'dnamphenoage' (Levine DNAm PhenoAge)."
        ),
        "clocks_run": clocks_used,
        "grimage2": {
            "included": "grimage2" in clocks_used,
            "skip_reason": grimage_skip_reason,
        },
        "our_clock": {
            "label": OUR_CLOCK_LABEL,
            "model_path": str(model),
        },
        "split": {
            "gse40279_test_size": 0.2,
            "gse40279_random_state": 42,
            "n_train_expected": 524,
            "n_test_expected": 132,
            "never_score_train_with_our_model": True,
        },
        "cpg_restriction": {
            "n_probes": len(restrict),
            "note": "Columns restricted to the union of Illumina features required by selected clocks.",
        },
        "cohort_notes": {
            gse87571.name: gse87571.notes,
            gse40279.name: gse40279.notes,
        },
        "pyaging_uns": uns_by_cohort,
        "training_overlap": overlap,
        "outputs": {
            "metrics": str(metrics_path.relative_to(REPO_ROOT)),
            "acceleration": str(accel_path.relative_to(REPO_ROOT)),
            "sex_models": str(sex_path.relative_to(REPO_ROOT)),
            "acceleration_correlation": str(corr_path.relative_to(REPO_ROOT)),
            "per_sample": str(sample_path.relative_to(REPO_ROOT)),
        },
        "release_check": release_checks,
    }
    provenance_path = out_dir / "provenance.json"
    provenance_path.write_text(
        json.dumps(provenance, indent=2, default=str) + "\n", encoding="utf-8"
    )
    typer.echo(f"Wrote {provenance_path}")
    for check in release_checks:
        typer.echo(f"release-check {check['path']}: {check['status']}")
    if any(c["status"] != "passed" for c in release_checks):
        raise typer.Exit(code=1)
    typer.echo("Done.")


if __name__ == "__main__":
    typer.run(main)
