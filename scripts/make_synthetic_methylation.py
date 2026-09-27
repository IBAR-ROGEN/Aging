#!/usr/bin/env python3
"""Generate a synthetic methylation beta matrix and run the existing clock harness.

Synthetic-only: does not download GEO, call GEOparse/GEOquery, or fetch a series
matrix. CpG identifiers are read from the harness clock artifact
(``models/fixtures/fixture_clock_cg_test.pkl`` or the joblib fallback), never from
a published coefficient table.

Example:
    uv run python scripts/make_synthetic_methylation.py
"""

from __future__ import annotations

import importlib.util
import json
import shutil
import tempfile
import traceback
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
import pandas as pd
import typer
from sklearn.impute import SimpleImputer

from rogen_aging.clock.evaluate import build_feature_matrix
from rogen_aging.config import find_repo_root

REPO_ROOT = find_repo_root()
HARNESS_PATH = REPO_ROOT / "scripts" / "clock" / "evaluate_methylation_clock.py"
MODEL_PKL = REPO_ROOT / "models" / "fixtures" / "fixture_clock_cg_test.pkl"
MODEL_JOBLIB = REPO_ROOT / "models" / "methylation_clock_v1.joblib"
SYNTHETIC_DIR = REPO_ROOT / "synthetic"
FILE_PREFIX = "SYNTHETIC — "
TITLE_PREFIX = "SYNTHETIC — "
METRIC_DISCLAIMER = "synthetic data, not a performance estimate"

N_SAMPLES = 200
AGE_MIN = 20.0
AGE_MAX = 90.0
DROP_CPG_FRACTION = 0.05
NAN_VALUE_FRACTION = 0.02
N_ORPHAN_SAMPLES = 20
DEFAULT_SEED = 20260830
# Residual injected before clipping; not tuned for a target correlation.
AGE_RESIDUAL_SD_YEARS = 8.0
BETA_BASELINE_SD = 0.12
BETA_CLIP_LO = 0.02
BETA_CLIP_HI = 0.98


@dataclass
class ModeResult:
    """Pass/fail record for one injected failure mode."""

    name: str
    outcome: str
    missing_cpg_behaviour: str
    detail: str
    error: str | None = None
    n_samples_scored: int | None = None
    n_imputed_missing_cpgs: int | None = None
    imputed_cpg_ids: list[str] = field(default_factory=list)
    warning_messages: list[str] = field(default_factory=list)


def _load_harness() -> ModuleType:
    """Import ``scripts/clock/evaluate_methylation_clock.py`` by path."""
    if not HARNESS_PATH.is_file():
        raise FileNotFoundError(
            "Clock evaluation harness not found. Searched for "
            f"{HARNESS_PATH.relative_to(REPO_ROOT)} (beta matrix → age predictions). "
            "Also considered scripts/clock/validate_clock.py and "
            "rogen_aging.clock.evaluate.evaluate_clock. Halting — will not write a replacement clock."
        )
    name = "evaluate_methylation_clock_synthetic"
    spec = importlib.util.spec_from_file_location(name, HARNESS_PATH)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load harness from {HARNESS_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def resolve_coefficient_file() -> Path:
    """Return the harness clock artifact that stores CpG coefficients.

    Returns:
        Path to an existing pickle or joblib clock.

    Raises:
        FileNotFoundError: If neither configured artifact exists.
    """
    if MODEL_PKL.is_file():
        return MODEL_PKL
    if MODEL_JOBLIB.is_file():
        return MODEL_JOBLIB
    raise FileNotFoundError(
        "Harness coefficient file is missing. Searched:\n"
        f"  - {MODEL_PKL}\n"
        f"  - {MODEL_JOBLIB}\n"
        "Halting — will not substitute a published clock's coefficients."
    )


def load_clock_coefficients(model_path: Path, harness: ModuleType) -> tuple[Any, pd.Series]:
    """Read CpG names and ElasticNet weights from the harness artifact.

    Args:
        model_path: Pickle or joblib path used by the evaluator.
        harness: Loaded ``evaluate_methylation_clock`` module.

    Returns:
        ``(model, coefficient_series)`` indexed by probe ID.
    """
    model = harness.load_elasticnet_clock(model_path)
    weights = harness.extract_cpg_coefficients(model)
    if weights.empty:
        raise ValueError(f"Coefficient set from {model_path} is empty.")
    return model, weights


def describe_missing_cpg_behaviour(model: Any) -> str:
    """Name the harness imputation policy for CpGs absent from the matrix.

    Args:
        model: Loaded clock estimator or pipeline.

    Returns:
        One of ``mean-imputed (training SimpleImputer)``,
        ``mean-imputed (test-set / global mean)``, or ``unknown``.
    """
    if hasattr(model, "named_steps"):
        for step in model.named_steps.values():
            if isinstance(step, SimpleImputer) and hasattr(step, "statistics_"):
                strategy = str(getattr(step, "strategy", "mean"))
                return (
                    "mean-imputed (training-set column means via sklearn "
                    f"SimpleImputer(strategy={strategy!r}); absent probes are left "
                    "NaN in build_feature_matrix and filled at predict time from "
                    "imputer.statistics_ — not dropped, not test-set mean)"
                )
    return (
        "mean-imputed (test-set column mean, else global mean of present CpGs, "
        "in rogen_aging.clock.evaluate.build_feature_matrix; not dropped)"
    )


def simulate_beta_matrix(
    ages: np.ndarray,
    coefficients: pd.Series,
    *,
    intercept: float,
    train_means: np.ndarray,
    rng: np.random.Generator,
    residual_sd: float = AGE_RESIDUAL_SD_YEARS,
    baseline_sd: float = BETA_BASELINE_SD,
) -> pd.DataFrame:
    """Simulate betas whose linear clock tracks age with a fixed residual.

    Generation is not tuned: one residual SD and clipping are written to
    GENERATION_PARAMS and left unchanged.

    Args:
        ages: Chronological ages in years.
        coefficients: Clock weights indexed by CpG ID.
        intercept: Clock intercept.
        train_means: Per-CpG training means (imputer statistics when present).
        rng: NumPy generator.
        residual_sd: Gaussian noise added to the target age (years).
        baseline_sd: SD of the unstructured beta draw around training means.

    Returns:
        Sample-by-CpG DataFrame of clipped beta values in ``(0, 1)``.
    """
    n = int(ages.shape[0])
    names = [str(c) for c in coefficients.index]
    weights = coefficients.to_numpy(dtype=float)
    means = np.asarray(train_means, dtype=float)
    if means.shape[0] != len(names):
        means = np.full(len(names), 0.5, dtype=float)

    betas = rng.normal(loc=means, scale=baseline_sd, size=(n, len(names)))
    betas = np.clip(betas, BETA_CLIP_LO, BETA_CLIP_HI)

    target = ages + rng.normal(0.0, residual_sd, size=n)
    current = intercept + betas @ weights
    denom = float(weights @ weights)
    if denom > 0.0:
        delta = target - current
        betas = betas + np.outer(delta, weights / denom)
    betas = np.clip(betas, BETA_CLIP_LO, BETA_CLIP_HI)
    return pd.DataFrame(betas, columns=names)


def inject_failure_modes(
    meth: pd.DataFrame,
    meta: pd.DataFrame,
    *,
    rng: np.random.Generator,
    drop_fraction: float,
    nan_fraction: float,
    n_orphan: int,
    apply_drop_cpgs: bool,
    apply_nans: bool,
    apply_orphans: bool,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Copy the complete cohort and inject selected failure modes.

    Args:
        meth: Complete ``sample_id`` + ``cg*`` matrix.
        meta: Complete ``sample_id`` + ``chronological_age`` table.
        rng: NumPy generator.
        drop_fraction: Fraction of clock CpGs to remove entirely.
        nan_fraction: Fraction of remaining beta cells to set to NaN.
        n_orphan: Number of matrix sample IDs to omit from metadata.
        apply_drop_cpgs: If True, drop CpG columns.
        apply_nans: If True, insert NaN cells.
        apply_orphans: If True, drop orphan IDs from metadata.

    Returns:
        ``(matrix, metadata, injection_record)``.
    """
    out_meth = meth.copy()
    out_meta = meta.copy()
    cpg_cols = [c for c in out_meth.columns if str(c).startswith("cg")]
    record: dict[str, Any] = {
        "dropped_cpgs": [],
        "n_nan_cells": 0,
        "orphan_sample_ids": [],
    }

    if apply_drop_cpgs and cpg_cols:
        n_drop = int(round(drop_fraction * len(cpg_cols)))
        if drop_fraction > 0.0 and n_drop < 1:
            n_drop = 1
        n_drop = min(n_drop, max(len(cpg_cols) - 1, 0))
        dropped = sorted(
            rng.choice(np.array(cpg_cols, dtype=object), size=n_drop, replace=False).tolist()
        )
        out_meth = out_meth.drop(columns=dropped)
        record["dropped_cpgs"] = [str(x) for x in dropped]
        cpg_cols = [c for c in out_meth.columns if str(c).startswith("cg")]

    if apply_nans and cpg_cols:
        n_cells = int(len(out_meth) * len(cpg_cols))
        n_nan = int(round(nan_fraction * n_cells))
        if nan_fraction > 0.0 and n_nan < 1:
            n_nan = 1
        n_nan = min(n_nan, n_cells)
        row_idx = rng.integers(0, len(out_meth), size=n_nan)
        col_idx = rng.integers(0, len(cpg_cols), size=n_nan)
        for r, c in zip(row_idx, col_idx, strict=False):
            out_meth.iat[int(r), out_meth.columns.get_loc(cpg_cols[int(c)])] = np.nan
        record["n_nan_cells"] = n_nan

    if apply_orphans and n_orphan > 0:
        ids = out_meth["sample_id"].astype(str).tolist()
        n_keep = max(len(ids) - n_orphan, 0)
        keep_ids = ids[:n_keep]
        orphan_ids = ids[n_keep:]
        out_meta = out_meta.loc[out_meta["sample_id"].astype(str).isin(keep_ids)].copy()
        record["orphan_sample_ids"] = [str(x) for x in orphan_ids]

    return out_meth, out_meta, record


def _prefix_figure_titles(harness: ModuleType) -> None:
    """Prefix panel titles and label on-figure metrics as synthetic."""
    orig_a = harness.plot_panel_a
    orig_b = harness.plot_panel_b
    orig_c = harness.plot_panel_c

    def plot_panel_a(
        ax: Any,
        chronological_age: np.ndarray,
        predicted_age: np.ndarray,
        metrics: dict[str, Any],
    ) -> None:
        orig_a(ax, chronological_age, predicted_age, metrics)
        ax.set_title(f"{TITLE_PREFIX}{ax.get_title()}")
        for text in ax.texts:
            body = text.get_text()
            if "MAE" in body and METRIC_DISCLAIMER not in body:
                text.set_text(f"{body}\n({METRIC_DISCLAIMER})")

    def plot_panel_b(ax: Any, chronological_age: np.ndarray, residual: np.ndarray) -> None:
        orig_b(ax, chronological_age, residual)
        ax.set_title(f"{TITLE_PREFIX}{ax.get_title()}")

    def plot_panel_c(ax: Any, weights: pd.Series, gene_map: dict[str, str], top_n: int) -> None:
        orig_c(ax, weights, gene_map, top_n)
        ax.set_title(f"{TITLE_PREFIX}{ax.get_title()}")

    harness.plot_panel_a = plot_panel_a
    harness.plot_panel_b = plot_panel_b
    harness.plot_panel_c = plot_panel_c


def _write_cohort(meth: pd.DataFrame, meta: pd.DataFrame, dest: Path) -> tuple[Path, Path]:
    """Write parquet + CSV under ``dest`` and return those paths."""
    dest.mkdir(parents=True, exist_ok=True)
    meth_path = dest / "meth.parquet"
    meta_path = dest / "meta.csv"
    meth.to_parquet(meth_path, index=False)
    meta.to_csv(meta_path, index=False)
    return meth_path, meta_path


def run_harness_on_cohort(
    harness: ModuleType,
    *,
    meth_path: Path,
    meta_path: Path,
    model_path: Path,
    metrics_path: Path,
    figure_stem: Path,
    mode_name: str,
    missing_cpg_behaviour: str,
) -> ModeResult:
    """Run the existing evaluator and record handled vs crashed.

    Args:
        harness: Loaded evaluation module.
        meth_path: Beta matrix parquet.
        meta_path: Phenotype CSV.
        model_path: Clock artifact.
        metrics_path: Destination JSON.
        figure_stem: Figure path stem.
        mode_name: Failure-mode label for the report.
        missing_cpg_behaviour: Precomputed imputation description.

    Returns:
        Outcome record (never raises to the caller).
    """
    warning_messages: list[str] = []
    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            result = harness.run_validation(
                methylation_path=meth_path,
                meta_path=meta_path,
                model_path=model_path,
                metrics_path=metrics_path,
                figure_stem=figure_stem,
                annotation_path=None,
                top_n_cpgs=int(getattr(harness, "TOP_N_CPGS", 25)),
                skip_manifest_check=True,
                allow_positional_align=False,
            )
            warning_messages = [str(w.message) for w in caught]
        n_imputed = int(result.get("n_imputed_missing_cpgs", 0))
        n_scored = int(result.get("n_samples", 0))
        imputed_ids = [
            line.split("'")[1]
            for line in warning_messages
            if "expected by the model is absent" in line and "'" in line
        ]
        if n_imputed > 0:
            behaviour = "mean-imputed"
            detail = (
                f"Scored {n_scored} samples; mean-imputed {n_imputed} absent clock "
                f"CpG(s) {imputed_ids or '(ids in warnings)'} "
                "(training SimpleImputer column means; not dropped, not errored)."
            )
        elif "(b)" in mode_name:
            behaviour = "n/a (none missing)"
            detail = (
                f"Scored {n_scored} samples. No clock CpG was dropped. "
                "NaN cells were mean-imputed at predict time by the training "
                "SimpleImputer (not dropped, not errored)."
            )
        elif "(c)" in mode_name:
            behaviour = "n/a (none missing)"
            detail = (
                f"Inner join on sample_id dropped matrix rows absent from "
                f"metadata; scored {n_scored} of {N_SAMPLES} samples. "
                "No clock CpG was missing."
            )
        else:
            behaviour = "n/a (none missing)"
            detail = (
                f"Scored {n_scored} samples; every clock CpG was present. "
                "Cell-level NaNs, if any, were left for the training-mean "
                "SimpleImputer. Unmatched sample IDs, if any, were dropped by "
                "inner join on sample_id."
            )
        return ModeResult(
            name=mode_name,
            outcome="handled",
            missing_cpg_behaviour=behaviour,
            detail=detail,
            n_samples_scored=n_scored,
            n_imputed_missing_cpgs=n_imputed,
            imputed_cpg_ids=imputed_ids,
            warning_messages=warning_messages,
        )
    except Exception as exc:
        return ModeResult(
            name=mode_name,
            outcome="crashed",
            missing_cpg_behaviour="errored",
            detail=f"{type(exc).__name__}: {exc}",
            error="".join(traceback.format_exception(type(exc), exc, exc.__traceback__)),
            warning_messages=warning_messages,
        )


def inspect_alignment(
    harness: ModuleType,
    meth_path: Path,
    meta_path: Path,
    model: Any,
) -> dict[str, Any]:
    """Join the cohort and align features the same way the harness does.

    Args:
        harness: Evaluation module.
        meth_path: Beta matrix parquet.
        meta_path: Phenotype CSV.
        model: Loaded clock.

    Returns:
        Alignment diagnostics (sample overlap, imputed CpGs, NaN cells).
    """
    wide = harness.load_validation_cohort(meth_path, meta_path, allow_positional_align=False)
    x, imputed = build_feature_matrix(wide, model)
    n_nan = int(x.isna().to_numpy().sum())
    return {
        "n_rows_after_id_join": int(len(wide)),
        "n_cg_columns_in_joined_table": int(sum(str(c).startswith("cg") for c in wide.columns)),
        "n_model_features": int(x.shape[1]),
        "imputed_missing_cpgs": [str(c) for c in imputed],
        "n_nan_cells_in_aligned_X": n_nan,
    }


def write_generation_params(path: Path, payload: dict[str, Any]) -> None:
    """Write GENERATION_PARAMS.json."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def write_run_report(
    path: Path, *, results: list[ModeResult], alignment: dict[str, Any], params_path: Path
) -> None:
    """Write the pass/fail table. No unlabeled accuracy metrics."""
    lines = [
        f"# {FILE_PREFIX}clock harness synthetic preflight",
        "",
        "This file is a pass/fail record of how the existing evaluator treats",
        "failure modes that the real external methylation matrix is expected to",
        "have. It is not a clock validation result.",
        "",
        "## Failure modes",
        "",
        "| Failure mode | Outcome | Missing-CpG behaviour | What the harness did |",
        "|--------------|---------|----------------------|----------------------|",
    ]
    for row in results:
        what = row.detail.replace("|", "/")
        behaviour = row.missing_cpg_behaviour.replace("|", "/")
        lines.append(f"| {row.name} | {row.outcome} | {behaviour} | {what} |")
    lines.extend(
        [
            "",
            "## Missing CpG imputation (will bias the real run)",
            "",
            "The evaluator does **not drop** clock CpGs that are absent from the",
            "beta matrix, and it does **not error**. It **mean-imputes** them.",
            "",
            "Named behaviour: **training-column-mean imputation** via the fitted",
            "`sklearn.impute.SimpleImputer(strategy='mean')` step on the clock",
            "pipeline. `rogen_aging.clock.evaluate.build_feature_matrix` inserts a",
            "column of NaN for each missing probe and leaves those NaNs for the",
            "pipeline imputer, which fills them from training `statistics_`",
            "(the per-CpG mean on the training matrix). This is **not** a test-set",
            "mean and **not** a row mean. On a real external dataset, every missing",
            "clock probe is replaced by the training-cohort mean methylation for",
            "that site, which shrinks predicted age toward the training intercept.",
            "",
            "Cell-level NaNs in CpGs that *are* present are handled the same way:",
            "left as NaN and filled by the same training-mean SimpleImputer.",
            "",
            "Sample IDs present in the matrix but absent from metadata are dropped",
            "by an inner join on `sample_id` (`load_validation_cohort`). They are",
            "not scored and they do not crash the run when any IDs still overlap.",
            "",
            "## Alignment diagnostics (synthetic data, not a performance estimate)",
            "",
            f"- Rows after ID join: {alignment.get('n_rows_after_id_join')}",
            f"- Clock CpGs mean-imputed because they were absent: "
            f"{alignment.get('imputed_missing_cpgs')}",
            f"- NaN cells left in the aligned matrix for the training imputer: "
            f"{alignment.get('n_nan_cells_in_aligned_X')}",
            "",
            f"Generating parameters: `{params_path.as_posix()}`.",
            "",
        ]
    )
    crashed = [r for r in results if r.outcome == "crashed"]
    if crashed:
        lines.append("## Crash traces")
        lines.append("")
        for row in crashed:
            lines.append(f"### {row.name}")
            lines.append("")
            lines.append("```")
            lines.append((row.error or row.detail).rstrip())
            lines.append("```")
            lines.append("")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def _imputer_means(model: Any, n_features: int) -> np.ndarray:
    """Return training imputer means, or 0.5 if the model has none."""
    if hasattr(model, "named_steps"):
        for step in model.named_steps.values():
            stats = getattr(step, "statistics_", None)
            if stats is not None:
                arr = np.asarray(stats, dtype=float)
                if arr.shape[0] == n_features:
                    return arr
    return np.full(n_features, 0.5, dtype=float)


def _intercept(model: Any) -> float:
    """Return the ElasticNet intercept, or 0.0 if missing."""
    est = model
    if hasattr(model, "named_steps"):
        est = list(model.named_steps.values())[-1]
    value = getattr(est, "intercept_", 0.0)
    return float(np.ravel(value)[0])


def main(
    seed: int = typer.Option(
        DEFAULT_SEED, "--seed", help="RNG seed for ages, betas, and failure-mode injection."
    ),
) -> None:
    """Write a synthetic beta matrix, run the clock harness, and record pass/fail."""
    if not HARNESS_PATH.is_file():
        typer.echo(
            "STOP: no script that takes a beta matrix and produces age predictions.\n"
            f"Searched: {HARNESS_PATH}\n"
            "Also searched: scripts/clock/validate_clock.py, "
            "src/rogen_aging/clock/evaluate.py (evaluate_clock). "
            "Will not write a replacement clock.",
            err=True,
        )
        raise typer.Exit(code=1)

    coefficient_path = resolve_coefficient_file()
    harness = _load_harness()
    _prefix_figure_titles(harness)
    model, coefficients = load_clock_coefficients(coefficient_path, harness)
    cpg_ids = [str(c) for c in coefficients.index]
    intercept = _intercept(model)
    train_means = _imputer_means(model, len(cpg_ids))
    missing_cpg_behaviour = describe_missing_cpg_behaviour(model)

    rng = np.random.default_rng(seed)
    ages = rng.uniform(AGE_MIN, AGE_MAX, size=N_SAMPLES)
    sample_ids = [f"SYN{i:04d}" for i in range(N_SAMPLES)]
    betas = simulate_beta_matrix(
        ages,
        coefficients,
        intercept=intercept,
        train_means=train_means,
        rng=rng,
    )
    complete_meth = pd.concat(
        [pd.DataFrame({"sample_id": sample_ids}), betas.reset_index(drop=True)],
        axis=1,
    )
    complete_meta = pd.DataFrame({"sample_id": sample_ids, "chronological_age": ages})

    # Independent RNGs so isolated probes do not consume the combined-run stream.
    drop_rng = np.random.default_rng(seed + 1)
    nan_rng = np.random.default_rng(seed + 2)
    orphan_rng = np.random.default_rng(seed + 3)
    combined_rng = np.random.default_rng(seed + 4)

    SYNTHETIC_DIR.mkdir(parents=True, exist_ok=True)
    leftover_scratch = SYNTHETIC_DIR / "_scratch"
    if leftover_scratch.exists():
        shutil.rmtree(leftover_scratch)

    isolated_specs: list[tuple[str, bool, bool, bool, np.random.Generator]] = [
        ("(a) drop 5% of clock CpGs entirely", True, False, False, drop_rng),
        ("(b) set 2% of values to NaN", False, True, False, nan_rng),
        ("(c) 20 sample IDs absent from metadata", False, False, True, orphan_rng),
    ]

    results: list[ModeResult] = []
    isolated_records: dict[str, Any] = {}
    with tempfile.TemporaryDirectory(prefix="synthetic_clock_") as tmp:
        tmp_path = Path(tmp)
        for name, drop, nans, orphans, mode_rng in isolated_specs:
            meth_i, meta_i, rec = inject_failure_modes(
                complete_meth,
                complete_meta,
                rng=mode_rng,
                drop_fraction=DROP_CPG_FRACTION,
                nan_fraction=NAN_VALUE_FRACTION,
                n_orphan=N_ORPHAN_SAMPLES,
                apply_drop_cpgs=drop,
                apply_nans=nans,
                apply_orphans=orphans,
            )
            isolated_records[name] = rec
            slug = name[1:4].strip(") ").replace(" ", "_")
            dest = tmp_path / slug
            meth_path, meta_path = _write_cohort(meth_i, meta_i, dest)
            results.append(
                run_harness_on_cohort(
                    harness,
                    meth_path=meth_path,
                    meta_path=meta_path,
                    model_path=coefficient_path,
                    metrics_path=dest / "metrics.json",
                    figure_stem=dest / "figure",
                    mode_name=name,
                    missing_cpg_behaviour=missing_cpg_behaviour,
                )
            )

    combined_meth, combined_meta, combined_record = inject_failure_modes(
        complete_meth,
        complete_meta,
        rng=combined_rng,
        drop_fraction=DROP_CPG_FRACTION,
        nan_fraction=NAN_VALUE_FRACTION,
        n_orphan=N_ORPHAN_SAMPLES,
        apply_drop_cpgs=True,
        apply_nans=True,
        apply_orphans=True,
    )

    meth_out = SYNTHETIC_DIR / f"{FILE_PREFIX}beta_matrix.parquet"
    meta_out = SYNTHETIC_DIR / f"{FILE_PREFIX}metadata.csv"
    metrics_out = SYNTHETIC_DIR / f"{FILE_PREFIX}clock_metrics.json"
    figure_stem = SYNTHETIC_DIR / f"{FILE_PREFIX}Figure_Epigenetic_Clock_Panels"
    params_out = SYNTHETIC_DIR / "GENERATION_PARAMS.json"
    report_out = SYNTHETIC_DIR / "SYNTHETIC_RUN_REPORT.md"

    combined_meth.to_parquet(meth_out, index=False)
    combined_meta.to_csv(meta_out, index=False)

    params = {
        "disclaimer": METRIC_DISCLAIMER,
        "seed": seed,
        "n_samples": N_SAMPLES,
        "age_min": AGE_MIN,
        "age_max": AGE_MAX,
        "age_distribution": "uniform on [age_min, age_max)",
        "coefficient_file": str(coefficient_path.relative_to(REPO_ROOT)),
        "n_clock_cpgs": len(cpg_ids),
        "cpg_ids": cpg_ids,
        "coefficients": {k: float(v) for k, v in coefficients.items()},
        "intercept": intercept,
        "imputer_training_means": {k: float(v) for k, v in zip(cpg_ids, train_means, strict=True)},
        "beta_generation": {
            "method": "unstructured_normal_around_training_means then rank-1 adjustment along the coefficient vector so the linear clock matches age plus fixed Gaussian residual, then clip to (0.02, 0.98)",
            "matrix_layout": "samples as rows; one column per CpG in the clock coefficient set",
            "age_residual_sd_years": AGE_RESIDUAL_SD_YEARS,
            "beta_baseline_sd": BETA_BASELINE_SD,
            "clip_lo": BETA_CLIP_LO,
            "clip_hi": BETA_CLIP_HI,
            "tuned_to_improve_fit": False,
        },
        "failure_modes": {
            "drop_cpg_fraction": DROP_CPG_FRACTION,
            "nan_value_fraction": NAN_VALUE_FRACTION,
            "n_samples_absent_from_metadata": N_ORPHAN_SAMPLES,
            "combined": combined_record,
            "isolated": isolated_records,
        },
        "outputs": {
            "beta_matrix": str(meth_out.relative_to(REPO_ROOT)),
            "metadata": str(meta_out.relative_to(REPO_ROOT)),
            "metrics": str(metrics_out.relative_to(REPO_ROOT)),
            "figure_stem": str(figure_stem.relative_to(REPO_ROOT)),
            "report": str(report_out.relative_to(REPO_ROOT)),
        },
    }
    write_generation_params(params_out, params)

    combined_result = run_harness_on_cohort(
        harness,
        meth_path=meth_out,
        meta_path=meta_out,
        model_path=coefficient_path,
        metrics_path=metrics_out,
        figure_stem=figure_stem,
        mode_name="combined (a)+(b)+(c)",
        missing_cpg_behaviour=missing_cpg_behaviour,
    )
    results.append(combined_result)

    if metrics_out.is_file():
        payload = json.loads(metrics_out.read_text(encoding="utf-8"))
        payload["disclaimer"] = METRIC_DISCLAIMER
        metrics_out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    alignment: dict[str, Any]
    try:
        alignment = inspect_alignment(harness, meth_out, meta_out, model)
    except Exception as exc:
        alignment = {"error": f"{type(exc).__name__}: {exc}"}

    write_run_report(report_out, results=results, alignment=alignment, params_path=params_out)

    typer.echo(f"Coefficient file: {coefficient_path}")
    typer.echo(f"n clock CpGs: {len(cpg_ids)}")
    typer.echo(f"Matrix: {meth_out}")
    typer.echo(f"Metadata: {meta_out}")
    typer.echo(f"Report: {report_out}")
    for row in results:
        typer.echo(f"  {row.name}: {row.outcome}")

    if any(r.outcome == "crashed" for r in results):
        raise typer.Exit(code=1)


if __name__ == "__main__":
    typer.run(main)
