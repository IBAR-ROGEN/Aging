"""Train the Hannum GSE40279 clock from the prepared wide matrix.

Probe filters that use methylation or age are fit on the training split only.
Thresholds below are fixed before any look at GSE87571 performance.
"""

from __future__ import annotations

import hashlib
import json
import logging
import subprocess
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import polars as pl
import sklearn
from sklearn.model_selection import train_test_split

from rogen_aging.clock.external_data import (
    MATRIX1_FILENAME,
    MATRIX2_FILENAME,
    SUPP_MATRIX1_URL,
    SUPP_MATRIX2_URL,
    _ensure_supplementary_matrix,
)
from rogen_aging.clock.train import train_clock
from rogen_aging.config import get_config

logger = logging.getLogger(__name__)

# Fixed before external evaluation. Not adjusted to improve GSE87571 MAE.
MISSING_FRACTION_MAX = 0.05
AGE_ABS_CORRELATION_MIN = 0.2
MAX_PROBES_FOR_ELASTICNET = 8_000
CV_FOLDS_WHEN_WIDE = 5
_BATCH = 2_000
_FIXTURE_PREFIX = "cg_test_"


@dataclass(frozen=True)
class ProbeFilterReport:
    """Probe counts after each pre-specified filter."""

    n_probes_input: int
    n_probes_overlap_gse87571: int
    n_probes_after_missingness: int
    n_probes_after_age_correlation: int
    n_probes_trained: int
    n_train_samples: int
    n_test_samples: int
    columns: list[str]
    cv: int
    cv_note: str


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _git_commit(repo_root: Path) -> str:
    result = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=repo_root,
        check=False,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        return "unknown"
    return result.stdout.strip()


def gse87571_cpg_ids_from_parquet(path: Path) -> list[str] | None:
    """Return cg columns when the parquet is a real 450K matrix, else None.

    The demo file written by ``write_clock_fixtures`` uses ``cg_test_*`` names.
    That header is not the GSE87571 probe set.
    """
    if not path.is_file():
        return None
    names = pl.scan_parquet(path).collect_schema().names()
    cg = [name for name in names if name.startswith("cg")]
    if not cg or all(name.startswith(_FIXTURE_PREFIX) for name in cg):
        return None
    return cg


def _iter_id_ref(path: Path) -> list[str]:
    """Read the first TSV field of each data row. Skips the header."""
    import gzip

    ids: list[str] = []

    def consume(handle: Any) -> None:
        header = handle.readline()
        if not header:
            raise ValueError(f"Empty GSE87571 matrix: {path}")
        first = str(header).split("\t", 1)[0].strip().strip('"')
        if first != "ID_REF":
            raise ValueError(f"Expected ID_REF header in {path}, got {first!r}.")
        for line in handle:
            probe = str(line).split("\t", 1)[0].strip().strip('"')
            if probe.startswith("cg"):
                ids.append(probe)

    if path.suffix.lower() == ".gz":
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
            consume(handle)
    else:
        with path.open("rt", encoding="utf-8", errors="replace") as handle:
            consume(handle)
    return ids


def gse87571_cpg_ids_from_geo(geo_cache_dir: Path) -> list[str]:
    """Read CpG IDs from the GSE87571 supplementary matrices (ID_REF only)."""
    found: set[str] = set()
    for url, filename in (
        (SUPP_MATRIX1_URL, MATRIX1_FILENAME),
        (SUPP_MATRIX2_URL, MATRIX2_FILENAME),
    ):
        path = _ensure_supplementary_matrix(url, filename, geo_cache_dir)
        found.update(_iter_id_ref(path))
    if not found:
        raise ValueError(f"No cg probe IDs found in GSE87571 matrices under {geo_cache_dir}.")
    return sorted(found)


def resolve_gse87571_cpg_ids(processed_parquet: Path, geo_cache_dir: Path) -> list[str]:
    """Use a real processed header when present, otherwise the GEO matrices."""
    from_parquet = gse87571_cpg_ids_from_parquet(processed_parquet)
    if from_parquet is not None:
        return from_parquet
    return gse87571_cpg_ids_from_geo(geo_cache_dir)


def _batch_train_stats(
    values: np.ndarray,
    age: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return missing fraction and Pearson r with age. ``values`` is ``(n, p)``."""
    x = values.astype(np.float64, copy=False)
    y = age.astype(np.float64, copy=False)
    valid = np.isfinite(x) & np.isfinite(y)[:, None]
    n = valid.sum(axis=0).astype(np.float64)
    n_rows = float(x.shape[0])
    missing = 1.0 - (n / n_rows)
    x0 = np.where(valid, x, 0.0)
    y0 = np.where(valid, y[:, None], 0.0)
    sum_x = x0.sum(axis=0)
    sum_y = y0.sum(axis=0)
    sum_xx = np.square(x0).sum(axis=0)
    sum_yy = np.square(y0).sum(axis=0)
    sum_xy = (x0 * y[:, None]).sum(axis=0)
    numer = n * sum_xy - sum_x * sum_y
    denom = np.sqrt((n * sum_xx - np.square(sum_x)) * (n * sum_yy - np.square(sum_y)))
    with np.errstate(invalid="ignore", divide="ignore"):
        correlation = numer / denom
    correlation[n < 10] = np.nan
    return missing, correlation


def select_training_probes(
    frame: pd.DataFrame,
    gse87571_cpgs: set[str],
    *,
    test_size: float,
    random_state: int,
) -> ProbeFilterReport:
    """Intersect, drop sparse probes, then keep age-correlated probes on the train split."""
    cg_columns = [column for column in frame.columns if str(column).startswith("cg")]
    overlap = [column for column in cg_columns if column in gse87571_cpgs]
    if not overlap:
        raise ValueError("No GSE40279 cg columns are present in the GSE87571 probe set.")
    ages = pd.to_numeric(frame["chronological_age"], errors="coerce")
    if ages.isna().any():
        raise ValueError("chronological_age contains missing values; refusing to drop samples.")
    indices = np.arange(len(frame))
    train_idx, test_idx = train_test_split(
        indices,
        test_size=test_size,
        random_state=random_state,
    )
    train_age = ages.to_numpy(dtype=np.float64)[train_idx]
    kept_after_missing: list[str] = []
    correlations: dict[str, float] = {}
    train_frame = frame.iloc[train_idx]
    for start in range(0, len(overlap), _BATCH):
        batch_cols = overlap[start : start + _BATCH]
        values = (
            train_frame.loc[:, batch_cols]
            .apply(pd.to_numeric, errors="coerce")
            .to_numpy(dtype=np.float64)
        )
        missing, correlation = _batch_train_stats(values, train_age)
        for column, miss, corr in zip(batch_cols, missing, correlation, strict=True):
            if float(miss) > MISSING_FRACTION_MAX:
                continue
            kept_after_missing.append(column)
            if np.isfinite(corr) and abs(float(corr)) >= AGE_ABS_CORRELATION_MIN:
                correlations[column] = abs(float(corr))
    correlated = sorted(correlations, key=correlations.__getitem__, reverse=True)
    trained = correlated[:MAX_PROBES_FOR_ELASTICNET]
    if not trained:
        raise ValueError(
            "Age-correlation filter on the training split kept zero probes "
            f"at |r| >= {AGE_ABS_CORRELATION_MIN}."
        )
    cfg_cv = int(get_config().clock.elasticnet.cv)
    if len(trained) > 2_000 and cfg_cv > CV_FOLDS_WHEN_WIDE:
        cv = CV_FOLDS_WHEN_WIDE
        cv_note = (
            f"cv lowered from {cfg_cv} to {cv} because {len(trained)} training probes "
            "would make ElasticNetCV exceed about 2 hours on a laptop. "
            "No other hyperparameter was changed."
        )
    else:
        cv = cfg_cv
        cv_note = f"cv left at the config value {cfg_cv}."
    return ProbeFilterReport(
        n_probes_input=len(cg_columns),
        n_probes_overlap_gse87571=len(overlap),
        n_probes_after_missingness=len(kept_after_missing),
        n_probes_after_age_correlation=len(correlated),
        n_probes_trained=len(trained),
        n_train_samples=int(len(train_idx)),
        n_test_samples=int(len(test_idx)),
        columns=trained,
        cv=cv,
        cv_note=cv_note,
    )


def write_training_matrix(
    source: Path,
    destination: Path,
    columns: list[str],
) -> None:
    """Write the sample-complete matrix restricted to ``columns`` (row order unchanged)."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    frame = pl.read_parquet(source, columns=["sample_id", "chronological_age", *columns])
    frame.write_parquet(destination, compression="zstd")


def select_probes_from_parquet(
    source: Path,
    gse87571_cpgs: set[str],
    *,
    test_size: float,
    random_state: int,
) -> ProbeFilterReport:
    """Filter probes in column batches so the full 450K table is not materialized."""
    schema_names = pl.scan_parquet(source).collect_schema().names()
    cg_columns = [name for name in schema_names if name.startswith("cg")]
    overlap = [name for name in cg_columns if name in gse87571_cpgs]
    if not overlap:
        raise ValueError("No GSE40279 cg columns are present in the GSE87571 probe set.")
    ages = pl.read_parquet(source, columns=["chronological_age"]).get_column("chronological_age")
    if ages.null_count() or not ages.dtype.is_numeric():
        raise ValueError("chronological_age contains missing values; refusing to drop samples.")
    age_np = ages.to_numpy().astype(np.float64, copy=False)
    indices = np.arange(age_np.shape[0])
    train_idx, test_idx = train_test_split(
        indices,
        test_size=test_size,
        random_state=random_state,
    )
    train_age = age_np[train_idx]
    kept_after_missing: list[str] = []
    correlations: dict[str, float] = {}
    logger.info("Overlap with GSE87571: %s probes", len(overlap))
    for start in range(0, len(overlap), _BATCH):
        batch_cols = overlap[start : start + _BATCH]
        logger.info("Probe filter batch %s / %s", start, len(overlap))
        batch = pl.read_parquet(source, columns=batch_cols).to_numpy()
        values = batch[train_idx]
        missing, correlation = _batch_train_stats(values, train_age)
        for column, miss, corr in zip(batch_cols, missing, correlation, strict=True):
            if float(miss) > MISSING_FRACTION_MAX:
                continue
            kept_after_missing.append(column)
            if np.isfinite(corr) and abs(float(corr)) >= AGE_ABS_CORRELATION_MIN:
                correlations[column] = abs(float(corr))
    correlated = sorted(correlations, key=correlations.__getitem__, reverse=True)
    trained = correlated[:MAX_PROBES_FOR_ELASTICNET]
    if not trained:
        raise ValueError(
            "Age-correlation filter on the training split kept zero probes "
            f"at |r| >= {AGE_ABS_CORRELATION_MIN}."
        )
    cfg_cv = int(get_config().clock.elasticnet.cv)
    if len(trained) > 2_000 and cfg_cv > CV_FOLDS_WHEN_WIDE:
        cv = CV_FOLDS_WHEN_WIDE
        cv_note = (
            f"cv lowered from {cfg_cv} to {cv} because {len(trained)} training probes "
            "would make ElasticNetCV exceed about 2 hours on a laptop. "
            "No other hyperparameter was changed."
        )
    else:
        cv = cfg_cv
        cv_note = f"cv left at the config value {cfg_cv}."
    return ProbeFilterReport(
        n_probes_input=len(cg_columns),
        n_probes_overlap_gse87571=len(overlap),
        n_probes_after_missingness=len(kept_after_missing),
        n_probes_after_age_correlation=len(correlated),
        n_probes_trained=len(trained),
        n_train_samples=int(len(train_idx)),
        n_test_samples=int(len(test_idx)),
        columns=trained,
        cv=cv,
        cv_note=cv_note,
    )


def write_coefficient_table(model_path: Path, csv_path: Path) -> int:
    """Write non-zero CpG weights. Returns the non-zero count."""
    import joblib

    model = joblib.load(model_path)
    step = model.named_steps["elasticnet"]
    names = [str(name) for name in model.feature_names_in_]
    coef = np.ravel(step.coef_)
    nonzero = np.flatnonzero(coef)
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    pl.DataFrame(
        {
            "cpg": [names[index] for index in nonzero],
            "weight": coef[nonzero].astype(np.float64),
        }
    ).write_csv(csv_path)
    return int(nonzero.size)


def train_published_gse40279(
    input_parquet: Path,
    *,
    gse87571_processed: Path,
    geo_cache_dir: Path,
    output_model: Path,
    output_metrics: Path,
    coefficient_csv: Path,
    provenance_path: Path,
    training_matrix: Path,
    repo_root: Path,
) -> dict[str, Any]:
    """Filter on the training split, then fit with :func:`train_clock`."""
    cfg = get_config().clock
    test_size = float(cfg.test_size)
    random_state = int(cfg.random_state)
    cpg_ids = resolve_gse87571_cpg_ids(gse87571_processed, geo_cache_dir)
    report = select_probes_from_parquet(
        input_parquet,
        set(cpg_ids),
        test_size=test_size,
        random_state=random_state,
    )
    write_training_matrix(input_parquet, training_matrix, report.columns)
    metrics = train_clock(
        training_matrix,
        output_model,
        output_metrics,
        test_size=test_size,
        random_state=random_state,
        cv=report.cv,
    )
    write_coefficient_table(output_model, coefficient_csv)
    write_provenance(
        provenance_path,
        input_parquet=input_parquet,
        report=report,
        metrics=metrics,
        coefficient_csv=coefficient_csv,
        model_path=output_model,
        repo_root=repo_root,
    )
    return {"metrics": metrics, "probe_filter": report}


def write_provenance(
    path: Path,
    *,
    input_parquet: Path,
    report: ProbeFilterReport,
    metrics: dict[str, Any],
    coefficient_csv: Path,
    model_path: Path,
    repo_root: Path,
) -> dict[str, Any]:
    """Write the sidecar JSON. Metrics values are copied from ``train_clock`` output."""
    payload: dict[str, Any] = {
        "geo_accession": "GSE40279",
        "model_path": str(model_path),
        "input_parquet": str(input_parquet),
        "input_parquet_sha256": _sha256(input_parquet),
        "n_cpgs_features": report.n_probes_trained,
        "n_train_samples": report.n_train_samples,
        "n_test_samples": report.n_test_samples,
        "probe_counts": {
            "input_cg": report.n_probes_input,
            "overlap_gse87571": report.n_probes_overlap_gse87571,
            "after_missingness_le_5pct_on_train": report.n_probes_after_missingness,
            "after_abs_age_correlation_on_train": report.n_probes_after_age_correlation,
            "trained": report.n_probes_trained,
        },
        "filters": {
            "missing_fraction_max": MISSING_FRACTION_MAX,
            "age_abs_correlation_min": AGE_ABS_CORRELATION_MIN,
            "max_probes": MAX_PROBES_FOR_ELASTICNET,
            "correlation_split": "training split only; same test_size and random_state as train_clock",
        },
        "n_nonzero_coefficients": int(metrics["n_cpgs_selected_nonzero"]),
        "nonzero_cpg_csv": str(coefficient_csv),
        "hyperparameters": {
            "alpha": metrics["alpha"],
            "l1_ratio": metrics["l1_ratio"],
            "cv": report.cv,
            "cv_note": report.cv_note,
        },
        "random_state": metrics["random_state"],
        "sklearn_version": sklearn.__version__,
        "git_commit": _git_commit(repo_root),
        "timestamp_utc": datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "holdout_metrics_from_train_clock": {
            "test_mae": metrics["test_mae"],
            "test_rmse": metrics["test_rmse"],
            "test_pearson_r": metrics["test_pearson_r"],
        },
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return payload
