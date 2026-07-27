"""Evaluate a saved epigenetic clock on held-out methylation + age data."""

from __future__ import annotations

import json
import pickle
import warnings
from pathlib import Path
from typing import Any, cast

import joblib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import pearsonr
from sklearn.base import BaseEstimator
from sklearn.metrics import mean_absolute_error

from rogen_aging.clock.data import load_wide_table


def load_model(model_path: Path) -> Any:
    """Load a trained clock from ``.joblib``, ``.pkl``, or legacy pickle paths.

    Args:
        model_path: Path to a serialized sklearn estimator or pipeline.

    Returns:
        The deserialized model object.

    Raises:
        FileNotFoundError: If ``model_path`` is not a file.
    """
    if not model_path.is_file():
        raise FileNotFoundError(f"Model file not found: {model_path}")
    suffix = model_path.suffix.lower()
    if suffix == ".joblib":
        return joblib.load(model_path)
    if suffix in (".pkl", ".pickle"):
        with model_path.open("rb") as handle:
            return pickle.load(handle)
    try:
        return joblib.load(model_path)
    except Exception:
        with model_path.open("rb") as handle:
            return pickle.load(handle)


def _cg_feature_columns(df: pd.DataFrame) -> list[str]:
    return [c for c in df.columns if str(c).startswith("cg")]


def _extract_feature_names_in(model: Any) -> list[str] | None:
    """Return feature names if the fitted estimator recorded them.

    Args:
        model: Fitted estimator or pipeline that may expose
            ``feature_names_in_``.

    Returns:
        Feature name strings, or ``None`` if unavailable.
    """
    if hasattr(model, "feature_names_in_"):
        names = getattr(model, "feature_names_in_", None)
        if names is not None:
            return [str(x) for x in names]
    if hasattr(model, "named_steps"):
        steps = getattr(model, "named_steps", {})
        for _name, step in reversed(list(steps.items())):
            if hasattr(step, "feature_names_in_"):
                names = getattr(step, "feature_names_in_", None)
                if names is not None:
                    return [str(x) for x in names]
    return None


def _n_features_in(model: Any) -> int | None:
    if hasattr(model, "n_features_in_"):
        n = getattr(model, "n_features_in_", None)
        if n is not None:
            return int(n)
    if hasattr(model, "named_steps"):
        for _name, step in reversed(list(getattr(model, "named_steps").items())):
            if hasattr(step, "n_features_in_"):
                n = getattr(step, "n_features_in_", None)
                if n is not None:
                    return int(n)
    return None


def _has_fitted_imputer(model: Any) -> bool:
    """Return True when ``model`` is a Pipeline with a fitted imputer step.

    Args:
        model: Fitted estimator or sklearn Pipeline.

    Returns:
        Whether a step exposes ``statistics_`` (e.g. ``SimpleImputer``).
    """
    if not hasattr(model, "named_steps"):
        return False
    for step in getattr(model, "named_steps").values():
        if hasattr(step, "statistics_"):
            return True
    return False


def _imputer_statistics(model: Any) -> np.ndarray | None:
    """Return fitted imputer ``statistics_`` aligned to model features, if any.

    Args:
        model: Fitted estimator or sklearn Pipeline.

    Returns:
        1-D array of fill values, or ``None``.
    """
    if not hasattr(model, "named_steps"):
        return None
    for step in getattr(model, "named_steps").values():
        stats = getattr(step, "statistics_", None)
        if stats is not None:
            return np.asarray(stats, dtype=float)
    return None


def build_feature_matrix(
    df: pd.DataFrame,
    model: Any,
) -> tuple[pd.DataFrame, list[str]]:
    """Align test CpGs to training features.

    When ``model`` is a Pipeline with a fitted imputer, missing values are left
    as NaN so ``model.predict`` applies training imputation statistics. Bare
    estimators without an imputer are filled from training ``statistics_`` when
    available, otherwise from test-set column / global means.

    Args:
        df: Wide test table with ``chronological_age`` and ``cg*`` columns.
        model: Fitted clock estimator or pipeline (preferably with
            ``feature_names_in_``).

    Returns:
        A pair ``(X, imputed_names)`` where ``X`` is the aligned feature
        matrix and ``imputed_names`` lists training sites absent from ``df``.

    Raises:
        ValueError: If ``chronological_age`` or ``cg*`` columns are missing,
            or feature count cannot be reconciled without names.
    """
    if "chronological_age" not in df.columns:
        raise ValueError("Test data must include a 'chronological_age' column.")

    cg_cols = _cg_feature_columns(df)
    if not cg_cols:
        raise ValueError("No feature columns starting with 'cg' were found in the test data.")

    expected = _extract_feature_names_in(model)
    imputed: list[str] = []
    defer_to_pipeline = _has_fitted_imputer(model)
    train_stats = _imputer_statistics(model)

    if expected is not None:
        x = pd.DataFrame(index=df.index)
        present_cg = df.reindex(columns=cg_cols).apply(pd.to_numeric, errors="coerce")
        flat_mean = float(np.nanmean(present_cg.to_numpy(dtype=float))) if cg_cols else 0.5
        if not np.isfinite(flat_mean):
            flat_mean = 0.5

        for idx, name in enumerate(expected):
            train_fill = (
                float(train_stats[idx])
                if train_stats is not None and idx < len(train_stats) and np.isfinite(train_stats[idx])
                else None
            )
            if name in df.columns:
                col = pd.to_numeric(df[name], errors="coerce")
                if defer_to_pipeline:
                    x[name] = col
                elif train_fill is not None:
                    x[name] = col.fillna(train_fill)
                else:
                    fill = float(np.nanmean(col.to_numpy())) if col.notna().any() else flat_mean
                    if not np.isfinite(fill):
                        fill = flat_mean
                    x[name] = col.fillna(fill)
            else:
                fill = train_fill if train_fill is not None else flat_mean
                warnings.warn(
                    f"CpG '{name}' expected by the model is absent from test data; "
                    + (
                        "leaving NaN for pipeline imputer."
                        if defer_to_pipeline
                        else f"filling with training/global mean ({fill:.6g})."
                    ),
                    stacklevel=2,
                )
                imputed.append(name)
                x[name] = np.nan if defer_to_pipeline else fill

        return x, imputed

    n_feat = _n_features_in(model)
    if n_feat is not None and len(cg_cols) != n_feat:
        raise ValueError(
            f"Model expects {n_feat} features but found {len(cg_cols)} 'cg*' columns, "
            "and the model has no feature_names_in_ to align or impute. "
            "Export the training feature list or refit with a DataFrame so names are stored."
        )

    x = df.reindex(columns=cg_cols).apply(pd.to_numeric, errors="coerce")
    if defer_to_pipeline:
        return x, imputed
    row_mean = x.mean(axis=1)
    x = x.T.fillna(row_mean).T
    col_mean = x.mean(axis=0)
    x = x.fillna(col_mean)
    x = x.fillna(0.5)
    return x, imputed


def assign_age_decade(ages: pd.Series) -> pd.Series:
    """Bin chronological ages into labeled decade intervals.

    Args:
        ages: Chronological ages in years.

    Returns:
        Categorical decade labels (``<20``, ``20-29``, …, ``90+``).
    """
    bins = [-np.inf, 20, 30, 40, 50, 60, 70, 80, 90, np.inf]
    labels = ["<20", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80-89", "90+"]
    return pd.cut(ages.astype(float), bins=bins, labels=labels, right=False)


def plot_residuals(age: np.ndarray, residual: np.ndarray, out_path: Path) -> None:
    """Save a scatter of residuals vs chronological age.

    Args:
        age: Chronological ages (years).
        residual: Predicted age minus chronological age (years).
        out_path: Destination PNG path.
    """
    fig, ax = plt.subplots(figsize=(7.0, 5.5))
    ax.scatter(age, residual, alpha=0.75, edgecolors="black", linewidths=0.25, s=36)
    ax.axhline(0.0, color="crimson", linewidth=1.2, linestyle="--", label="Zero residual")
    ax.set_xlabel("Chronological age (years)")
    ax.set_ylabel("Residual (predicted − chronological, years)")
    ax.set_title("Clock residuals vs chronological age")
    ax.grid(True, linestyle=":", linewidth=0.6, alpha=0.75)
    ax.legend(loc="best", frameon=True)
    fig.tight_layout()
    fig.savefig(out_path, dpi=160)
    plt.close(fig)


def plot_mae_by_decade(decade_df: pd.DataFrame, out_path: Path) -> None:
    """Save a bar chart of MAE by age decade.

    Args:
        decade_df: Table with ``decade`` and ``mae`` columns.
        out_path: Destination PNG path.
    """
    label_order = ["<20", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80-89", "90+"]
    seen = set(decade_df["decade"].astype(str))
    plot_order = [lab for lab in label_order if lab in seen]
    fig, ax = plt.subplots(figsize=(8.0, 5.0))
    sns.barplot(
        data=decade_df,
        x="decade",
        y="mae",
        order=plot_order,
        ax=ax,
        color="#4C72B0",
        edgecolor="black",
    )
    ax.set_xlabel("Age decade")
    ax.set_ylabel("Mean absolute error (years)")
    ax.set_title("MAE by chronological age decade")
    ax.grid(True, axis="y", linestyle=":", linewidth=0.6, alpha=0.75)
    fig.tight_layout()
    fig.savefig(out_path, dpi=160)
    plt.close(fig)


def evaluate_clock(
    model_path: Path,
    test_data: Path,
    output_dir: Path,
) -> dict[str, Any]:
    """Run held-out evaluation; write figures and ``validation_metrics.json``.

    Args:
        model_path: Path to a saved clock (``.joblib`` / pickle).
        test_data: Path to a wide test table with ages and ``cg*`` columns.
        output_dir: Directory for metrics JSON and residual/MAE figures.

    Returns:
        Metrics dictionary including overall MAE, Pearson r, decade MAE, and
        ``metrics_path`` pointing at the written JSON file.

    Raises:
        FileNotFoundError: If the model or test data path is missing.
        ValueError: If ages are entirely non-numeric/missing or features
            cannot be aligned.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    model = load_model(model_path)
    df = load_wide_table(test_data)

    y = pd.to_numeric(df["chronological_age"], errors="coerce")
    if y.isna().all():
        raise ValueError("chronological_age is entirely non-numeric or missing.")
    valid_mask = y.notna()
    if not bool(valid_mask.all()):
        dropped = int((~valid_mask).sum())
        warnings.warn(f"Dropping {dropped} rows with invalid chronological_age.", stacklevel=2)
    df = df.loc[valid_mask].copy()
    y = y.loc[valid_mask]

    x, imputed_names = build_feature_matrix(df, model)
    estimator = cast(BaseEstimator, model)
    y_pred = np.asarray(estimator.predict(x), dtype=float)

    mae_overall = float(mean_absolute_error(y.to_numpy(), y_pred))
    r_value, r_p = pearsonr(y.to_numpy(dtype=float), y_pred.astype(float))

    residual = y_pred - y.to_numpy(dtype=float)
    decades = assign_age_decade(y)
    eval_df = pd.DataFrame(
        {
            "chronological_age": y.to_numpy(),
            "predicted_age": y_pred,
            "residual": residual,
            "decade": decades,
        }
    )
    decade_mae = eval_df.groupby("decade", observed=True)["residual"].apply(
        lambda s: float(np.mean(np.abs(s.to_numpy(dtype=float))))
    )
    decade_table = decade_mae.rename("mae").reset_index()

    metrics: dict[str, Any] = {
        "mae_overall": mae_overall,
        "pearson_r": float(r_value),
        "pearson_p": float(r_p),
        "n_samples": int(len(y)),
        "n_features_used": int(x.shape[1]),
        "imputed_missing_cpgs": imputed_names,
        "mae_by_decade": decade_mae.astype(float).to_dict(),
    }
    metrics_path = output_dir / "validation_metrics.json"
    metrics_path.write_text(json.dumps(metrics, indent=2), encoding="utf-8")

    plot_residuals(y.to_numpy(dtype=float), residual, output_dir / "Fig_Clock_Residuals.png")
    plot_mae_by_decade(decade_table, output_dir / "Fig_Clock_MAE_by_decade.png")

    return {**metrics, "metrics_path": str(metrics_path)}
