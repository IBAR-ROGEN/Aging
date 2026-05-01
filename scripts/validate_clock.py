#!/usr/bin/env python3
"""Validate a pre-trained epigenetic clock (e.g. ElasticNet) on held-out bedMethyl-style data.

Loads a pickled/joblib model, aligns CpG features (with optional mean imputation for
missing sites), reports MAE and Pearson r, stratifies MAE by age decade, and saves
residual and decade MAE figures under ``--output_dir``.
"""

from __future__ import annotations

import argparse
import json
import pickle
import sys
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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Evaluate a pre-trained epigenetic clock on held-out methylation + age data.",
    )
    parser.add_argument(
        "--model_path",
        type=Path,
        required=True,
        help="Path to a trained model (.pkl or .joblib).",
    )
    parser.add_argument(
        "--test_data",
        type=Path,
        required=True,
        help="Path to test table (.parquet or .csv) with chronological_age and cg* CpG columns.",
    )
    parser.add_argument(
        "--output_dir",
        type=Path,
        required=True,
        help="Directory for figures and metrics JSON (created if missing).",
    )
    return parser.parse_args()


def load_model(model_path: Path) -> Any:
    if not model_path.is_file():
        raise FileNotFoundError(f"Model file not found: {model_path}")
    suffix = model_path.suffix.lower()
    if suffix == ".joblib":
        return joblib.load(model_path)
    if suffix == ".pkl" or suffix == ".pickle":
        with model_path.open("rb") as handle:
            return pickle.load(handle)
    # Best-effort: joblib often works for extensionless or legacy paths
    try:
        return joblib.load(model_path)
    except Exception:
        with model_path.open("rb") as handle:
            return pickle.load(handle)


def load_test_table(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"Test data file not found: {path}")
    suffix = path.suffix.lower()
    if suffix == ".parquet":
        return pd.read_parquet(path)
    if suffix in (".csv", ".tsv"):
        sep = "\t" if suffix == ".tsv" else ","
        return pd.read_csv(path, sep=sep)
    raise ValueError(f"Unsupported test data extension: {path.suffix} (use .parquet or .csv)")


def _cg_feature_columns(df: pd.DataFrame) -> list[str]:
    return [c for c in df.columns if str(c).startswith("cg")]


def _extract_feature_names_in(model: Any) -> list[str] | None:
    """Return feature names if the fitted estimator recorded them (e.g. fit on DataFrame)."""
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


def build_feature_matrix(
    df: pd.DataFrame,
    model: Any,
) -> tuple[pd.DataFrame, list[str]]:
    """Return X aligned to training features and the list of model CpGs filled by imputation."""
    if "chronological_age" not in df.columns:
        raise ValueError("Test data must include a 'chronological_age' column.")

    cg_cols = _cg_feature_columns(df)
    if not cg_cols:
        raise ValueError("No feature columns starting with 'cg' were found in the test data.")

    expected = _extract_feature_names_in(model)
    imputed: list[str] = []

    if expected is not None:
        X = pd.DataFrame(index=df.index)
        present_cg = df.reindex(columns=cg_cols).apply(pd.to_numeric, errors="coerce")
        flat_mean = float(np.nanmean(present_cg.to_numpy(dtype=float))) if cg_cols else 0.5
        if not np.isfinite(flat_mean):
            flat_mean = 0.5

        for name in expected:
            if name in df.columns:
                col = pd.to_numeric(df[name], errors="coerce")
                fill = float(np.nanmean(col.to_numpy())) if col.notna().any() else flat_mean
                if not np.isfinite(fill):
                    fill = flat_mean
                X[name] = col.fillna(fill)
            else:
                warnings.warn(
                    f"CpG '{name}' expected by the model is absent from test data; "
                    f"filling with global mean imputation ({flat_mean:.6g}).",
                    stacklevel=2,
                )
                imputed.append(name)
                X[name] = flat_mean

        return X, imputed

    n_feat = _n_features_in(model)
    if n_feat is not None and len(cg_cols) != n_feat:
        raise ValueError(
            f"Model expects {n_feat} features but found {len(cg_cols)} 'cg*' columns, "
            "and the model has no feature_names_in_ to align or impute. "
            "Export the training feature list or refit with a DataFrame so names are stored."
        )

    X = df.reindex(columns=cg_cols).apply(pd.to_numeric, errors="coerce")
    row_mean = X.mean(axis=1)
    X = X.T.fillna(row_mean).T
    col_mean = X.mean(axis=0)
    X = X.fillna(col_mean)
    X = X.fillna(0.5)
    return X, imputed


def assign_age_decade(ages: pd.Series) -> pd.Series:
    bins = [-np.inf, 20, 30, 40, 50, 60, 70, 80, 90, np.inf]
    labels = ["<20", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80-89", "90+"]
    return pd.cut(ages.astype(float), bins=bins, labels=labels, right=False)


def plot_residuals(age: np.ndarray, residual: np.ndarray, out_path: Path) -> None:
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


def main() -> None:
    args = parse_args()
    output_dir: Path = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    model = load_model(cast(Path, args.model_path))
    df = load_test_table(cast(Path, args.test_data))

    y = pd.to_numeric(df["chronological_age"], errors="coerce")
    if y.isna().all():
        raise ValueError("chronological_age is entirely non-numeric or missing.")
    valid_mask = y.notna()
    if not bool(valid_mask.all()):
        dropped = int((~valid_mask).sum())
        warnings.warn(f"Dropping {dropped} rows with invalid chronological_age.", stacklevel=2)
    df = df.loc[valid_mask].copy()
    y = y.loc[valid_mask]

    X, imputed_names = build_feature_matrix(df, model)
    estimator = cast(BaseEstimator, model)
    # Pass DataFrame when the model was fit with feature names so sklearn validates columns.
    y_pred = np.asarray(estimator.predict(X), dtype=float)

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

    metrics = {
        "mae_overall": mae_overall,
        "pearson_r": float(r_value),
        "pearson_p": float(r_p),
        "n_samples": int(len(y)),
        "n_features_used": int(X.shape[1]),
        "imputed_missing_cpgs": imputed_names,
        "mae_by_decade": decade_mae.astype(float).to_dict(),
    }
    metrics_path = output_dir / "validation_metrics.json"
    metrics_path.write_text(json.dumps(metrics, indent=2), encoding="utf-8")

    plot_residuals(y.to_numpy(dtype=float), residual, output_dir / "Fig_Clock_Residuals.png")
    plot_mae_by_decade(decade_table, output_dir / "Fig_Clock_MAE_by_decade.png")

    print(json.dumps({"metrics_path": str(metrics_path), **{k: v for k, v in metrics.items() if k != "imputed_missing_cpgs"}}, indent=2))
    if imputed_names:
        print(f"Imputed {len(imputed_names)} missing model CpGs (see metrics JSON for names).", file=sys.stderr)


if __name__ == "__main__":
    main()
