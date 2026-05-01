#!/usr/bin/env python3
"""Train an ElasticNetCV epigenetic clock from a sample-by-CpG table plus chronological age.

GSE40279 (Hannum 2013, Illumina 450K whole blood) must be obtained from GEO and converted
to a table with rows=samples, CpG columns whose names start with ``cg``, and a
``chronological_age`` column. That download and conversion are not implemented here.
"""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
from typing import Any

import joblib
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from sklearn.impute import SimpleImputer
from sklearn.linear_model import ElasticNetCV
from sklearn.metrics import mean_absolute_error, root_mean_squared_error
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline

logger = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    """CLI: input table, output paths, split and seed."""
    p = argparse.ArgumentParser(description="Train ElasticNetCV clock on methylation + age table.")
    p.add_argument("--input_data", type=Path, required=True, help="Parquet or CSV: samples x (cg* + chronological_age).")
    p.add_argument("--output_model", type=Path, required=True, help="Path to write trained model (.pkl).")
    p.add_argument("--output_metrics", type=Path, required=True, help="Path to write training metrics JSON.")
    p.add_argument("--test_size", type=float, default=0.2, help="Fraction held out for test (default 0.2).")
    p.add_argument("--random_state", type=int, default=42, help="Random seed (default 42).")
    return p.parse_args()


def load_table(path: Path) -> pd.DataFrame:
    """Load Parquet or CSV/TSV into a DataFrame."""
    if not path.is_file():
        raise FileNotFoundError(path)
    suf = path.suffix.lower()
    if suf == ".parquet":
        return pd.read_parquet(path)
    if suf in (".csv", ".tsv"):
        return pd.read_csv(path, sep="\t" if suf == ".tsv" else ",")
    raise ValueError(f"Unsupported extension {path.suffix}; use .parquet or .csv")


def features_and_target(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.Series]:
    """Return CpG feature matrix (``cg*`` columns) and age target."""
    if "chronological_age" not in df.columns:
        raise KeyError("Expected column 'chronological_age'")
    cpg_cols = [c for c in df.columns if c.startswith("cg")]
    if not cpg_cols:
        raise ValueError("No feature columns starting with 'cg'")
    return df[cpg_cols].copy(), df["chronological_age"].copy()


def main() -> None:
    """Train/test split, fit pipeline, write model + metrics JSON, print summary."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    args = parse_args()
    X, y = features_and_target(load_table(args.input_data))
    y = pd.to_numeric(y, errors="coerce")
    if y.isna().any():
        raise ValueError("chronological_age contains non-numeric or missing values; drop or fix rows first.")
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=args.test_size, random_state=args.random_state
    )
    enet = ElasticNetCV(
        l1_ratio=[0.1, 0.5, 0.7, 0.9, 0.95, 0.99, 1.0],
        n_alphas=20,
        cv=10,
        random_state=args.random_state,
        max_iter=5000,
    )
    pipe: Pipeline = Pipeline([("imputer", SimpleImputer(strategy="mean")), ("elasticnet", enet)])
    logger.info("Fitting Pipeline(imputer, ElasticNetCV) …")
    pipe.fit(X_train, y_train)
    y_pred = pipe.predict(X_test)
    mae = float(mean_absolute_error(y_test, y_pred))
    rmse = float(root_mean_squared_error(y_test, y_pred))
    r, _ = pearsonr(y_test, y_pred)
    coef = np.ravel(pipe.named_steps["elasticnet"].coef_)
    names = list(X_train.columns)
    selected = [names[i] for i in np.flatnonzero(coef)]
    args.output_model.parent.mkdir(parents=True, exist_ok=True)
    joblib.dump(pipe, args.output_model)
    metrics: dict[str, Any] = {
        "test_mae": mae,
        "test_rmse": rmse,
        "test_pearson_r": float(r),
        "alpha": float(pipe.named_steps["elasticnet"].alpha_),
        "l1_ratio": float(pipe.named_steps["elasticnet"].l1_ratio_),
        "n_cpgs_features": len(names),
        "n_cpgs_selected_nonzero": len(selected),
        "selected_cpgs": selected,
        "test_size": args.test_size,
        "random_state": args.random_state,
    }
    args.output_metrics.parent.mkdir(parents=True, exist_ok=True)
    args.output_metrics.write_text(json.dumps(metrics, indent=2), encoding="utf-8")
    logger.info(
        "Done: n_cpgs=%d alpha=%.6g l1_ratio=%.4g test_MAE=%.4f test_r=%.4f",
        len(names),
        metrics["alpha"],
        metrics["l1_ratio"],
        mae,
        r,
    )
    print(
        f"# CpGs used: {len(names)} | alpha: {metrics['alpha']:.6g} | l1_ratio: {metrics['l1_ratio']:.4g} | "
        f"test MAE: {mae:.4f} | test r: {r:.4f}"
    )


if __name__ == "__main__":
    main()
