"""Regression: refactored clock package matches pre-refactor GSE40279 training metrics."""

from __future__ import annotations

import json
from pathlib import Path

import joblib
import pytest
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from sklearn.impute import SimpleImputer
from sklearn.linear_model import ElasticNetCV
from sklearn.metrics import mean_absolute_error, root_mean_squared_error
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline

from rogen_aging.clock.data import load_wide_table, split_features_target
from rogen_aging.clock.train import train_clock

REPO_ROOT = Path(__file__).resolve().parent.parent
MOCK_WIDE = REPO_ROOT / "test_data" / "mock_clock_wide.csv"


def _legacy_train_metrics(
    input_data: Path,
    *,
    test_size: float = 0.2,
    random_state: int = 42,
) -> dict[str, float]:
    """Pre-refactor training logic from ``scripts/train_clock_on_gse40279.py``."""
    df = load_wide_table(input_data)
    x, y = split_features_target(df)
    y = pd.to_numeric(y, errors="coerce")
    if y.isna().any():
        raise ValueError("chronological_age contains non-numeric or missing values")

    x_train, x_test, y_train, y_test = train_test_split(
        x,
        y,
        test_size=test_size,
        random_state=random_state,
    )
    enet = ElasticNetCV(
        l1_ratio=[0.1, 0.5, 0.7, 0.9, 0.95, 0.99, 1.0],
        n_alphas=20,
        cv=10,
        random_state=random_state,
        max_iter=5000,
    )
    pipe: Pipeline = Pipeline([("imputer", SimpleImputer(strategy="mean")), ("elasticnet", enet)])
    pipe.fit(x_train, y_train)
    y_pred = pipe.predict(x_test)
    mae = float(mean_absolute_error(y_test, y_pred))
    rmse = float(root_mean_squared_error(y_test, y_pred))
    r, _ = pearsonr(y_test, y_pred)
    return {
        "test_mae": mae,
        "test_rmse": rmse,
        "test_pearson_r": float(r),
        "alpha": float(pipe.named_steps["elasticnet"].alpha_),
        "l1_ratio": float(pipe.named_steps["elasticnet"].l1_ratio_),
    }


def test_train_clock_matches_legacy_metrics(tmp_path: Path) -> None:
    assert MOCK_WIDE.is_file(), f"Missing fixture: {MOCK_WIDE}"
    legacy = _legacy_train_metrics(MOCK_WIDE, test_size=0.2, random_state=42)

    model_path = tmp_path / "clock.pkl"
    metrics_path = tmp_path / "train_metrics.json"
    refactored = train_clock(
        MOCK_WIDE,
        model_path,
        metrics_path,
        test_size=0.2,
        random_state=42,
    )

    assert model_path.is_file()
    assert metrics_path.is_file()
    joblib.load(model_path)

    for key in ("test_mae", "test_rmse", "test_pearson_r", "alpha", "l1_ratio"):
        assert refactored[key] == pytest.approx(legacy[key], rel=0.0, abs=1e-9)

    on_disk = json.loads(metrics_path.read_text(encoding="utf-8"))
    assert on_disk["test_mae"] == pytest.approx(legacy["test_mae"], rel=0.0, abs=1e-9)
