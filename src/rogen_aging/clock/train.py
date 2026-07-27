"""Train an ElasticNet epigenetic clock and persist model + metrics."""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any

import joblib
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from sklearn.metrics import mean_absolute_error, root_mean_squared_error
from sklearn.model_selection import train_test_split

from rogen_aging.clock.data import load_wide_table, split_features_target
from rogen_aging.clock.model import make_clock_pipeline

logger = logging.getLogger(__name__)


def train_clock(
    input_data: Path,
    output_model: Path,
    output_metrics: Path,
    *,
    test_size: float = 0.2,
    random_state: int = 42,
) -> dict[str, Any]:
    """Fit a clock pipeline on a wide table and write model + metrics JSON.

    Returns the metrics dictionary written to ``output_metrics``.
    """
    df = load_wide_table(input_data)
    x, y = split_features_target(df)
    y = pd.to_numeric(y, errors="coerce")
    if y.isna().any():
        raise ValueError("chronological_age contains non-numeric or missing values; drop or fix rows first.")

    x_train, x_test, y_train, y_test = train_test_split(
        x,
        y,
        test_size=test_size,
        random_state=random_state,
    )
    pipe = make_clock_pipeline(random_state=random_state)
    logger.info("Fitting Pipeline(imputer, ElasticNetCV) …")
    pipe.fit(x_train, y_train)
    y_pred = pipe.predict(x_test)

    mae = float(mean_absolute_error(y_test, y_pred))
    rmse = float(root_mean_squared_error(y_test, y_pred))
    r, _ = pearsonr(y_test, y_pred)
    coef = np.ravel(pipe.named_steps["elasticnet"].coef_)
    names = list(x_train.columns)
    selected = [names[i] for i in np.flatnonzero(coef)]

    output_model.parent.mkdir(parents=True, exist_ok=True)
    joblib.dump(pipe, output_model)

    metrics: dict[str, Any] = {
        "test_mae": mae,
        "test_rmse": rmse,
        "test_pearson_r": float(r),
        "alpha": float(pipe.named_steps["elasticnet"].alpha_),
        "l1_ratio": float(pipe.named_steps["elasticnet"].l1_ratio_),
        "n_cpgs_features": len(names),
        "n_cpgs_selected_nonzero": len(selected),
        "selected_cpgs": selected,
        "test_size": test_size,
        "random_state": random_state,
    }
    output_metrics.parent.mkdir(parents=True, exist_ok=True)
    output_metrics.write_text(json.dumps(metrics, indent=2), encoding="utf-8")
    logger.info(
        "Done: n_cpgs=%d alpha=%.6g l1_ratio=%.4g test_MAE=%.4f test_r=%.4f",
        len(names),
        metrics["alpha"],
        metrics["l1_ratio"],
        mae,
        r,
    )
    return metrics
