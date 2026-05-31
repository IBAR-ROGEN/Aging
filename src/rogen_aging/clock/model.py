"""Sklearn pipeline factory for epigenetic clock training."""

from __future__ import annotations

from sklearn.impute import SimpleImputer
from sklearn.linear_model import ElasticNetCV
from sklearn.pipeline import Pipeline


def make_clock_pipeline(
    *,
    random_state: int = 42,
    l1_ratio: list[float] | None = None,
    n_alphas: int = 20,
    cv: int = 10,
    max_iter: int = 5000,
) -> Pipeline:
    """Build a mean-imputation + ElasticNetCV pipeline for wide CpG tables.

    The imputer is fit on the training split only when used inside
    :func:`rogen_aging.clock.train.train_clock`.
    """
    ratios = l1_ratio if l1_ratio is not None else [0.1, 0.5, 0.7, 0.9, 0.95, 0.99, 1.0]
    enet = ElasticNetCV(
        l1_ratio=ratios,
        n_alphas=n_alphas,
        cv=cv,
        random_state=random_state,
        max_iter=max_iter,
    )
    return Pipeline([("imputer", SimpleImputer(strategy="mean")), ("elasticnet", enet)])
