"""Sklearn pipeline factory for epigenetic clock training."""

from __future__ import annotations

from sklearn.impute import SimpleImputer
from sklearn.linear_model import ElasticNetCV
from sklearn.pipeline import Pipeline

from rogen_aging.config import get_config


def make_clock_pipeline(
    *,
    random_state: int | None = None,
    l1_ratio: list[float] | None = None,
    alphas: int | None = None,
    cv: int | None = None,
    max_iter: int | None = None,
) -> Pipeline:
    """Build a mean-imputation + ElasticNetCV pipeline for wide CpG tables.

    The imputer is fit on the training split only when used inside
    :func:`rogen_aging.clock.train.train_clock`. Hyperparameter defaults are
    read from the active pipeline config (``clock.elasticnet``).

    Args:
        random_state: Seed passed to ``ElasticNetCV``.
        l1_ratio: Candidate L1 ratios for ElasticNetCV.
        alphas: Number of alpha values along the regularization path.
        cv: Number of cross-validation folds.
        max_iter: Maximum coordinate-descent iterations per fit.

    Returns:
        An unfitted ``sklearn.pipeline.Pipeline`` with ``imputer`` and
        ``elasticnet`` steps.
    """
    cfg = get_config().clock
    enet_cfg = cfg.elasticnet
    ratios = list(l1_ratio) if l1_ratio is not None else [float(x) for x in enet_cfg.l1_ratio]
    enet = ElasticNetCV(
        l1_ratio=ratios,
        alphas=int(alphas if alphas is not None else enet_cfg.alphas),
        cv=int(cv if cv is not None else enet_cfg.cv),
        random_state=int(random_state if random_state is not None else cfg.random_state),
        max_iter=int(max_iter if max_iter is not None else enet_cfg.max_iter),
    )
    return Pipeline([("imputer", SimpleImputer(strategy="mean")), ("elasticnet", enet)])
