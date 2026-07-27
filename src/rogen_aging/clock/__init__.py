"""Epigenetic clock training, evaluation, and external validation loaders."""

from rogen_aging.clock.data import (
    load_gse87571,
    load_romanian_cohort,
    load_wide_table,
    save_as_parquet,
    split_features_target,
    write_mock_romanian_cohort,
)
from rogen_aging.clock.evaluate import evaluate_clock
from rogen_aging.clock.model import make_clock_pipeline
from rogen_aging.clock.train import train_clock

__all__ = [
    "evaluate_clock",
    "load_gse87571",
    "load_romanian_cohort",
    "load_wide_table",
    "make_clock_pipeline",
    "save_as_parquet",
    "split_features_target",
    "train_clock",
    "write_mock_romanian_cohort",
]
