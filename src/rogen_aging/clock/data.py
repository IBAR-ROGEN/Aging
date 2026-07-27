"""Data loaders for epigenetic clock training and external validation."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl

from rogen_aging.clock.external_data import load_gse87571, save_as_parquet

__all__ = [
    "load_gse87571",
    "load_romanian_cohort",
    "load_wide_table",
    "save_as_parquet",
    "split_features_target",
    "write_mock_romanian_cohort",
]


def load_wide_table(path: Path) -> pd.DataFrame:
    """Load a GSE40279-style wide table (Parquet, CSV, or TSV)."""
    if not path.is_file():
        raise FileNotFoundError(path)
    suffix = path.suffix.lower()
    if suffix == ".parquet":
        return pd.read_parquet(path)
    if suffix in (".csv", ".tsv"):
        sep = "\t" if suffix == ".tsv" else ","
        return pd.read_csv(path, sep=sep)
    raise ValueError(f"Unsupported extension {path.suffix}; use .parquet, .csv, or .tsv")


def split_features_target(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.Series]:
    """Return CpG feature matrix (``cg*`` columns) and chronological age target."""
    if "chronological_age" not in df.columns:
        raise KeyError("Expected column 'chronological_age'")
    cpg_cols = [c for c in df.columns if str(c).startswith("cg")]
    if not cpg_cols:
        raise ValueError("No feature columns starting with 'cg'")
    return df[cpg_cols].copy(), df["chronological_age"].copy()


def write_mock_romanian_cohort(
    data_dir: Path,
    *,
    n_samples: int = 120,
    n_cpgs: int = 80,
    random_state: int = 42,
) -> None:
    """Create deterministic mock methylation + metadata CSVs for development."""
    data_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(random_state)
    sample_ids = [f"ROM{i:04d}" for i in range(1, n_samples + 1)]
    ages = rng.normal(52.0, 14.0, size=n_samples).clip(25.0, 92.0)

    beta = rng.uniform(0.15, 0.85, size=(n_samples, n_cpgs)).astype(np.float64)
    signal_idx = rng.choice(n_cpgs, size=max(12, n_cpgs // 5), replace=False)
    for j in signal_idx:
        beta[:, j] += 0.018 * (ages - ages.mean())
    beta = np.clip(beta, 0.02, 0.98)

    cpg_names = [f"cg_mock_{k:05d}" for k in range(1, n_cpgs + 1)]
    meth = pl.DataFrame(beta, schema=cpg_names).with_columns(pl.Series("sample_id", sample_ids))
    meta = pl.DataFrame({"sample_id": sample_ids, "chronological_age": ages.astype(np.float64)})
    meth.write_csv(data_dir / "methylation_matrix.csv")
    meta.write_csv(data_dir / "metadata.csv")


def load_romanian_cohort(
    data_dir: Path,
    *,
    regenerate_mock: bool = False,
    random_state: int = 42,
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """Load Romanian-style two-file cohort; return X, y, and sample IDs."""
    meth_path = data_dir / "methylation_matrix.csv"
    meta_path = data_dir / "metadata.csv"
    if regenerate_mock or not meth_path.is_file() or not meta_path.is_file():
        write_mock_romanian_cohort(data_dir, random_state=random_state)

    meth = pl.read_csv(meth_path)
    meta = pl.read_csv(meta_path)
    if "sample_id" not in meth.columns or "sample_id" not in meta.columns:
        raise ValueError("Both tables must include a 'sample_id' column.")
    if "chronological_age" not in meta.columns:
        raise ValueError("Metadata must include 'chronological_age'.")

    joined = meta.join(meth, on="sample_id", how="inner").sort("sample_id")
    feature_cols = [c for c in joined.columns if c not in ("sample_id", "chronological_age")]
    if not feature_cols:
        raise ValueError("No CpG / feature columns found after join.")

    x_df = joined.select(feature_cols)
    x = x_df.to_numpy()
    y = joined["chronological_age"].to_numpy()
    sample_ids = joined["sample_id"].to_list()
    return x, y, sample_ids
