"""Parquet loading and synthetic cohort generation for the EDA dashboard."""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl
import streamlit as st

from rogen_aging.eda_dashboard.schema import ensure_epigenetic_age_acceleration, normalize_column_names

REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_MERGED_PATH = (REPO_ROOT / "data" / "merged_cohort.parquet").resolve()


def default_merged_parquet_path() -> Path:
    env = os.environ.get("ROGEN_MERGED_COHORT_PARQUET")
    if env:
        return Path(env).expanduser().resolve()
    return DEFAULT_MERGED_PATH


@st.cache_data(show_spinner="Loading merged cohort (Parquet)…")
def load_merged_parquet(path_str: str) -> pd.DataFrame:
    """Load the integration pipeline merged table from Parquet via Polars."""
    df = pl.read_parquet(path_str).to_pandas()
    df = normalize_column_names(df)
    return ensure_epigenetic_age_acceleration(df)


@st.cache_data(show_spinner="Synthesizing in-memory cohort…")
def load_synthetic_cohort(*, n_samples: int = 320, random_seed: int = 42) -> pd.DataFrame:
    """Build a reproducible mock multi-omics-style cohort for offline exploration."""
    rng = np.random.default_rng(random_seed)
    n = int(n_samples)
    age = rng.normal(58.0, 14.0, n).clip(22.0, 92.0)
    sex = rng.choice(np.array(["Female", "Male"], dtype=object), size=n, replace=True)
    disease = rng.choice(
        np.array(["Control", "Case", "Prodromal"], dtype=object),
        size=n,
        replace=True,
        p=np.array([0.55, 0.35, 0.10]),
    )
    noise = rng.normal(0.0, 4.0, n)
    epi_age = age + noise
    hdl = rng.normal(52.0, 12.0, n).clip(25.0, 120.0)
    ldl = rng.normal(120.0, 35.0, n).clip(40.0, 220.0)
    bmi = rng.normal(27.0, 4.5, n).clip(18.0, 48.0)

    def sample_genotype(maf: float, size: int) -> np.ndarray:
        p0 = (1.0 - maf) ** 2
        p1 = 2.0 * maf * (1.0 - maf)
        p2 = maf**2
        return rng.choice(np.array([0, 1, 2], dtype=np.int8), size=size, replace=True, p=np.array([p0, p1, p2]))

    geno_5882 = sample_genotype(0.28, n)
    geno_7412 = sample_genotype(0.12, n)

    df = pl.DataFrame(
        {
            "Sample_ID": [f"ROGEN-{i:04d}" for i in range(n)],
            "Chronological_Age": age.astype(np.float64),
            "Sex": sex,
            "Disease_Status": disease,
            "Epigenetic_Age": epi_age.astype(np.float64),
            "HDL_Cholesterol": hdl.astype(np.float64),
            "LDL_Cholesterol": ldl.astype(np.float64),
            "BMI": bmi.astype(np.float64),
            "Phenotype_Score": rng.beta(2.0, 5.0, n).astype(np.float64),
            "rs5882_CETP": geno_5882,
            "rs7412_APOE": geno_7412,
        }
    )
    pdf = df.to_pandas()
    pdf = ensure_epigenetic_age_acceleration(pdf)
    return pdf
