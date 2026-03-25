"""Column resolution and cohort schema helpers for the EDA dashboard."""

from __future__ import annotations

import re
from typing import Iterable

import pandas as pd

# Canonical names expected by dashboard tabs after normalization
AGE_CANDIDATES: tuple[str, ...] = (
    "Chronological_Age",
    "chronological_age",
    "Age",
    "age",
)
SEX_CANDIDATES: tuple[str, ...] = ("Sex", "sex", "SEX", "gender", "Gender")
DISEASE_CANDIDATES: tuple[str, ...] = (
    "Disease_Status",
    "disease_status",
    "Disease",
    "disease",
    "Case_Control",
    "case_control",
)
EPIGEN_AGE_CANDIDATES: tuple[str, ...] = (
    "Epigenetic_Age",
    "epigenetic_age",
    "DNAmAge",
    "dnam_age",
)
SAMPLE_ID_CANDIDATES: tuple[str, ...] = ("Sample_ID", "sample_id", "IID", "participant_id")


def _first_present(columns: Iterable[str], candidates: tuple[str, ...]) -> str | None:
    cols = set(columns)
    for c in candidates:
        if c in cols:
            return c
    return None


def resolve_column(df: pd.DataFrame, candidates: tuple[str, ...], *, label: str) -> str:
    found = _first_present(df.columns, candidates)
    if found is None:
        raise KeyError(f"Merged cohort is missing a {label} column. Tried: {candidates}")
    return found


def normalize_column_names(df: pd.DataFrame) -> pd.DataFrame:
    """Rename common aliases to canonical dashboard names when unambiguous."""
    rename_map: dict[str, str] = {}
    if _first_present(df.columns, AGE_CANDIDATES) and "Chronological_Age" not in df.columns:
        src = _first_present(df.columns, AGE_CANDIDATES)
        if src is not None:
            rename_map[src] = "Chronological_Age"
    if _first_present(df.columns, SEX_CANDIDATES) and "Sex" not in df.columns:
        src = _first_present(df.columns, SEX_CANDIDATES)
        if src is not None:
            rename_map[src] = "Sex"
    if _first_present(df.columns, DISEASE_CANDIDATES) and "Disease_Status" not in df.columns:
        src = _first_present(df.columns, DISEASE_CANDIDATES)
        if src is not None:
            rename_map[src] = "Disease_Status"
    if _first_present(df.columns, EPIGEN_AGE_CANDIDATES) and "Epigenetic_Age" not in df.columns:
        src = _first_present(df.columns, EPIGEN_AGE_CANDIDATES)
        if src is not None:
            rename_map[src] = "Epigenetic_Age"
    if _first_present(df.columns, SAMPLE_ID_CANDIDATES) and "Sample_ID" not in df.columns:
        src = _first_present(df.columns, SAMPLE_ID_CANDIDATES)
        if src is not None:
            rename_map[src] = "Sample_ID"
    if rename_map:
        return df.rename(columns=rename_map)
    return df


def list_snp_columns(columns: Iterable[str]) -> list[str]:
    """Columns that look like SNP / variant IDs (e.g. rs5882_CETP, CETP_rs5882)."""
    out: list[str] = []
    for c in columns:
        if re.search(r"rs\d+", str(c), flags=re.IGNORECASE):
            out.append(str(c))
    return sorted(out)


def ensure_epigenetic_age_acceleration(df: pd.DataFrame) -> pd.DataFrame:
    """Add Epigenetic_Age_Acceleration if Epigenetic_Age and Chronological_Age exist."""
    if "Epigenetic_Age_Acceleration" in df.columns:
        return df
    if "Epigenetic_Age" in df.columns and "Chronological_Age" in df.columns:
        out = df.copy()
        out["Epigenetic_Age_Acceleration"] = out["Epigenetic_Age"] - out["Chronological_Age"]
        return out
    return df


def clinical_numeric_columns_for_heatmap(df: pd.DataFrame) -> list[str]:
    """Continuous columns suitable for a clinical correlation heatmap."""
    skip_substrings = ("sample", "id", "IID")
    snp_cols = set(list_snp_columns(df.columns))
    numeric = df.select_dtypes(include=["number"]).columns.tolist()
    selected: list[str] = []
    for c in numeric:
        cl = c.lower()
        if c in snp_cols:
            continue
        if any(s in cl for s in skip_substrings):
            continue
        if df[c].notna().sum() < 3:
            continue
        if df[c].std(skipna=True) == 0 or (df[c].std(skipna=True) != df[c].std(skipna=True)):
            continue
        selected.append(c)
    return selected
