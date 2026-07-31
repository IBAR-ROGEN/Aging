"""Typed pandas DataFrame schemas used across ``rogen_aging``.

Mypy cannot see DataFrame column names, so we document expected structures with
``TypedDict`` definitions and enforce them at runtime with the helpers below.
"""

from __future__ import annotations

from typing import Final, NotRequired, TypedDict

import pandas as pd

# ---------------------------------------------------------------------------
# Gene–LA-SNP overlap table (input to genomics validation)
# ---------------------------------------------------------------------------

OverlapTableRow = TypedDict(
    "OverlapTableRow",
    {
        "Gene Symbol": str,
        "SNP Identifier": str,
        "SNP Association": str,
        "Gene Location": NotRequired[str],
        "SNP Location": NotRequired[str],
        "Start": NotRequired[int],
        "End": NotRequired[int],
        "SNP PubMed ID": NotRequired[int],
    },
)

OVERLAP_REQUIRED_COLUMNS: Final[tuple[str, ...]] = (
    "Gene Symbol",
    "SNP Identifier",
    "SNP Association",
    "Gene Location",
    "SNP Location",
    "Start",
    "End",
)

OVERLAP_ASSOCIATION_VALUES: Final[frozenset[str]] = frozenset({"significant", "non-significant"})


# ---------------------------------------------------------------------------
# Epigenetic clock wide tables
# ---------------------------------------------------------------------------


class ClockWideRow(TypedDict):
    """Wide methylation table row: chronological age plus CpG beta columns."""

    chronological_age: float


CLOCK_REQUIRED_COLUMNS: Final[tuple[str, ...]] = ("chronological_age",)


# ---------------------------------------------------------------------------
# EDA merged cohort (canonical names after normalization)
# ---------------------------------------------------------------------------


class EdaCohortRow(TypedDict):
    """Merged clinical / epigenetic cohort after dashboard normalization."""

    Sample_ID: str
    Chronological_Age: float
    Sex: NotRequired[str]
    Disease_Status: NotRequired[str]
    Epigenetic_Age: NotRequired[float]


EDA_REQUIRED_COLUMNS: Final[tuple[str, ...]] = (
    "Sample_ID",
    "Chronological_Age",
)


def missing_columns(df: pd.DataFrame, required: tuple[str, ...]) -> list[str]:
    """Return required column names absent from ``df`` (stable order)."""
    present = set(df.columns)
    return [name for name in required if name not in present]


def assert_columns(
    df: pd.DataFrame,
    required: tuple[str, ...],
    *,
    label: str,
) -> pd.DataFrame:
    """Raise ``ValueError`` when ``df`` is missing any required columns.

    Args:
        df: Frame to validate.
        required: Column names that must be present.
        label: Human-readable schema name used in the error message.

    Returns:
        The same ``df`` when validation succeeds.
    """
    absent = missing_columns(df, required)
    if absent:
        raise ValueError(f"{label} is missing required columns: {absent}")
    return df


def assert_overlap_table_schema(df: pd.DataFrame) -> pd.DataFrame:
    """Validate the gene–LA-SNP overlap table schema used by genomics checks.

    Args:
        df: Candidate overlap table.

    Returns:
        ``df`` unchanged when the schema is valid.

    Raises:
        ValueError: If required columns are missing, the frame is empty, or
            ``SNP Association`` contains unexpected values.
    """
    assert_columns(df, OVERLAP_REQUIRED_COLUMNS, label="Overlap table")
    if df.empty:
        raise ValueError("Overlap table has no rows")

    associations = {str(value) for value in df["SNP Association"].dropna().unique().tolist()}
    unexpected = sorted(associations - OVERLAP_ASSOCIATION_VALUES)
    if unexpected:
        raise ValueError(
            "Overlap table has unexpected SNP Association values: "
            f"{unexpected}; expected subset of {sorted(OVERLAP_ASSOCIATION_VALUES)}"
        )

    # Lightweight dtype sanity for genomic coordinates when present.
    for col in ("Start", "End"):
        series = pd.to_numeric(df[col], errors="coerce")
        if series.notna().any() and (series.dropna() < 0).any():
            raise ValueError(f"Overlap table column {col!r} contains negative coordinates")

    return df


def assert_clock_wide_schema(df: pd.DataFrame) -> pd.DataFrame:
    """Validate a GSE40279-style wide methylation table.

    Args:
        df: Wide table with ``chronological_age`` and ``cg*`` feature columns.

    Returns:
        ``df`` unchanged when the schema is valid.
    """
    assert_columns(df, CLOCK_REQUIRED_COLUMNS, label="Clock wide table")
    cpg_cols = [name for name in df.columns if str(name).startswith("cg")]
    if not cpg_cols:
        raise ValueError("Clock wide table has no CpG columns starting with 'cg'")
    return df


def assert_eda_cohort_schema(df: pd.DataFrame) -> pd.DataFrame:
    """Validate a dashboard cohort after column normalization.

    Args:
        df: Merged cohort frame.

    Returns:
        ``df`` unchanged when the schema is valid.
    """
    return assert_columns(df, EDA_REQUIRED_COLUMNS, label="EDA cohort")
