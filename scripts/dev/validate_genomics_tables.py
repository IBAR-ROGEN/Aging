#!/usr/bin/env python3
"""Lightweight offline genomics overlap-table schema check for pre-commit.

The full audit lives at ``analysis/validate_genomics_tables/validate_genomics_tables.py``
and hits network APIs. This script only validates the pandas DataFrame schema
(required columns, association labels, coordinate sanity) against a minimal
fixture so silent schema drift fails before commit.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import typer

from rogen_aging.pandas_schemas import (
    OVERLAP_REQUIRED_COLUMNS,
    assert_overlap_table_schema,
)

app = typer.Typer(add_completion=False, help=__doc__)

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_FIXTURE = REPO_ROOT / "test_data" / "genomics_overlap_minimal.csv"


def load_overlap_table(path: Path) -> pd.DataFrame:
    """Load an overlap table from CSV or Excel.

    Args:
        path: Path to a ``.csv``, ``.tsv``, ``.xlsx``, or ``.xls`` file.

    Returns:
        Loaded DataFrame.
    """
    if not path.is_file():
        raise FileNotFoundError(path)
    suffix = path.suffix.lower()
    if suffix == ".csv":
        return pd.read_csv(path)
    if suffix == ".tsv":
        return pd.read_csv(path, sep="\t")
    if suffix in {".xlsx", ".xls"}:
        return pd.read_excel(path, engine="openpyxl")
    raise ValueError(f"Unsupported overlap table extension: {path.suffix}")


def validate_overlap_frame(df: pd.DataFrame) -> None:
    """Run schema checks and emit a one-line summary.

    Args:
        df: Candidate overlap table.
    """
    assert_overlap_table_schema(df)
    n_genes = df["Gene Symbol"].nunique(dropna=True)
    n_snps = df["SNP Identifier"].nunique(dropna=True)
    typer.echo(
        "genomics schema OK: "
        f"rows={len(df)} genes={n_genes} snps={n_snps} "
        f"required_columns={len(OVERLAP_REQUIRED_COLUMNS)}"
    )


@app.command()
def main(
    fixture: Path = typer.Option(
        DEFAULT_FIXTURE,
        "--fixture",
        help="Minimal overlap-table fixture (CSV/TSV/XLSX) for offline schema checks",
    ),
    input_path: Path | None = typer.Option(
        None,
        "--input",
        help="Optional real overlap table to validate instead of the fixture",
    ),
) -> None:
    """Validate overlap-table schema without network calls."""
    path = input_path if input_path is not None else fixture
    typer.echo(f"Validating genomics overlap schema: {path}")
    frame = load_overlap_table(path)
    validate_overlap_frame(frame)


if __name__ == "__main__":
    app()
