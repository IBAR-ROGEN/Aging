#!/usr/bin/env python3
"""Build an empty longevity-allele curation worksheet for the 47-variant panel.

Prefills ``rsid`` and ``gene`` from ``variants_47_input_v2.csv``. All other
columns from the existing ``longevity_allele_direction.csv`` header stay empty.
No direction, allele, PMID, population, or note is invented.

Example:
    uv run python scripts/curation/make_curation_worksheet.py
"""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import polars as pl
import typer

from rogen_aging.config import find_repo_root
from rogen_aging.longevity.allele_harmonize import is_palindromic

REPO_ROOT = find_repo_root()
DEFAULT_VARIANTS = REPO_ROOT / "data" / "processed" / "variants_47_input_v2.csv"
DEFAULT_SCHEMA = REPO_ROOT / "data" / "curation" / "longevity_allele_direction.csv"
DEFAULT_PENDING = REPO_ROOT / "docs" / "PENDING_ALLELES.md"
DEFAULT_OUT = REPO_ROOT / "data" / "curation" / "longevity_allele_direction.csv"
DEFAULT_PROVENANCE = REPO_ROOT / "data" / "curation" / "make_curation_worksheet_provenance.json"
IDENTIFIER_COLUMNS = ("rsid", "gene")
PALINDROME_FLAG = "orientation cannot be inferred from frequency"
EXPECTED_N = 47


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _git_commit() -> str:
    completed = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=REPO_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    return completed.stdout.strip()


def _repo_relative(path: Path) -> str:
    resolved = path.resolve()
    root = REPO_ROOT.resolve()
    try:
        return resolved.relative_to(root).as_posix()
    except ValueError:
        return str(path)


def read_schema_columns(schema_csv: Path) -> list[str]:
    """Return the header of the existing curation CSV. Add no columns."""
    assert schema_csv.is_file(), f"MISSING INPUT: {schema_csv}"
    with schema_csv.open("r", encoding="utf-8", newline="") as handle:
        header = handle.readline().strip("\n").strip("\r")
    columns = [part.strip() for part in header.split(",")]
    assert columns, f"STOP: empty header in {schema_csv}"
    for required in IDENTIFIER_COLUMNS:
        assert required in columns, f"STOP: schema missing identifier column {required!r}"
    return columns


def load_panel(path: Path) -> pl.DataFrame:
    """Load the 47-variant panel and assert the ground-truth row count."""
    assert path.is_file(), f"MISSING INPUT: {path}"
    frame = pl.read_csv(
        path,
        schema_overrides={
            "chrom": pl.Utf8,
            "pos": pl.Utf8,
            "ref": pl.Utf8,
            "alt": pl.Utf8,
            "rsid": pl.Utf8,
            "gene_symbol": pl.Utf8,
            "alt_status": pl.Utf8,
        },
    )
    required = {"rsid", "gene_symbol", "ref", "alt"}
    missing = sorted(required - set(frame.columns))
    assert not missing, f"STOP: panel missing columns: {missing}"
    assert (
        frame.height == EXPECTED_N
    ), f"ROW COUNT DISAGREES: computed={frame.height} ground_truth={EXPECTED_N}"
    return frame


def build_worksheet(panel: pl.DataFrame, columns: list[str]) -> pl.DataFrame:
    """Prefill identifiers only; leave every curation column empty."""
    empty_cols = [name for name in columns if name not in IDENTIFIER_COLUMNS]
    data: dict[str, list[str]] = {
        "rsid": [str(value) for value in panel["rsid"].to_list()],
        "gene": [str(value) for value in panel["gene_symbol"].to_list()],
    }
    for name in empty_cols:
        data[name] = [""] * panel.height
    worksheet = pl.DataFrame(data).select(columns)
    assert (
        worksheet.height == EXPECTED_N
    ), f"ROW COUNT DISAGREES: computed={worksheet.height} ground_truth={EXPECTED_N}"
    for name in empty_cols:
        nonempty = worksheet.filter(pl.col(name).cast(pl.Utf8).str.strip_chars() != "")
        assert (
            nonempty.height == 0
        ), f"STOP: curation column {name!r} was filled; worksheet must leave it empty"
    return worksheet


def palindrome_rows(panel: pl.DataFrame) -> list[dict[str, str]]:
    """Return panel rows whose ref/alt are A/T or C/G."""
    rows: list[dict[str, str]] = []
    for row in panel.iter_rows(named=True):
        ref = str(row["ref"] or "")
        alt = str(row["alt"] or "")
        if not is_palindromic(ref, alt):
            continue
        rows.append(
            {
                "rsid": str(row["rsid"]),
                "gene": str(row["gene_symbol"]),
                "ref": ref.upper(),
                "alt": alt.upper(),
                "flag": PALINDROME_FLAG,
            }
        )
    return rows


def write_provenance(
    path: Path,
    *,
    input_hashes: dict[str, dict[str, str]],
    out: Path,
    palindromes: list[dict[str, str]],
) -> None:
    """Write tool versions, input hashes, git commit, and timestamp next to output."""
    payload: dict[str, Any] = {
        "generated_at_utc": datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "git_commit": _git_commit(),
        "python": sys.version.split()[0],
        "polars": pl.__version__,
        "typer": typer.__version__,
        "expected_n_rows": EXPECTED_N,
        "inputs": input_hashes,
        "output": {
            "path": _repo_relative(out),
            "sha256": _sha256(out),
            "n_rows": EXPECTED_N,
        },
        "palindromic_variants": palindromes,
        "palindrome_flag": PALINDROME_FLAG,
        "note": (
            "Identifiers prefilled from the panel. Curation columns left empty. "
            "No allele, direction, PMID, population, or note was filled by this script."
        ),
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def main(
    variants: Path = typer.Option(
        DEFAULT_VARIANTS,
        "--variants",
        help="Panel CSV (variants_47_input_v2.csv).",
    ),
    schema_csv: Path = typer.Option(
        DEFAULT_SCHEMA,
        "--schema-csv",
        help="Existing longevity_allele_direction.csv; header only is the schema.",
    ),
    pending: Path = typer.Option(
        DEFAULT_PENDING,
        "--pending",
        help="PENDING_ALLELES.md (referenced for provenance; not used to fill cells).",
    ),
    out: Path = typer.Option(
        DEFAULT_OUT,
        "--out",
        help="Worksheet CSV to write (same schema; curation columns empty).",
    ),
    provenance: Path = typer.Option(
        DEFAULT_PROVENANCE,
        "--provenance",
        help="Provenance JSON written next to the worksheet run.",
    ),
) -> None:
    """Write a 47-row curation worksheet with identifiers only."""
    assert pending.is_file(), f"MISSING INPUT: {pending}"
    columns = read_schema_columns(schema_csv)
    # Hash schema before overwrite when --out equals --schema-csv.
    input_hashes = {
        "variants": {
            "path": _repo_relative(variants),
            "sha256": _sha256(variants),
        },
        "schema_csv": {
            "path": _repo_relative(schema_csv),
            "sha256": _sha256(schema_csv),
        },
        "pending_alleles": {
            "path": _repo_relative(pending),
            "sha256": _sha256(pending),
            "role": "referenced; worksheet does not copy allele decisions",
        },
    }
    panel = load_panel(variants)
    worksheet = build_worksheet(panel, columns)
    palindromes = palindrome_rows(panel)

    out.parent.mkdir(parents=True, exist_ok=True)
    worksheet.write_csv(out)
    assert out.is_file(), f"STOP: failed to write {out}"
    written = pl.read_csv(out, schema_overrides={name: pl.Utf8 for name in columns})
    assert (
        written.height == EXPECTED_N
    ), f"ROW COUNT DISAGREES: computed={written.height} ground_truth={EXPECTED_N}"
    assert (
        list(written.columns) == columns
    ), f"STOP: output columns {written.columns} != schema {columns}"

    write_provenance(
        provenance,
        input_hashes=input_hashes,
        out=out,
        palindromes=palindromes,
    )

    typer.echo(f"wrote {_repo_relative(out)} rows={written.height}")
    typer.echo(f"schema_columns={columns}")
    typer.echo(f"provenance={_repo_relative(provenance)}")
    if palindromes:
        typer.echo(f"palindromic_variants={len(palindromes)} ({PALINDROME_FLAG})")
        for row in palindromes:
            typer.echo(
                f"  {row['rsid']} {row['gene']} {row['ref']}/{row['alt']}: {PALINDROME_FLAG}"
            )
    else:
        typer.echo("palindromic_variants=0")


if __name__ == "__main__":
    typer.run(main)
