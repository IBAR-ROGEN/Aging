#!/usr/bin/env python3
"""Panel gene overlap with Broad Human MitoCarta3.0 (descriptive only).

Reads gene membership and pathway text only from the downloaded MitoCarta
workbook. Never invents a mitochondrial gene list.

One output row per distinct panel gene from
``data/processed/variants_47_input_v2.csv`` (column ``gene_symbol``), joined to:

- MitoCarta Symbol / MitoPathways from ``data/external/Human.MitoCarta3.0.xlsx``
- ``evidence_labels`` from ``analysis/longevity/integrated_evidence.csv``
- minimum ``q_bh`` from ``analysis/longevity/longevity_gwas_lookup.csv``

``results/eqtl_associations_full.csv`` is a required input and is recorded in
provenance (integrated evidence already carries eQTL labels).

Descriptive overlap counts only. No enrichment, hypergeometric, or background
test. The panel is literature-curated, so no valid background exists.

Example:
    uv run python scripts/mito/mitocarta_panel_overlap.py
"""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import pandas as pd
import polars as pl
import typer

from rogen_aging.config import find_repo_root

REPO_ROOT = find_repo_root()
DEFAULT_MITOCARTA = REPO_ROOT / "data" / "external" / "Human.MitoCarta3.0.xlsx"
DEFAULT_PANEL = REPO_ROOT / "data" / "processed" / "variants_47_input_v2.csv"
DEFAULT_GWAS = REPO_ROOT / "analysis" / "longevity" / "longevity_gwas_lookup.csv"
DEFAULT_EQTL = REPO_ROOT / "results" / "eqtl_associations_full.csv"
DEFAULT_INTEGRATED = REPO_ROOT / "analysis" / "longevity" / "integrated_evidence.csv"
DEFAULT_OUT_DIR = REPO_ROOT / "analysis" / "mito"
MITOCARTA_SHEET = "A Human MitoCarta3.0"
MITOCARTA_SYMBOL_COL = "Symbol"
MITOCARTA_PATHWAY_COL = "MitoPathways"
LABEL_ORDER = ("GWAS support", "eQTL support", "population difference")
EXPECTED_PANEL_ROWS = 47

app = typer.Typer(add_completion=False, no_args_is_help=False)


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


def _cell(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and pd.isna(value):
        return ""
    text = str(value).strip()
    if text.lower() in {"nan", "none", "null"}:
        return ""
    return text


def _tool_versions() -> dict[str, str]:
    return {
        "python": sys.version.split()[0],
        "polars": pl.__version__,
        "pandas": pd.__version__,
        "openpyxl": _package_version("openpyxl"),
    }


def _package_version(name: str) -> str:
    from importlib.metadata import PackageNotFoundError, version

    try:
        return version(name)
    except PackageNotFoundError as exc:
        raise AssertionError(f"MISSING PACKAGE: {name}") from exc


def load_panel_genes(path: Path) -> pl.DataFrame:
    """Distinct panel genes from ``gene_symbol`` (one row per gene)."""
    assert path.is_file(), f"MISSING INPUT: {_repo_relative(path)}"
    frame = pl.read_csv(
        path,
        schema_overrides={
            "rsid": pl.Utf8,
            "gene_symbol": pl.Utf8,
            "alt_status": pl.Utf8,
        },
    )
    required = {"rsid", "gene_symbol", "alt_status"}
    missing = sorted(required - set(frame.columns))
    assert not missing, f"STOP: panel missing columns: {missing}"
    assert frame.height == EXPECTED_PANEL_ROWS, (
        f"ROW COUNT DISAGREES: computed={frame.height} " f"ground_truth={EXPECTED_PANEL_ROWS}"
    )
    n_resolved = frame.filter(pl.col("alt_status") == "resolved").height
    assert n_resolved == EXPECTED_PANEL_ROWS, (
        f"RESOLVED COUNT DISAGREES: computed={n_resolved} " f"ground_truth={EXPECTED_PANEL_ROWS}"
    )
    genes = (
        frame.select(pl.col("gene_symbol").str.strip_chars().alias("gene")).unique().sort("gene")
    )
    assert genes.filter(pl.col("gene") == "").height == 0, "STOP: empty gene_symbol"
    return genes


def load_mitocarta(path: Path) -> pl.DataFrame:
    """Load Symbol and MitoPathways from the downloaded Broad workbook only."""
    assert path.is_file(), (
        f"MISSING INPUT: {_repo_relative(path)}. "
        "Download Broad Human MitoCarta3.0 and place it at this path "
        "(convert .xls to .xlsx if needed). "
        "This script never writes a mitochondrial gene list from memory."
    )
    workbook = pd.ExcelFile(path, engine="openpyxl")
    assert MITOCARTA_SHEET in workbook.sheet_names, (
        f"STOP: expected sheet {MITOCARTA_SHEET!r}; " f"found {workbook.sheet_names!r}"
    )
    raw = pd.read_excel(workbook, sheet_name=MITOCARTA_SHEET, dtype=str)
    columns = list(raw.columns)
    for required in (MITOCARTA_SYMBOL_COL, MITOCARTA_PATHWAY_COL):
        assert required in columns, (
            f"STOP: MitoCarta sheet missing column {required!r}; " f"found columns={columns!r}"
        )
    frame = pl.from_pandas(raw[[MITOCARTA_SYMBOL_COL, MITOCARTA_PATHWAY_COL]].copy()).rename(
        {
            MITOCARTA_SYMBOL_COL: "gene",
            MITOCARTA_PATHWAY_COL: "mitocarta_pathway",
        }
    )
    frame = frame.with_columns(
        pl.col("gene").map_elements(_cell, return_dtype=pl.Utf8),
        pl.col("mitocarta_pathway").map_elements(_cell, return_dtype=pl.Utf8),
    )
    frame = frame.filter(pl.col("gene") != "")
    assert frame.height > 0, "STOP: MitoCarta sheet has no Symbol rows"
    # One pathway string per Symbol as given in the file (first row if duplicated).
    frame = (
        frame.unique(subset=["gene"], keep="first")
        .with_columns(pl.lit(True).alias("_in_mitocarta"))
        .sort("gene")
    )
    return frame


def load_integrated_evidence(path: Path) -> pl.DataFrame:
    """Gene-level evidence labels (union across rsids for the same gene)."""
    assert path.is_file(), f"MISSING INPUT: {_repo_relative(path)}"
    frame = pl.read_csv(
        path,
        schema_overrides={
            "rsid": pl.Utf8,
            "gene": pl.Utf8,
            "evidence_labels": pl.Utf8,
        },
    )
    required = {"rsid", "gene", "evidence_labels"}
    missing = sorted(required - set(frame.columns))
    assert not missing, f"STOP: integrated evidence missing columns: {missing}"

    by_gene: dict[str, set[str]] = {}
    for row in frame.iter_rows(named=True):
        gene = _cell(row["gene"])
        assert gene, "STOP: empty gene in integrated evidence"
        labels = {part.strip() for part in _cell(row["evidence_labels"]).split("|") if part.strip()}
        by_gene.setdefault(gene, set()).update(labels)

    order_index = {name: index for index, name in enumerate(LABEL_ORDER)}

    def _join(labels: set[str]) -> str:
        return "|".join(
            sorted(labels, key=lambda name: (order_index.get(name, len(LABEL_ORDER)), name))
        )

    return pl.DataFrame(
        {
            "gene": sorted(by_gene),
            "evidence_labels": [_join(by_gene[gene]) for gene in sorted(by_gene)],
        }
    )


def load_min_gwas_q(path: Path) -> pl.DataFrame:
    """Minimum BH q per gene among GWAS rows with a finite ``q_bh``."""
    assert path.is_file(), f"MISSING INPUT: {_repo_relative(path)}"
    frame = pl.read_csv(
        path,
        schema_overrides={
            "rsid": pl.Utf8,
            "gene_symbol": pl.Utf8,
            "q_bh": pl.Float64,
            "status": pl.Utf8,
        },
    )
    required = {"rsid", "gene_symbol", "q_bh", "status"}
    missing = sorted(required - set(frame.columns))
    assert not missing, f"STOP: GWAS lookup missing columns: {missing}"
    present = frame.filter(pl.col("q_bh").is_not_null()).select(
        pl.col("gene_symbol").str.strip_chars().alias("gene"),
        pl.col("q_bh"),
    )
    return present.group_by("gene").agg(pl.col("q_bh").min().alias("min_gwas_q_bh")).sort("gene")


def require_eqtl(path: Path) -> dict[str, Any]:
    """Assert the eQTL associations table is present; return provenance fields."""
    assert path.is_file(), f"MISSING INPUT: {_repo_relative(path)}"
    frame = pl.read_csv(path, n_rows=5)
    required = {"canonical_rsid", "egene_symbol", "fdr_bh", "dataset_id"}
    missing = sorted(required - set(frame.columns))
    assert not missing, f"STOP: eQTL associations missing columns: {missing}"
    n_rows = sum(1 for _ in path.open(encoding="utf-8")) - 1
    return {
        "path": _repo_relative(path),
        "sha256": _sha256(path),
        "n_rows": n_rows,
        "columns_checked": sorted(required),
    }


def build_table(
    panel_genes: pl.DataFrame,
    mitocarta: pl.DataFrame,
    evidence: pl.DataFrame,
    min_q: pl.DataFrame,
) -> pl.DataFrame:
    """One row per panel gene with Y/N membership and joined annotations."""
    joined = (
        panel_genes.join(mitocarta, on="gene", how="left")
        .join(evidence, on="gene", how="left")
        .join(min_q, on="gene", how="left")
        .with_columns(
            pl.when(pl.col("_in_mitocarta").fill_null(False))
            .then(pl.lit("Y"))
            .otherwise(pl.lit("N"))
            .alias("in_MitoCarta"),
            pl.col("mitocarta_pathway").fill_null(""),
            pl.col("evidence_labels").fill_null(""),
        )
        .select(
            "gene",
            "in_MitoCarta",
            "mitocarta_pathway",
            "evidence_labels",
            "min_gwas_q_bh",
        )
        .sort("gene")
    )
    assert joined.height == panel_genes.height
    return joined


def descriptive_counts(table: pl.DataFrame) -> dict[str, Any]:
    """Counts only. No background comparison."""
    in_yes = table.filter(pl.col("in_MitoCarta") == "Y")
    in_no = table.filter(pl.col("in_MitoCarta") == "N")
    return {
        "n_panel_genes": table.height,
        "n_in_MitoCarta": in_yes.height,
        "n_not_in_MitoCarta": in_no.height,
        "genes_in_MitoCarta": in_yes["gene"].to_list(),
        "genes_not_in_MitoCarta": in_no["gene"].to_list(),
    }


@app.command()
def main(
    mitocarta: Path = typer.Option(DEFAULT_MITOCARTA, help="Broad Human.MitoCarta3.0.xlsx"),
    panel: Path = typer.Option(DEFAULT_PANEL, help="Panel variants CSV"),
    gwas: Path = typer.Option(DEFAULT_GWAS, help="Longevity GWAS lookup CSV"),
    eqtl: Path = typer.Option(DEFAULT_EQTL, help="eQTL associations CSV"),
    integrated: Path = typer.Option(DEFAULT_INTEGRATED, help="Integrated evidence CSV"),
    out_dir: Path = typer.Option(DEFAULT_OUT_DIR, help="Output directory"),
) -> None:
    """Write panel × MitoCarta descriptive overlap table and provenance."""
    panel_genes = load_panel_genes(panel)
    mito = load_mitocarta(mitocarta)
    evidence = load_integrated_evidence(integrated)
    min_q = load_min_gwas_q(gwas)
    eqtl_meta = require_eqtl(eqtl)

    table = build_table(panel_genes, mito, evidence, min_q)
    counts = descriptive_counts(table)

    out_dir.mkdir(parents=True, exist_ok=True)
    out_csv = out_dir / "mitocarta_panel_overlap.csv"
    out_prov = out_dir / "mitocarta_panel_overlap_provenance.json"
    table.write_csv(out_csv)

    timestamp = datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")
    provenance: dict[str, Any] = {
        "generated_at_utc": timestamp,
        "git_commit": _git_commit(),
        "tool_versions": _tool_versions(),
        "inputs": {
            "mitocarta": {
                "path": _repo_relative(mitocarta),
                "sha256": _sha256(mitocarta),
                "sheet": MITOCARTA_SHEET,
                "symbol_column": MITOCARTA_SYMBOL_COL,
                "pathway_column": MITOCARTA_PATHWAY_COL,
                "n_symbols": mito.height,
            },
            "panel": {
                "path": _repo_relative(panel),
                "sha256": _sha256(panel),
                "gene_column": "gene_symbol",
                "n_rows": EXPECTED_PANEL_ROWS,
                "n_unique_genes": panel_genes.height,
            },
            "gwas_lookup": {
                "path": _repo_relative(gwas),
                "sha256": _sha256(gwas),
                "q_column": "q_bh",
            },
            "eqtl_associations": eqtl_meta,
            "integrated_evidence": {
                "path": _repo_relative(integrated),
                "sha256": _sha256(integrated),
                "label_column": "evidence_labels",
                "gene_column": "gene",
            },
        },
        "outputs": {
            "table": _repo_relative(out_csv),
            "table_sha256": _sha256(out_csv),
            "provenance": _repo_relative(out_prov),
        },
        "counts": counts,
        "method_notes": [
            "Membership and pathway text come only from the MitoCarta workbook.",
            "No enrichment or hypergeometric test was run.",
            "The panel is literature-curated; no valid background set is defined.",
        ],
    }
    out_prov.write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")

    typer.echo(
        f"Panel genes: {counts['n_panel_genes']}. "
        f"In MitoCarta: {counts['n_in_MitoCarta']}. "
        f"Not in MitoCarta: {counts['n_not_in_MitoCarta']}."
    )
    if counts["genes_in_MitoCarta"]:
        typer.echo("Genes in MitoCarta: " + ", ".join(counts["genes_in_MitoCarta"]))
    typer.echo(f"Wrote {_repo_relative(out_csv)}")
    typer.echo(f"Wrote {_repo_relative(out_prov)}")


if __name__ == "__main__":
    app()
