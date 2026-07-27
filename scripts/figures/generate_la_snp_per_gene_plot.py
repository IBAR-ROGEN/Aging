#!/usr/bin/env python3
"""Supplementary figure: unique LA-SNPs per gene (horizontal bar chart).

Expects columns ``Gene`` and ``SNP_rsID`` by default. For spreadsheets that use
other headers (e.g. ``Gene Symbol`` / ``SNP Identifier``), pass ``--gene-column``
and ``--snp-column``.
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import typer

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
DEFAULT_INPUT = REPO_ROOT / "overlapping_genes_with_snps.xlsx"
DEFAULT_OUTPUT = REPO_ROOT / "figures" / "Fig_LA_SNPs_per_gene.png"
HIGHLIGHT_THRESHOLD = 3
COLOR_DEFAULT = "#5B7C99"
COLOR_HIGHLIGHT = "#C45C3E"

app = typer.Typer(add_completion=False, help=__doc__)


@app.command()
def main(
    input: Path = typer.Option(
        DEFAULT_INPUT,
        "--input",
        help="Excel path (.xlsx)",
        path_type=Path,
    ),
    output: Path = typer.Option(
        DEFAULT_OUTPUT,
        "--output",
        help="PNG path (300 DPI); companion PDF written alongside",
        path_type=Path,
    ),
    gene_column: str = typer.Option("Gene", "--gene-column", help="Gene name column"),
    snp_column: str = typer.Option("SNP_rsID", "--snp-column", help="SNP identifier column"),
) -> None:
    """Write a horizontal bar chart of unique LA-SNPs per gene."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    df = pd.read_excel(input)
    gcol, scol = gene_column, snp_column
    missing = {gcol, scol} - set(df.columns)
    if missing:
        logging.error("Missing columns %s (available: %s)", sorted(missing), list(df.columns))
        raise SystemExit(1)
    counts = df.groupby(gcol, sort=False)[scol].nunique().sort_values(ascending=False)
    counts = counts.iloc[::-1]  # barh: largest count at top
    genes, vals = counts.index.astype(str).tolist(), counts.to_numpy()
    colors = [COLOR_HIGHLIGHT if v >= HIGHLIGHT_THRESHOLD else COLOR_DEFAULT for v in vals]
    logging.info("Genes: %d, unique SNPs (table-wide): %d", df[gcol].nunique(), df[scol].nunique())

    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "axes.facecolor": "white",
            "figure.facecolor": "white",
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )
    fig, ax = plt.subplots(figsize=(6, max(3.0, 0.18 * len(genes))))
    bars = ax.barh(genes, vals, color=colors, edgecolor="none", height=0.72)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlabel("Number of Longevity-Associated SNPs")
    ax.set_ylabel("")
    ax.set_title("LA-SNP distribution across the 41-gene AD/PD overlap set")
    ax.bar_label(bars, labels=[str(int(v)) for v in vals], padding=3, fontsize=8)
    ax.margins(x=0.12)
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    output_pdf = output.with_suffix(".pdf")
    fig.savefig(output, dpi=300, bbox_inches="tight")
    fig.savefig(output_pdf, bbox_inches="tight")
    plt.close(fig)
    logging.info("Wrote %s", output)
    logging.info("Wrote %s", output_pdf)


if __name__ == "__main__":
    app()
