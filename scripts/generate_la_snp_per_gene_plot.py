#!/usr/bin/env python3
"""Supplementary figure: unique LA-SNPs per gene (horizontal bar chart).

Expects columns ``Gene`` and ``SNP_rsID`` by default. For spreadsheets that use
other headers (e.g. ``Gene Symbol`` / ``SNP Identifier``), pass ``--gene-column``
and ``--snp-column``.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_INPUT = REPO_ROOT / "overlapping_genes_with_snps.xlsx"
DEFAULT_OUTPUT = REPO_ROOT / "analysis" / "Fig_LA_SNPs_per_gene.png"
HIGHLIGHT_THRESHOLD = 3
COLOR_DEFAULT = "#5B7C99"
COLOR_HIGHLIGHT = "#C45C3E"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", type=Path, default=DEFAULT_INPUT, help="Excel path (.xlsx)")
    p.add_argument("--output", type=Path, default=DEFAULT_OUTPUT, help="PNG path (300 DPI)")
    p.add_argument("--gene-column", default="Gene", help="Gene name column")
    p.add_argument("--snp-column", default="SNP_rsID", help="SNP identifier column")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    df = pd.read_excel(args.input)
    gcol, scol = args.gene_column, args.snp_column
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
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=300, bbox_inches="tight")
    plt.close(fig)
    logging.info("Wrote %s", args.output)


if __name__ == "__main__":
    main()
