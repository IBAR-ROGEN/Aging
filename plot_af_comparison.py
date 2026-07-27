#!/usr/bin/env python3
"""Publication-quality allele-frequency comparison for longevity-associated SNPs.

Scatter plot: 1000 Genomes EUR MAF vs gnomAD v4 NFE MAF, with top |ΔAF| loci
highlighted and labelled. Saves PNG (300 DPI) and vector PDF.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

# ---------------------------------------------------------------------------
# Editable constants
# ---------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parent

INPUT_CSV = REPO_ROOT / "analysis" / "la_snp_maf_comparison.csv"
OUTPUT_STEM = REPO_ROOT / "figures" / "af_1kg_vs_gnomad_comparison"

N_LABEL = 8
FIGURE_DPI = 300
MAF_MARGIN = 0.02

REQUIRED_COLUMNS = ("rsid", "gene", "maf_1000g_eur", "maf_gnomad_nfe")

COLORS = {
    "scatter": "#b0b0b0",
    "highlight": "#c45c3e",
    "diagonal": "#666666",
    "label": "#333333",
    "grid": "#dddddd",
}

FONT_SIZES = {
    "base": 11,
    "axis_label": 12,
    "tick": 10,
    "annotation": 9,
}


def _load_data(path: Path) -> pd.DataFrame:
    if not path.is_file():
        print(f"Error: input CSV not found: {path}", file=sys.stderr)
        sys.exit(1)

    df = pd.read_csv(path)
    missing = [col for col in REQUIRED_COLUMNS if col not in df.columns]
    if missing:
        print(
            f"Error: input CSV is missing required column(s): {', '.join(missing)}",
            file=sys.stderr,
        )
        print(f"  Expected columns: {', '.join(REQUIRED_COLUMNS)}", file=sys.stderr)
        print(f"  Found columns: {', '.join(df.columns.astype(str))}", file=sys.stderr)
        sys.exit(1)

    work = df.loc[:, list(REQUIRED_COLUMNS)].copy()
    work["maf_1000g_eur"] = pd.to_numeric(work["maf_1000g_eur"], errors="coerce")
    work["maf_gnomad_nfe"] = pd.to_numeric(work["maf_gnomad_nfe"], errors="coerce")
    work["delta_af"] = work["maf_gnomad_nfe"] - work["maf_1000g_eur"]
    work["abs_delta"] = work["delta_af"].abs()
    return work


def _format_label(gene: object, rsid: object) -> str:
    gene_text = str(gene).strip()
    rsid_text = str(rsid).strip()
    if gene_text and gene_text.lower() != "nan":
        return f"{gene_text} {rsid_text}"
    return rsid_text


def _print_top_table(df: pd.DataFrame, n: int) -> None:
    top = df.nlargest(n, "abs_delta")
    display_cols = ["gene", "rsid", "maf_1000g_eur", "maf_gnomad_nfe", "delta_af", "abs_delta"]
    print(f"\nTop {n} loci by |ΔAF| (gnomAD NFE − 1000G EUR):\n")
    print(top[display_cols].to_string(index=False, float_format=lambda x: f"{x:.4f}"))
    print()


def _annotate_top_loci(ax: plt.Axes, top: pd.DataFrame) -> None:
    try:
        from adjustText import adjust_text
    except ImportError:
        adjust_text = None

    offsets = [(6, 6), (-6, 6), (6, -6), (-6, -6), (10, 0), (-10, 0), (0, 10), (0, -10)]
    texts: list[plt.Text] = []

    for idx, row in enumerate(top.itertuples(index=False)):
        label = _format_label(row.gene, row.rsid)
        if adjust_text is not None:
            text = ax.text(
                row.maf_1000g_eur,
                row.maf_gnomad_nfe,
                label,
                fontsize=FONT_SIZES["annotation"],
                color=COLORS["label"],
            )
            texts.append(text)
        else:
            dx, dy = offsets[idx % len(offsets)]
            ax.annotate(
                label,
                xy=(row.maf_1000g_eur, row.maf_gnomad_nfe),
                xytext=(dx, dy),
                textcoords="offset points",
                fontsize=FONT_SIZES["annotation"],
                color=COLORS["label"],
                ha="left" if dx >= 0 else "right",
                va="bottom" if dy >= 0 else "top",
            )

    if adjust_text is not None and texts:
        adjust_text(
            texts,
            ax=ax,
            arrowprops={"arrowstyle": "-", "color": COLORS["label"], "lw": 0.6},
        )


def _style_axes(ax: plt.Axes) -> None:
    ax.set_xlabel("MAF — 1000G European", fontsize=FONT_SIZES["axis_label"])
    ax.set_ylabel("MAF — gnomAD v4 NFE", fontsize=FONT_SIZES["axis_label"])
    ax.tick_params(labelsize=FONT_SIZES["tick"])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, color=COLORS["grid"], linewidth=0.8, linestyle="-", alpha=0.7)
    ax.set_axisbelow(True)


def plot_comparison(df: pd.DataFrame, output_stem: Path, n_label: int) -> None:
    paired = df.dropna(subset=["maf_1000g_eur", "maf_gnomad_nfe"])
    if paired.empty:
        print("Error: no SNPs with MAF values in both columns.", file=sys.stderr)
        sys.exit(1)

    _print_top_table(paired, n_label)

    plt.rcParams.update(
        {
            "font.size": FONT_SIZES["base"],
            "font.family": "sans-serif",
            "axes.facecolor": "white",
            "figure.facecolor": "white",
        }
    )

    fig, ax = plt.subplots(figsize=(6.5, 6.5), layout="constrained")

    ax.scatter(
        paired["maf_1000g_eur"],
        paired["maf_gnomad_nfe"],
        s=22,
        color=COLORS["scatter"],
        alpha=0.75,
        linewidths=0,
        zorder=1,
    )

    max_maf = float(max(paired["maf_1000g_eur"].max(), paired["maf_gnomad_nfe"].max()))
    upper = min(1.0, max_maf + MAF_MARGIN)
    ax.plot(
        [0, upper],
        [0, upper],
        linestyle="--",
        color=COLORS["diagonal"],
        linewidth=1.0,
        zorder=2,
    )
    ax.set_xlim(0, upper)
    ax.set_ylim(0, upper)
    ax.set_aspect("equal", adjustable="box")

    top = paired.nlargest(n_label, "abs_delta")
    ax.scatter(
        top["maf_1000g_eur"],
        top["maf_gnomad_nfe"],
        s=64,
        color=COLORS["highlight"],
        alpha=0.95,
        linewidths=0.6,
        edgecolors="white",
        zorder=3,
    )
    _annotate_top_loci(ax, top)

    _style_axes(ax)

    output_stem = Path(output_stem)
    output_stem.parent.mkdir(parents=True, exist_ok=True)
    png_path = output_stem.with_suffix(".png")
    pdf_path = output_stem.with_suffix(".pdf")
    fig.savefig(png_path, dpi=FIGURE_DPI, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved: {png_path}")
    print(f"Saved: {pdf_path}")


def main() -> None:
    df = _load_data(Path(INPUT_CSV))
    plot_comparison(df, Path(OUTPUT_STEM), N_LABEL)


if __name__ == "__main__":
    main()
