#!/usr/bin/env python3
"""Publication figure: 1000 Genomes EUR vs gnomAD v4 NFE allele-frequency comparison.

Panel A — scatter of AF_1kg vs AF_gnomad_nfe with identity line and top-discrepancy labels.
Panel B — horizontal bar chart of the largest |ΔAF| loci.

Gene names are merged from the UKB LA-SNP manifest (``Gene`` + ``SNP_rsID``).
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Configurable constants
# ---------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parent.parent.parent

INPUT_CSV = REPO_ROOT / "analysis" / "la_snp_af_1kg_vs_gnomad.csv"
MANIFEST_CSV = REPO_ROOT / "analysis" / "ukb_snp_manifest_v0.1.csv"
OUTPUT_DIR = REPO_ROOT / "figures"
FIG_BASENAME = "af_1kg_vs_gnomad_comparison"

DIFF_THRESHOLD = 0.05
FIGURE_DPI = 300
TOP_N_LABELS = 12

COLOR_CONCORDANT = "#888888"
COLOR_DISCREPANT = "#C45C3E"
FONT_SIZE = 11


def normalize_rsid(value: object) -> str:
    """Return a normalized rsID string (``rs`` prefix, stripped)."""
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return ""
    return text if text.lower().startswith("rs") else f"rs{text}"


def load_gene_map(manifest_path: Path) -> pd.Series:
    """Map rsID -> gene name(s) from the LA-SNP manifest."""
    manifest = pd.read_csv(manifest_path)
    required = {"Gene", "SNP_rsID"}
    missing = required - set(manifest.columns)
    if missing:
        raise ValueError(f"Manifest CSV missing columns: {sorted(missing)}")

    work = manifest[["Gene", "SNP_rsID"]].copy()
    work["rsID"] = work["SNP_rsID"].map(normalize_rsid)
    work["Gene"] = work["Gene"].astype(str).str.strip()
    work = work.loc[work["rsID"].ne("") & work["Gene"].ne("") & work["Gene"].ne("nan")]

    grouped = (
        work.groupby("rsID", sort=False)["Gene"]
        .apply(lambda genes: ", ".join(dict.fromkeys(genes.astype(str))))
        .rename("Gene")
    )
    return grouped


def load_comparison_table(path: Path) -> pd.DataFrame:
    """Load the 1KG vs gnomAD comparison CSV."""
    df = pd.read_csv(path)
    required = {"rsID", "AF_1kg", "AF_gnomad_nfe"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Comparison CSV missing columns: {sorted(missing)}")

    df = df.copy()
    df["rsID"] = df["rsID"].map(normalize_rsid)
    return df


def merge_genes(comparison: pd.DataFrame, gene_map: pd.Series) -> pd.DataFrame:
    """Attach gene names from the manifest; keep an existing ``Gene`` column when present."""
    out = comparison.copy()
    if "Gene" not in out.columns:
        out["Gene"] = out["rsID"].map(gene_map)
    else:
        manifest_genes = out["rsID"].map(gene_map)
        out["Gene"] = out["Gene"].fillna(manifest_genes)
    return out


def ensure_diff_columns(df: pd.DataFrame, threshold: float) -> pd.DataFrame:
    """Ensure ``abs_diff`` and ``large_diff`` exist, recomputing when needed."""
    out = df.copy()
    paired = out["AF_1kg"].notna() & out["AF_gnomad_nfe"].notna()

    if "abs_diff" not in out.columns:
        out["abs_diff"] = np.nan
    out.loc[paired, "abs_diff"] = (out.loc[paired, "AF_1kg"] - out.loc[paired, "AF_gnomad_nfe"]).abs()

    if "large_diff" not in out.columns:
        out["large_diff"] = False
    out.loc[paired, "large_diff"] = out.loc[paired, "abs_diff"] > threshold
    return out


def format_locus_label(row: pd.Series) -> str:
    """Build a scatter annotation label from rsID and optional gene name."""
    rsid = str(row["rsID"])
    gene = row.get("Gene")
    if pd.notna(gene) and str(gene).strip() and str(gene).strip().lower() != "nan":
        return f"{rsid} ({gene})"
    return rsid


def print_coverage_summary(df: pd.DataFrame, threshold: float) -> None:
    """Print headline counts for the prioritized SNP set."""
    total = len(df)
    paired_mask = df["AF_1kg"].notna() & df["AF_gnomad_nfe"].notna()
    paired = int(paired_mask.sum())
    flagged = int(df.loc[paired_mask, "large_diff"].sum())

    print("Allele-frequency comparison coverage")
    print(f"  Total SNPs in prioritized set: {total}")
    print(f"  SNPs with usable AF in both sources: {paired}")
    print(f"  SNPs flagged (|ΔAF| > {threshold:g}): {flagged}")


def plot_scatter_panel(ax: plt.Axes, paired: pd.DataFrame, top_n: int) -> None:
    """Panel A: AF scatter with identity line and top-discrepancy labels."""
    concordant = paired.loc[~paired["large_diff"]]
    discrepant = paired.loc[paired["large_diff"]]

    ax.scatter(
        concordant["AF_1kg"],
        concordant["AF_gnomad_nfe"],
        s=36,
        color=COLOR_CONCORDANT,
        alpha=0.85,
        linewidths=0,
        label=f"|ΔAF| ≤ {DIFF_THRESHOLD:g}",
    )
    if not discrepant.empty:
        ax.scatter(
            discrepant["AF_1kg"],
            discrepant["AF_gnomad_nfe"],
            s=48,
            color=COLOR_DISCREPANT,
            alpha=0.95,
            linewidths=0,
            label=f"|ΔAF| > {DIFF_THRESHOLD:g}",
        )

    max_af = float(max(paired["AF_1kg"].max(), paired["AF_gnomad_nfe"].max()))
    upper = min(1.0, max_af + 0.02)
    ax.plot([0, upper], [0, upper], linestyle="--", color="#666666", linewidth=1.0, label="y = x")
    ax.set_xlim(0, upper)
    ax.set_ylim(0, upper)
    ax.set_xlabel("1000 Genomes EUR AF")
    ax.set_ylabel("gnomAD v4 NFE AF")
    ax.set_title("Allele frequency: 1000 Genomes EUR vs gnomAD v4 NFE")
    ax.set_aspect("equal", adjustable="box")
    ax.legend(loc="lower right", frameon=False, fontsize=FONT_SIZE - 2)

    top = paired.sort_values("abs_diff", ascending=False).head(top_n)
    offsets = [(6, 6), (-6, 6), (6, -6), (-6, -6), (10, 0), (-10, 0), (0, 10), (0, -10)]
    for idx, (_, row) in enumerate(top.iterrows()):
        dx, dy = offsets[idx % len(offsets)]
        ax.annotate(
            format_locus_label(row),
            xy=(row["AF_1kg"], row["AF_gnomad_nfe"]),
            xytext=(dx, dy),
            textcoords="offset points",
            fontsize=FONT_SIZE - 2,
            color="#333333",
            ha="left" if dx >= 0 else "right",
            va="bottom" if dy >= 0 else "top",
        )


def plot_ranked_panel(ax: plt.Axes, paired: pd.DataFrame, top_n: int) -> None:
    """Panel B: horizontal bars of the top |ΔAF| loci."""
    top = paired.sort_values("abs_diff", ascending=False).head(top_n).iloc[::-1]
    y_pos = np.arange(len(top))
    colors = [COLOR_DISCREPANT if flag else COLOR_CONCORDANT for flag in top["large_diff"]]

    ax.barh(y_pos, top["abs_diff"].to_numpy(), color=colors, edgecolor="none", height=0.72)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(top["rsID"].astype(str), fontsize=FONT_SIZE - 1)
    ax.set_xlabel("|ΔAF| (1000 Genomes − gnomAD v4 NFE)")
    ax.set_title("Largest allele-frequency discrepancies")
    ax.axvline(DIFF_THRESHOLD, color="#666666", linestyle="--", linewidth=1.0, label=f"Threshold ({DIFF_THRESHOLD:g})")
    ax.legend(loc="lower right", frameon=False, fontsize=FONT_SIZE - 2)
    ax.margins(x=0.08)


def main() -> None:
    plt.rcParams.update(
        {
            "font.size": FONT_SIZE,
            "font.family": "sans-serif",
            "axes.facecolor": "white",
            "figure.facecolor": "white",
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )

    comparison = load_comparison_table(INPUT_CSV)
    gene_map = load_gene_map(MANIFEST_CSV)
    comparison = merge_genes(comparison, gene_map)
    comparison = ensure_diff_columns(comparison, DIFF_THRESHOLD)

    print_coverage_summary(comparison, DIFF_THRESHOLD)

    paired = comparison.dropna(subset=["AF_1kg", "AF_gnomad_nfe"]).copy()
    if paired.empty:
        raise ValueError("No SNPs with allele frequencies in both 1000 Genomes and gnomAD.")

    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.8), constrained_layout=True)
    plot_scatter_panel(axes[0], paired, TOP_N_LABELS)
    plot_ranked_panel(axes[1], paired, TOP_N_LABELS)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    png_path = OUTPUT_DIR / f"{FIG_BASENAME}.png"
    pdf_path = OUTPUT_DIR / f"{FIG_BASENAME}.pdf"
    fig.savefig(png_path, dpi=FIGURE_DPI, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)

    print(f"PNG: {png_path}")
    print(f"PDF: {pdf_path}")


if __name__ == "__main__":
    main()
