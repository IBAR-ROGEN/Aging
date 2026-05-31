#!/usr/bin/env python3
"""Activity 2.1.7.1 — hub-and-spoke network of LA genes clustered by functional pathway.

Reads ``overlapping_genes_with_snps.xlsx`` (columns ``Gene``, ``SNP_rsID``). Pathway
assignments come from ``--pathway-map`` CSV (``Gene``, ``Pathway``) or fall back to four
hardcoded pathway groups below. Each pathway is a hub node; genes are spokes scaled by
unique LA-SNP count. Writes PNG (300 DPI) and PDF to ``analysis/``.
"""

from __future__ import annotations

import argparse
import logging
import math
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
DEFAULT_INPUT = REPO_ROOT / "overlapping_genes_with_snps.xlsx"
DEFAULT_OUTPUT = REPO_ROOT / "analysis" / "Fig_LA_SNP_network.png"

# Fallback pathway → gene lists (Activity 2.1.7.1). Extend to cover all 41 genes.
DEFAULT_PATHWAY_GENES: dict[str, tuple[str, ...]] = {
    "Proteostasis": (
        "HSPA1A",
        "HSPA1B",
        "HSPA1L",
    ),
    "Lipid Metabolism": (
        "CETP",
        "APOC1",
    ),
    "Mitochondrial Integrity": (
        "NDUFS1",
    ),
    "Immune Regulation": (
        "HLA-DQB1",
        "NLRC5",
        "SDC4",
    ),
}

PATHWAY_PALETTE: tuple[str, ...] = (
    "#5B7C99",
    "#5B8C5A",
    "#C4A35A",
    "#C45C3E",
    "#7B6FA6",
    "#6B7280",
)
UNASSIGNED_PATHWAY = "Unassigned"
UNASSIGNED_COLOR = "#94a3b8"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", type=Path, default=DEFAULT_INPUT, help="Excel path (.xlsx)")
    p.add_argument(
        "--pathway-map",
        type=Path,
        default=None,
        help="CSV with Gene, Pathway columns (overrides hardcoded groups)",
    )
    p.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help="PNG output path (companion PDF written alongside)",
    )
    p.add_argument("--gene-column", default="Gene", help="Gene name column")
    p.add_argument("--snp-column", default="SNP_rsID", help="SNP identifier column")
    return p.parse_args()


def _load_snp_counts(
    input_path: Path,
    *,
    gene_column: str,
    snp_column: str,
) -> pd.Series:
    df = pd.read_excel(input_path)
    missing = {gene_column, snp_column} - set(df.columns)
    if missing:
        msg = f"Missing columns {sorted(missing)} (available: {list(df.columns)})"
        raise ValueError(msg)
    counts = df.groupby(gene_column, sort=False)[snp_column].nunique()
    return counts.sort_values(ascending=False)


def _load_pathway_map_csv(path: Path) -> dict[str, str]:
    df = pd.read_csv(path)
    missing = {"Gene", "Pathway"} - set(df.columns)
    if missing:
        msg = f"Pathway map missing columns {sorted(missing)} (available: {list(df.columns)})"
        raise ValueError(msg)
    gene_to_pathway: dict[str, str] = {}
    for gene, pathway in zip(df["Gene"], df["Pathway"], strict=True):
        gene_to_pathway[str(gene).strip()] = str(pathway).strip()
    return gene_to_pathway


def _default_gene_to_pathway() -> dict[str, str]:
    mapping: dict[str, str] = {}
    for pathway, genes in DEFAULT_PATHWAY_GENES.items():
        for gene in genes:
            mapping[gene] = pathway
    return mapping


def _pathway_color_map(pathways: list[str]) -> dict[str, str]:
    colors: dict[str, str] = {}
    palette_idx = 0
    for pathway in sorted(pathways):
        if pathway == UNASSIGNED_PATHWAY:
            colors[pathway] = UNASSIGNED_COLOR
            continue
        colors[pathway] = PATHWAY_PALETTE[palette_idx % len(PATHWAY_PALETTE)]
        palette_idx += 1
    return colors


def _group_genes_by_pathway(
    snp_counts: pd.Series,
    gene_to_pathway: dict[str, str],
) -> tuple[dict[str, list[str]], list[str]]:
    unmapped: list[str] = []
    grouped: dict[str, list[str]] = {}

    for gene in snp_counts.index.astype(str):
        pathway = gene_to_pathway.get(gene)
        if pathway is None:
            unmapped.append(gene)
            pathway = UNASSIGNED_PATHWAY
        grouped.setdefault(pathway, []).append(gene)

    for pathway in grouped:
        grouped[pathway].sort(key=lambda g: (-int(snp_counts[g]), g))
    return grouped, unmapped


def _hub_spoke_positions(
    grouped: dict[str, list[str]],
    *,
    r_hub: float = 2.05,
    r_gene: float = 3.85,
    gene_spread_rad: float = 0.13,
) -> dict[str, tuple[float, float]]:
    """Place pathway hubs on a circle; genes fan outward along each spoke."""
    pos: dict[str, tuple[float, float]] = {}
    pathways = sorted(grouped.keys())
    n_pathways = len(pathways)

    for i, pathway in enumerate(pathways):
        theta_c = math.pi / 2 - (2 * math.pi * i) / n_pathways
        hub_id = f"hub::{pathway}"
        pos[hub_id] = (r_hub * math.cos(theta_c), r_hub * math.sin(theta_c))

        genes = grouped[pathway]
        k = len(genes)
        if k == 1:
            offsets = (0.0,)
        else:
            half = (k - 1) / 2
            offsets = tuple((j - half) * gene_spread_rad for j in range(k))

        for off, gene in zip(offsets, genes, strict=True):
            theta_g = theta_c + off
            pos[gene] = (r_gene * math.cos(theta_g), r_gene * math.sin(theta_g))

    return pos


def _gene_node_size(n_snps: int, *, snp_min: int, snp_max: int) -> float:
    lo, hi = 480.0, 4000.0
    if snp_max <= snp_min:
        return (lo + hi) / 2
    t = (n_snps - snp_min) / (snp_max - snp_min)
    return lo + t * (hi - lo)


def _build_graph(
    grouped: dict[str, list[str]],
    snp_counts: pd.Series,
    pathway_colors: dict[str, str],
) -> nx.Graph:
    g = nx.Graph()
    for pathway, genes in grouped.items():
        hub_id = f"hub::{pathway}"
        g.add_node(
            hub_id,
            kind="pathway",
            label=pathway.replace(" ", "\n"),
            color=pathway_colors[pathway],
        )
        for gene in genes:
            n_snps = int(snp_counts[gene])
            g.add_node(
                gene,
                kind="gene",
                label=f"{gene}\n(n={n_snps})",
                n_snps=n_snps,
                color=pathway_colors[pathway],
            )
            g.add_edge(hub_id, gene)
    return g


def _configure_matplotlib() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica Neue", "Helvetica", "DejaVu Sans"],
            "font.size": 9.0,
            "axes.facecolor": "white",
            "figure.facecolor": "white",
            "axes.linewidth": 0.0,
        }
    )


def render_network(
    snp_counts: pd.Series,
    grouped: dict[str, list[str]],
    *,
    output_png: Path,
    dpi: int = 300,
    figsize_in: tuple[float, float] = (8.0, 8.0),
) -> tuple[Path, Path]:
    _configure_matplotlib()
    output_png = output_png.resolve()
    output_png.parent.mkdir(parents=True, exist_ok=True)
    output_pdf = output_png.with_suffix(".pdf")

    pathway_colors = _pathway_color_map(list(grouped.keys()))
    graph = _build_graph(grouped, snp_counts, pathway_colors)
    pos = _hub_spoke_positions(grouped)

    gene_snps = [int(snp_counts[gene]) for gene in snp_counts.index.astype(str)]
    snp_min, snp_max = min(gene_snps), max(gene_snps)

    fig, ax = plt.subplots(figsize=figsize_in, dpi=dpi, facecolor="white")
    ax.set_facecolor("white")
    ax.axis("off")
    ax.set_aspect("equal", adjustable="datalim")

    nx.draw_networkx_edges(
        graph,
        pos,
        ax=ax,
        width=1.1,
        edge_color="#cbd5e1",
        alpha=0.95,
    )

    pathway_nodes = [n for n, d in graph.nodes(data=True) if d.get("kind") == "pathway"]
    gene_nodes = [n for n, d in graph.nodes(data=True) if d.get("kind") == "gene"]

    nx.draw_networkx_nodes(
        graph,
        pos,
        nodelist=pathway_nodes,
        ax=ax,
        node_size=4600,
        node_color=[graph.nodes[n]["color"] for n in pathway_nodes],
        linewidths=1.0,
        edgecolors="#334155",
        alpha=0.93,
    )
    nx.draw_networkx_nodes(
        graph,
        pos,
        nodelist=gene_nodes,
        ax=ax,
        node_size=[
            _gene_node_size(graph.nodes[n]["n_snps"], snp_min=snp_min, snp_max=snp_max)
            for n in gene_nodes
        ],
        node_color=[graph.nodes[n]["color"] for n in gene_nodes],
        linewidths=0.9,
        edgecolors="#475569",
        alpha=0.95,
    )

    for node_id, data in graph.nodes(data=True):
        x, y = pos[node_id]
        if data.get("kind") == "pathway":
            color, size, weight = "#0f172a", 8.2, "700"
        else:
            color, size, weight = "#0f172a", 7.4, "600"
        ax.text(
            x,
            y,
            data["label"],
            ha="center",
            va="center",
            fontsize=size,
            fontweight=weight,
            color=color,
            zorder=6,
            linespacing=0.95,
        )

    n_genes = len(gene_nodes)
    n_snps = int(snp_counts.sum())
    ax.set_title(
        f"LA-SNP network: {n_genes} genes, {n_snps} unique SNPs by functional pathway",
        fontsize=10.5,
        fontweight="600",
        color="#0f172a",
        pad=12,
    )

    fig.tight_layout(pad=0.15)
    for path, fmt in ((output_png, "png"), (output_pdf, "pdf")):
        fig.savefig(
            path,
            format=fmt,
            dpi=dpi if fmt == "png" else None,
            bbox_inches="tight",
            pad_inches=0.03,
            facecolor="white",
            edgecolor="none",
        )
    plt.close(fig)
    return output_png, output_pdf


def main() -> None:
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    snp_counts = _load_snp_counts(
        args.input,
        gene_column=args.gene_column,
        snp_column=args.snp_column,
    )
    logging.info(
        "Loaded %d genes, %d unique SNPs from %s",
        len(snp_counts),
        int(snp_counts.sum()),
        args.input,
    )

    if args.pathway_map is not None:
        gene_to_pathway = _load_pathway_map_csv(args.pathway_map)
        logging.info("Pathway assignments from %s (%d entries)", args.pathway_map, len(gene_to_pathway))
    else:
        gene_to_pathway = _default_gene_to_pathway()
        logging.info("Using hardcoded pathway groups (%d gene assignments)", len(gene_to_pathway))

    grouped, unmapped = _group_genes_by_pathway(snp_counts, gene_to_pathway)
    if unmapped:
        logging.warning(
            "%d gene(s) not in pathway map → %s: %s",
            len(unmapped),
            UNASSIGNED_PATHWAY,
            ", ".join(unmapped),
        )

    png_path, pdf_path = render_network(snp_counts, grouped, output_png=args.output)
    logging.info("Wrote %s", png_path)
    logging.info("Wrote %s", pdf_path)


if __name__ == "__main__":
    main()
