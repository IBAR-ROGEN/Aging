#!/usr/bin/env python3
"""Figure 1, Panel C: publication-ready static network of LA-SNP mechanisms (networkx + matplotlib).

Layout: hub-and-spoke. Central hub connects to four pathway category nodes; each category
connects to its gene nodes. Gene node area scales with LA-SNP count.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import networkx as nx
import seaborn as sns
import typer
from adjustText import adjust_text
from matplotlib.colors import to_hex

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_OUT_DIR = REPO_ROOT / "analysis"

# --- Hardcoded biological structure (LA-SNPs in neurodegeneration context) ---

HUB_ID = "hub"
HUB_LABEL = "Genetic Resilience\nagainst Neurodegeneration"

# One Seaborn pass → stable hex fills (colorblind, slightly desaturated for print).
_BRANCH_FILL_COLORS: tuple[str, ...] = tuple(
    to_hex(c) for c in sns.color_palette("colorblind", n_colors=4, desat=0.85)
)


@dataclass(frozen=True)
class GeneSpec:
    """One gene leaf: stable id, display label, LA-SNP count."""

    node_id: str
    label: str
    n_snps: int


@dataclass(frozen=True)
class BranchSpec:
    """Primary pathway branch: color name, category node id/label, genes, mechanistic annotation."""

    branch_id: str
    category_label: str
    color_hex: str
    genes: tuple[GeneSpec, ...]
    annotation: str


BRANCHES: tuple[BranchSpec, ...] = (
    BranchSpec(
        branch_id="lipid",
        category_label="Lipid Metabolism",
        color_hex=_BRANCH_FILL_COLORS[0],
        genes=(
            GeneSpec("CETP", "CETP\n(n=10)", 10),
            GeneSpec("APOC1", "APOC1\n(n=2)", 2),
        ),
        annotation="Altered kinetics /\nHigher HDL-C",
    ),
    BranchSpec(
        branch_id="proteostasis",
        category_label="Protein Homeostasis",
        color_hex=_BRANCH_FILL_COLORS[1],
        genes=(
            GeneSpec("HSPA1A", "HSPA1A\n(n=3)", 3),
            GeneSpec("HSPA1B", "HSPA1B\n(n=3)", 3),
            GeneSpec("HSPA1L", "HSPA1L\n(n=2)", 2),
        ),
        annotation="Enhanced folding\nefficiency",
    ),
    BranchSpec(
        branch_id="mito",
        category_label="Mitochondrial\nFunction",
        color_hex=_BRANCH_FILL_COLORS[2],
        genes=(GeneSpec("NDUFS1", "NDUFS1\n(n=5)", 5),),
        annotation="Optimized electron transfer /\nReduced ROS",
    ),
    BranchSpec(
        branch_id="immune",
        category_label="Immune Regulation",
        color_hex=_BRANCH_FILL_COLORS[3],
        genes=(
            GeneSpec("HLA_DQB1", "HLA-DQB1\n(n=2)", 2),
            GeneSpec("NLRC5", "NLRC5\n(n=1)", 1),
            GeneSpec("SDC4", "SDC4\n(n=2)", 2),
        ),
        annotation="Balanced response /\nAttenuated inflammaging",
    ),
)


def _hub_and_spoke_positions(
    *,
    r_pathway: float,
    r_gene: float,
    gene_spread_rad: float,
) -> dict[str, tuple[float, float]]:
    """Place hub at origin; four pathways on a circle; genes fan slightly around each spoke.

    Branch index 0..3 maps to angles starting at the top (pi/2) and proceeding clockwise
    so Lipid sits at 12 o'clock, Protein at 9, Mitochondria at 6, Immune at 3.
    """
    pos: dict[str, tuple[float, float]] = {HUB_ID: (0.0, 0.0)}
    n_branches = len(BRANCHES)

    for i, br in enumerate(BRANCHES):
        # Clockwise from top: theta = pi/2 - 2*pi*i/n
        theta_c = math.pi / 2 - (2 * math.pi * i) / n_branches
        px = r_pathway * math.cos(theta_c)
        py = r_pathway * math.sin(theta_c)
        cat_node = f"cat_{br.branch_id}"
        pos[cat_node] = (px, py)

        k = len(br.genes)
        if k == 1:
            offsets = (0.0,)
        else:
            half = (k - 1) / 2
            offsets = tuple((j - half) * gene_spread_rad for j in range(k))

        for off, g in zip(offsets, br.genes, strict=True):
            theta_g = theta_c + off
            gx = r_gene * math.cos(theta_g)
            gy = r_gene * math.sin(theta_g)
            pos[g.node_id] = (gx, gy)

    return pos


def _build_graph() -> nx.Graph:
    """Undirected conceptual diagram; edges represent mechanistic grouping."""
    g = nx.Graph()
    g.add_node(HUB_ID, kind="hub", label=HUB_LABEL)

    for br in BRANCHES:
        cat = f"cat_{br.branch_id}"
        g.add_node(
            cat,
            kind="category",
            label=br.category_label,
            color=br.color_hex,
        )
        g.add_edge(HUB_ID, cat)
        for gene in br.genes:
            g.add_node(
                gene.node_id,
                kind="gene",
                label=gene.label,
                n_snps=gene.n_snps,
                color=br.color_hex,
            )
            g.add_edge(cat, gene.node_id)
    return g


def _gene_node_size(n_snps: int, *, snp_min: int, snp_max: int) -> float:
    """Map SNP count to matplotlib scatter *area* (node_size), monotonic and visible range."""
    lo, hi = 520.0, 4200.0
    if snp_max <= snp_min:
        return (lo + hi) / 2
    t = (n_snps - snp_min) / (snp_max - snp_min)
    return lo + t * (hi - lo)


def _configure_matplotlib_fonts() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica Neue", "Helvetica", "DejaVu Sans"],
            "font.size": 9.0,
            "axes.linewidth": 0.0,
        }
    )


def render_figure1c_mechanisms(
    out_dir: Path,
    *,
    dpi: int = 300,
    figsize_in: tuple[float, float] = (7.5, 7.5),
) -> tuple[Path, Path]:
    """Draw the network and write PNG (raster) + PDF (vector) with tight bounding boxes."""
    _configure_matplotlib_fonts()
    out_dir = out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    png_path = out_dir / "Figure1C_Mechanisms.png"
    pdf_path = out_dir / "Figure1C_Mechanisms.pdf"

    g = _build_graph()
    pos = _hub_and_spoke_positions(r_pathway=1.85, r_gene=3.55, gene_spread_rad=0.14)

    all_snps = [d["n_snps"] for _, d in g.nodes(data=True) if d.get("kind") == "gene"]
    snp_min, snp_max = min(all_snps), max(all_snps)

    fig, ax = plt.subplots(figsize=figsize_in, dpi=dpi, facecolor="white")
    ax.set_facecolor("white")
    ax.axis("off")
    ax.set_aspect("equal", adjustable="datalim")

    # Edges: thin neutral lines behind nodes
    nx.draw_networkx_edges(
        g,
        pos,
        ax=ax,
        width=1.15,
        edge_color="#cbd5e1",
        alpha=0.95,
    )

    # Nodes by kind
    hub_nodes = [n for n, d in g.nodes(data=True) if d.get("kind") == "hub"]
    cat_nodes = [n for n, d in g.nodes(data=True) if d.get("kind") == "category"]
    gene_nodes = [n for n, d in g.nodes(data=True) if d.get("kind") == "gene"]

    nx.draw_networkx_nodes(
        g,
        pos,
        nodelist=hub_nodes,
        ax=ax,
        node_size=9800,
        node_color="#1e293b",
        linewidths=1.2,
        edgecolors="#0f172a",
    )
    cat_colors = [g.nodes[n]["color"] for n in cat_nodes]
    nx.draw_networkx_nodes(
        g,
        pos,
        nodelist=cat_nodes,
        ax=ax,
        node_size=4200,
        node_color=cat_colors,
        linewidths=1.0,
        edgecolors="#334155",
        alpha=0.92,
    )
    gene_colors = [g.nodes[n]["color"] for n in gene_nodes]
    gene_sizes = [_gene_node_size(g.nodes[n]["n_snps"], snp_min=snp_min, snp_max=snp_max) for n in gene_nodes]
    nx.draw_networkx_nodes(
        g,
        pos,
        nodelist=gene_nodes,
        ax=ax,
        node_size=gene_sizes,
        node_color=gene_colors,
        linewidths=0.9,
        edgecolors="#475569",
        alpha=0.95,
    )

    labels: dict[str, str] = {n: d["label"] for n, d in g.nodes(data=True) if "label" in d}
    for nid, text in labels.items():
        kind = g.nodes[nid].get("kind")
        if kind == "hub":
            color, size, weight = "#f8fafc", 8.6, "700"
        else:
            color, size, weight = "#0f172a", 7.9, "600"
        x, y = pos[nid]
        ax.text(
            x,
            y,
            text,
            ha="center",
            va="center",
            fontsize=size,
            fontweight=weight,
            color=color,
            zorder=6,
            linespacing=0.95,
        )

    # Mechanistic annotations near each pathway cluster (initial guess positions; adjust_text nudges them)
    annotation_texts: list[Any] = []
    n_branches = len(BRANCHES)
    r_ann = 2.65
    for i, br in enumerate(BRANCHES):
        theta_c = math.pi / 2 - (2 * math.pi * i) / n_branches
        # Slight outward bias so text sits between pathway and gene ring
        ax_off = 0.42 * math.cos(theta_c)
        ay_off = 0.42 * math.sin(theta_c)
        tx = r_ann * math.cos(theta_c) + ax_off
        ty = r_ann * math.sin(theta_c) + ay_off
        t = ax.text(
            tx,
            ty,
            br.annotation,
            ha="center",
            va="center",
            fontsize=7.8,
            color="#334155",
            linespacing=1.15,
            zorder=5,
        )
        annotation_texts.append(t)

    # Reposition only (no leader lines): avoids FancyArrowPatch warnings with equal aspect + tight bbox.
    adjust_text(
        annotation_texts,
        ax=ax,
        expand=(1.35, 1.45),
        force_text=(0.45, 0.55),
        force_points=(0.2, 0.25),
        lim=520,
        arrowprops=None,
    )

    fig.tight_layout(pad=0.15)
    for path, fmt in ((png_path, "png"), (pdf_path, "pdf")):
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
    return png_path, pdf_path


app = typer.Typer(add_completion=False, no_args_is_help=True)


@app.command()
def main(
    out_dir: Path = typer.Option(
        DEFAULT_OUT_DIR,
        "--out-dir",
        "-o",
        help="Directory for Figure1C_Mechanisms.png and .pdf.",
    ),
    dpi: int = typer.Option(300, "--dpi", min=72, max=600, help="PNG resolution (DPI)."),
) -> None:
    png, pdf = render_figure1c_mechanisms(out_dir, dpi=dpi)
    typer.echo(f"Wrote {png}")
    typer.echo(f"Wrote {pdf}")


if __name__ == "__main__":
    app()
