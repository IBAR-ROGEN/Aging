#!/usr/bin/env python3
"""Render the longevity conceptual network diagram to PNG (matplotlib).

Mirrors the structure and color semantics of ``LongevityNetworkDiagram.tsx``.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import typer
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_OUT = REPO_ROOT / "figures" / "longevity_network_diagram.png"

COLORS = {
    "core": {"fill": "#eff6ff", "accent": "#3b82f6", "text": "#1e293b"},
    "pathophysiology": {"fill": "#eff6ff", "accent": "#3b82f6", "text": "#1e293b"},
    "intersection": {"fill": "#eff6ff", "accent": "#3b82f6", "text": "#1e293b"},
    "mechanism": {"fill": "#fffbeb", "accent": "#f59e0b", "text": "#1e293b"},
    "axis": {"fill": "#ecfdf5", "accent": "#10b981", "text": "#1e293b"},
    "genes": {"fill": "#faf5ff", "accent": "#a855f7", "text": "#1e293b"},
    "outcome": {"fill": "#ecfdf5", "accent": "#10b981", "text": "#1e293b"},
}


@dataclass(frozen=True)
class NodeSpec:
    nid: str
    cx: float
    cy: float
    w: float
    h: float
    lines: tuple[str, ...]
    category: str


@dataclass(frozen=True)
class EdgeSpec:
    src: str
    tgt: str
    dashed: bool
    thick: bool


NODES: tuple[NodeSpec, ...] = (
    NodeSpec(
        "core",
        7.15,
        10.35,
        3.5,
        0.72,
        ("Exceptional Longevity:", "Active Genetic Resilience"),
        "core",
    ),
    NodeSpec(
        "adpd",
        7.0,
        9.35,
        2.6,
        0.55,
        ("AD & PD Pathophysiology",),
        "pathophysiology",
    ),
    NodeSpec(
        "intersection",
        6.2,
        8.2,
        3.4,
        0.62,
        ("Dysregulated Network: 41 Genes, 70 LA-SNPs",),
        "intersection",
    ),
    NodeSpec(
        "qualitative",
        4.0,
        6.85,
        3.0,
        0.78,
        (
            "Qualitative Effects",
            "(Protein Efficiency & Stability)",
        ),
        "mechanism",
    ),
    NodeSpec(
        "eqtl",
        9.9,
        6.85,
        2.4,
        0.62,
        ("eQTL Effects (Gene Expression)",),
        "mechanism",
    ),
    NodeSpec("axis1", 1.35, 5.35, 2.0, 0.5, ("Protein Homeostasis",), "axis"),
    NodeSpec(
        "genes1",
        1.35,
        4.25,
        2.0,
        0.55,
        ("HSPA1A, HSPA1B, HSPA1L",),
        "genes",
    ),
    NodeSpec(
        "outcome1",
        1.35,
        3.1,
        2.0,
        0.55,
        ("Proteotoxic stress resilience",),
        "outcome",
    ),
    NodeSpec(
        "axis2",
        3.95,
        5.35,
        2.35,
        0.55,
        ("Lipid Metabolism & Mitochondria",),
        "axis",
    ),
    NodeSpec(
        "genes2",
        3.95,
        4.25,
        2.35,
        0.55,
        ("CETP (rs5882), NDUFS1",),
        "genes",
    ),
    NodeSpec(
        "outcome2",
        3.95,
        3.1,
        2.35,
        0.55,
        ("Altered HDL & limits ROS",),
        "outcome",
    ),
    NodeSpec("axis3", 6.55, 5.35, 1.85, 0.5, ("Immune Regulation",), "axis"),
    NodeSpec(
        "genes3",
        6.55,
        4.25,
        1.85,
        0.55,
        ("NLRC5, HLA-DQB1",),
        "genes",
    ),
    NodeSpec(
        "outcome3",
        6.55,
        3.1,
        1.85,
        0.62,
        ("Attenuates 'Inflammaging'",),
        "outcome",
    ),
    NodeSpec("axis4", 9.15, 5.35, 1.95, 0.5, ("Antioxidant Factors",), "axis"),
    NodeSpec(
        "genes4",
        9.15,
        4.25,
        1.95,
        0.5,
        ("HMOX1, GPX1",),
        "genes",
    ),
    NodeSpec(
        "outcome4",
        9.15,
        3.1,
        1.95,
        0.55,
        ("Preserves mitochondrial structure",),
        "outcome",
    ),
)

EDGES: tuple[EdgeSpec, ...] = (
    EdgeSpec("core", "adpd", False, False),
    EdgeSpec("adpd", "intersection", False, False),
    EdgeSpec("intersection", "qualitative", False, True),
    EdgeSpec("intersection", "eqtl", True, False),
    EdgeSpec("qualitative", "axis1", False, False),
    EdgeSpec("qualitative", "axis2", False, False),
    EdgeSpec("qualitative", "axis3", False, False),
    EdgeSpec("qualitative", "axis4", False, False),
    EdgeSpec("axis1", "genes1", False, False),
    EdgeSpec("genes1", "outcome1", False, False),
    EdgeSpec("axis2", "genes2", False, False),
    EdgeSpec("genes2", "outcome2", False, False),
    EdgeSpec("axis3", "genes3", False, False),
    EdgeSpec("genes3", "outcome3", False, False),
    EdgeSpec("axis4", "genes4", False, False),
    EdgeSpec("genes4", "outcome4", False, False),
)


def _node_map() -> dict[str, NodeSpec]:
    return {n.nid: n for n in NODES}


def _box_bounds(n: NodeSpec) -> tuple[float, float, float, float]:
    return (n.cx - n.w / 2, n.cy - n.h / 2, n.w, n.h)


def render_longevity_network_diagram(out_path: Path, *, dpi: int = 200) -> Path:
    out_path = out_path.resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)

    nm = _node_map()
    fig_w, fig_h = 14.0, 11.5
    fig, ax = plt.subplots(
        figsize=(fig_w, fig_h),
        dpi=dpi,
        facecolor="#f8fafc",
    )
    ax.set_xlim(0, fig_w)
    ax.set_ylim(0, fig_h)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.patch.set_facecolor("#f1f5f9")

    for n in NODES:
        pal = COLORS[n.category]
        x, y, w, h = _box_bounds(n)
        pad = 0.02
        box = FancyBboxPatch(
            (x + pad, y + pad),
            w - 2 * pad,
            h - 2 * pad,
            boxstyle="round,pad=0.02,rounding_size=0.08",
            facecolor=pal["fill"],
            edgecolor="#e2e8f0",
            linewidth=1.0,
            mutation_aspect=0.35,
        )
        ax.add_patch(box)
        accent = FancyBboxPatch(
            (x + pad, y + pad),
            0.07,
            h - 2 * pad,
            boxstyle="round,pad=0,rounding_size=0.02",
            facecolor=pal["accent"],
            edgecolor="none",
        )
        ax.add_patch(accent)

        fs = 9.2 if len(n.lines) > 1 or max(len(s) for s in n.lines) > 34 else 9.8
        line_h = 0.115
        start_y = n.cy + (len(n.lines) - 1) * line_h / 2
        for i, line in enumerate(n.lines):
            ax.text(
                n.cx + 0.06,
                start_y - i * line_h,
                line,
                ha="center",
                va="center",
                fontsize=fs,
                color=pal["text"],
                fontweight="600" if i == 0 else "500",
            )

    for e in EDGES:
        a, b = nm[e.src], nm[e.tgt]
        x1, y1 = a.cx, a.cy - a.h / 2
        x2, y2 = b.cx, b.cy + b.h / 2
        lw = 2.2 if e.thick else 1.35
        ls = (0, (5, 4)) if e.dashed else "solid"
        arr = FancyArrowPatch(
            (x1, y1),
            (x2, y2),
            arrowstyle="-|>",
            mutation_scale=12,
            linewidth=lw,
            color="#64748b",
            linestyle=ls,
            shrinkA=2,
            shrinkB=2,
            clip_on=False,
        )
        ax.add_patch(arr)

    legend_y = 0.55
    ax.text(
        0.35,
        legend_y,
        "Core / intersection",
        fontsize=8,
        color="#475569",
        va="center",
    )
    ax.add_patch(
        plt.Rectangle((0.12, legend_y - 0.08), 0.14, 0.16, color="#3b82f6", zorder=5)
    )
    ax.text(
        2.15,
        legend_y,
        "Mechanism",
        fontsize=8,
        color="#475569",
        va="center",
    )
    ax.add_patch(
        plt.Rectangle((1.92, legend_y - 0.08), 0.14, 0.16, color="#f59e0b", zorder=5)
    )
    ax.text(
        3.45,
        legend_y,
        "Axes / outcomes",
        fontsize=8,
        color="#475569",
        va="center",
    )
    ax.add_patch(
        plt.Rectangle((3.22, legend_y - 0.08), 0.14, 0.16, color="#10b981", zorder=5)
    )
    ax.text(
        5.05,
        legend_y,
        "Target genes",
        fontsize=8,
        color="#475569",
        va="center",
    )
    ax.add_patch(
        plt.Rectangle((4.82, legend_y - 0.08), 0.14, 0.16, color="#a855f7", zorder=5)
    )
    ax.text(
        6.55,
        legend_y,
        "Solid: qualitative path",
        fontsize=8,
        color="#475569",
        va="center",
        fontweight="600",
    )
    ax.plot(
        [8.05, 8.85],
        [legend_y, legend_y],
        color="#64748b",
        linewidth=1.5,
        linestyle=(0, (4, 3)),
    )
    ax.text(
        9.0,
        legend_y,
        "eQTL path",
        fontsize=8,
        color="#475569",
        va="center",
    )

    ax.set_title(
        "Longevity resilience, LA-SNPs, and neurodegenerative pathophysiology",
        fontsize=12.5,
        fontweight="600",
        color="#0f172a",
        pad=14,
    )

    fig.tight_layout(pad=0.4)
    fig.savefig(out_path, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    return out_path


app = typer.Typer(add_completion=False, no_args_is_help=True)


@app.command()
def main(
    out: Path = typer.Option(
        DEFAULT_OUT,
        "--out",
        "-o",
        help="Output PNG path.",
    ),
    dpi: int = typer.Option(200, "--dpi", min=72, max=600, help="Raster resolution."),
) -> None:
    path = render_longevity_network_diagram(out, dpi=dpi)
    typer.echo(f"Wrote {path}")


if __name__ == "__main__":
    app()
