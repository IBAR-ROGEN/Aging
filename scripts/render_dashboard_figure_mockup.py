"""Render a static PNG of the ROGEN EDA dashboard figure mockup (matplotlib).

Mirrors the layout of ``components/DashboardFigureMockup.tsx`` for manuscripts when
Node/Playwright is unavailable. For pixel-perfect React output, use
``components/dashboard-figure-render`` (``npm run capture``).

Example:

    uv run python scripts/render_dashboard_figure_mockup.py
"""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import typer
from matplotlib import patches
from matplotlib.patches import FancyBboxPatch

app = typer.Typer(add_completion=False)

REPO_ROOT = Path(__file__).resolve().parent.parent


def _mock_scatter_points() -> tuple[np.ndarray, np.ndarray]:
    """Same LCG recipe as ``DashboardFigureMockup.tsx`` MOCK_SCATTER_DATA."""
    xs: list[float] = []
    ys: list[float] = []
    s = 42_069
    for i in range(50):
        s = (s * 1_103_515_245 + 12_345) & 0xFFFF_FFFF
        u = s / 4_294_967_295
        t = i / 49
        x = 40 + 60 * t + (u - 0.5) * 3.5
        noise = math.sin(i * 1.7 + u * 6) * 3.2 + (u - 0.48) * 4
        y = x + noise
        y = max(38.0, min(102.0, y))
        xs.append(round(x, 2))
        ys.append(round(y, 2))
    return np.array(xs), np.array(ys)


def _draw_round_rect(
    fig: plt.Figure,
    xy: tuple[float, float],
    w: float,
    h: float,
    *,
    facecolor: str,
    edgecolor: str | None = None,
    linewidth: float = 0,
    radius: float = 0.008,
) -> None:
    box = FancyBboxPatch(
        xy,
        w,
        h,
        boxstyle=f"round,pad=0,rounding_size={radius}",
        transform=fig.transFigure,
        facecolor=facecolor,
        edgecolor=edgecolor,
        linewidth=linewidth,
        clip_on=False,
    )
    fig.add_artist(box)


@app.command()
def main(
    output: Path = typer.Option(
        REPO_ROOT / "analysis" / "dashboard_figure_mockup.png",
        help="Output PNG path.",
        path_type=Path,
    ),
    dpi: int = typer.Option(220, help="Figure resolution."),
) -> None:
    """Write a wide dashboard mockup figure to ``output``."""
    output = output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)

    # Academic palette (slate / sky / white)
    C = {
        "bg": "#f1f5f9",
        "white": "#ffffff",
        "side": "#f8fafc",
        "border": "#cbd5e1",
        "slate4": "#94a3b8",
        "slate5": "#64748b",
        "slate6": "#475569",
        "slate7": "#334155",
        "slate9": "#0f172a",
        "sky2": "#bae6fd",
        "sky6": "#0284c7",
        "sky8": "#075985",
        "emerald": "#10b981",
        "grid": "#e2e8f0",
    }

    fig_w, fig_h = 14.0, 8.2
    fig = plt.figure(figsize=(fig_w, fig_h), dpi=dpi, facecolor=C["bg"])

    # Outer "application window"
    _draw_round_rect(
        fig,
        (0.03, 0.04),
        0.94,
        0.90,
        facecolor=C["white"],
        edgecolor=C["border"],
        linewidth=1.2,
        radius=0.012,
    )

    # Sidebar column (~20%)
    _draw_round_rect(fig, (0.045, 0.055), 0.185, 0.87, facecolor=C["side"], radius=0.006)
    fig.add_artist(
        patches.Rectangle(
            (0.045, 0.055),
            0.185,
            0.87,
            transform=fig.transFigure,
            facecolor="none",
            edgecolor=C["border"],
            linewidth=0.8,
        )
    )

    sx = 0.055
    sy = 0.88
    fig.text(sx, sy, "ROGEN EDA Settings", fontsize=11, fontweight="bold", color=C["slate9"])
    sy -= 0.045
    fig.text(sx, sy, "Merged cohort path", fontsize=7.5, color=C["slate6"])
    sy -= 0.028
    _draw_round_rect(fig, (sx, sy - 0.022), 0.165, 0.028, facecolor=C["white"], edgecolor=C["border"], linewidth=0.6)
    fig.text(sx + 0.008, sy - 0.012, "data/merged_cohort.parquet", fontsize=7, color=C["slate7"], family="monospace")
    sy -= 0.055
    _draw_round_rect(fig, (sx, sy - 0.032), 0.165, 0.036, facecolor="#f8fafc", edgecolor=C["border"], linewidth=0.5)
    fig.text(sx + 0.01, sy - 0.012, "Use Synthetic Mock Cohort", fontsize=8, color=C["slate7"], va="center")
    _draw_round_rect(fig, (sx + 0.12, sy - 0.026), 0.038, 0.022, facecolor=C["emerald"], radius=0.003)
    circ = plt.Circle((sx + 0.148, sy - 0.015), 0.008, transform=fig.transFigure, color=C["white"], zorder=5)
    fig.add_artist(circ)

    sy -= 0.075
    fig.text(sx, sy, "GLOBAL FILTERS", fontsize=7, fontweight="bold", color=C["slate5"])
    sy -= 0.035
    fig.text(sx, sy, "Age range", fontsize=7.5, color=C["slate6"])
    fig.text(sx + 0.12, sy, "40 — 100 years", fontsize=7.5, color=C["slate5"], ha="right")
    sy -= 0.018
    fig.add_artist(
        patches.FancyBboxPatch(
            (sx, sy - 0.008),
            0.165,
            0.012,
            boxstyle="round,pad=0,rounding_size=0.003",
            transform=fig.transFigure,
            facecolor=C["grid"],
            edgecolor="none",
        )
    )
    fig.add_artist(
        patches.FancyBboxPatch(
            (sx + 0.006, sy - 0.008),
            0.135,
            0.012,
            boxstyle="round,pad=0,rounding_size=0.003",
            transform=fig.transFigure,
            facecolor=C["sky2"],
            edgecolor="none",
            alpha=0.85,
        )
    )
    fig.add_artist(
        plt.Circle((sx + 0.02, sy - 0.002), 0.007, transform=fig.transFigure, facecolor=C["white"], edgecolor=C["sky6"], linewidth=1.2)
    )
    fig.add_artist(
        plt.Circle((sx + 0.145, sy - 0.002), 0.007, transform=fig.transFigure, facecolor=C["white"], edgecolor=C["sky6"], linewidth=1.2)
    )

    def chips_block(y_top: float, label: str, chips: tuple[str, ...]) -> None:
        fig.text(sx, y_top, label, fontsize=7.5, color=C["slate6"])
        yb = y_top - 0.038
        _draw_round_rect(fig, (sx, yb), 0.165, 0.036, facecolor=C["white"], edgecolor=C["border"], linewidth=0.6)
        cx = sx + 0.012
        cy = yb + 0.022
        step = 0.044 if len(chips) <= 2 else 0.038
        for ch in chips:
            fig.text(
                cx,
                cy,
                ch,
                fontsize=6,
                color=C["sky8"],
                va="center",
                bbox=dict(boxstyle="round,pad=0.2", facecolor="#e0f2fe", edgecolor="#7dd3fc", linewidth=0.5),
            )
            cx += step

    chips_block(sy - 0.045, "Sex", ("Female", "Male"))
    chips_block(sy - 0.14, "Disease status", ("Control", "Case", "Prodromal"))

    # Main content
    mx = 0.255
    mw = 0.70
    fig.text(
        mx,
        0.905,
        "ROGEN Multi-Omics Exploratory Data Analysis",
        fontsize=15,
        fontweight="bold",
        color=C["slate9"],
    )
    fig.text(
        mx,
        0.865,
        "Integrating whole-genome sequence-derived variants, DNA methylation–based age estimates,",
        fontsize=8.5,
        color=C["slate6"],
    )
    fig.text(
        mx,
        0.845,
        "and structured clinical phenotypes to characterize aging-related biology in the merged cohort.",
        fontsize=8.5,
        color=C["slate6"],
    )

    # Tabs
    ty = 0.798
    tw = mw / 3
    tab_labels = ("Clinical & Phenotypic", "Epigenetic Clock Validation", "LA-SNPs")
    for i, lab in enumerate(tab_labels):
        tx = mx + i * tw
        active = i == 1
        fc = C["white"] if active else "#f8fafc"
        ec = C["sky6"] if active else "none"
        lw = 2.0 if active else 0
        _draw_round_rect(fig, (tx + 0.005, ty - 0.038), tw - 0.02, 0.04, facecolor=fc, edgecolor=ec if active else C["border"], linewidth=lw if active else 0.4)
        fig.text(
            tx + tw / 2,
            ty - 0.018,
            lab,
            ha="center",
            va="center",
            fontsize=8 if active else 7.5,
            fontweight="bold" if active else "normal",
            color=C["sky8"] if active else C["slate5"],
        )

    fig.text(mx, 0.745, "Epigenetic clock validation", fontsize=11, fontweight="bold", color=C["slate9"])
    fig.text(
        mx,
        0.718,
        "Chronological age versus predicted DNAm age with ordinary least squares fit; metrics summarize agreement for the filtered cohort.",
        fontsize=8,
        color=C["slate6"],
    )

    # Metric cards
    card_y = 0.615
    card_w = (mw - 0.02) / 3
    metrics = (
        ("MEAN ABSOLUTE ERROR (MAE)", "3.42 years"),
        ("PEARSON CORRELATION (r)", "0.89"),
        ("COHORT SIZE (N)", "842"),
    )
    for i, (title, val) in enumerate(metrics):
        cx = mx + i * (card_w + 0.008)
        _draw_round_rect(fig, (cx, card_y - 0.065), card_w, 0.07, facecolor=C["white"], edgecolor=C["border"], linewidth=0.7)
        fig.text(cx + 0.015, card_y - 0.018, title, fontsize=6.5, color=C["slate5"], fontweight="bold")
        fig.text(cx + 0.015, card_y - 0.048, val, fontsize=14, color=C["slate9"], fontweight="bold")

    # Scatter axes region
    ax_left, ax_bottom, ax_w, ax_h = mx, 0.09, mw, 0.48
    ax = fig.add_axes([ax_left, ax_bottom, ax_w, ax_h])
    ax.set_facecolor(C["white"])
    for spine in ax.spines.values():
        spine.set_color(C["border"])
        spine.set_linewidth(0.8)
    ax.grid(True, linestyle="--", linewidth=0.6, color=C["grid"], alpha=0.9)
    ax.set_xlim(40, 100)
    ax.set_ylim(40, 100)
    ax.set_xticks(np.arange(40, 101, 10))
    ax.set_yticks(np.arange(40, 101, 10))
    ax.tick_params(labelsize=8, colors=C["slate6"])
    ax.set_xlabel("Chronological Age (years)", fontsize=9, color=C["slate7"], fontweight="medium", labelpad=6)
    ax.set_ylabel("Predicted DNAm Age (years)", fontsize=9, color=C["slate7"], fontweight="medium", labelpad=6)
    ax.set_title("Chronological age vs. predicted DNAm age", fontsize=9.5, color=C["slate7"], pad=8, fontweight="medium")

    ax.plot([40, 100], [40, 100], color=C["slate5"], linewidth=2.0, linestyle=(0, (6, 5)), zorder=1, label="OLS (y = x)")
    xs, ys = _mock_scatter_points()
    ax.scatter(xs, ys, s=38, c=C["sky6"], alpha=0.88, edgecolors="white", linewidths=0.4, zorder=3)

    fig.text(
        mx + mw / 2,
        0.065,
        "Dashed line: identity / OLS reference (y = x). Points: independent samples (mock data for illustration).",
        ha="center",
        fontsize=7,
        color=C["slate5"],
    )

    fig.savefig(output, bbox_inches="tight", facecolor=C["bg"], edgecolor="none", pad_inches=0.15)
    plt.close(fig)
    typer.echo(f"Wrote {output}")


if __name__ == "__main__":
    app()
