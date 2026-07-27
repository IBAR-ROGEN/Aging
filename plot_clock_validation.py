#!/usr/bin/env python3
"""Predicted vs chronological age scatter for methylation clock validation.

Reads a per-sample CSV (chronological_age, predicted_age) and saves PNG + PDF.
Default cohort label targets the GSE40279 proxy validation set.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_absolute_error, r2_score

# ---------------------------------------------------------------------------
# Editable constants
# ---------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parent

INPUT_CSV = REPO_ROOT / "figures" / "validation_gse40279" / "per_sample_predictions.csv"
OUTPUT_STEM = REPO_ROOT / "figures" / "validation_gse40279" / "clock_validation_gse40279"

DATASET_LABEL = "GSE40279"

SCATTER_COLOR = "#404040"
IDENTITY_LINE_COLOR = "#666666"
REGRESSION_LINE_COLOR = "#d95f02"

DRAW_REGRESSION_LINE = True

FIGURE_DPI = 300
SCATTER_SIZE = 28
SCATTER_ALPHA = 0.75

FONT_SIZES = {
    "base": 11,
    "axis_label": 12,
    "tick": 10,
    "title": 12,
    "metrics": 11,
    "legend": 10,
}

REQUIRED_COLUMNS = ("chronological_age", "predicted_age")


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
    work["chronological_age"] = pd.to_numeric(work["chronological_age"], errors="coerce")
    work["predicted_age"] = pd.to_numeric(work["predicted_age"], errors="coerce")
    valid = work["chronological_age"].notna() & work["predicted_age"].notna()
    dropped = int((~valid).sum())
    if dropped:
        print(
            f"Warning: dropping {dropped} row(s) with non-numeric or missing age values.",
            file=sys.stderr,
        )
    work = work.loc[valid]
    if work.empty:
        print("Error: no samples with valid chronological_age and predicted_age.", file=sys.stderr)
        sys.exit(1)
    return work


def _axis_limits(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    lo = float(min(x.min(), y.min()))
    hi = float(max(x.max(), y.max()))
    pad = 0.05 * (hi - lo) if hi > lo else 1.0
    return lo - pad, hi + pad


def plot_validation(
    df: pd.DataFrame,
    output_stem: Path,
    *,
    dataset_label: str,
    draw_regression_line: bool,
) -> None:
    x = df["chronological_age"].to_numpy(dtype=float)
    y = df["predicted_age"].to_numpy(dtype=float)

    mae = float(mean_absolute_error(x, y))
    r2 = float(r2_score(x, y))
    n = int(len(x))

    lim_lo, lim_hi = _axis_limits(x, y)

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
        x,
        y,
        s=SCATTER_SIZE,
        alpha=SCATTER_ALPHA,
        color=SCATTER_COLOR,
        edgecolors="white",
        linewidths=0.4,
        zorder=2,
    )

    ax.plot(
        [lim_lo, lim_hi],
        [lim_lo, lim_hi],
        linestyle="--",
        color=IDENTITY_LINE_COLOR,
        linewidth=1.2,
        label="y = x",
        zorder=1,
    )

    if draw_regression_line:
        reg = LinearRegression()
        reg.fit(x.reshape(-1, 1), y)
        reg_x = np.array([lim_lo, lim_hi])
        reg_y = reg.predict(reg_x.reshape(-1, 1))
        ax.plot(
            reg_x,
            reg_y,
            color=REGRESSION_LINE_COLOR,
            linewidth=1.2,
            alpha=0.85,
            label="Linear fit",
            zorder=1,
        )

    ax.set_xlim(lim_lo, lim_hi)
    ax.set_ylim(lim_lo, lim_hi)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Chronological age (years)", fontsize=FONT_SIZES["axis_label"])
    ax.set_ylabel("Predicted age (years)", fontsize=FONT_SIZES["axis_label"])
    ax.set_title(dataset_label, fontsize=FONT_SIZES["title"])
    ax.tick_params(labelsize=FONT_SIZES["tick"])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    if draw_regression_line:
        ax.legend(loc="lower right", frameon=True, fontsize=FONT_SIZES["legend"])

    ax.text(
        0.03,
        0.97,
        f"MAE = {mae:.2f} yr\nR² = {r2:.3f}\nn = {n}",
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=FONT_SIZES["metrics"],
        bbox={
            "boxstyle": "round,pad=0.35",
            "facecolor": "white",
            "edgecolor": "0.8",
            "alpha": 0.95,
        },
    )

    output_stem = Path(output_stem)
    output_stem.parent.mkdir(parents=True, exist_ok=True)
    png_path = output_stem.with_suffix(".png")
    pdf_path = output_stem.with_suffix(".pdf")
    fig.savefig(png_path, dpi=FIGURE_DPI, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)

    print(f"MAE: {mae:.4f} years")
    print(f"R²: {r2:.4f}")
    print(f"n samples: {n}")
    print(f"Saved: {png_path}")
    print(f"Saved: {pdf_path}")


def main() -> None:
    df = _load_data(Path(INPUT_CSV))
    plot_validation(
        df,
        Path(OUTPUT_STEM),
        dataset_label=DATASET_LABEL,
        draw_regression_line=DRAW_REGRESSION_LINE,
    )


if __name__ == "__main__":
    main()
