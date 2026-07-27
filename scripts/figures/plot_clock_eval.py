"""Publication-style figure for methylation clock external validation (GSE87571).

Panel A: predicted vs chronological age scatter.
Panel B: top CpG sites by absolute ElasticNet weight.

If EVAL_CSV is missing, predictions are recomputed from MODEL_PATH + TEST_DATA_PATH
using the same feature alignment logic as ``rogen_aging.clock.evaluate``.
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import typer
from scipy.stats import linregress, pearsonr

from rogen_aging.clock.data import load_wide_table
from rogen_aging.clock.evaluate import build_feature_matrix, load_model

# ---------------------------------------------------------------------------
# Configurable constants
# ---------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parent.parent.parent

# Optional per-sample table from ``rogen-clock evaluate`` (not written by default).
# Expected columns: chronological_age, predicted_age; sample_id optional.
EVAL_CSV: Path | None = REPO_ROOT / "figures" / "validation_gse87571" / "per_sample_predictions.csv"

# Used when EVAL_CSV is absent or does not exist.
MODEL_PATH = REPO_ROOT / "analysis" / "gse40279_elasticnet_clock.pkl"
TEST_DATA_PATH = REPO_ROOT / "data" / "gse87571.parquet"

OUTPUT_DIR = REPO_ROOT / "figures" / "validation_gse87571"
FIG_BASENAME = "clock_eval_gse87571"

TOP_N_CPgs = 25
FIGURE_DPI = 300
FONT_SIZE = 11
POSITIVE_COLOR = "#2166ac"
NEGATIVE_COLOR = "#b2182b"
SCATTER_COLOR = "#404040"

app = typer.Typer(add_completion=False, help=__doc__)


def _feature_names(model: Any) -> list[str] | None:
    """Return ``feature_names_in_`` from a Pipeline or bare estimator."""
    if hasattr(model, "feature_names_in_"):
        return [str(x) for x in model.feature_names_in_]
    if hasattr(model, "named_steps"):
        for step in reversed(list(model.named_steps.values())):
            if hasattr(step, "feature_names_in_"):
                return [str(x) for x in step.feature_names_in_]
    return None


def predict_ages(model: Any, x: pd.DataFrame) -> np.ndarray:
    """Predict ages via ``model.predict`` (Pipeline or bare ElasticNet).

    Args:
        model: Fitted sklearn Pipeline or ElasticNet.
        x: Feature matrix aligned to training CpG columns.

    Returns:
        1-D array of predicted ages.
    """
    return np.asarray(model.predict(x), dtype=float)


def load_or_compute_eval_table(
    *,
    eval_csv: Path | None,
    model_path: Path,
    test_data_path: Path,
) -> pd.DataFrame:
    """Load per-sample ages from CSV, or compute them from model + test cohort.

    Args:
        eval_csv: Optional per-sample predictions CSV; used when the file exists.
        model_path: Serialized clock model (used when recomputing).
        test_data_path: Wide methylation table with ``chronological_age``.

    Returns:
        DataFrame with ``chronological_age``, ``predicted_age``, and optional
        ``sample_id``.

    Raises:
        ValueError: If ``eval_csv`` exists but lacks required columns.
        FileNotFoundError: If neither a usable eval CSV nor test data is available.
    """
    if eval_csv is not None and eval_csv.is_file():
        df = pd.read_csv(eval_csv)
        required = {"chronological_age", "predicted_age"}
        missing = required - set(df.columns)
        if missing:
            raise ValueError(f"EVAL_CSV missing columns: {sorted(missing)}")
        return df

    if not test_data_path.is_file():
        raise FileNotFoundError(
            f"No EVAL_CSV at {eval_csv} and test data not found at {test_data_path}. "
            "Run training/evaluation first or set EVAL_CSV to a predictions table."
        )

    model = load_model(model_path)
    wide = load_wide_table(test_data_path)
    y = pd.to_numeric(wide["chronological_age"], errors="coerce")
    valid = y.notna()
    if not bool(valid.all()):
        warnings.warn(f"Dropping {(~valid).sum()} rows with invalid chronological_age.", stacklevel=2)
    wide = wide.loc[valid].copy()
    y = y.loc[valid]

    x, _imputed = build_feature_matrix(wide, model)
    y_pred = predict_ages(model, x)

    out = pd.DataFrame(
        {
            "chronological_age": y.to_numpy(dtype=float),
            "predicted_age": y_pred,
        },
        index=wide.index,
    )
    if out.index.name is not None or not isinstance(out.index, pd.RangeIndex):
        out.insert(0, "sample_id", out.index.astype(str))
    return out.reset_index(drop=True)


def extract_cpg_weights(model_path: Path) -> pd.Series:
    """Return probe ID -> ElasticNet coefficient for all model features.

    Args:
        model_path: Path to a serialized Pipeline (``elasticnet`` step) or bare
            ElasticNet.

    Returns:
        Series of coefficients indexed by CpG probe IDs.

    Raises:
        ValueError: If coefficients or feature names cannot be recovered.
    """
    model = load_model(model_path)
    if hasattr(model, "named_steps"):
        named_steps = model.named_steps
        if "elasticnet" not in named_steps:
            raise ValueError("Pipeline model missing 'elasticnet' step.")
        enet = named_steps["elasticnet"]
    else:
        enet = model
    if not hasattr(enet, "coef_"):
        raise ValueError("Model has no coef_; expected ElasticNet or Pipeline with elasticnet.")
    coef = np.ravel(enet.coef_)
    names = _feature_names(model)
    if names is None:
        raise ValueError("Model has no feature_names_in_; cannot label CpG coefficients.")
    return pd.Series(coef, index=names, name="coefficient")


def plot_predicted_vs_chronological(ax: plt.Axes, eval_df: pd.DataFrame) -> tuple[float, float, int]:
    """Scatter with identity and regression lines; return MAE, r, n.

    Args:
        ax: Matplotlib axes for the scatter panel.
        eval_df: Table with ``chronological_age`` and ``predicted_age``.

    Returns:
        Tuple of ``(mae, pearson_r, n_samples)``.
    """
    x = eval_df["chronological_age"].to_numpy(dtype=float)
    y = eval_df["predicted_age"].to_numpy(dtype=float)
    n = int(len(x))

    mae = float(np.mean(np.abs(y - x)))
    r_value, _ = pearsonr(x, y)

    ax.scatter(x, y, s=28, alpha=0.75, color=SCATTER_COLOR, edgecolors="white", linewidths=0.4)

    lo = float(min(x.min(), y.min()))
    hi = float(max(x.max(), y.max()))
    pad = 0.05 * (hi - lo) if hi > lo else 1.0
    lim_lo, lim_hi = lo - pad, hi + pad
    ax.plot([lim_lo, lim_hi], [lim_lo, lim_hi], linestyle="--", color="0.45", linewidth=1.2, label="y = x")

    slope, intercept, _, _, _ = linregress(x, y)
    reg_x = np.array([lim_lo, lim_hi])
    ax.plot(reg_x, slope * reg_x + intercept, color="#d95f02", linewidth=1.4, label="Linear fit")

    ax.set_xlim(lim_lo, lim_hi)
    ax.set_ylim(lim_lo, lim_hi)
    ax.set_xlabel("Chronological age (years)")
    ax.set_ylabel("Predicted age (years)")
    ax.set_title("Methylation clock — external validation (GSE87571)")
    ax.set_aspect("equal", adjustable="box")
    ax.legend(loc="lower right", frameon=True, fontsize=FONT_SIZE - 1)

    ax.text(
        0.03,
        0.97,
        f"MAE = {mae:.2f} yr\nPearson r = {r_value:.3f}\nn = {n}",
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=FONT_SIZE,
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "edgecolor": "0.8", "alpha": 0.95},
    )
    return mae, float(r_value), n


def plot_top_cpgs(ax: plt.Axes, weights: pd.Series, top_n: int) -> None:
    """Horizontal bar chart of the top |weight| CpG probes.

    Args:
        ax: Matplotlib axes for the bar panel.
        weights: ElasticNet coefficients indexed by CpG ID.
        top_n: Number of largest-|coefficient| sites to show.
    """
    top = weights.reindex(weights.abs().sort_values(ascending=False).head(top_n).index)
    top = top.sort_values()

    colors = [POSITIVE_COLOR if v >= 0 else NEGATIVE_COLOR for v in top.to_numpy()]
    y_pos = np.arange(len(top))
    ax.barh(y_pos, top.to_numpy(), color=colors, edgecolor="0.2", linewidth=0.4)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(top.index, fontsize=FONT_SIZE - 1)
    ax.axvline(0.0, color="0.3", linewidth=0.8)
    ax.set_xlabel("ElasticNet coefficient")
    ax.set_title("Top CpG sites by model weight")


@app.command()
def main(
    model_path: Path = typer.Option(
        MODEL_PATH,
        "--model-path",
        help="Serialized clock model (.pkl / .joblib)",
        path_type=Path,
    ),
    test_data: Path = typer.Option(
        TEST_DATA_PATH,
        "--test-data",
        help="Wide test table used when eval CSV is missing",
        path_type=Path,
    ),
    eval_csv: Path | None = typer.Option(
        EVAL_CSV,
        "--eval-csv",
        help="Optional per-sample predictions CSV (chronological_age, predicted_age)",
        path_type=Path,
    ),
    output_dir: Path = typer.Option(
        OUTPUT_DIR,
        "--output-dir",
        help="Directory for PNG and PDF outputs",
        path_type=Path,
    ),
    basename: str = typer.Option(
        FIG_BASENAME,
        "--basename",
        help="Output filename stem (without extension)",
    ),
    top_n: int = typer.Option(
        TOP_N_CPgs,
        "--top-n",
        help="Number of top |weight| CpG sites for panel B",
    ),
) -> None:
    """Write the clock validation figure (PNG + PDF)."""
    plt.rcParams.update({"font.size": FONT_SIZE})

    eval_df = load_or_compute_eval_table(
        eval_csv=eval_csv,
        model_path=model_path,
        test_data_path=test_data,
    )

    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.8), constrained_layout=True)
    mae, r_value, n = plot_predicted_vs_chronological(axes[0], eval_df)

    weights = extract_cpg_weights(model_path)
    plot_top_cpgs(axes[1], weights, top_n)

    output_dir.mkdir(parents=True, exist_ok=True)
    png_path = output_dir / f"{basename}.png"
    pdf_path = output_dir / f"{basename}.pdf"
    fig.savefig(png_path, dpi=FIGURE_DPI, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)

    print(f"MAE: {mae:.4f} years")
    print(f"Pearson r: {r_value:.4f}")
    print(f"n samples: {n}")
    print(f"PNG: {png_path}")
    print(f"PDF: {pdf_path}")


if __name__ == "__main__":
    app()
