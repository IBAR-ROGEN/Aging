"""Publication-style figure for methylation clock external validation (GSE87571).

Panel A: predicted vs chronological age scatter.
Panel B: top CpG sites by absolute ElasticNet weight.

If EVAL_CSV is missing, predictions are recomputed from MODEL_PATH + TEST_DATA_PATH
using the same feature alignment logic as ``rogen_aging.clock.evaluate``.
"""

from __future__ import annotations

import pickle
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import linregress, pearsonr

# ---------------------------------------------------------------------------
# Configurable constants
# ---------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parent

# Optional per-sample table from ``rogen-clock evaluate`` (not written by default).
# Expected columns: chronological_age, predicted_age; sample_id optional.
EVAL_CSV: Path | None = REPO_ROOT / "analysis" / "validation_gse87571" / "per_sample_predictions.csv"

# Used when EVAL_CSV is absent or does not exist.
MODEL_PATH = REPO_ROOT / "analysis" / "gse40279_elasticnet_clock.pkl"
TEST_DATA_PATH = REPO_ROOT / "data" / "gse87571.parquet"

OUTPUT_DIR = REPO_ROOT / "analysis" / "validation_gse87571" / "figures"
FIG_BASENAME = "clock_eval_gse87571"

TOP_N_CPgs = 25
FIGURE_DPI = 300
FONT_SIZE = 11
POSITIVE_COLOR = "#2166ac"
NEGATIVE_COLOR = "#b2182b"
SCATTER_COLOR = "#404040"


def load_model(model_path: Path) -> object:
    """Load a sklearn Pipeline saved with ``joblib.dump`` (``train_clock`` default)."""
    if not model_path.is_file():
        raise FileNotFoundError(f"Model file not found: {model_path}")
    try:
        import joblib

        return joblib.load(model_path)
    except Exception:
        with model_path.open("rb") as handle:
            return pickle.load(handle)


def load_wide_table(path: Path) -> pd.DataFrame:
    """Load a wide methylation table (Parquet, CSV, or TSV)."""
    suffix = path.suffix.lower()
    if suffix == ".parquet":
        return pd.read_parquet(path)
    if suffix in {".csv", ".tsv"}:
        sep = "\t" if suffix == ".tsv" else ","
        return pd.read_csv(path, sep=sep)
    raise ValueError(f"Unsupported test data format: {path}")


def _cg_columns(df: pd.DataFrame) -> list[str]:
    return [c for c in df.columns if str(c).startswith("cg")]


def _feature_names(model: object) -> list[str] | None:
    if hasattr(model, "feature_names_in_"):
        return [str(x) for x in model.feature_names_in_]
    if hasattr(model, "named_steps"):
        for step in reversed(list(model.named_steps.values())):
            if hasattr(step, "feature_names_in_"):
                return [str(x) for x in step.feature_names_in_]
    return None


def build_feature_matrix(df: pd.DataFrame, model: object) -> pd.DataFrame:
    """Align test CpGs to training features; mean-impute missing model sites."""
    if "chronological_age" not in df.columns:
        raise ValueError("Test data must include a 'chronological_age' column.")

    cg_cols = _cg_columns(df)
    if not cg_cols:
        raise ValueError("No feature columns starting with 'cg' were found in the test data.")

    expected = _feature_names(model)
    if expected is None:
        raise ValueError("Model has no feature_names_in_; cannot align CpG columns.")

    present_cg = df.reindex(columns=cg_cols).apply(pd.to_numeric, errors="coerce")
    flat_mean = float(np.nanmean(present_cg.to_numpy(dtype=float))) if cg_cols else 0.5
    if not np.isfinite(flat_mean):
        flat_mean = 0.5

    x = pd.DataFrame(index=df.index)
    for name in expected:
        if name in df.columns:
            col = pd.to_numeric(df[name], errors="coerce")
            fill = float(np.nanmean(col.to_numpy())) if col.notna().any() else flat_mean
            if not np.isfinite(fill):
                fill = flat_mean
            x[name] = col.fillna(fill)
        else:
            warnings.warn(
                f"CpG '{name}' expected by the model is absent from test data; "
                f"filling with global mean ({flat_mean:.6g}).",
                stacklevel=2,
            )
            x[name] = flat_mean
    return x


def predict_ages(model: object, x: pd.DataFrame) -> np.ndarray:
    """Run imputer + ElasticNet prediction without importing sklearn."""
    imputer = model.named_steps["imputer"]
    enet = model.named_steps["elasticnet"]
    x_imputed = x.to_numpy(dtype=float)
    nan_mask = np.isnan(x_imputed)
    if nan_mask.any():
        x_imputed = np.where(nan_mask, imputer.statistics_, x_imputed)
    return np.asarray(enet.predict(x_imputed), dtype=float)


def load_or_compute_eval_table() -> pd.DataFrame:
    """Load per-sample ages from CSV, or compute them from model + test cohort."""
    if EVAL_CSV is not None and EVAL_CSV.is_file():
        df = pd.read_csv(EVAL_CSV)
        required = {"chronological_age", "predicted_age"}
        missing = required - set(df.columns)
        if missing:
            raise ValueError(f"EVAL_CSV missing columns: {sorted(missing)}")
        return df

    if not TEST_DATA_PATH.is_file():
        raise FileNotFoundError(
            f"No EVAL_CSV at {EVAL_CSV} and test data not found at {TEST_DATA_PATH}. "
            "Run training/evaluation first or set EVAL_CSV to a predictions table."
        )

    model = load_model(MODEL_PATH)
    wide = load_wide_table(TEST_DATA_PATH)
    y = pd.to_numeric(wide["chronological_age"], errors="coerce")
    valid = y.notna()
    if not bool(valid.all()):
        warnings.warn(f"Dropping {(~valid).sum()} rows with invalid chronological_age.", stacklevel=2)
    wide = wide.loc[valid].copy()
    y = y.loc[valid]

    x = build_feature_matrix(wide, model)
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
    """Return probe ID -> ElasticNet coefficient for all model features."""
    model = load_model(model_path)
    enet = model.named_steps["elasticnet"]
    coef = np.ravel(enet.coef_)
    names = _feature_names(model)
    if names is None:
        raise ValueError("Model has no feature_names_in_; cannot label CpG coefficients.")
    return pd.Series(coef, index=names, name="coefficient")


def plot_predicted_vs_chronological(ax: plt.Axes, eval_df: pd.DataFrame) -> tuple[float, float, int]:
    """Scatter with identity and regression lines; return MAE, r, n."""
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
    """Horizontal bar chart of the top |weight| CpG probes."""
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


def main() -> None:
    plt.rcParams.update({"font.size": FONT_SIZE})

    eval_df = load_or_compute_eval_table()
    mae, r_value, n = 0.0, 0.0, 0

    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.8), constrained_layout=True)
    mae, r_value, n = plot_predicted_vs_chronological(axes[0], eval_df)

    weights = extract_cpg_weights(MODEL_PATH)
    plot_top_cpgs(axes[1], weights, TOP_N_CPgs)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    png_path = OUTPUT_DIR / f"{FIG_BASENAME}.png"
    pdf_path = OUTPUT_DIR / f"{FIG_BASENAME}.pdf"
    fig.savefig(png_path, dpi=FIGURE_DPI, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)

    print(f"MAE: {mae:.4f} years")
    print(f"Pearson r: {r_value:.4f}")
    print(f"n samples: {n}")
    print(f"PNG: {png_path}")
    print(f"PDF: {pdf_path}")


if __name__ == "__main__":
    main()
