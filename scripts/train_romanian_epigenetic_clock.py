#!/usr/bin/env python3
"""Train a mock Elastic Net epigenetic aging clock for a Romanian cohort.

Loads a synthetic methylation matrix and sample metadata (chronological age),
fits ElasticNetCV with 5-fold cross-validation over alpha and l1_ratio,
reports MAE and Pearson r on a held-out test set, and writes a scatter plot.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import polars as pl
from typer import Option, Typer
from scipy.stats import pearsonr
from sklearn.linear_model import ElasticNetCV
from sklearn.metrics import mean_absolute_error
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

app = Typer(add_completion=False, no_args_is_help=True)

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_DATA_DIR = REPO_ROOT / "data" / "mock_romanian_cohort"


def write_mock_romanian_cohort(data_dir: Path, n_samples: int = 120, n_cpgs: int = 80) -> None:
    """Create deterministic mock methylation + metadata for development."""
    data_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(42)
    sample_ids = [f"ROM{i:04d}" for i in range(1, n_samples + 1)]
    ages = rng.normal(52.0, 14.0, size=n_samples).clip(25.0, 92.0)

    # A subset of CpGs carry a synthetic age signal; the rest is noise (realistic sparsity).
    beta = rng.uniform(0.15, 0.85, size=(n_samples, n_cpgs)).astype(np.float64)
    signal_idx = rng.choice(n_cpgs, size=max(12, n_cpgs // 5), replace=False)
    for j in signal_idx:
        beta[:, j] += 0.018 * (ages - ages.mean())
    beta = np.clip(beta, 0.02, 0.98)

    cpg_names = [f"cg_mock_{k:05d}" for k in range(1, n_cpgs + 1)]
    meth = pl.DataFrame(beta, schema=cpg_names).with_columns(pl.Series("sample_id", sample_ids))
    meta = pl.DataFrame({"sample_id": sample_ids, "chronological_age": ages.astype(np.float64)})
    meth.write_csv(data_dir / "methylation_matrix.csv")
    meta.write_csv(data_dir / "metadata.csv")


def load_cohort(data_dir: Path, *, regenerate_mock: bool = False) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """Return feature matrix X, target y (chronological age), and ordered sample IDs."""
    meth_path = data_dir / "methylation_matrix.csv"
    meta_path = data_dir / "metadata.csv"
    if regenerate_mock or not meth_path.is_file() or not meta_path.is_file():
        write_mock_romanian_cohort(data_dir)

    meth = pl.read_csv(meth_path)
    meta = pl.read_csv(meta_path)
    if "sample_id" not in meth.columns or "sample_id" not in meta.columns:
        raise ValueError("Both tables must include a 'sample_id' column.")
    if "chronological_age" not in meta.columns:
        raise ValueError("Metadata must include 'chronological_age'.")

    joined = meta.join(meth, on="sample_id", how="inner").sort("sample_id")
    feature_cols = [c for c in joined.columns if c not in ("sample_id", "chronological_age")]
    if not feature_cols:
        raise ValueError("No CpG / feature columns found after join.")

    x_df = joined.select(feature_cols)
    x = x_df.to_numpy()
    y = joined["chronological_age"].to_numpy()
    sample_ids = joined["sample_id"].to_list()
    return x, y, sample_ids


@app.command()
def main(
    data_dir: Path = Option(
        DEFAULT_DATA_DIR,
        help="Directory containing methylation_matrix.csv and metadata.csv.",
    ),
    test_size: float = Option(0.2, min=0.05, max=0.5, help="Fraction held out for evaluation."),
    random_state: int = Option(42, help="Random seed for the train/test split."),
    output_plot: Path = Option(
        REPO_ROOT / "figures" / "romanian_mock_epigenetic_clock_scatter.png",
        help="Where to save the scatter plot (PNG).",
    ),
    regenerate_mock: bool = Option(
        False,
        "--regenerate-mock",
        help="Overwrite mock CSVs in data_dir with freshly simulated values.",
    ),
) -> None:
    x, y, _sample_ids = load_cohort(data_dir, regenerate_mock=regenerate_mock)
    x_train, x_test, y_train, y_test = train_test_split(
        x,
        y,
        test_size=test_size,
        random_state=random_state,
    )

    l1_ratios = [0.1, 0.3, 0.5, 0.7, 0.9, 0.95, 0.99, 1.0]
    model = Pipeline(
        steps=[
            ("scale", StandardScaler()),
            (
                "enet",
                ElasticNetCV(
                    l1_ratio=l1_ratios,
                    cv=5,
                    random_state=random_state,
                    max_iter=20000,
                    n_jobs=-1,
                ),
            ),
        ]
    )
    model.fit(x_train, y_train)
    y_pred = model.predict(x_test)

    mae = float(mean_absolute_error(y_test, y_pred))
    r_value, _p = pearsonr(y_test, y_pred)

    enet: ElasticNetCV = model.named_steps["enet"]
    print("Romanian cohort mock epigenetic clock (Elastic Net)")
    print(f"  Training samples: {len(y_train)}, test samples: {len(y_test)}")
    print(f"  Selected l1_ratio (L1 / L1+L2 mix): {enet.l1_ratio_:.4f}")
    print(f"  Selected alpha (regularization strength): {enet.alpha_:.6f}")
    print(f"  Test MAE (years): {mae:.3f}")
    print(f"  Test Pearson r (chronological vs predicted): {r_value:.4f}")

    slope, intercept = np.polyfit(y_test, y_pred, 1)
    x_line = np.linspace(float(np.min(y_test)), float(np.max(y_test)), 200)
    y_line = slope * x_line + intercept

    output_plot.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.0, 6.5))
    ax.scatter(y_test, y_pred, alpha=0.75, edgecolors="black", linewidths=0.3, s=42)
    ax.plot(x_line, y_line, color="crimson", linewidth=2.0, label="Line of best fit")
    ax.plot([y_test.min(), y_test.max()], [y_test.min(), y_test.max()], "k--", linewidth=1.0, label="Identity (y = x)")
    ax.set_xlabel("Chronological age (years)")
    ax.set_ylabel("Predicted epigenetic age (years)")
    ax.set_title("Mock Romanian cohort: chronological vs predicted age")
    ax.grid(True, linestyle=":", linewidth=0.6, alpha=0.7)
    ax.legend(loc="upper left", frameon=True)
    text = f"MAE = {mae:.2f} y\nr = {r_value:.3f}"
    ax.text(
        0.98,
        0.02,
        text,
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=11,
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.88, "edgecolor": "#333333"},
    )
    fig.tight_layout()
    fig.savefig(output_plot, dpi=160)
    plt.close(fig)
    print(f"  Saved plot: {output_plot}")


if __name__ == "__main__":
    app()
