#!/usr/bin/env python3
"""Train a mock Elastic Net epigenetic aging clock for a Romanian cohort.

Loads a synthetic methylation matrix and sample metadata (chronological age),
fits ElasticNetCV with 5-fold cross-validation over alpha and l1_ratio,
reports MAE and Pearson r on a held-out test set, and writes a scatter plot.

Romanian mock cohort I/O lives in :mod:`rogen_aging.clock.data`; this script keeps the
StandardScaler + ElasticNetCV training path used for the mock Romanian demo.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr
from sklearn.linear_model import ElasticNetCV
from sklearn.metrics import mean_absolute_error
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from typer import Option, Typer

from rogen_aging.clock.data import load_romanian_cohort
from rogen_aging.config import cfg_path, find_repo_root, get_config
from rogen_aging.config.cli import config_option, load_cli_config

app = Typer(add_completion=False, no_args_is_help=True)


@app.command()
def main(
    config: Path | None = config_option(),
    data_dir: Path | None = Option(
        None,
        help="Directory containing methylation_matrix.csv and metadata.csv. Default: from config.",
    ),
    test_size: float | None = Option(
        None, min=0.05, max=0.5, help="Fraction held out for evaluation. Default: from config."
    ),
    random_state: int | None = Option(None, help="Random seed. Default: from config."),
    output_plot: Path | None = Option(
        None,
        help="Where to save the scatter plot (PNG). Default: from config.",
    ),
    regenerate_mock: bool = Option(
        False,
        "--regenerate-mock",
        help="Overwrite mock CSVs in data_dir with freshly simulated values.",
    ),
) -> None:
    cfg = load_cli_config(config)
    repo_root = find_repo_root()
    resolved_data_dir = data_dir or cfg_path(cfg, "paths", "methylation", "romanian_mock_dir")
    # Historical bug used scripts/ as root; keep resolving against real repo root.
    if not resolved_data_dir.is_absolute():
        resolved_data_dir = repo_root / resolved_data_dir
    resolved_test_size = float(cfg.clock.test_size if test_size is None else test_size)
    resolved_random_state = int(cfg.clock.random_state if random_state is None else random_state)
    resolved_plot = output_plot or cfg_path(cfg, "paths", "clock_eval", "romanian_plot")
    romanian = cfg.clock.romanian
    l1_ratios = [float(x) for x in romanian.l1_ratio]

    x, y, _sample_ids = load_romanian_cohort(
        resolved_data_dir,
        regenerate_mock=regenerate_mock,
        random_state=resolved_random_state,
    )
    x_train, x_test, y_train, y_test = train_test_split(
        x,
        y,
        test_size=resolved_test_size,
        random_state=resolved_random_state,
    )

    model = Pipeline(
        steps=[
            ("scale", StandardScaler()),
            (
                "enet",
                ElasticNetCV(
                    l1_ratio=l1_ratios,
                    cv=int(romanian.cv),
                    random_state=resolved_random_state,
                    max_iter=int(romanian.max_iter),
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

    resolved_plot.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.0, 6.5))
    ax.scatter(y_test, y_pred, alpha=0.75, edgecolors="black", linewidths=0.3, s=42)
    ax.plot(x_line, y_line, color="crimson", linewidth=2.0, label="Line of best fit")
    ax.plot(
        [y_test.min(), y_test.max()],
        [y_test.min(), y_test.max()],
        "k--",
        linewidth=1.0,
        label="Identity (y = x)",
    )
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
    fig.savefig(resolved_plot, dpi=int(get_config().clock.figure_dpi))
    plt.close(fig)
    print(f"  Saved plot: {resolved_plot}")


if __name__ == "__main__":
    app()
