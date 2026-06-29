# Clock external-validation figure (`plot_clock_eval.py`)

**Project:** IBAR-ROGEN Aging  
**Activity:** 2.1.10.1 — methylation aging clock (GSE40279 train, GSE87571 validate)  
**Script:** [`plot_clock_eval.py`](../plot_clock_eval.py) (repo root)  
**Related:** [CLOCK_LIBRARY.md](CLOCK_LIBRARY.md) · [GSE40279_CLOCK_TRAINING.md](GSE40279_CLOCK_TRAINING.md) · [ACTIVITIES.md](ACTIVITIES.md#21101--methylation-aging-clock)

## Purpose

`plot_clock_eval.py` produces a **publication-style two-panel figure** from a trained ElasticNet methylation clock and its held-out external validation on GSE87571:

| Panel | Content |
|-------|---------|
| **A** | Predicted vs chronological age scatter (one point per sample), y = x reference, linear fit, MAE and Pearson r computed from the data |
| **B** | Top N CpG probes by absolute ElasticNet coefficient (horizontal bar chart; positive and negative weights in distinct colors) |

Outputs are written as **PNG (300 dpi)** and **PDF** under a configurable output directory.

This complements the built-in figures from **`rogen-clock evaluate`** (`Fig_Clock_Residuals.png`, `Fig_Clock_MAE_by_decade.png`), which focus on residuals and decade-stratified MAE rather than a predicted-vs-chronological scatter with model weights.

## Prerequisites

Complete the train → external cohort → evaluate workflow first:

```bash
uv sync

# 1. Train on GSE40279-style wide table (cg* + chronological_age)
uv run rogen-clock train \
  --input_data data/gse40279_beta_age.parquet \
  --output_model analysis/gse40279_elasticnet_clock.pkl \
  --output_metrics analysis/gse40279_train_metrics.json

# 2. Build GSE87571 evaluation cohort
uv run python -m rogen_aging.clock.external_data \
  --output data/gse87571.parquet \
  --geo-cache-dir data/geo

# 3. Run library evaluation (residual / decade figures + metrics JSON)
uv run rogen-clock evaluate \
  --model_path analysis/gse40279_elasticnet_clock.pkl \
  --test_data data/gse87571.parquet \
  --output_dir analysis/validation_gse87571
```

See [GSE40279_CLOCK_TRAINING.md](GSE40279_CLOCK_TRAINING.md) for input formats, GEO download notes, and `restrict_to_cpgs` options.

## Run the figure script

```bash
uv run python plot_clock_eval.py
```

Stdout reports **MAE**, **Pearson r**, **n samples**, and the paths to the saved PNG and PDF.

Default outputs:

| File | Description |
|------|-------------|
| `analysis/validation_gse87571/figures/clock_eval_gse87571.png` | Raster figure (300 dpi) |
| `analysis/validation_gse87571/figures/clock_eval_gse87571.pdf` | Vector figure |

## Configuration

All paths and styling are **constants at the top of** `plot_clock_eval.py` (no CLI flags). Edit these before running:

| Constant | Default | Role |
|----------|---------|------|
| `EVAL_CSV` | `analysis/validation_gse87571/per_sample_predictions.csv` | Optional per-sample table (see below) |
| `MODEL_PATH` | `analysis/gse40279_elasticnet_clock.pkl` | Fitted `Pipeline` from `rogen-clock train` |
| `TEST_DATA_PATH` | `data/gse87571.parquet` | External validation wide table |
| `OUTPUT_DIR` | `analysis/validation_gse87571/figures` | Directory for PNG/PDF |
| `FIG_BASENAME` | `clock_eval_gse87571` | Output filename stem |
| `TOP_N_CPgs` | `25` | Number of probes in panel B |
| `FIGURE_DPI` | `300` | PNG resolution |
| `FONT_SIZE` | `11` | Base matplotlib font size |
| `POSITIVE_COLOR` / `NEGATIVE_COLOR` | blue / red | Bar colors in panel B |

## Input data: per-sample ages

### Option A — Recompute from model + test table (default)

If `EVAL_CSV` is **missing**, the script loads `MODEL_PATH` and `TEST_DATA_PATH`, aligns CpG columns to the model’s `feature_names_in_` (mean-imputing absent probes, matching `rogen_aging.clock.evaluate`), and predicts ages for panel A.

**Column mapping for panel A:**

| Column | Axis / meaning |
|--------|----------------|
| `chronological_age` | x — true age in years |
| `predicted_age` | y — model prediction in years |

Sample identifiers come from the wide table index when present (e.g. GSM IDs from `load_gse87571`).

### Option B — Optional predictions CSV

`rogen-clock evaluate` does **not** write a per-sample CSV by default; it only saves aggregate `validation_metrics.json` and residual/decade plots. To drive panel A from a file, export a CSV with at least:

```text
chronological_age,predicted_age
58.0,56.3
...
```

Optional column: `sample_id`. Point `EVAL_CSV` at that file; if it exists, the script skips recomputation.

## Model coefficients (panel B)

Panel B reads **ElasticNet weights** from the saved pipeline:

```python
import joblib

model = joblib.load("analysis/gse40279_elasticnet_clock.pkl")
coef = model.named_steps["elasticnet"].coef_
cpg_ids = list(model.feature_names_in_)
```

The script selects the top `TOP_N_CPgs` probes by `abs(coefficient)` and sorts them for the horizontal bar chart.

Training metrics JSON (`analysis/gse40279_train_metrics.json`) lists **`selected_cpgs`** (non-zero probes) but not numeric weights; use the pickled model for panel B.

## Dependencies

Plotting and statistics use **pandas**, **numpy**, **matplotlib**, and **scipy** only. Loading the trained pipeline requires **joblib** (same format as `train_clock`); sklearn must be available in the environment because the saved object is a fitted sklearn `Pipeline`.

## Distinction from other clock figures

| Tool | Data | Output |
|------|------|--------|
| **`plot_clock_eval.py`** | Real model + GSE87571 (or EVAL CSV) | Two-panel scatter + top CpGs (PNG/PDF) |
| **`rogen-clock evaluate`** | Same test table | `validation_metrics.json`, residuals, MAE-by-decade PNGs |
| **`scripts/figures/generate_clock_validation.py`** | Simulated demo data via `methylation_visualizations` | `analysis/Fig3_Clock_Validation.png` (layout mock, not real cohort) |

## Related documentation

- [CLOCK_LIBRARY.md](CLOCK_LIBRARY.md) — package API and CLI  
- [GSE40279_CLOCK_TRAINING.md](GSE40279_CLOCK_TRAINING.md) — training and GSE87571 external validation  
- [FIGURES.md](FIGURES.md) — full figure script index  
- [ACTIVITIES.md](ACTIVITIES.md#21101--methylation-aging-clock) — activity 2.1.10.1 map
