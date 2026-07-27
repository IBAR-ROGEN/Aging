# Epigenetic clock library (`rogen_aging.clock`)

**Project:** IBAR-ROGEN Aging  
**Activity:** 2.1.10.1 — methylation aging clock (GSE40279 training, GSE87571 external validation)  
**Package:** `src/rogen_aging/clock/`  
**Canonical CLI:** `uv run rogen-clock` or `scripts/clock/run_clock.py`  
**Activity map:** [ACTIVITIES.md](ACTIVITIES.md#21101--methylation-aging-clock)

## Overview

Clock training and evaluation logic lives in the installable **`rogen_aging.clock`** package.

| Entry | Role |
|-------|------|
| **`rogen-clock train\|evaluate`** | **Recommended** unified CLI |
| `scripts/clock/run_clock.py` | Same as console entry |
| `scripts/clock/train_clock_on_gse40279.py` | **Deprecated** GSE40279 wrapper |
| `scripts/clock/validate_clock.py` | **Deprecated** evaluation wrapper |
| `scripts/clock/train_romanian_epigenetic_clock.py` | Romanian mock demo (StandardScaler + ElasticNetCV; not the main library path) |

Flat paths `scripts/run_clock.py`, `scripts/validate_clock.py`, etc. are deprecation shims.

## Package modules

| Module | Purpose |
|--------|---------|
| `data.py` | `load_wide_table`, Romanian mock cohort I/O |
| `model.py` | `make_clock_pipeline()` — `SimpleImputer` + `ElasticNetCV` (`alphas=20`, `cv=10`, fixed `l1_ratio` grid) |
| `train.py` | `train_clock()` |
| `evaluate.py` | `evaluate_clock()` |
| `external_data.py` | GSE87571 → Parquet for external validation |

```python
from rogen_aging.clock import train_clock, evaluate_clock, load_wide_table, load_gse87571
```

## CLI examples

```bash
uv run rogen-clock train \
  --input_data data/gse40279_beta_age.parquet \
  --output_model analysis/gse40279_elasticnet_clock.pkl \
  --output_metrics analysis/gse40279_train_metrics.json

uv run rogen-clock evaluate \
  --model_path analysis/gse40279_elasticnet_clock.pkl \
  --test_data data/gse87571.parquet \
  --output_dir figures/validation_gse87571
```

Wide-table convention: `cg*` probe columns + `chronological_age` target.

## Publication figure (external validation)

After train + evaluate on GSE87571, generate a two-panel scatter (predicted vs chronological age) and top-CpG weight bar chart:

```bash
uv run python scripts/figures/plot_clock_eval.py
# → figures/validation_gse87571/clock_eval_gse87571.png + .pdf
```

See **[CLOCK_EVAL_FIGURES.md](CLOCK_EVAL_FIGURES.md)** for configuration constants, input options, and how this relates to `rogen-clock evaluate` outputs.

For final validation metrics (MAE, median AE, Pearson *r*, age-stratified MAE) and a three-panel publication figure written under `outputs/`:

```bash
uv run python scripts/clock/evaluate_methylation_clock.py
# → outputs/clock_metrics.json
# → outputs/figures/Figure_Epigenetic_Clock_Panels.png/.pdf
```

Required inputs are listed in [`INPUT_MANIFEST.md`](../INPUT_MANIFEST.md) (bare ElasticNet pickle at `models/ro_clock_elasticnet_gse40279.pkl`). Override paths or pass `--skip-manifest-check` for ad-hoc runs.

See **[METHYLATION_CLOCK_VALIDATION.md](METHYLATION_CLOCK_VALIDATION.md)**.

## Tests

```bash
uv run pytest tests/test_clock_regression.py tests/test_evaluate_methylation_clock.py tests/test_package_imports.py -q
```

`test_clock_regression.py` asserts that `train_clock()` reproduces the pre-refactor GSE40279 training metrics on `test_data/mock_clock_wide.csv`. Both code paths call the real `sklearn.linear_model.ElasticNetCV` with **`alphas=20`** (20-point alpha grid per `l1_ratio`). Use `alphas`, not the removed `n_alphas` keyword, so tests pass on scikit-learn **1.9+** in CI.

`test_evaluate_methylation_clock.py` covers manifest preflight, bare-ElasticNet loading (Pipeline / ElasticNetCV rejection), age-stratum metrics, cohort alignment edge cases, and end-to-end metric/figure writers.

## Related documentation

- [GSE40279 Clock Training](GSE40279_CLOCK_TRAINING.md)
- [Clock eval figure (`plot_clock_eval.py`)](CLOCK_EVAL_FIGURES.md)
- [Methylation clock validation (`evaluate_methylation_clock.py`)](METHYLATION_CLOCK_VALIDATION.md)
- [Romanian Epigenetic Clock](ROMANIAN_EPIGENETIC_CLOCK.md)
- [WORKFLOWS.md](WORKFLOWS.md)

---

**Last updated:** July 27, 2026
