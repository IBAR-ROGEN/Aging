# Epigenetic clock library (`rogen_aging.clock`)

**Project:** IBAR-ROGEN Aging  
**Package:** `src/rogen_aging/clock/`  
**Canonical CLI:** `uv run rogen-clock` or `scripts/clock/run_clock.py`

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
| `model.py` | `make_clock_pipeline()` — `SimpleImputer` + `ElasticNetCV` |
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
  --output_dir analysis/clock_validation/
```

Wide-table convention: `cg*` probe columns + `chronological_age` target.

## Tests

```bash
uv run pytest tests/test_clock_regression.py tests/test_package_imports.py -q
```

## Related documentation

- [GSE40279 Clock Training](GSE40279_CLOCK_TRAINING.md)
- [Romanian Epigenetic Clock](ROMANIAN_EPIGENETIC_CLOCK.md)
- [WORKFLOWS.md](WORKFLOWS.md)

---

**Last updated:** May 31, 2026
