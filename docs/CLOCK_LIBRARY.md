# Epigenetic clock library (`rogen_aging.clock`)

**Project:** IBAR-ROGEN Aging  
**Package:** `src/rogen_aging/clock/`  
**Unified CLI:** `scripts/run_clock.py`

## Overview

Clock training and evaluation logic lives in the installable **`rogen_aging.clock`** package. Legacy scripts remain as thin wrappers for backward compatibility:

| Script | Role |
|--------|------|
| `scripts/run_clock.py` | Unified **`train`** / **`evaluate`** subcommands |
| `scripts/train_clock_on_gse40279.py` | GSE40279-style wide-table training wrapper |
| `scripts/train_romanian_epigenetic_clock.py` | Romanian mock cohort demo (StandardScaler + ElasticNetCV) |
| `scripts/validate_clock.py` | Held-out evaluation wrapper |

## Package modules

| Module | Purpose |
|--------|---------|
| `data.py` | `load_wide_table`, `split_features_target`, Romanian mock cohort I/O; re-exports GSE87571 loader |
| `model.py` | `make_clock_pipeline()` — `SimpleImputer` + `ElasticNetCV` |
| `train.py` | `train_clock()` — fit, metrics JSON, joblib model |
| `evaluate.py` | `evaluate_clock()` — MAE, r, decade MAE, scatter/residual figures |
| `external_data.py` | GSE87571 GEO download / merge → Parquet for external validation |

Public API (from `rogen_aging.clock`):

```python
from rogen_aging.clock import train_clock, evaluate_clock, load_wide_table, load_gse87571
```

## Unified CLI

### Train

```bash
uv run python scripts/run_clock.py train \
  --input_data data/gse40279_beta_age.parquet \
  --output_model analysis/gse40279_elasticnet_clock.pkl \
  --output_metrics analysis/gse40279_train_metrics.json
```

### Evaluate

```bash
uv run python scripts/run_clock.py evaluate \
  --model_path analysis/gse40279_elasticnet_clock.pkl \
  --test_data data/gse87571.parquet \
  --output_dir analysis/clock_validation/
```

Same wide-table convention: `cg*` probe columns + `chronological_age` target.

## Tests

```bash
uv run pytest tests/test_clock_regression.py tests/test_package_imports.py -q
```

`test_clock_regression.py` asserts refactored `train_clock()` matches the pre-refactor GSE40279 training metrics on `test_data/mock_clock_wide.csv`.

## Related documentation

- [GSE40279 Clock Training](GSE40279_CLOCK_TRAINING.md)
- [Romanian Epigenetic Clock](ROMANIAN_EPIGENETIC_CLOCK.md)
- [Code Modules Reference](CODE_MODULES_REFERENCE.md)

---

**Last updated:** May 31, 2026
