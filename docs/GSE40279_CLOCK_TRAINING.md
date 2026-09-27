# GSE40279 (Hannum 2013) Elastic Net clock training

**Project:** IBAR-ROGEN Aging  
**Canonical CLI:** `uv run rogen-clock prepare-gse40279`, `train`, or `evaluate`  
**Library:** `src/rogen_aging/clock/` — see **[docs/CLOCK_LIBRARY.md](CLOCK_LIBRARY.md)**  
**Legacy wrappers:** `scripts/clock/train_clock_on_gse40279.py`, `scripts/clock/validate_clock.py` (deprecated)

## Overview

[GSE40279](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40279) is a public Illumina HumanMethylation450 whole-blood dataset (Hannum et al., 2013; on the order of hundreds of samples). This repository provides a **training-only** CLI that fits a **`Pipeline(SimpleImputer, ElasticNetCV)`** on a **wide** table: one row per sample, one column per CpG probe (names starting with `cg`), plus **`chronological_age`**.

`rogen-clock prepare-gse40279` downloads the GEO series matrix and writes the wide beta table. It does not train a model.

```bash
uv run rogen-clock prepare-gse40279 \
  --output data/methylation/GSE40279_processed.parquet
```

| Path | Role |
|------|------|
| `data/methylation/GSE40279_processed.parquet` | Wide betas (`sample_id`, `chronological_age`, `cg*` float32) |
| `models/gse40279_hannum450k_elasticnet.joblib` | Production clock. [CLOCK_ARTIFACT_AUDIT.md](CLOCK_ARTIFACT_AUDIT.md) flags it CANDIDATE REAL |
| `models/gse40279_hannum450k_elasticnet.provenance.json` | Training provenance |
| `models/gse40279_hannum450k_elasticnet_cpgs.csv` | Non-zero CpGs and weights |
| `models/fixtures/fixture_clock_cg_test.pkl` | Demo clock. Probe names are `cg_test_*` |

`rogen-clock evaluate` refuses a model whose features are all `cg_test_*` unless `--demo` is passed. The fixture path above is that demo clock. The production clock is `models/gse40279_hannum450k_elasticnet.joblib`.

Which serializations are fixtures is the audit report, not a filename. Regenerate it with `uv run python scripts/audit_clock_artifacts.py`.

### Provenance fields

`models/gse40279_hannum450k_elasticnet.provenance.json` records `geo_accession`, `model_path`, `input_parquet`, `input_parquet_sha256`, `n_cpgs_features`, `n_train_samples`, `n_test_samples`, `probe_counts` (`input_cg`, `overlap_gse87571`, `after_missingness_le_5pct_on_train`, `after_abs_age_correlation_on_train`, `trained`), `filters`, `n_nonzero_coefficients`, `nonzero_cpg_csv`, `hyperparameters` (`alpha`, `l1_ratio`, `cv`, `cv_note`), `random_state`, `sklearn_version`, `git_commit`, `timestamp_utc`, and `holdout_metrics_from_train_clock`.

From that file: `n_train_samples` 524, `n_test_samples` 132, `n_cpgs_features` 8000, `probe_counts.input_cg` 470043, `probe_counts.overlap_gse87571` 470043, `probe_counts.after_missingness_le_5pct_on_train` 470043, `probe_counts.after_abs_age_correlation_on_train` 51165, `probe_counts.trained` 8000, `n_nonzero_coefficients` 2314. `hyperparameters.cv` is 5. The `cv_note` says cv was lowered from 10 to 5 because 8000 training probes would make ElasticNetCV exceed about 2 hours on a laptop, and that no other hyperparameter was changed.

Hold-out metrics in `holdout_metrics_from_train_clock` of `models/gse40279_hannum450k_elasticnet.provenance.json` (the same three values are in `outputs/gse40279_hannum450k_train_metrics.json`): `test_mae` 3.5248637199401855, `test_rmse` 4.727581024169922, `test_pearson_r` 0.9566776752471924.

GSE87571 evaluation of this joblib: not yet run. `outputs/clock_metrics.json` is an older file whose `model_path` is `models/ro_clock_elasticnet_gse40279.pkl`, not the Hannum joblib.

## Expected input format

| Requirement | Detail |
|-------------|--------|
| Layout | Rows = samples, columns = features + target |
| CpG columns | Any column whose name **starts with** `cg` (Illumina-style probe IDs) |
| Target | Single column **`chronological_age`** (years), numeric with no missing values |
| Missing β | Allowed in CpG columns; **training split** is mean-imputed per feature (imputer fit on train only) |
| File types | **`.parquet`**, **`.csv`**, or **`.tsv`** |

## Model and evaluation

1. **Train/test split:** `sklearn.model_selection.train_test_split` with `--test_size` (default `0.2`) and `--random_state` (default `42`).
2. **Pipeline:** `SimpleImputer(strategy="mean")` → **`ElasticNetCV`** with `cv=10`, `l1_ratio=[0.1, 0.5, 0.7, 0.9, 0.95, 0.99, 1.0]`, `alphas=20`, `max_iter=5000`.
3. **Held-out metrics:** MAE, RMSE (`root_mean_squared_error`), Pearson **r** (`scipy.stats.pearsonr`) on the test split.

`alphas=20` means **20 automatically generated regularization strengths** on the log-spaced path (one grid per `l1_ratio` value). This is the scikit-learn **1.9+** spelling; older releases used the deprecated `n_alphas=20` alias for the same behaviour.

The saved object is the **fitted `Pipeline`**, written with **`joblib.dump`**. It exposes **`feature_names_in_`**, so **`rogen-clock evaluate`** can align probes and mean-impute missing expected CpGs on new cohorts.

### scikit-learn compatibility

| Topic | Detail |
|-------|--------|
| Minimum version | `scikit-learn>=1.6.0` in `pyproject.toml` |
| Hyperparameter | Use **`alphas=20`** in `make_clock_pipeline()` and tests — required from **1.9** onward (`n_alphas` was removed) |
| Regression test | `tests/test_clock_regression.py` compares refactored `train_clock()` against inlined legacy training logic on `test_data/mock_clock_wide.csv` |
| CI | GitHub Actions runs `uv sync --extra dev` then `uv run pytest -q` on `ubuntu-latest` (see `.github/workflows/ci.yml`) |

## CLI

| Argument | Description |
|----------|-------------|
| `--input_data` | Path to Parquet or CSV/TSV wide table |
| `--output_model` | Path for the pickled pipeline (`.pkl` / `.joblib` extension as you prefer) |
| `--output_metrics` | Path for training metrics JSON |
| `--test_size` | Test fraction (default `0.2`) |
| `--random_state` | Random seed (default `42`) |

## Outputs

1. **`--output_model`** — Fitted sklearn `Pipeline` (imputer + `ElasticNetCV`).
2. **`--output_metrics`** — JSON including `test_mae`, `test_rmse`, `test_pearson_r`, chosen `alpha`, `l1_ratio`, total `cg*` count, **`selected_cpgs`** (probes with non-zero elastic-net coefficients), and split metadata.

Stdout prints a one-line summary: number of CpGs, `alpha`, `l1_ratio`, test MAE, test **r**.

## Usage

```bash
uv sync
uv run rogen-clock train --help
```

Example (paths are illustrative):

```bash
uv run rogen-clock train \
  --input_data data/methylation/GSE40279_processed.parquet \
  --output_model models/gse40279_hannum450k_elasticnet.joblib \
  --output_metrics outputs/gse40279_hannum450k_train_metrics.json
```

Equivalent script path:

```bash
uv run python scripts/clock/run_clock.py train \
  --input_data data/methylation/GSE40279_processed.parquet \
  --output_model models/gse40279_hannum450k_elasticnet.joblib \
  --output_metrics outputs/gse40279_hannum450k_train_metrics.json
```

## Validating the saved model

Use **`uv run rogen-clock evaluate`** with a held-out table in the same wide format. See **[docs/ROMANIAN_EPIGENETIC_CLOCK.md](ROMANIAN_EPIGENETIC_CLOCK.md#held-out-validation-validate_clockpy)** for CLI arguments, figures, and `validation_metrics.json`.

### External validation (GSE87571)

[GSE87571](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE87571) is a public Illumina 450K whole-blood cohort (hundreds of samples, wide age range) suitable for **external** checks on clocks trained on GSE40279-style inputs.

Use the installable module **`rogen_aging.clock.external_data`** to build an evaluation-compatible Parquet (or call **`load_gse87571`** then **`save_as_parquet`** in Python):

| Item | Detail |
|------|--------|
| Module | `src/rogen_aging/clock/external_data.py` |
| API | `load_gse87571(local_path=None, geo_cache_dir='./data/geo', restrict_to_cpgs=None)` → `pandas.DataFrame` with index = GSM, columns = `chronological_age` + `cg*` |
| CLI | `uv run python -m rogen_aging.clock.external_data --output path/to/gse87571.parquet` |
| GEO layout | The series matrix carries metadata only; probe β-values are read from supplementary **`GSE87571_matrix1of2.txt.gz`** and **`GSE87571_matrix2of2.txt.gz`**, merged automatically after download into `geo_cache_dir`. |
| `restrict_to_cpgs` | Optional sequence of probe IDs; keeps only the intersection with the cohort (for example training JSON **`selected_cpgs`**). |
| Manual fetch | If **GEOparse** or HTTPS download fails (timeouts, TLS, or corporate proxies), download **`GSE87571_series_matrix.txt.gz`** and the two supplementary matrix files from GEO and pass **`local_path`** to the series matrix; cache the `.gz` files under the same **`geo_cache_dir`** so the loader can find them. |

Example end-to-end (paths illustrative):

```bash
uv run python -m rogen_aging.clock.external_data \
  --output data/gse87571.parquet \
  --restrict-cpgs-file data/my_clock_selected_cpgs.txt

uv run rogen-clock evaluate \
  --model_path models/gse40279_hannum450k_elasticnet.joblib \
  --test_data data/gse87571.parquet \
  --output_dir figures/validation_gse87571

# Optional: publication two-panel figure (scatter + top CpG weights)
uv run python scripts/figures/plot_clock_eval.py
```

See **[docs/CLOCK_EVAL_FIGURES.md](CLOCK_EVAL_FIGURES.md)** for `plot_clock_eval.py` configuration, per-sample input options, and output paths.

## Related documentation

- **[docs/CLOCK_EVAL_FIGURES.md](CLOCK_EVAL_FIGURES.md)** — External-validation scatter + top-CpG figure (`plot_clock_eval.py`).
- **[docs/ROMANIAN_EPIGENETIC_CLOCK.md](ROMANIAN_EPIGENETIC_CLOCK.md)** — Romanian-style mock trainer and shared validation notes.
- **[docs/WORKFLOWS.md](WORKFLOWS.md)** — Clock workflow index.
- **`notebooks/02_methylation_pipeline/`** — Broader methylation and clock context.
