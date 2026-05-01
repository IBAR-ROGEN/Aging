# GSE40279 (Hannum 2013) Elastic Net clock training

**Project:** IBAR-ROGEN Aging  
**Script:** `scripts/train_clock_on_gse40279.py`  
**Validation:** `scripts/validate_clock.py` (same `cg*` + `chronological_age` convention as test data)

## Overview

[GSE40279](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40279) is a public Illumina HumanMethylation450 whole-blood dataset (Hannum et al., 2013; on the order of hundreds of samples). This repository provides a **training-only** CLI that fits a **`Pipeline(SimpleImputer, ElasticNetCV)`** on a **wide** table: one row per sample, one column per CpG probe (names starting with `cg`), plus **`chronological_age`**.

The script does **not** download GEO Series Matrix files, parse IDATs, or map probe annotations. You must obtain the data from NCBI GEO and convert it to the expected column layout yourself (for example in R with `minfi`/`GEOquery`, or in Python after exporting β-values).

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
2. **Pipeline:** `SimpleImputer(strategy="mean")` → **`ElasticNetCV`** with `cv=10`, `l1_ratio=[0.1, 0.5, 0.7, 0.9, 0.95, 0.99, 1.0]`, `n_alphas=20`, `max_iter=5000`.
3. **Held-out metrics:** MAE, RMSE (`root_mean_squared_error`), Pearson **r** (`scipy.stats.pearsonr`) on the test split.

The saved object is the **fitted `Pipeline`**, written with **`joblib.dump`**. It exposes **`feature_names_in_`**, so **`validate_clock.py`** can align probes and mean-impute missing expected CpGs on new cohorts.

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
uv run python scripts/train_clock_on_gse40279.py --help
```

Example (paths are illustrative):

```bash
uv run python scripts/train_clock_on_gse40279.py \
  --input_data data/gse40279_beta_age.parquet \
  --output_model analysis/gse40279_elasticnet_clock.pkl \
  --output_metrics analysis/gse40279_train_metrics.json
```

## Validating the saved model

Use **`scripts/validate_clock.py`** with a held-out table in the same wide format. See **[docs/ROMANIAN_EPIGENETIC_CLOCK.md](ROMANIAN_EPIGENETIC_CLOCK.md#held-out-validation-validate_clockpy)** for CLI arguments, figures, and `validation_metrics.json`.

### External validation (GSE87571)

[GSE87571](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE87571) is a public Illumina 450K whole-blood cohort (hundreds of samples, wide age range) suitable for **external** checks on clocks trained on GSE40279-style inputs.

Use the installable module **`rogen_aging.clock.external_data`** to build a **`validate_clock.py`**-compatible Parquet (or call **`load_gse87571`** then **`save_as_parquet`** in Python):

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

uv run python scripts/validate_clock.py \
  --model_path analysis/gse40279_elasticnet_clock.pkl \
  --test_data data/gse87571.parquet \
  --output_dir analysis/validation_gse87571
```

## Related documentation

- **[docs/ROMANIAN_EPIGENETIC_CLOCK.md](ROMANIAN_EPIGENETIC_CLOCK.md)** — Romanian-style mock trainer (`train_romanian_epigenetic_clock.py`) and shared validation notes for `validate_clock.py`.
- **`notebooks/02_methylation_pipeline/`** — Broader methylation and clock context.
