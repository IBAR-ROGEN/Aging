# Romanian cohort mock epigenetic clock (Elastic Net)

**Project:** IBAR-ROGEN Aging  
**Scripts:** `scripts/train_romanian_epigenetic_clock.py` (train), `scripts/validate_clock.py` (evaluate a saved model on held-out data). For a **wide-table** trainer aimed at public **GSE40279** (Hannum 2013) style inputs, see **`scripts/train_clock_on_gse40279.py`** and **[docs/GSE40279_CLOCK_TRAINING.md](GSE40279_CLOCK_TRAINING.md)**.  
**Repository:** [IBAR-ROGEN/Aging](https://github.com/IBAR-ROGEN/Aging)

## Overview

`train_romanian_epigenetic_clock.py` trains a **custom epigenetic aging clock** using **Elastic Net regression** with **`ElasticNetCV`** from scikit-learn. It is wired for a **Romanian-style cohort** naming convention (`ROM0001`, …) and uses **synthetic methylation** when real array data is not yet available.

The workflow is suitable for:

- Prototyping clock training before 450K/EPIC beta matrices are finalized
- Verifying sklearn pipelines (scaling, cross-validation, metrics, plots)
- Teaching or benchmarking MAE and Pearson correlation against chronological age

**Important:** Default inputs are **mock** (simulated β-values with a sparse age-related signal). For production clocks, replace the CSVs under `data/` with real preprocessed methylation and validated ages, and use rigorous nested cross-validation or external validation cohorts.

## Input files

Place files under a directory (default: `data/mock_romanian_cohort/`, git-ignored). The script **creates** both CSVs automatically if they are missing, or you can supply your own.

### `methylation_matrix.csv`

- **Rows:** one per sample  
- **Columns:** `sample_id` plus one column per probe (e.g. Illumina probe IDs or arbitrary feature names)  
- **Values:** methylation levels (typically β in \([0,1]\) for array data)

### `metadata.csv`

- **`sample_id`** — must match `methylation_matrix.csv`  
- **`chronological_age`** — age in years (regression target)

Samples are aligned by an **inner join** on `sample_id`; only overlapping IDs are used.

## Model and tuning

1. **`StandardScaler`** on all CpG (feature) columns (recommended for penalized regression on many correlated probes).  
2. **`ElasticNetCV`** with **`cv=5`**:  
   - Searches over a grid of **`l1_ratio`** values (L1 vs L2 mix in the elastic net penalty).  
   - For each `l1_ratio`, selects **`alpha`** by cross-validation.  
   - Chooses the `(l1_ratio, alpha)` pair with the best CV score.

Hyperparameters selected during CV are printed after fitting.

## Evaluation split

The script holds out a **random test fraction** (default **20%**) for reporting **MAE** (years) and **Pearson r** between chronological and predicted age. The internal CV in `ElasticNetCV` only uses the **training** portion passed into `fit` for tuning; the reported MAE and **r** refer to the **held-out test set** (not the full cohort), which avoids optimistically scoring on training data.

For publication-grade clocks, plan additional **external** validation and consider **nested CV** when comparing many preprocessing or feature sets.

## Outputs

1. **Stdout** — training/test sizes, chosen `l1_ratio` and `alpha`, test **MAE**, test **Pearson r**.  
2. **PNG figure** (default: `figures/romanian_mock_epigenetic_clock_scatter.png`) — scatter of chronological age (x) vs predicted age (y), **line of best fit**, **identity line** (y = x), and an annotation box with **MAE** and **r**.

The `figures/` directory is git-ignored (generated artifact); regenerate locally after cloning.

## Usage

```bash
uv sync
uv run python scripts/train_romanian_epigenetic_clock.py --help
```

**Default run** (creates mock CSVs under `data/mock_romanian_cohort/` if needed):

```bash
uv run python scripts/train_romanian_epigenetic_clock.py
```

**Regenerate mock data** after changing simulation defaults in the script (overwrites CSVs in the chosen `data_dir`):

```bash
uv run python scripts/train_romanian_epigenetic_clock.py --regenerate-mock
```

**Custom paths** (paths are whatever you pass; both flags accept any valid filesystem location):

```bash
uv run python scripts/train_romanian_epigenetic_clock.py \
  --data-dir cohort_data \
  --output-plot artifacts/clock_scatter.png \
  --test-size 0.25
```

## Dependencies

Declared in `pyproject.toml`: **polars**, **numpy**, **scikit-learn**, **scipy**, **matplotlib**, **typer**.

## Held-out validation (`validate_clock.py`)

**Script:** `scripts/validate_clock.py`  
Use this after you have a **saved** fitted estimator (for example from a training notebook or a Typer/CLI trainer that writes `joblib.dump` / `pickle`).

### Purpose

- Load a pre-trained clock and score a **held-out** cohort in one table (bedMethyl-style or wide β matrix).
- Report **MAE** (years) and **Pearson r** (chronological vs predicted age).
- Stratify **MAE by chronological age decade** (`<20`, `20-29`, …, `90+`) using `pandas.cut`.
- Save **figures** and a **JSON** summary for manuscripts or QC.

### CLI arguments

| Argument | Description |
|----------|-------------|
| `--model_path` | Path to `.joblib` or `.pkl` fitted model (pipeline or estimator). |
| `--test_data` | `.parquet` or `.csv` (`.tsv` tab-separated supported). |
| `--output_dir` | Directory created if needed; receives plots and metrics. |

### Test data expectations

- A column **`chronological_age`** (years).
- Feature columns whose names **start with `cg`** (CpG probe IDs or site labels).

### Outputs (under `--output_dir`)

| File | Content |
|------|---------|
| `validation_metrics.json` | Overall MAE, Pearson r and p-value, sample and feature counts, list of **imputed** model CpGs (if any), MAE per decade. |
| `Fig_Clock_Residuals.png` | Scatter: residual (predicted − chronological) vs chronological age; horizontal line at 0. |
| `Fig_Clock_MAE_by_decade.png` | Bar chart: MAE by age decade. |

Stdout prints the same metrics (JSON); stderr notes how many CpGs were imputed when applicable.

### Feature alignment and missing CpGs

- If the fitted model exposes **`feature_names_in_`** (typical when `fit` was called on a **pandas `DataFrame`**), the script builds `X` in that column order. Any **expected probe missing** from the test file is filled with the **mean of all present `cg*` values** in the test set (fallback `0.5` if undefined), and a **warning** is emitted.
- If **no** feature names are stored (for example the model was fit on a **NumPy** array only), the script uses **all `cg*` columns** in file order and imputes **NaNs** via row/column means, then `0.5`. If **`n_features_in_`** is known and does not match the number of `cg*` columns, the script **raises** a clear error (cannot guess probe identity). For that case, refit on a `DataFrame` with named columns or export the training feature list and preprocess test data to match.

### Example

```bash
uv sync
uv run python scripts/validate_clock.py \
  --model_path artifacts/romanian_clock.joblib \
  --test_data data/romanian_holdout.parquet \
  --output_dir figures/clock_validation_run1
```

### External cohort (GSE87571)

To validate a GSE40279-trained clock on a second public 450K whole-blood series, use **`rogen_aging.clock.external_data`** (`load_gse87571` / `python -m rogen_aging.clock.external_data`) to produce a Parquet with **`chronological_age`** and **`cg*`** columns, then point **`validate_clock.py`** at that file. See **[docs/GSE40279_CLOCK_TRAINING.md](GSE40279_CLOCK_TRAINING.md#external-validation-gse87571)**.

## Related material

- **[docs/GSE40279_CLOCK_TRAINING.md](GSE40279_CLOCK_TRAINING.md)** — Train an Elastic Net clock from a sample × CpG table (+ `chronological_age`) for GEO GSE40279-style data (`train_clock_on_gse40279.py`), and notes on **GSE87571** external validation (`rogen_aging.clock.external_data`).  
- **[docs/EDA_MOCK_INTEGRATION.md](EDA_MOCK_INTEGRATION.md)** — EDA on mock clinical/epigenetic age tables (`test_data/mock_epigenetic_clinical.csv`).  
- **Notebook:** `notebooks/02_methylation_pipeline/MethylationClocks.ipynb` — broader epigenetic clock context.  
- **Synthetic Romanian VCF:** [docs/SYNTHETIC_ROMANIAN_VCF_GENERATOR.md](SYNTHETIC_ROMANIAN_VCF_GENERATOR.md) (genotype mock data; separate from this methylation clock script).
