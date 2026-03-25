# Romanian cohort mock epigenetic clock (Elastic Net)

**Project:** IBAR-ROGEN Aging  
**Script:** `scripts/train_romanian_epigenetic_clock.py`  
**Repository:** [IBAR-ROGEN/Aging](https://github.com/IBAR-ROGEN/Aging)

## Overview

This script trains a **custom epigenetic aging clock** using **Elastic Net regression** with **`ElasticNetCV`** from scikit-learn. It is wired for a **Romanian-style cohort** naming convention (`ROM0001`, …) and uses **synthetic methylation** when real array data is not yet available.

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

## Related material

- **[docs/EDA_MOCK_INTEGRATION.md](EDA_MOCK_INTEGRATION.md)** — EDA on mock clinical/epigenetic age tables (`test_data/mock_epigenetic_clinical.csv`).  
- **Notebook:** `notebooks/02_methylation_pipeline/MethylationClocks.ipynb` — broader epigenetic clock context.  
- **Synthetic Romanian VCF:** [docs/SYNTHETIC_ROMANIAN_VCF_GENERATOR.md](SYNTHETIC_ROMANIAN_VCF_GENERATOR.md) (genotype mock data; separate from this methylation clock script).
