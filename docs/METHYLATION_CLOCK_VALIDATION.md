# Methylation clock validation (`evaluate_methylation_clock.py`)

**Project:** IBAR-ROGEN Aging  
**Activity:** 2.1.10.1 — methylation aging clock (GSE40279 train, GSE87571 validate)  
**Script:** [`evaluate_methylation_clock.py`](../evaluate_methylation_clock.py)  
**Library:** `rogen_aging.clock` (`load_model`, `build_feature_matrix`)  
**Related:** [CLOCK_LIBRARY.md](CLOCK_LIBRARY.md) · [CLOCK_EVAL_FIGURES.md](CLOCK_EVAL_FIGURES.md) · [GSE40279_CLOCK_TRAINING.md](GSE40279_CLOCK_TRAINING.md)

## Purpose

Publication-oriented **external validation** of a trained ElasticNet DNAm-age clock on an independent GSE87571 cohort. The script:

1. Joins a processed methylation beta matrix with phenotype ages.
2. Aligns CpGs to the model’s training features (mean-imputing missing probes).
3. Computes MAE, median absolute error, Pearson *r*, and age-stratified MAE (`<30`, `30–60`, `>60`).
4. Writes a three-panel figure: predicted-vs-chronological scatter (with 95% CI), residual plot, and top-25 CpG weights labeled by nearest gene.

This complements `rogen-clock evaluate` (decade MAE / residual PNGs) and `plot_clock_eval.py` (two-panel figure under `figures/`).

## Required inputs

| Input | Default path | Description |
|-------|--------------|-------------|
| Methylation matrix | `data/methylation/GSE87571_processed.parquet` | Wide beta table (`sample_id` + `cg*` columns; probes×samples also accepted) |
| Phenotype metadata | `data/methylation/GSE87571_meta.csv` | `sample_id` (or GEO accession) + `chronological_age` / `age` |
| Trained model | `models/ro_clock_elasticnet_gse40279.pkl` | Fitted bare `sklearn.linear_model.ElasticNet` (Pipelines / ElasticNetCV rejected) |
| Probe→gene annotation (optional) | `data/methylation/HM450_probe_annotation.csv` | `IlmnID` + `UCSC_RefGene_Name` (falls back to Horvath S3 table) |

Preflight: the script reads [`INPUT_MANIFEST.md`](../INPUT_MANIFEST.md) and aborts if any required file is missing.

## Outputs

| Output | Path |
|--------|------|
| Metrics JSON | `outputs/clock_metrics.json` |
| Figure (raster) | `outputs/figures/Figure_Epigenetic_Clock_Panels.png` (300 dpi) |
| Figure (vector) | `outputs/figures/Figure_Epigenetic_Clock_Panels.pdf` |
| Markdown summary | Printed to stdout |

## CLI usage

```bash
uv sync

uv run python evaluate_methylation_clock.py

# Explicit paths / options
uv run python evaluate_methylation_clock.py \
  --methylation data/methylation/GSE87571_processed.parquet \
  --meta data/methylation/GSE87571_meta.csv \
  --model models/ro_clock_elasticnet_gse40279.pkl \
  --metrics-out outputs/clock_metrics.json \
  --figure-stem outputs/figures/Figure_Epigenetic_Clock_Panels \
  --annotation data/methylation/HM450_probe_annotation.csv \
  --top-n 25

# Skip INPUT_MANIFEST.md preflight when overriding paths in tests/CI
uv run python evaluate_methylation_clock.py \
  --model /tmp/ro_clock.pkl \
  --methylation /tmp/meth.parquet \
  --meta /tmp/meta.csv \
  --skip-manifest-check
```

## Technical notes

- **Estimator contract:** the pickle must unpickle to exactly `sklearn.linear_model.ElasticNet`. `Pipeline`, `ElasticNetCV`, and other regressors raise `TypeError` (no silent retrain).
- **Feature alignment:** uses `rogen_aging.clock.evaluate.build_feature_matrix` to reorder CpGs to `feature_names_in_` and mean-impute probes missing from GSE87571.
- **Strata:** middle bin is closed `[30, 60]`; empty strata report `MAE = null` in JSON / `NA` in the markdown summary.
- **Tests:** `uv run pytest tests/test_evaluate_methylation_clock.py -q`

## Metrics definition

| Metric | Definition |
|--------|------------|
| MAE | Mean \|chronological − predicted\| (years) |
| Median Absolute Error | Median \|chronological − predicted\| |
| Pearson *r* | Correlation of chronological vs predicted age |
| Age-stratified MAE | MAE within `<30`, `30–60` (inclusive), and `>60` years |
| Residual (panel B) | chronological − predicted |

## See also

- [CLOCK_LIBRARY.md](CLOCK_LIBRARY.md) — package API and `rogen-clock` CLI
- [CLOCK_EVAL_FIGURES.md](CLOCK_EVAL_FIGURES.md) — two-panel `plot_clock_eval.py` figure
- [ACTIVITIES.md](ACTIVITIES.md#21101--methylation-aging-clock) — activity 2.1.10.1 map

---

**Last updated:** July 27, 2026
