# Published clocks vs retrained Hannum

Benchmark of pyaging DNA-methylation clocks against the ElasticNet Hannum model
retrained on GSE40279 (`models/gse40279_hannum450k_elasticnet.joblib`).

Numbers below are taken only from the files under `analysis/clocks_benchmark/`
after the repair rerun.

## Cohorts

| Cohort | Samples | Notes |
|--------|---------|--------|
| GSE87571 | 729 with age + sex | Fair external comparison for published Hannum |
| GSE40279 held-out | 132 (`test_size=0.2`, `random_state=42`) | Betas from the **full GEO series matrix** (not the training-prepared parquet). Retrained model never scores train rows; published Hannum is **in-sample** here |

Sex is parsed from GEO `gender:` / `sex:` characteristics. Ethnicity (GSE40279 only)
is reported descriptively; no formal ethnicity test.

## Clocks

Target pyaging clocks: `horvath2013`, `hannum`, `dnamphenoage`, `dunedinpace`.
(`phenoage` in pyaging is clinical Phenotypic Age, not DNAm PhenoAge.)

`grimage2` is included. `female` (0/1) and chronological `age` are **data columns**
passed to `df_to_adata` and are **not** in `metadata_cols`. After predict, provenance
asserts that neither appears in `grimage2` `missing_features`.

Retrained model label: `rogen_hannum_retrain`. Held-out MAE on GSE40279 must stay
≈ 3.5249 (n=132).

## How to run

```bash
uv run python scripts/clock/benchmark_published_clocks.py
```

Inputs:

- `data/geo/GSE87571_series_matrix.txt.gz` + supplementary beta matrices
- `data/geo/GSE40279_series_matrix.txt.gz` (benchmark matrix source)
- `models/gse40279_hannum450k_elasticnet.joblib`

CpG columns are restricted to the union of probes required by the selected clocks
plus the retrained Hannum feature list. pyaging uses `df_to_adata(..., imputer_strategy="knn")`.
CpGs absent from the GEO file itself stay missing until that KNN step
(`provenance.json` → `gse40279_geo_probe_report`).

## Metrics

Per clock × cohort:

- MAE vs chronological age (**empty for DunedinPACE**)
- Pearson *r*
- Mean signed error (pred − age; **empty for DunedinPACE**)
- Age acceleration = residual of clock ~ chronological age (within cohort); kept for DunedinPACE
- Sex model: `accel ~ female + age`; female coefficient, 95% CI, p; Bonferroni across clocks within cohort
- Acceleration correlation matrix across clocks (per cohort)
- GSE40279 ethnicity: counts / mean accel only

## Outputs

Under `analysis/clocks_benchmark/`:

| File | Contents |
|------|----------|
| `metrics_by_clock_cohort.csv` | Aggregate accuracy |
| `acceleration_by_clock_cohort.csv` | Mean/SD acceleration (overall and by sex) |
| `sex_acceleration_models.csv` | Sex coefficients + Bonferroni |
| `acceleration_correlation.csv` | Accel correlations |
| `ethnicity_accel_descriptive.csv` | GSE40279 ethnicity summaries |
| `predicted_vs_age_*.png` | Pred vs age panels |
| `acceleration_by_sex_*.png` | Accel by sex |
| `provenance.json` | pyaging version, clock list, split recipe, missing features before/after, overlap notes, release checks (repo-relative paths) |
| `release_checks/*.json` | `rogen-release-check --table-kind summary` |

Per-sample predictions and `run.log` are written for inspection but gitignored.

## Training-cohort overlap

See `provenance.json` → `training_overlap`. In short:

- Published **hannum**: trained on GSE40279 → held-out metrics are in-sample; GSE87571 is the fair comparison.
- **rogen_hannum_retrain**: fit on the GSE40279 train split only.
- Other listed clocks: no documented GSE40279 / GSE87571 training overlap flagged from pyaging metadata + paper knowledge.

## Observed results

Recorded from the repair rerun (pyaging **0.5.7**; clocks
`horvath2013`, `hannum`, `dnamphenoage`, `dunedinpace`, `grimage2`,
`rogen_hannum_retrain`). Release checks on aggregate tables: **passed**.
Retrained held-out MAE = **3.524864** (n=132).

### Accuracy (MAE / Pearson *r* / mean signed error)

**GSE87571 (n=729)**

| Clock | MAE | *r* | Mean signed error |
|-------|-----|-----|-------------------|
| horvath2013 | 3.90 | 0.971 | +0.33 |
| hannum | 5.13 | 0.975 | +4.01 |
| dnamphenoage | 9.34 | 0.968 | −8.90 |
| dunedinpace | — | 0.593 | — |
| grimage2 | 6.58 | 0.972 | +5.41 |
| rogen_hannum_retrain | 3.76 | 0.987 | +3.10 |

**GSE40279 held-out (n=132)** — published `hannum` is in-sample

| Clock | MAE | *r* | Mean signed error |
|-------|-----|-----|-------------------|
| horvath2013 | 4.92 | 0.933 | −2.21 |
| hannum | 5.55 | 0.950 | +4.59 |
| dnamphenoage | 9.90 | 0.843 | −8.77 |
| dunedinpace | — | 0.009 | — |
| grimage2 | 6.69 | 0.905 | +4.09 |
| rogen_hannum_retrain | 3.52 | 0.957 | −0.05 |

### Sex acceleration (`accel ~ female + age`)

Female−male coefficient; Bonferroni across 6 clocks within cohort.

**GSE87571**

| Clock | Sex coef | 95% CI | Bonferroni *p* |
|-------|----------|--------|----------------|
| horvath2013 | −1.00 | −1.73 to −0.27 | 0.044 |
| hannum | −2.09 | −2.73 to −1.45 | 1.6×10⁻⁹ |
| dnamphenoage | +0.22 | −0.63 to +1.07 | 1 |
| dunedinpace | −0.010 | −0.026 to +0.006 | 1 |
| grimage2 | −0.82 | −1.43 to −0.22 | 0.047 |
| rogen_hannum_retrain | −0.24 | −0.75 to +0.28 | 1 |

**GSE40279 held-out**

| Clock | Sex coef | 95% CI | Bonferroni *p* |
|-------|----------|--------|----------------|
| horvath2013 | −1.93 | −3.64 to −0.21 | 0.167 |
| hannum | −3.20 | −4.65 to −1.75 | 1.6×10⁻⁴ |
| dnamphenoage | −1.48 | −4.06 to +1.10 | 1 |
| dunedinpace | −0.017 | −0.069 to +0.035 | 1 |
| grimage2 | −2.36 | −4.06 to −0.65 | 0.043 |
| rogen_hannum_retrain | −0.73 | −2.26 to +0.80 | 1 |

### Ethnicity (GSE40279 held-out, descriptive)

97 Caucasian − European; 35 Hispanic − Mexican. Mean accel by ethnicity is in
`ethnicity_accel_descriptive.csv` (no formal test).

### Missing features (before vs after pyaging)

Counts from `provenance.json` (`missing_features_before` / `missing_features_after`).
GSE40279 GEO file lacked **88** requested Illumina probes
(`gse40279_geo_probe_report.n_probes_absent_from_geo`); those stay missing until KNN.

| Cohort | Clock | Before | After |
|--------|-------|--------|-------|
| GSE87571 | horvath2013 | 0 | 0 |
| GSE87571 | hannum | 0 | 0 |
| GSE87571 | dnamphenoage | 0 | 0 |
| GSE87571 | dunedinpace | 0 | 0 |
| GSE87571 | grimage2 | 0 | 0 |
| GSE40279_heldout | horvath2013 | 0 | 0 |
| GSE40279_heldout | hannum | 0 | 0 |
| GSE40279_heldout | dnamphenoage | 0 | 0 |
| GSE40279_heldout | dunedinpace | 16 | 16 |
| GSE40279_heldout | grimage2 | 72 | 72 |

`grimage2` after-lists are CpG IDs only (no `female` / `age`).
