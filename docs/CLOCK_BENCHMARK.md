# Published clocks vs retrained Hannum

Benchmark of pyaging DNA-methylation clocks against the ElasticNet Hannum model
retrained on GSE40279 (`models/gse40279_hannum450k_elasticnet.joblib`).

## Cohorts

| Cohort | Samples | Notes |
|--------|---------|--------|
| GSE87571 | 729 with age + sex | Fair external comparison for published Hannum |
| GSE40279 held-out | 132 (`test_size=0.2`, `random_state=42`) | Retrained model never scores train rows; published Hannum is **in-sample** here |

Sex is parsed from GEO `gender:` / `sex:` characteristics. Ethnicity (GSE40279 only)
is reported descriptively; no formal ethnicity test.

## Clocks

Target pyaging clocks: `horvath2013`, `hannum`, `dnamphenoage`, `dunedinpace`.
(`phenoage` in pyaging is clinical Phenotypic Age, not DNAm PhenoAge.)

`grimage2` ran successfully with `female` + chronological `age` metadata
(pyaging 0.5.7).

Retrained model label: `rogen_hannum_retrain`.

## How to run

```bash
uv run python scripts/clock/benchmark_published_clocks.py
```

Inputs (already in-repo / GEO cache):

- `data/geo/GSE87571_series_matrix.txt.gz` + supplementary beta matrices
- `data/methylation/GSE40279_processed.parquet` + `data/geo/GSE40279_series_matrix.txt.gz`
- `models/gse40279_hannum450k_elasticnet.joblib`

CpG columns are restricted to the union of probes required by the selected clocks
plus the retrained Hannum feature list (28 724 Illumina probes in the recorded run).
pyaging uses `df_to_adata(..., imputer_strategy="knn")`.

## Metrics

Per clock × cohort:

- MAE vs chronological age (**omitted for DunedinPACE**)
- Pearson *r*, mean signed error (pred − age)
- Age acceleration = residual of clock ~ chronological age (within cohort)
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
| `provenance.json` | pyaging version, clock list, split recipe, percent_na / missing features, overlap notes, release checks |
| `release_checks/*.json` | `rogen-release-check --table-kind summary` |

Per-sample predictions and `run.log` are written for inspection but gitignored.

## Training-cohort overlap

See `provenance.json` → `training_overlap`. In short:

- Published **hannum**: trained on GSE40279 → held-out metrics are in-sample; GSE87571 is the fair comparison.
- **rogen_hannum_retrain**: fit on the GSE40279 train split only.
- Other listed clocks: no documented GSE40279 / GSE87571 training overlap flagged from pyaging metadata + paper knowledge.

## Observed results

Recorded from the end-to-end run (pyaging **0.5.7**; clocks
`horvath2013`, `hannum`, `dnamphenoage`, `dunedinpace`, `grimage2`,
`rogen_hannum_retrain`). Release checks on aggregate tables: **passed**.

### Accuracy (MAE / Pearson *r* / mean signed error)

**GSE87571 (n=729)**

| Clock | MAE | *r* | Mean signed error |
|-------|-----|-----|-------------------|
| horvath2013 | 3.90 | 0.971 | +0.33 |
| hannum | 5.13 | 0.975 | +4.01 |
| dnamphenoage | 9.34 | 0.968 | −8.90 |
| dunedinpace | — (rate) | 0.593 | −46.42 (pace − age; not age units) |
| grimage2 | 18.82 | 0.754 | +15.20 |
| rogen_hannum_retrain | 3.76 | 0.987 | +3.10 |

**GSE40279 held-out (n=132)** — published `hannum` is in-sample

| Clock | MAE | *r* | Mean signed error |
|-------|-----|-----|-------------------|
| horvath2013 | 4.92 | 0.933 | −2.21 |
| hannum | 6.11 | 0.952 | +5.40 |
| dnamphenoage | 9.90 | 0.843 | −8.77 |
| dunedinpace | — (rate) | 0.009 | −63.58 (pace − age; not age units) |
| grimage2 | 13.41 | 0.206 | +3.83 |
| rogen_hannum_retrain | 3.52 | 0.957 | −0.05 |

### Sex acceleration (`accel ~ female + age`)

Female−male coefficient; Bonferroni across 6 clocks within cohort.

**GSE87571:** published `hannum` sex coef −2.09 (95% CI −2.73 to −1.45; Bonferroni *p* = 1.7×10⁻⁹);
`horvath2013` −1.00 (−1.73 to −0.27; Bonferroni *p* = 0.044). Other clocks NS after Bonferroni
(`rogen_hannum_retrain` −0.24, *p*_Bonf = 1).

**GSE40279 held-out:** published `hannum` −3.13 (−4.54 to −1.72; Bonferroni *p* = 1.4×10⁻⁴).
Others NS after Bonferroni.

### Ethnicity (GSE40279 held-out, descriptive)

97 Caucasian − European; 35 Hispanic − Mexican. Mean accel by ethnicity is in
`ethnicity_accel_descriptive.csv` (no formal test).

### Missing features (pyaging `*_percent_na`)

GSE87571: GrimAge2 0.19% NA (2 missing features); others 0%.
GSE40279 held-out: hannum 2.8% (2), dunedinpace 0.08% (16), grimage2 7.2% (74).
