# Clock sex check (chrX / chrY sensitivity)

Non-zero coefficients in `models/gse40279_hannum450k_elasticnet.joblib` mapped with the Zhou HM450 hg38 manifest (`data/dryrun/reference/HM450.hg38.manifest.tsv.gz`).

| Item | Value |
|------|-------|
| Non-zero coefficients | 2314 |
| Non-zero on chrX | 0 |
| Non-zero on chrY | 0 |
| Unmapped non-zero | 0 |
| Sum \|coef\| on chrX+chrY | 0 |
| Sum \|coef\| all non-zero | 1352.22 |
| chrX+chrY share of sum \|coef\| | 0 |

The noXY sensitivity used the same GSE40279 data, split (`test_size=0.2`, `random_state=42`), missingness and age-correlation filters, probe cap (8000), and ElasticNetCV settings as the main provenance, except chrX and chrY probes are removed after the missingness filter and before the age-correlation filter. Removing chrX and chrY changed two candidate probes, neither in the top 8,000, so the retrained model was byte-identical and is not stored separately. Provenance only: `models/gse40279_hannum450k_elasticnet_noXY.provenance.json`.

## Accuracy and sex model

| Model | Cohort | n | MAE | r | Mean signed error | Sex coef (F−M) | Sex p |
|-------|--------|---|-----|---|-------------------|----------------|-------|
| main | GSE40279_heldout | 132 | 3.52486 | 0.956678 | -0.0548695 | -0.729104 | 0.347549 |
| main | GSE87571 | 729 | 3.76019 | 0.986745 | 3.10323 | -0.235038 | 0.367368 |
| noXY | GSE40279_heldout | 132 | 3.52486 | 0.956678 | -0.0548695 | -0.729104 | 0.347549 |
| noXY | GSE87571 | 729 | 3.76019 | 0.986745 | 3.10323 | -0.235038 | 0.367368 |

Neither the main model nor the noXY model shows a sex difference (accel ~ female + age, p < 0.05) on GSE40279 held-out or GSE87571.
