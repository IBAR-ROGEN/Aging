# Input Manifest — Methylation Clock Validation

Required inputs for [`evaluate_methylation_clock.py`](evaluate_methylation_clock.py)
(activity 2.1.10.1: GSE40279-trained ElasticNet clock, GSE87571 external validation).

| Path | Role | Required |
|------|------|----------|
| `data/methylation/GSE87571_processed.parquet` | Independent validation beta matrix (`sample_id` + `cg*` columns) | yes |
| `data/methylation/GSE87571_meta.csv` | Phenotype metadata (`sample_id` + `chronological_age`) | yes |
| `models/ro_clock_elasticnet_gse40279.pkl` | Fitted `sklearn.linear_model.ElasticNet` (not a Pipeline / ElasticNetCV) | yes |
| `data/methylation/HM450_probe_annotation.csv` | Probe → nearest gene (`IlmnID`, `UCSC_RefGene_Name`) for panel C labels | no |

## Outputs (written by the evaluation script)

| Path | Description |
|------|-------------|
| `outputs/clock_metrics.json` | MAE, median AE, Pearson *r*, age-stratified MAE |
| `outputs/figures/Figure_Epigenetic_Clock_Panels.pdf` | Vector three-panel figure |
| `outputs/figures/Figure_Epigenetic_Clock_Panels.png` | 300 DPI raster figure |

## Notes

- The model file must unpickle to **exactly** `sklearn.linear_model.ElasticNet`.
  Pipelines and `ElasticNetCV` objects are rejected; do not retrain at evaluation time.
- Optional annotation falls back to Horvath S3 (`test_data/gb-2013-14-10-r115-S3.csv`) when absent.
