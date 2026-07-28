# Input Manifest — ROGEN Aging Pipelines

## Activity A.2.1.8.1 — July prioritized-variant functional annotation

Required inputs for [`scripts/ukb/run_july_annotation_pipeline.py`](scripts/ukb/run_july_annotation_pipeline.py)
(production run: 47 prioritized GRCh38 variants → Supplementary Table 1).

| Path | Role | Required |
|------|------|----------|
| `data/processed/variants_47_input.csv` | 47 prioritized variants (`chrom`, `pos`, `ref`, `alt`, `rsid`, `gene_symbol`) | yes |
| `data/scores/alphagenome_raw.parquet` | Pre-computed AlphaGenome score matrix | yes |
| `data/scores/alphamissense_raw.parquet` | Pre-computed AlphaMissense score matrix | yes |

### Outputs (written by the annotation script)

| Path | Description |
|------|-------------|
| `outputs/Supplementary_Table_1_Annotated_Variants.xlsx` | Three-sheet workbook: Combined_Master, High_Impact_Functional, GTEx_eQTL_Summary |

### Notes

- Default CLI paths target **production** (`outputs/`), not `outputs/demo/`.
- Pass `--demo` only for offline fixture smoke tests (writes under `outputs/demo/`).
- GTEx Portal API uses `datasetId=gtex_v8`; Ensembl VEP uses GRCh38 REST (`hgvs=1`).

---

## Methylation Clock Validation

Required inputs for [`scripts/clock/evaluate_methylation_clock.py`](scripts/clock/evaluate_methylation_clock.py)
(activity 2.1.10.1: GSE40279-trained ElasticNet clock, GSE87571 external validation).

| Path | Role | Required |
|------|------|----------|
| `data/methylation/GSE87571_processed.parquet` | Independent validation beta matrix (`sample_id` + `cg*` columns) | yes |
| `data/methylation/GSE87571_meta.csv` | Phenotype metadata (`sample_id` + `chronological_age`) | yes |
| `models/ro_clock_elasticnet_gse40279.pkl` | Fitted ElasticNet **or** Pipeline ending in ElasticNet/ElasticNetCV | yes |
| `data/methylation/HM450_probe_annotation.csv` | Probe → nearest gene (`IlmnID`, `UCSC_RefGene_Name`) for panel C labels | no |

## Outputs (written by the evaluation script)

| Path | Description |
|------|-------------|
| `outputs/clock_metrics.json` | MAE, median AE, Pearson *r*, age-stratified MAE |
| `outputs/figures/Figure_Epigenetic_Clock_Panels.pdf` | Vector three-panel figure |
| `outputs/figures/Figure_Epigenetic_Clock_Panels.png` | 300 DPI raster figure |

## Notes

- Preferred model path is `models/ro_clock_elasticnet_gse40279.pkl`.
- If that pickle is missing, the evaluator accepts
  `models/methylation_clock_v1.joblib` (Pipeline) and/or can materialize the
  pickle via `uv run python scripts/dev/write_pipeline_fixtures.py` or
  `uv run python scripts/clock/evaluate_methylation_clock.py --demo`.
- Optional annotation falls back to Horvath S3 (`test_data/gb-2013-14-10-r115-S3.csv`) when absent.

---

## Activity A.2.1.11.1 — Integrative multi-omics (variant × tissue × phenotype)

Required upstream artefact for [`scripts/integrative/run_pipeline.py`](scripts/integrative/run_pipeline.py)
(production defaults; not `--demo`).

| Path | Role | Required |
|------|------|----------|
| `outputs/Supplementary_Table_1_Annotated_Variants.xlsx` | July workbook (`Combined_Master` + `GTEx_eQTL_Summary`) | yes (or parquet siblings below) |
| `outputs/Supplementary_Table_1_Combined_Master.parquet` | Parsed Combined_Master (auto-written from Excel if missing) | preferred |
| `outputs/Supplementary_Table_1_GTEx_eQTL_Summary.parquet` | Parsed long eQTL table (`slope`→`nes`) | preferred |

### Outputs

| Path | Description |
|------|-------------|
| `analysis/integrative/results/annotated_variants.parquet` | Variants remapped to tissue eQTL summaries |
| `analysis/integrative/results/eqtl_summary.parquet` | Per-rsID GTEx summary |
| `analysis/integrative/results/variant_risks.parquet` | Channel scores + `composite_risk` |
