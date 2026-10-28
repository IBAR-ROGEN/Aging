# Input Manifest — ROGEN Aging Pipelines

## Activity A.2.1.8.1 — July prioritized-variant functional annotation

Required inputs for [`scripts/ukb/run_july_annotation_pipeline.py`](scripts/ukb/run_july_annotation_pipeline.py)
(production run: 47 prioritized GRCh38 variants → Supplementary Table 1).

| Path | Role | Required |
|------|------|----------|
| `data/processed/variants_47_input.csv` | Original 47-variant allele assignments. Kept for provenance. Not the panel input. | yes |
| `data/processed/variants_47_input_v2.csv` | Panel and `paths.july.variants` input. Eleven corrected alts; `rs5882` remains unresolved. | yes |
| `data/scores/alphagenome_raw.parquet` | Pre-computed AlphaGenome score matrix | yes |
| `data/scores/alphamissense_raw.parquet` | Pre-computed AlphaMissense score matrix | yes |

### Outputs (written by the annotation script)

| Path | Description |
|------|-------------|
| `outputs/Supplementary_Table_1_Annotated_Variants.xlsx` | Three-sheet workbook: Combined_Master, High_Impact_Functional, GTEx_eQTL_Summary. STALE relative to `variants_47_input_v2.csv` (not rebuilt; the July pipeline joins AlphaGenome). Affected rsids: rs28366003, rs1800774, rs41383, rs9916344, rs7207422, rs1801318, rs1981429, rs139170, rs155979, rs524533, rs882696. |

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
| `data/methylation/GSE40279_processed.parquet` | GSE40279 wide betas (`sample_id`, `chronological_age`, `cg*` float32) from `rogen-clock prepare-gse40279` | yes |
| `data/methylation/GSE87571_processed.parquet` | Independent validation beta matrix (`sample_id` + `cg*` columns) | yes |
| `data/methylation/GSE87571_meta.csv` | Phenotype metadata (`sample_id` + `chronological_age`) | yes |
| `models/gse40279_hannum450k_elasticnet.joblib` | GSE40279 ElasticNetCV pipeline (`rogen-clock train`) | yes |
| `models/gse40279_hannum450k_elasticnet.provenance.json` | Training provenance (checksum, probe counts, CV hyperparameters, git commit) | yes |
| `models/gse40279_hannum450k_elasticnet_cpgs.csv` | Non-zero CpGs and weights | yes |
| `models/fixtures/fixture_clock_cg_test.pkl` | Demo clock (`cg_test_*` only). Not a GSE40279 result | no |
| `data/methylation/HM450_probe_annotation.csv` | Probe → nearest gene (`IlmnID`, `UCSC_RefGene_Name`) for panel C labels | no |

## Outputs (written by the evaluation script)

| Path | Description |
|------|-------------|
| `outputs/clock_metrics.json` | MAE, median AE, Pearson *r*, age-stratified MAE |
| `outputs/figures/Figure_Epigenetic_Clock_Panels.pdf` | Vector three-panel figure |
| `outputs/figures/Figure_Epigenetic_Clock_Panels.png` | 300 DPI raster figure |

## Notes

- The trained clock path is `models/gse40279_hannum450k_elasticnet.joblib`.
- `models/fixtures/fixture_clock_cg_test.pkl` is the old demo clock (`cg_test_*`). `rogen-clock evaluate` refuses it unless `--demo` is passed.
- `uv run python scripts/audit_clock_artifacts.py` writes [`docs/CLOCK_ARTIFACT_AUDIT.md`](docs/CLOCK_ARTIFACT_AUDIT.md). Read that report for which serializations are fixtures. Do not copy metrics out of it by hand.
- `uv run python scripts/dev/write_pipeline_fixtures.py` or
  `uv run python scripts/clock/evaluate_methylation_clock.py --demo`
  materializes the demo fixture only.
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

---

## September 2026: secure intake and dry runs

From the provenance files currently in the tree:

| Artefact | What it shows |
|----------|---------------|
| `analysis/panel/ro_longevity_panel_provenance.json` | Panel input `data/processed/variants_47_input_v2.csv` (`n_rows` 47). gnomAD v4 NFE from `data/geo/gnomad_r4_nfe_cache.json` (`cache_hits` 47, `new_queries` 0). Ensembl REST release 116, host `https://rest.ensembl.org`, populations TSI/IBS/CEU/GBR/FIN (`cache_hits` 47, `new_queries` 0). TSI vs NFE count test `n_tests` 47. Both panel CSV release checks `status` `passed`, return code 0. `sanity_check.flag_count` 0. `generated_at_utc` 2026-10-28T07:56:18Z. |
| `analysis/panel/release_checks/ro_longevity_panel.json` | `release_blocked` false for the panel frequency table (summary kind). |
| `analysis/panel/release_checks/ro_longevity_power_grid.json` | `release_blocked` false for the power grid (summary kind). |

Intake (`rogen-intake`) is exercised on synthetic fixtures (`tests/test_intake.py`) and on the public HG002 short-read dry run below. No Romanian sequencing delivery is in the tree.

**HG002 short-read dry run (done, 2026-10-28).** Driver: `scripts/dryrun/hg002_shortread_dryrun.py`. Region `chr20:10000000-20000000`, Illumina BAM downsample `-s 0.1`, minimap2, untuned bcftools, `rtg vcfeval` vs GIAB NISTv4.2.1. Results: [`docs/HG002_DRYRUN.md`](docs/HG002_DRYRUN.md) and `analysis/dryrun/hg002_shortread_metrics.json`. Large intermediates stay under gitignored `data/dryrun/`.

**ONT vs bisulfite HG002 dry run (done, 2026-10-28).** Region `chr20:10000000-20000000`, ONT CRAM via HTTPS + `.crai` (never `all.bam`), `rogen-intake`, `modkit` 0.6.4 `--cpg --combine-strands`. Results: [`docs/METHYLATION_CLOCK_VALIDATION.md`](docs/METHYLATION_CLOCK_VALIDATION.md) and `analysis/dryrun/hg002_ont_methylation_metrics.json`. Large intermediates stay under gitignored `data/dryrun/`.

