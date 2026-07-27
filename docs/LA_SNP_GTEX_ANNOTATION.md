# LA-SNP GTEx eQTL Annotation

**Project:** IBAR-ROGEN Aging  
**Activity:** 2.1.7.1 (manuscript supplementary table)  
**Script:** `scripts/ukb/annotate_la_snps_gtex.py` (deprecated shim: `annotate_la_snps_gtex.py` at repo root)  
**Related:** [LA_SNP_VEP_ANNOTATION.md](LA_SNP_VEP_ANNOTATION.md), [ALPHAGENOME_ANALYSIS_EXPLANATION.md](ALPHAGENOME_ANALYSIS_EXPLANATION.md)

## Overview

Annotate the ~70 longevity-associated SNP (LA-SNP) set with **significant single-tissue cis-eQTL** evidence from the [GTEx Portal API v2](https://gtexportal.org/api/v2/redoc). The script:

1. Merges the curated AlphaGenome LA-SNP table with GRCh38 coordinates from the VEP annotation output.
2. Resolves each rsID to a GTEx `variantId` via `GET /api/v2/dataset/variant`.
3. Queries `GET /api/v2/association/singleTissueEqtl` for **brain** and **whole-blood** tissues (neurodegeneration angle).
4. Caches raw JSON under `analysis/gtex_cache/` so re-runs skip live queries.
5. Writes a long-format table (one row per SNP–gene–tissue eQTL) as CSV and Excel.

This workflow uses **GTEx Portal API endpoints only** — no bulk download and no UK Biobank data.

## Input files

The script merges two local tables at runtime and validates that all 58 unique rsIDs align:

| File | Column | Rows | Unique rsIDs | Role |
|------|--------|------|--------------|------|
| `analysis/alphagenome/alphagenome_impact_analysis.csv` | `snp` | 70 | 58 | Curated LA-SNP list (AlphaGenome impact table) |
| `analysis/vep_annotation/la_snp_vep_annotations.xlsx` | `rsID` | 58 | 58 | GRCh38 coordinates from Ensembl VEP ([LA_SNP_VEP_ANNOTATION.md](LA_SNP_VEP_ANNOTATION.md)) |

Run the VEP annotation step first if `la_snp_vep_annotations.xlsx` is missing:

```bash
uv run python scripts/ukb/annotate_la_snps_vep.py
```

To switch inputs or tissues, edit the **CONSTANTS** block at the top of `scripts/ukb/annotate_la_snps_gtex.py`.

## Target tissues

Default `TARGET_TISSUES` (13 brain regions + whole blood):

- `Brain_Amygdala`, `Brain_Anterior_cingulate_cortex_BA24`, `Brain_Caudate_basal_ganglia`, `Brain_Cerebellar_Hemisphere`, `Brain_Cerebellum`, `Brain_Cortex`, `Brain_Frontal_Cortex_BA9`, `Brain_Hippocampus`, `Brain_Hypothalamus`, `Brain_Nucleus_accumbens_basal_ganglia`, `Brain_Putamen_basal_ganglia`, `Brain_Spinal_cord_cervical_c-1`, `Brain_Substantia_nigra`, `Whole_Blood`

## GTEx API endpoints

| Step | Endpoint | Key query params | Fields parsed |
|------|----------|------------------|---------------|
| Variant lookup | `GET /api/v2/dataset/variant` | `snpId`, `datasetId` (`gtex_v10`) | `variantId`, `snpId`, `chromosome`, `pos` |
| Location fallback | same | `chromosome`, `pos`, `datasetId` | same (when rsID lookup returns no hit) |
| Single-tissue eQTLs | `GET /api/v2/association/singleTissueEqtl` | `variantId`, repeated `tissueSiteDetailId`, `datasetId`, `page` | `snpId`, `variantId`, `geneSymbol`, `tissueSiteDetailId`, `nes`, `pValue` |

**Note:** The eQTL endpoint requires a GTEx `variantId` (e.g. `chr2_206129125_A_G_b38`); passing an rsID directly as `variantId` returns empty results.

## Prerequisites

```bash
uv sync
```

Dependencies: `requests`, `pandas`, `openpyxl` (Excel export). Network access is required on the first run; subsequent runs are mostly offline if the cache is warm.

## Run

From the repo root:

```bash
uv run python scripts/ukb/annotate_la_snps_gtex.py
```

Example summary (values vary by GTEx release):

```
--- Summary ---
LA-SNPs queried:              58
Variants resolved:            58
SNPs with >=1 target eQTL:    22
Total eQTL hits (long table): 303
Unresolved SNPs:              0
CSV output:                   .../analysis/gtex_annotation/la_snp_gtex_eqtls.csv
Excel output:                 .../analysis/gtex_annotation/la_snp_gtex_eqtls.xlsx
GTEx cache dir:               .../analysis/gtex_cache
```

## Output columns

| Column | Source |
|--------|--------|
| `rsID` | Input / GTEx `snpId` |
| `gtex_variant_id` | GTEx `variantId` |
| `gene_symbol` | GTEx `geneSymbol` |
| `tissue` | GTEx `tissueSiteDetailId` |
| `nes` | GTEx normalized effect size |
| `p_value` | GTEx nominal p-value |

SNPs that cannot be resolved to a GTEx variant are listed in `analysis/gtex_annotation/la_snp_gtex_unresolved.txt`.

## API etiquette and caching

| Setting | Default | Purpose |
|---------|---------|---------|
| `REQUEST_DELAY_SEC` | `0.5` | Minimum pause between live API requests |
| `MAX_RETRIES` | `4` | Retries on HTTP 429 / 503 with exponential backoff |
| `DATASET_ID` | `gtex_v10` | GTEx release queried |
| `CACHE_DIR` | `analysis/gtex_cache/` | Raw JSON per variant lookup, eQTL page, and consolidated eQTL result |

Cached responses are reused automatically. Delete individual cache files to refresh specific queries, or remove `analysis/gtex_cache/` to re-query all.

## Output paths (git-ignored)

| Path | Contents |
|------|----------|
| `analysis/gtex_annotation/la_snp_gtex_eqtls.csv` | Long-format eQTL table |
| `analysis/gtex_annotation/la_snp_gtex_eqtls.xlsx` | Same table, Excel format |
| `analysis/gtex_annotation/la_snp_gtex_unresolved.txt` | rsIDs with no GTEx variant match (if any) |
| `analysis/gtex_cache/*.json` | Raw GTEx API JSON |

## Manuscript workflow

1. Confirm `alphagenome_impact_analysis.csv` and `la_snp_vep_annotations.xlsx` cover the same 58 unique LA-SNPs (the script validates this at startup).
2. Run `uv run python scripts/ukb/annotate_la_snps_gtex.py`.
3. Import `analysis/gtex_annotation/la_snp_gtex_eqtls.xlsx` into the supplementary **GTEx eQTL evidence** table.
4. Cross-reference `gene_symbol` hits with the longevity gene list and VEP `gene_symbols` column where relevant.
