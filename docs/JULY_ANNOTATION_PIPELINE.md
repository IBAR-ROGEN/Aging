# July prioritized-variant annotation pipeline

**Script:** `run_july_annotation_pipeline.py` (repo root)  
**Related:** [LA_SNP_VEP_ANNOTATION.md](LA_SNP_VEP_ANNOTATION.md) · [LA_SNP_GTEX_ANNOTATION.md](LA_SNP_GTEX_ANNOTATION.md) · [GENOMICS_ANALYSIS.md](GENOMICS_ANALYSIS.md)

## Overview

Batch functional annotation for the prioritized GRCh38 variant set used in
manuscript supplementary tables. The pipeline:

1. Loads candidate variants from CSV.
2. Queries **GTEx Portal API v2** (`datasetId=gtex_v8`) for significant
   single-tissue cis-eQTLs in brain regions + whole blood (slope/`nes`, p-value).
3. Queries **Ensembl VEP REST** for transcript consequences, SIFT/PolyPhen, and
   HGVS nomenclature (optional local JSONL override).
4. Joins pre-computed **AlphaGenome** and **AlphaMissense** score matrices.
5. Writes a styled three-sheet Excel workbook via `openpyxl`.

API responses are cached under `data/cache/july_annotation/`. Progress uses
`tqdm`; structured logs use `loguru`.

## Inputs

| Path | Role | Required columns |
|------|------|------------------|
| `data/processed/prioritized_variants.csv` | Candidate variants | `chrom`, `pos`, `ref`, `alt`, `rsid`, `gene_symbol` |
| `data/scores/alphagenome_raw.parquet` | AlphaGenome scores | locus and/or `rsid` + `alphagenome_*` (aliases such as `ref_score` accepted) |
| `data/scores/alphamissense_raw.parquet` | AlphaMissense scores | locus and/or `rsid` + `alphamissense_score` / `alphamissense_class` |
| `data/processed/vep_local.jsonl` (optional) | Offline VEP payloads | JSONL keyed by `rsid` / `variant_key` |

## Output

`outputs/Supplementary_Table_1_Annotated_Variants.xlsx`

| Sheet | Contents |
|-------|----------|
| `Combined_Master` | One row per variant: VEP + Alpha scores + GTEx summary |
| `High_Impact_Functional` | AlphaMissense > 0.5 **or** VEP impact HIGH/MODERATE |
| `GTEx_eQTL_Summary` | Long table: SNP–gene–tissue eQTLs (`slope`, `p_value`) |

## CLI usage

```bash
uv sync
uv run python run_july_annotation_pipeline.py

# Warm-cache / offline mode (no live API on cache miss)
uv run python run_july_annotation_pipeline.py --cache-only

# Custom paths
uv run python run_july_annotation_pipeline.py \
  --variants data/processed/prioritized_variants.csv \
  --alphagenome data/scores/alphagenome_raw.parquet \
  --alphamissense data/scores/alphamissense_raw.parquet \
  --output outputs/Supplementary_Table_1_Annotated_Variants.xlsx \
  --verbose
```

## Dependencies

Declared in `pyproject.toml`: `polars`, `requests`, `typer`, `openpyxl`,
`loguru`, `tqdm`, `pyarrow` (parquet).
