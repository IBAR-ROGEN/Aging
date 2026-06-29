# LA-SNP VEP Functional Consequence Annotation

**Project:** IBAR-ROGEN Aging  
**Activity:** 2.1.7.1 (manuscript supplementary table)  
**Script:** `scripts/ukb/annotate_la_snps_vep.py` (deprecated shim: `annotate_la_snps_vep.py` at repo root)  
**Related:** [ALPHAGENOME_ANALYSIS_EXPLANATION.md](ALPHAGENOME_ANALYSIS_EXPLANATION.md), [LA_SNP_PUBLIC_FREQUENCY_PIPELINE.md](LA_SNP_PUBLIC_FREQUENCY_PIPELINE.md)

## Overview

Build the **inferred functional consequences** table for the ~70 longevity-associated SNP (LA-SNP) set using the [Ensembl VEP REST API](https://rest.ensembl.org/) on **GRCh38**. The script:

1. Reads rsIDs from a local file (no VCF, no bulk download).
2. Calls `GET /vep/human/id/{rsid}` once per variant.
3. Caches each raw JSON response under `analysis/vep_cache/` so re-runs skip live queries.
4. Writes a tidy one-row-per-rsID table as CSV and Excel.

This workflow uses **public Ensembl REST endpoints only** — no UK Biobank data.

## Input rsID sources

The script defaults to the curated AlphaGenome LA-SNP table:

| File | Column | Rows | Unique rsIDs | Notes |
|------|--------|------|--------------|-------|
| `analysis/alphagenome/alphagenome_impact_analysis.csv` | `snp` | 70 | 58 | **Default.** Gene–SNP pairs from `scripts/alphagenome/alphagenome_sequence_comparer.py` |
| `overlapping_genes_with_snps.xlsx` | `SNP Identifier` | 225 | 173 (rs-formatted) | Canonical longevity-map overlap table; requires export or column rename |
| Plain `.txt` | one rsID per line | — | — | Set `INPUT_RSID_FILE` to your list |

To switch inputs, edit the **CONSTANTS** block at the top of `scripts/ukb/annotate_la_snps_vep.py`:

```python
INPUT_RSID_FILE = Path("analysis/alphagenome/alphagenome_impact_analysis.csv")
RSID_COLUMN = "snp"
```

For the full Excel overlap table, export rs-formatted IDs first or point `INPUT_RSID_FILE` at a CSV with an `SNP_rsID` / `rsID` column.

## Prerequisites

```bash
uv sync
```

Dependencies used by the script: `requests`, `pandas`, `openpyxl` (Excel export via pandas). Network access is required on the first run; subsequent runs are mostly offline if the cache is warm.

## Run

From the repo root:

```bash
uv run python scripts/ukb/annotate_la_snps_vep.py
```

The script prints progress per rsID and a summary at the end:

```
--- Summary ---
rsIDs queried:   58
rsIDs annotated: 58
rsIDs not found: 0
CSV output:      .../analysis/vep_annotation/la_snp_vep_annotations.csv
Excel output:    .../analysis/vep_annotation/la_snp_vep_annotations.xlsx
VEP cache dir:   .../analysis/vep_cache
```

## Output columns

| Column | Source |
|--------|--------|
| `rsID` | Input |
| `chromosome` | VEP `seq_region_name` |
| `position_GRCh38` | VEP `start` |
| `ref_allele` / `alt_allele` | VEP `allele_string` (first allele = ref; remaining = alt, comma-separated) |
| `most_severe_consequence` | VEP top-level field |
| `gene_symbols` | Unique `gene_symbol` values from `transcript_consequences` (`;`-separated) |
| `SIFT` / `PolyPhen` | From the highest-impact transcript that reports predictions |

Variants with no VEP hit (HTTP 404 or empty payload) are listed in `analysis/vep_annotation/la_snp_vep_not_found.txt`.

## API etiquette and caching

| Setting | Default | Purpose |
|---------|---------|---------|
| `REQUEST_DELAY_SEC` | `0.34` | Minimum pause between live API requests |
| `MAX_RETRIES` | `4` | Retries on HTTP 429 / 503 with exponential backoff |
| `CACHE_DIR` | `analysis/vep_cache/` | One `{rsid}.json` file per variant |

Cached responses are reused automatically. Delete a single cache file to refresh one variant, or remove `analysis/vep_cache/` to re-query all.

## Output paths (git-ignored)

| Path | Contents |
|------|----------|
| `analysis/vep_annotation/la_snp_vep_annotations.csv` | Main supplementary table |
| `analysis/vep_annotation/la_snp_vep_annotations.xlsx` | Same table, Excel format |
| `analysis/vep_annotation/la_snp_vep_not_found.txt` | rsIDs with no VEP result (if any) |
| `analysis/vep_cache/*.json` | Raw Ensembl VEP JSON per rsID |

Add these under `analysis/` in `.gitignore` if you prefer not to commit cache artefacts; the script and this doc remain in version control.

## Manuscript workflow

1. Confirm the rsID input file matches the LA-SNP set used in the main text (~70 pairs / 58 unique rsIDs from AlphaGenome, or your chosen export from `overlapping_genes_with_snps.xlsx`).
2. Run `uv run python scripts/ukb/annotate_la_snps_vep.py`.
3. Import `analysis/vep_annotation/la_snp_vep_annotations.xlsx` into the supplementary **Inferred functional consequences** table.
4. Cross-check coordinates against `uv run rogen-ukb-manifest build` output when the UKB manifest is available.
