# Genomics analysis pipeline (GRCh38/hg38)

Independent validation, enrichment testing, and functional annotation for the
longevity-associated SNP (LA-SNP) overlap table used in the IBAR-ROGEN aging
manuscript.

**Genome build:** GRCh38/hg38 only. All API lookups, coordinates, and reports
use this assembly.

**Principles:**

- Never invent rsIDs, coordinates, alleles, scores, or p-values — use input
  files or authoritative APIs (NCBI dbSNP, Ensembl VEP, myvariant.info, GWAS
  Catalog).
- Cache every raw API response under `results/cache/` as JSON.
- When computed values disagree with manuscript-stated values, record both and
  propose corrections — do not silently overwrite.

## Pipeline overview

```text
overlapping_genes_with_snps.xlsx
        │
        ▼
┌───────────────────────────────────┐
│ 1. validate_genomics_tables       │  → validation_report.md, snps_validated.csv
└───────────────────────────────────┘
        │
        ├──────────────────────────────────────┐
        ▼                                      ▼
┌───────────────────────────┐    ┌────────────────────────────────────┐
│ 2. overlap_enrichment     │    │ 3. variant_functional_annotation   │
│    + Supplementary T3     │    │    VEP + AlphaMissense + GWAS      │
│    + longevitymap.sqlite  │    │    coding vs regulatory split      │
└───────────────────────────┘    └────────────────────────────────────┘
        │                                      │
        ▼                                      ▼
 overlap_enrichment.md              variant_functional_annotation.csv
                                    coding_vs_noncoding_summary.md
```

Related prior work (same LA-SNP set):

| Component | Script | Doc |
|-----------|--------|-----|
| AlphaGenome regulatory predictions | `scripts/alphagenome/` | [ALPHAGENOME_ANALYSIS_EXPLANATION.md](ALPHAGENOME_ANALYSIS_EXPLANATION.md) |
| VEP consequence table | `scripts/ukb/annotate_la_snps_vep.py` | [LA_SNP_VEP_ANNOTATION.md](LA_SNP_VEP_ANNOTATION.md) |
| GTEx eQTL lookup | `scripts/ukb/annotate_la_snps_gtex.py` | [LA_SNP_GTEX_ANNOTATION.md](LA_SNP_GTEX_ANNOTATION.md) |

## Required inputs

| File | Location | Used by |
|------|----------|---------|
| `overlapping_genes_with_snps.xlsx` | repo root | Step 1 |
| `Supplementary Table 3.xlsx` | `data/` | Step 2 (AD/PD cluster ENTREZ lists) |
| `longevitymap.sqlite` | `data/` | Step 2 (built from HAGR CSV if missing) |
| `snps_validated.csv` | `results/` | Steps 2–3 (output of step 1) |

Optional: `data/ad_pd_gpl_ids.txt` — one GPL accession per line for microarray
platform universe (b) in enrichment analysis.

## Run all three steps

From the repository root:

```bash
# 1. Recompute table counts and resolve legacy SNP aliases
uv run python analysis/validate_genomics_tables/validate_genomics_tables.py \
  --input overlapping_genes_with_snps.xlsx \
  --output-dir results

# 2. Test cluster ∩ LongevityMap overlap enrichment
uv run python analysis/overlap_enrichment/run_overlap_enrichment.py \
  --cluster-table "data/Supplementary Table 3.xlsx" \
  --longevity-db data/longevitymap.sqlite \
  --snps-validated results/snps_validated.csv \
  --output-dir results

# 3. Annotate variants: coding vs non-coding, AlphaMissense, GWAS
uv run python analysis/variant_functional_annotation/run_variant_functional_annotation.py \
  --input results/snps_validated.csv \
  --output-dir results
```

Each module has its own `requirements.txt` (pinned) and `README.md` under
`analysis/<module>/`.

## Outputs (`results/` — git-ignored except `.gitkeep`)

| File | Step | Description |
|------|------|-------------|
| `validation_report.md` | 1 | Stated vs computed metrics, flagged issues |
| `snps_validated.csv` | 1 | Canonical rsIDs, resolution provenance, coordinates |
| `overlap_enrichment.md` | 2 | Contingency tables, Fisher/hypergeometric tests |
| `overlap_enrichment_stats.csv` | 2 | Machine-readable enrichment results |
| `variant_functional_annotation.csv` | 3 | One row per unique rsID with VEP, AlphaMissense, GWAS |
| `coding_vs_noncoding_summary.md` | 3 | Class counts + AlphaGenome-ready variant table |
| `cache/` | 1–3 | Raw JSON API responses (auditable) |

## Key findings (current run)

### Table validation

- **41 significant unique genes** — matches manuscript.
- **48 unique canonical rsIDs** after GRCh38 resolution (manuscript cites 58 from
  AlphaGenome curated table; Excel has duplicate rows, legacy aliases, and
  gene-name placeholders).
- Flagged: CETP I405V (`rs1273184461`) ≠ rs5882 (different GRCh38 loci); unresolved
  HSPA1A/B/L legacy names; HLA and repeat entries excluded.

### Overlap enrichment

- |A ∩ B| = **41** — matches manuscript and `snps_validated.csv`.
- Under **honest background** (union of DE-tested cluster genes, universe c):
  overlap is **not significant** after BH correction.
- Significance appears only with all protein-coding background (~20k genes) —
  report primary result using universe (c).

### Variant functional annotation

| Class | Count |
|-------|------:|
| coding-missense | 7 |
| coding-synonymous | 2 |
| non-coding (UTR + intronic + regulatory) | 39 |

AlphaGenome “no expression change” is interpretable for the **39 non-coding**
variants only. **Coding missense** variants (7) use **AlphaMissense** for protein-level
claims.

## Module reference

| Path | Script |
|------|--------|
| `analysis/validate_genomics_tables/` | `validate_genomics_tables.py` |
| `analysis/overlap_enrichment/` | `run_overlap_enrichment.py` |
| `analysis/variant_functional_annotation/` | `run_variant_functional_annotation.py` |

See also [ACTIVITIES.md](ACTIVITIES.md) (activity IDs **2.1.7.2**–**2.1.7.4**).
