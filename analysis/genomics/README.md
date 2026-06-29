# Genomics analysis scripts (GRCh38/hg38)

Deterministic pipelines for manuscript table validation, cluster–LongevityMap
overlap enrichment, and variant functional annotation.

**Full guide:** [docs/GENOMICS_ANALYSIS.md](../../docs/GENOMICS_ANALYSIS.md)

| Module | Script | Purpose |
|--------|--------|---------|
| [validate_genomics_tables](../validate_genomics_tables/) | `validate_genomics_tables.py` | Recompute counts from `overlapping_genes_with_snps.xlsx`; resolve legacy aliases |
| [overlap_enrichment](../overlap_enrichment/) | `run_overlap_enrichment.py` | Fisher/hypergeometric enrichment vs LongevityMap |
| [variant_functional_annotation](../variant_functional_annotation/) | `run_variant_functional_annotation.py` | VEP + AlphaMissense + GWAS; coding vs regulatory split |

Outputs are written to `results/` (git-ignored). API JSON caches live under
`results/cache/`.
