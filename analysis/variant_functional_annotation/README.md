# Variant functional annotation

Annotates each unique resolved rsID in `results/snps_validated.csv` with Ensembl VEP
(GRCh38), AlphaMissense (via myvariant.info/dbNSFP), CADD, SIFT/PolyPhen, and optional
GWAS Catalog disease-trait associations. Splits variants into coding vs regulatory classes
so AlphaGenome "no expression change" claims apply only to the non-coding subset.

## Run

From the repo root:

```bash
uv run python analysis/variant_functional_annotation/run_variant_functional_annotation.py
```

Options:

- `--input` — validated SNP CSV (default: `results/snps_validated.csv`)
- `--output-dir` — write CSV and markdown here (default: `results/`)
- `--cache-dir` — raw JSON cache (default: `results/cache/variant_annotation/`)
- `--skip-gwas` — skip GWAS Catalog trait lookups (faster)

## Outputs

| File | Description |
|------|-------------|
| `results/variant_functional_annotation.csv` | One row per unique rsID |
| `results/coding_vs_noncoding_summary.md` | Class counts + AlphaGenome-ready table |
| `results/cache/variant_annotation/` | Cached VEP, myvariant, and GWAS JSON |

Existing VEP responses under `analysis/vep_cache/` are reused on first run when present.

See [docs/GENOMICS_ANALYSIS.md](../../docs/GENOMICS_ANALYSIS.md) for the full pipeline.
