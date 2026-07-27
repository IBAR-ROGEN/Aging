# Genomics table validation (GRCh38/hg38)

Independent recomputation of manuscript counts and coordinates from
`overlapping_genes_with_snps.xlsx` at the repository root.

**Genome build:** GRCh38/hg38 only. All coordinate checks use Ensembl Variation
REST (primary) with myvariant.info and NCBI dbSNP as fallbacks for legacy alias
resolution.

## Inputs

| File | Location | Description |
|------|----------|-------------|
| `overlapping_genes_with_snps.xlsx` | repo root | Gene–LA-SNP overlap table (225 rows) |
| `Supplementary Table 3.xlsx` | optional `data/` | Cluster gene lists (AD/PD up/down) |

## Outputs

| File | Description |
|------|-------------|
| `results/validation_report.md` | Stated vs computed table, flagged issues, verdicts |
| `results/snps_validated.csv` | De-duplicated rows with canonical rsID column |
| `results/cache/*.json` | Raw API responses (auditable) |

## Run

From the repository root:

```bash
uv run python analysis/validate_genomics_tables/validate_genomics_tables.py \
  --input overlapping_genes_with_snps.xlsx \
  --output-dir results
```

Optional cluster check (when `data/Supplementary Table 3.xlsx` is available):

```bash
uv run python analysis/validate_genomics_tables/validate_genomics_tables.py \
  --input overlapping_genes_with_snps.xlsx \
  --cluster-table data/Supplementary\ Table\ 3.xlsx \
  --output-dir results
```

## Rerun

1. Ensure dependencies: `uv sync` (or `uv pip install -r analysis/validate_genomics_tables/requirements.txt`).
2. Delete `results/cache/` to force fresh API calls, or delete individual JSON files to refresh one variant.
3. Re-run the command above. The script is deterministic given identical inputs and cached API responses.

See [docs/GENOMICS_ANALYSIS.md](../../docs/GENOMICS_ANALYSIS.md) for the full pipeline.
