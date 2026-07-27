# Cluster ∩ LongevityMap overlap enrichment (GRCh38/hg38)

Tests whether AD/PD cluster genes overlap LongevityMap significant longevity genes
more than expected by chance.

## Required inputs

| File | Location | Description |
|------|----------|-------------|
| `Supplementary Table 3.xlsx` | `data/` | Four sheets: AD/PD up/down cluster ENTREZ lists |
| `longevitymap.sqlite` | `data/` | LongevityMap DB (built automatically if missing) |
| `snps_validated.csv` | `results/` (fallback: `analysis/results/`) | Cross-check for expected overlap of 41 genes |
| `ad_pd_gpl_ids.txt` | `data/` (optional) | One GPL accession per line for universe (b) |

LongevityMap raw CSV source (if sqlite missing): [HAGR LongevityMap download](https://www.genomics.senescence.info/longevity/longevity_genes.zip)

## Outputs

| File | Description |
|------|-------------|
| `analysis/results/overlap_enrichment.md` | Contingency tables, statistics, interpretation |
| `results/overlap_enrichment.md` | Same (default `--output-dir results`) |
| `analysis/results/overlap_enrichment_stats.csv` | Machine-readable test results |
| `analysis/results/cache/` | Cached API/GEO responses |

## Run

From repository root:

```bash
uv run python analysis/overlap_enrichment/run_overlap_enrichment.py \
  --cluster-table "data/Supplementary Table 3.xlsx" \
  --longevity-db data/longevitymap.sqlite \
  --snps-validated results/snps_validated.csv \
  --output-dir results
```

Build LongevityMap sqlite only:

```bash
uv run python analysis/overlap_enrichment/run_overlap_enrichment.py build-sqlite \
  --longevity-db data/longevitymap.sqlite
```

## Universe definitions

1. **(a) All protein-coding genes** — GENCODE v46 primary assembly gene symbols (cached).
2. **(b) Microarray platform union** — genes annotated on all GPL platforms listed in `data/ad_pd_gpl_ids.txt` (honest background when platforms are documented).
3. **(c) Meta-analysis DE-tested genes** — union of all four cluster lists (primary honest background when platform file absent).

Primary interpretation uses **(b)** when GPL IDs are available, otherwise **(c)**.

See [docs/GENOMICS_ANALYSIS.md](../../docs/GENOMICS_ANALYSIS.md) for the full pipeline.
