# LA-SNP Public Allele-Frequency Pipeline

**Project:** IBAR-ROGEN Aging  
**Activity:** 2.1.8.1  
**CLI:** `uv run rogen-ukb-manifest`, `uv run rogen-compare-af-gnomad`  
**Scripts:** `scripts/ukb/la_snp_lookup.py`, `scripts/ukb/compare_af_gnomad.py`  
**Package:** `rogen_aging.ukb.manifest`, `rogen_aging.ukb.gnomad`  
**Related:** [Synthetic UKB-RAP Generator](SYNTHETIC_UKB_RAP_GENERATOR.md), [UKB Pre-Commit Hook](UKB_PRE_COMMIT_HOOK.md)

## Overview

Before spending UK Biobank extraction credits, validate the ~70 longevity-associated SNP (LA-SNP) set against **public** reference frequencies:

1. **Build** a GRCh38 manifest from the gene–SNP overlap Excel table (Ensembl REST).
2. **Extract** allele frequencies from indexed **1000 Genomes Project** GRCh38 VCFs (public proxy for expected UKB AFs).
3. **Compare** those AFs to **gnomAD v4** non-Finnish European (`nfe`) frequencies via the public GraphQL API.

This workflow uses **no UK Biobank participant data** and makes **no** DNAnexus, dx-toolkit, or dxFUSE calls.

## Prerequisites

- `overlapping_genes_with_snps.xlsx` at the repo root (or pass `--input`).
- For the extract step: 1000 Genomes GRCh38 VCFs, bgzipped and tabix-indexed (e.g. `ALL.chr*.vcf.gz` under a local `data/` path — git-ignored).
- Network access for the initial Ensembl and gnomAD queries; reruns can be mostly offline via caches.

Install dependencies:

```bash
uv sync
```

(`cyvcf2` is required for indexed VCF region queries.)

## Step 1 — Build manifest (v0.1)

Resolve rs IDs to GRCh38 coordinates and add UKB-oriented chunk labels:

```bash
uv run rogen-ukb-manifest build \
  --input overlapping_genes_with_snps.xlsx \
  --output analysis/ukb_snp_manifest_v0.1.csv
```

Equivalent script path:

```bash
uv run python scripts/ukb/la_snp_lookup.py build \
  --input overlapping_genes_with_snps.xlsx \
  --output analysis/ukb_snp_manifest_v0.1.csv
```

Sanity-check the CSV in **`notebooks/05_ukb_exploration/UKB_LA_SNP_FirstContact.ipynb`**.

## Step 2 — Extract 1KG frequencies (v0.2)

Pull the same loci from indexed 1KG VCFs by coordinate (no full-chromosome streaming):

```bash
uv run rogen-ukb-manifest extract \
  --manifest analysis/ukb_snp_manifest_v0.1.csv \
  --vcf-glob 'data/1kg/ALL.chr*.vcf.gz' \
  --output analysis/la_snp_1kg_frequencies.csv
```

**Output columns:** `rsID`, `chrom`, `pos`, `ref`, `alt`, `AF`, `N_called`.

SNPs missing from 1KG are logged and written with empty `AF` / alleles — the run does not crash.

## Step 3 — Compare to gnomAD v4 (NFE)

Join 1KG AFs to gnomAD v4 (`gnomad_r4`) NFE frequencies:

```bash
uv run rogen-compare-af-gnomad \
  --input analysis/la_snp_1kg_frequencies.csv \
  --output analysis/la_snp_af_1kg_vs_gnomad.csv \
  --scatter figures/af_1kg_vs_gnomad_scatter.png
```

### gnomAD lookup behaviour

- Primary: GRCh38 variant ID `chrom-pos-ref-alt` from the 1KG table.
- Fallback: single-base region query matched by rsID when alleles are missing.
- Population: **`nfe`** (non-Finnish European), aligned with 1KG EUR comparison.
- Batched GraphQL requests with pacing; responses cached at **`data/geo/gnomad_r4_nfe_cache.json`** (git-ignored). Use **`--refresh-cache`** to re-fetch.

### Comparison outputs

| Artifact | Description |
|----------|-------------|
| `analysis/la_snp_af_1kg_vs_gnomad.csv` | Per-SNP table: `rsID`, `AF_1kg`, `AF_gnomad_nfe`, `abs_diff`, `large_diff` (`\|diff\| > 0.05`) |
| `figures/af_1kg_vs_gnomad_scatter.png` | Scatter with identity line; large-diff points highlighted |
| Log / stderr | Lists rsIDs missing from gnomAD or lacking NFE AF |

## Step 4 — Summarize for reporting (optional)

Turn the comparison CSV into a manuscript-ready Markdown paragraph and a top-|ΔAF| table:

```bash
uv run rogen-compare-af-gnomad summarize \
  --input analysis/la_snp_af_1kg_vs_gnomad.csv \
  --output analysis/af_comparison_summary.md
```

Dev helper (same logic):

```bash
uv run python scripts/dev/summarize_af_comparison.py \
  --input analysis/la_snp_af_1kg_vs_gnomad.csv \
  --output analysis/af_comparison_summary.md
```

| Artifact | Description |
|----------|-------------|
| `analysis/af_comparison_summary.md` | Headline concordance stats (paired SNPs, concordant vs discordant at \|ΔAF\| &lt; 0.05) plus top five discordant loci |

**Tests:** `tests/test_af_comparison_summary.py`.

## Typical workflow diagram

```
overlapping_genes_with_snps.xlsx
        │
        ▼  rogen-ukb-manifest build (Ensembl)
analysis/ukb_snp_manifest_v0.1.csv
        │
        ▼  rogen-ukb-manifest extract (1KG VCFs)
analysis/la_snp_1kg_frequencies.csv
        │
        ▼  rogen-compare-af-gnomad (gnomAD GraphQL + cache)
analysis/la_snp_af_1kg_vs_gnomad.csv
figures/af_1kg_vs_gnomad_scatter.png
        │
        ▼  rogen-compare-af-gnomad summarize (optional)
analysis/af_comparison_summary.md
```

After validation, proceed to synthetic UKB-RAP generation or real UKB extraction planning — see [Synthetic UKB-RAP Generator](SYNTHETIC_UKB_RAP_GENERATOR.md).

## Security and compliance

- Output CSVs and caches under `analysis/` and `data/` are git-ignored (`.csv` and `data/` rules).
- `scripts/ukb/la_snp_lookup.py` and `src/rogen_aging/ukb/` are whitelisted in the pre-commit hook for intentional `UKB_` manifest column names; they must never contain participant IDs.
- Do not commit `.vcf` / `.vcf.gz` files; keep 1KG VCFs in `data/`.

## Related documentation

- [WORKFLOWS.md](WORKFLOWS.md) — UKB workflow index
- [UKB Pre-Commit Hook](UKB_PRE_COMMIT_HOOK.md)
- [Synthetic UKB-RAP Generator](SYNTHETIC_UKB_RAP_GENERATOR.md)

---

**Last updated:** May 31, 2026
