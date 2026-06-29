# LA-SNP allele-frequency comparison figure (`plot_af_comparison.py`)

**Project:** IBAR-ROGEN Aging  
**Activity:** 2.1.8.1 — LA-SNP manifest + public AF validation  
**Script:** [`scripts/figures/plot_af_comparison.py`](../scripts/figures/plot_af_comparison.py) (deprecated shim: [`plot_af_comparison.py`](../plot_af_comparison.py) at repo root)  
**Related:** [LA_SNP_PUBLIC_FREQUENCY_PIPELINE.md](LA_SNP_PUBLIC_FREQUENCY_PIPELINE.md) · [ACTIVITIES.md](ACTIVITIES.md#2181--la-snp-manifest--public-af-validation)

## Purpose

`plot_af_comparison.py` produces a **publication-style two-panel figure** from the output of `rogen-compare-af-gnomad`:

| Panel | Content |
|-------|---------|
| **A** | Scatter of 1000 Genomes EUR AF (`AF_1kg`) vs gnomAD v4 NFE AF (`AF_gnomad_nfe`), y = x reference, concordant points in grey and flagged discrepancies in red; top loci labelled as `rsID (Gene)` |
| **B** | Horizontal bar chart of the top N loci by \|ΔAF\|, y-axis labelled by rsID, threshold line at 0.05 |

Gene names are **merged from the LA-SNP manifest** (`analysis/ukb_snp_manifest_v0.1.csv`, columns `Gene` + `SNP_rsID`). When one rsID maps to multiple genes, labels use comma-separated gene symbols.

Stdout reports a **coverage summary**: total SNPs, SNPs with AF in both sources, and SNPs flagged at the discrepancy threshold.

This complements the single-panel scatter written by **`rogen-compare-af-gnomad --scatter`**, which is useful for quick QA but lacks ranked discrepancies, gene labels, and PDF export.

## Prerequisites

Complete the public AF pipeline through Step 3 first:

```bash
uv sync

uv run rogen-ukb-manifest build \
  --input overlapping_genes_with_snps.xlsx \
  --output analysis/ukb_snp_manifest_v0.1.csv

uv run rogen-ukb-manifest extract \
  --manifest analysis/ukb_snp_manifest_v0.1.csv \
  --vcf-glob 'data/1kg/ALL.chr*.vcf.gz' \
  --output analysis/la_snp_1kg_frequencies.csv

uv run rogen-compare-af-gnomad \
  --input analysis/la_snp_1kg_frequencies.csv \
  --output analysis/la_snp_af_1kg_vs_gnomad.csv \
  --scatter figures/af_1kg_vs_gnomad_scatter.png
```

See [LA_SNP_PUBLIC_FREQUENCY_PIPELINE.md](LA_SNP_PUBLIC_FREQUENCY_PIPELINE.md) for 1KG VCF layout, gnomAD cache behaviour, and optional Markdown summary (Step 4).

## Run the figure script

```bash
uv run python scripts/figures/plot_af_comparison.py
```

Example stdout:

```text
Allele-frequency comparison coverage
  Total SNPs in prioritized set: 70
  SNPs with usable AF in both sources: 68
  SNPs flagged (|ΔAF| > 0.05): 4
PNG: …/figures/af_1kg_vs_gnomad_comparison.png
PDF: …/figures/af_1kg_vs_gnomad_comparison.pdf
```

Default outputs:

| File | Description |
|------|-------------|
| `figures/af_1kg_vs_gnomad_comparison.png` | Raster figure (300 dpi) |
| `figures/af_1kg_vs_gnomad_comparison.pdf` | Vector figure |

## Configuration

All paths and styling are **constants at the top of** `plot_af_comparison.py` (no CLI flags). Edit these before running:

| Constant | Default | Role |
|----------|---------|------|
| `INPUT_CSV` | `analysis/la_snp_af_1kg_vs_gnomad.csv` | Comparison table from `rogen-compare-af-gnomad` |
| `MANIFEST_CSV` | `analysis/ukb_snp_manifest_v0.1.csv` | Gene names merged onto comparison rows |
| `OUTPUT_DIR` | `figures/` | Directory for PNG/PDF |
| `FIG_BASENAME` | `af_1kg_vs_gnomad_comparison` | Output filename stem |
| `DIFF_THRESHOLD` | `0.05` | Flag threshold for `large_diff` and panel styling (recomputed if column absent) |
| `TOP_N_LABELS` | `12` | Top \|ΔAF\| loci labelled in panel A and shown in panel B |
| `FIGURE_DPI` | `300` | PNG resolution |
| `FONT_SIZE` | `11` | Base matplotlib font size |

## Input table columns

The comparison CSV (**required**):

| Column | Meaning |
|--------|---------|
| `rsID` | dbSNP identifier |
| `AF_1kg` | 1000 Genomes EUR allele frequency |
| `AF_gnomad_nfe` | gnomAD v4 non-Finnish European AF |

**Optional** (recomputed when missing):

| Column | Meaning |
|--------|---------|
| `abs_diff` | \|AF_1kg − AF_gnomad_nfe\| |
| `large_diff` | Boolean flag; `True` when `abs_diff > DIFF_THRESHOLD` (default 0.05) |

The manifest CSV (**required for gene labels**):

| Column | Meaning |
|--------|---------|
| `SNP_rsID` | Join key (normalized to `rs…` form) |
| `Gene` | Gene symbol for scatter annotations |

## Gene merge behaviour

1. Load `MANIFEST_CSV` and build an rsID → gene map (multiple genes per SNP joined with `", "`).
2. Left-join onto the comparison table on `rsID`.
3. If the comparison CSV already contains a `Gene` column, manifest values fill only missing entries.

Panel A labels use `rs429358 (APOE)` when a gene is known; panel B y-axis shows rsID only.

## Dependencies

**pandas**, **numpy**, and **matplotlib** only (no scipy, no network access at plot time).

## Distinction from other AF outputs

| Tool | Output |
|------|--------|
| **`rogen-compare-af-gnomad`** | `analysis/la_snp_af_1kg_vs_gnomad.csv` + optional single-panel `--scatter` PNG |
| **`rogen-compare-af-gnomad summarize`** | `analysis/af_comparison_summary.md` with headline stats and top-five discordant table (generate first; not shipped in repo) |
| **`plot_af_comparison.py`** | Two-panel PNG/PDF with gene-labelled scatter and ranked \|ΔAF\| bars |

## Related documentation

- [LA_SNP_PUBLIC_FREQUENCY_PIPELINE.md](LA_SNP_PUBLIC_FREQUENCY_PIPELINE.md) — full manifest → 1KG → gnomAD workflow  
- [FIGURES.md](FIGURES.md) — manuscript figure script index  
- [ACTIVITIES.md](ACTIVITIES.md#2181--la-snp-manifest--public-af-validation) — activity 2.1.8.1 map
