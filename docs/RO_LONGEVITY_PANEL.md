# Illustrative detectability of the 47 prioritized variants

This note ranks the 47 prioritized GRCh38 variants by how readily an allelic case-control test could detect them **if** a Southeast European cohort were the size and odds ratio in the grid below. It is a power calculation on public reference frequencies. It is not a test of longevity in Romanians, and it does not use ROGEN genotypes.

Cohort sizes and odds ratios are illustrative. ROGEN sample size and design are not fixed. Frequencies are from public reference populations, not from Romanian data.

## Inputs

| File | Role |
|------|------|
| `data/processed/variants_47_input.csv` | The 47 prioritized variants (`chrom`, `pos`, `ref`, `alt`, `rsid`, `gene_symbol`). The script stops if this file is missing or does not have 47 data rows. It does not rebuild the list. |
| `manuscript/tables/41_gene_candidate_list.csv` | Gene-level `Longevity_Class` (Pro-Longevity, Anti-Longevity, Context-Dependent). Joined on `gene_symbol` for annotation only. |

`Longevity_Class` in that table is the LongevityMap class used in the manuscript. This script does not query LongevityMap.

`analysis/rogen_vs_gnomad_af.csv` is not an input. Its ROGEN allele-frequency column is synthetic.

## Frequencies

Two public sources are stored side by side. They are not averaged, and neither is called a Romanian frequency.

1. **gnomAD v4 non-Finnish European (NFE)** allele frequency of the alternate allele, from the existing cache `data/geo/gnomad_r4_nfe_cache.json`. The panel script reads that file and does not query gnomAD. An rsID missing from the cache keeps an empty NFE frequency and a note.
2. **1000 Genomes phase 3** EUR subpopulations **TSI, IBS, CEU, GBR, FIN**, from the Ensembl REST variation endpoint (GRCh38) through `src/rogen_aging/ensembl/`. The request includes `pops=1`, because the variation record otherwise omits population frequencies. The host is `config/default.yaml` key `apis.ensembl_rest_dated`, currently `https://rest.ensembl.org`. The June 2026 dated mirror returned HTTP 503 on 2026-09-27; the live host serves the same release. Before any variation query the script calls `GET /info/data` and stops if the reported release is not 116. TSI (Toscani in Italia) is the southern European group in that panel. Each subpopulation has its own column. If one rsID still fails after retries, that row keeps empty TSI, IBS, CEU, GBR, and FIN values and a note. The rest of the panel is still written. The empty cells are not filled from another source.

Minor-allele frequency (MAF) is the smaller of the alternate-allele frequency and its complement. A frequency that is missing in gnomAD or Ensembl stays empty, with a note on that row. It is not filled from another source.

For each variant the table also records:

- whether NFE MAF is below 0.01
- the minimum, maximum, and range of MAF across the EUR subpopulations that have a frequency (the range is empty when fewer than two subpopulations are present)
- `tsi_nfe_maf_abs_diff_gt_0_05_descriptive`: whether TSI MAF and NFE MAF differ by more than 0.05. This column is not a test.

The alternate allele, not the minor allele, is what is compared between 1000 Genomes phase 3 TSI and gnomAD v4 NFE. Each variant is a 2×2 table of alt-allele counts versus the other alleles. The test is two-sided Fisher exact, or Pearson chi-square when every expected count is at least 5. Benjamini-Hochberg q-values are computed across those tests. `tsi_differs_fdr05` is true when q < 0.05. A true flag means the alternate-allele frequency differs between public reference panels (1000G TSI vs gnomAD NFE). It is not a Romanian frequency and it is not a Southeast European frequency.

TSI allele number is the sum of Ensembl `allele_count` for `1000GENOMES:phase_3:TSI` when those counts are returned. If they are not, the allele number is 2 × 107 (the phase 3 TSI sample size), and that assumption is written in the provenance JSON. gnomAD NFE allele count and allele number are the cache fields `ac_nfe` and `an_nfe` for the panel alt. `af_gnomad_nfe` is not replaced. Columns `af_tsi`, `af_nfe`, and `af_diff` are that same alt-allele pair and their difference.

## Power

Design: 1:1 case-control (long-lived cases versus younger controls), two-sided allelic z-test, normal approximation. Allele counts are twice the number of people in each arm. The control frequency is the gnomAD v4 NFE MAF. Under allelic odds ratio \(R\) and control MAF \(p_0\), the case frequency is

\[
p_1 = \frac{R p_0}{1 - p_0 + R p_0}.
\]

The critical value uses the null standard error of the pooled frequency difference. The power uses the alternative standard error (`scipy.stats.norm`). Variants with missing, zero, or one NFE MAF have empty power. The control frequency is not imputed.

Grid:

| Axis | Values |
|------|--------|
| Cases per arm | 250, 500, 1000, 2500 |
| Allelic odds ratio | 1.1, 1.2, 1.5, 2.0 |
| Alpha | \(0.05/47\) (candidate-variant Bonferroni) and \(5 \times 10^{-8}\) |

The heatmap shows power at \(0.05/47\) only, with one panel per odds ratio. Rows are ordered by detectability rank.

**Rank.** Rank 1 is the highest illustrative power at 1000 cases per arm, allelic odds ratio 1.5, and alpha \(0.05/47\). Ties break toward a higher NFE MAF, then rsID. Missing power ranks last. The rank is about detectability under that one grid cell. It is not a LongevityMap priority and it is not evidence of a longevity association.

## Run

```bash
uv run python scripts/panel/ro_longevity_panel_power.py
```

Outputs, written under `analysis/panel/`:

| File | Contents |
|------|----------|
| `ro_longevity_panel.csv` | One row per variant: frequencies, flags, rank, notes |
| `ro_longevity_power_grid.csv` | One row per variant, cohort size, odds ratio, and alpha |
| `ro_longevity_power_heatmap.png` | Power heatmap at alpha \(0.05/47\) |
| `ro_longevity_panel_provenance.json` | gnomAD dataset, the Ensembl release and host reported by `/info/data`, query dates, cache hits versus new queries, and the release-check record |

After the two CSV tables are written, each is passed to `rogen-release-check --table PATH --table-kind summary --report PATH`. That kind runs the individual-level and direct-identifier checks and does not treat frequencies or power as cell counts. The return code, captured output, and JSON report path are stored in the provenance JSON. The scientific tables are not replaced by a sanitized copy.

The figure caption repeats the cohort-size disclaimer above.

## What this does not say

No variant is described as enriched, significant, or associated with longevity in Romanians. A high rank means the public-reference MAF and the chosen odds ratio would give a higher chance of rejecting the null in a study of that size. It does not mean the variant has that odds ratio, and it does not mean a ROGEN cohort of that size exists.
