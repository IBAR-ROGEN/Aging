# Longevity / lifespan GWAS lookup (candidate-gene replication)

Replication check of the prioritized longevity variants against large public
lifespan and longevity GWAS. **Most candidates will not replicate; that is a
valid result.** Effects are oriented to the panel alternate allele. This is not
a Romanian or ROGEN association analysis.

The Timmers 2019 statistics include UK Biobank parental lifespan data. Only
published summary statistics were used. No UK Biobank participant-level data
were accessed.

Some LongevityMap source studies may share cohorts with Deelen 2019 or
Timmers 2019. A hit is then not an independent replication. The result table
sets `possible_cohort_overlap = needs_curation` for every row until that is
reviewed by hand.

## Inputs

| Path | Role |
|------|------|
| `data/processed/variants_47_input_v2.csv` | Panel. Only `alt_status == "resolved"` rows are tested. Other statuses are written as `not_tested`. |
| `data/curation/longevity_allele_direction.csv` | Human curation of the longevity-associated allele (and optional PMID / population / notes). Prefill has `rsid` and `gene` only. Direction is computed only from this file. Empty allele → `direction_vs_reported_longevity_allele = needs_curation`. |
| `analysis/panel/allele_audit_47.csv` | Panel MAF (gnomAD NFE) for the palindrome guard. |
| `config/gwas_sources.yaml` | One block per summary-stat file: path, separator, columns, genome build, effect units, citation. Filled by inspecting each deposit README and header. |
| `data/external/gwas/*.gz` | Hand-downloaded summary statistics (git-ignored). |

`data/processed/allele_evidence_47.csv` column `longevitymap_alt_allele` is **not**
used for direction. That field was the LongevityMap alt allele used in September
for strand adjudication; it nearly always equals the panel alt and would only
restate the GWAS beta sign.

### Download (by hand)

1. **Timmers et al. 2019** parental lifespan — Edinburgh DataShare
   [doi:10.7488/ds/2463](https://doi.org/10.7488/ds/2463). File:
   `lifegen_phase2_bothpl_alldr_2017_09_18.tsv.gz`. Place under
   `data/external/gwas/`. README copy: [gwas/README_timmers2019_datashare.txt](gwas/README_timmers2019_datashare.txt).
2. **Deelen et al. 2019** longevity 90th and 99th percentile — GWAS Catalog
   [GCST008598](https://www.ebi.ac.uk/gwas/studies/GCST008598) /
   [GCST008599](https://www.ebi.ac.uk/gwas/studies/GCST008599). Files:
   `Results_90th_percentile.txt.gz`, `Results_99th_percentile.txt.gz`. README
   copies: [gwas/README_deelen2019_GCST008598.txt](gwas/README_deelen2019_GCST008598.txt),
   [gwas/README_deelen2019_GCST008599.txt](gwas/README_deelen2019_GCST008599.txt).

Deelen deposits have **no rsID column** in the table header (README mentions an
IMPUTE legend rsID). Those traits are matched on GRCh37 `Chr` + `Position`.
GRCh37 coordinates for panel rsIDs come from the Timmers file (`chr`, `pos`)
when present. rsIDs missing from Timmers are queried at
`https://grch37.rest.ensembl.org` (variation endpoint), cached under
`data/geo/ensembl_grch37_variation_cache/`, and recorded as
`grch37_source` = `timmers`, `ensembl_grch37`, or `not_found`.

### Effect units (from the source READMEs)

| Source | Units | Positive effect means |
|--------|-------|------------------------|
| Timmers 2019 | `beta1` = log hazard **protection** ratio (negation of log hazard ratio) | Longer parental lifespan for the effect allele |
| Deelen 2019 | `Beta` = LN(odds ratio) | Higher odds of surviving to the 90th/99th percentile |

## Run

```bash
uv run python scripts/longevity/longevity_gwas_lookup.py
```

The script scans each gzipped table lazily (Polars `scan_csv`) and keeps only
panel-matching rows. It never loads a full GWAS file into memory.

Per trait it prints: n present / absent / ambiguous / not_tested; observed
`p < 0.05` versus expected `0.05 × n_tests`; rsIDs with `q < 0.05`; rsIDs
passing Bonferroni. The same numbers are stored in the provenance JSON under
`multiple_testing.by_trait`.

### Allele harmonization

For each hit:

1. Same alleles as panel ref/alt → keep beta and EAF.
2. Swapped → flip beta, set EAF to `1 - eaf`.
3. Strand-flipped (complements match) → treat as same or swapped after complement.
4. A/T or C/G with panel MAF > 0.42 → `ambiguous_palindromic`; no flip.

Helpers live in `src/rogen_aging/longevity/allele_harmonize.py` (unit-tested).

### Direction versus reported longevity allele

Filled only from `data/curation/longevity_allele_direction.csv`:

- Empty `longevity_associated_allele` → `needs_curation`
- Allele not in `{panel_ref, panel_alt}` after complement strand check →
  `needs_curation` with a note in `direction_note`
- Otherwise `concordant` / `discordant` relative to the GWAS sign on the panel alt

Column name: `direction_vs_reported_longevity_allele`.

### Multiple testing

Over all successfully harmonized (`status=present`) variant × trait tests:

- Benjamini–Hochberg q-values
- Bonferroni threshold `0.05 / n_tests`
- Observed count of `p < 0.05` next to the null expectation `0.05 × n_tests`

## Outputs

| Path | Contents |
|------|----------|
| `analysis/longevity/longevity_gwas_lookup.csv` | One row per variant × trait: beta on panel alt, se, p, direction, status, curated direction, cohort-overlap flag, GRCh37 source, q, Bonferroni flag |
| `analysis/longevity/longevity_gwas_lookup_provenance.json` | File sha256, citations, GRCh37 sources, per-trait and overall multiple-testing counts, release-check record |
| `analysis/longevity/longevity_gwas_lookup_forest.png` | Forest of BH q < 0.05 hits, or top 10 by p if none |
| `analysis/longevity/release_checks/longevity_gwas_lookup.json` | `rogen-release-check` summary report |

Status values: `present`, `absent`, `ambiguous`, `not_tested`.

## Downstream

Per-variant public-evidence labels (GWAS × eQTL × TSI–NFE × power tier):

```bash
uv run python scripts/longevity/integrated_evidence.py
```

See [INTEGRATIVE_PIPELINE.md](INTEGRATIVE_PIPELINE.md#integrated-public-evidence-per-variant-labels).

## What this does not say

No variant is declared a longevity gene for Romanians or for ROGEN. A
non-significant or absent lookup is the expected outcome for most candidate
alleles and does not by itself refute prior candidate-gene literature.
