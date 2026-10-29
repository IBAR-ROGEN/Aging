# Technical report, 2026-10-30 (corrected)

Branch: `chore/sept-followup-to-main`. Numbers in §3 are produced by
`uv run python scripts/reports/build_results_table.py` from the cited files.
Analyses were not re-run for narrative sections except the integrated-evidence
refresh after the rs6435326 eQTL rule change.

The dated report `docs/reports/2026-09-27_technical_report.md` matches commit
`664881d` (`git diff 664881d -- docs/reports/2026-09-27_technical_report.md` is
empty).

## 1. Summary

September leftovers missing from `main` after the squash merge are on this
branch. Public HG002 short-read and ONT-versus-bisulfite dry runs finished on
the chr20 window. The published-clock benchmark ran (`d1b45e8`) and was
repaired (`98ebcf0`). Commit `8f13043` recorded that main’s tree already equals
`e0a2233` on this branch. The panel rebuild finished after Ensembl HTTP 500/503
(cache prefilled). Fix C recorded DT allele decisions: all 47
`alt_status=resolved`; GWAS lookup tests 47 variants. Integrated evidence was
refreshed; rs6435326’s significant palindromic eQTL hits count as eQTL support
with direction `ambiguous` (label counts now include that variant). Longevity
allele curation is still empty → direction `needs_curation`. ChrX/chrY
sensitivity produced a byte-identical clock; noXY joblib/CSVs were deleted,
provenance kept. No Romanian sequencing delivery is in the tree. PR to `main`
is not opened.

## 2. Sub-activities

### September leftovers on `chore/sept-followup-to-main`

An earlier merge attempt of `origin/feature/sept-secure-pipeline` into a branch
from `main` hit conflicts. Commit `8f13043` recorded that main’s tree already
equals `e0a2233` and is already on this branch. Commit `974539b` closed the
leftover checklist: sanity TSV + flag count in provenance, power-heatmap
caption, `analysis/eqtl/eqtl_datasets_used.tsv` (14 `dataset_id` rows),
`docs/PENDING_ALLELES.md`, and optional seq extra. `b6936fa` synced docs (that
commit also edited the September report; those edits are reverted so the
September file matches `664881d`). INPUT_MANIFEST later recorded the completed
HG002 dry runs (`339c814`, `11cc8b8`).

Status: DONE for the leftover doc/script items listed above.

### HG002 short-read dry run (public GIAB)

`scripts/dryrun/hg002_shortread_dryrun.py` downloads a chr20 window from a
pre-aligned Illumina BAM, runs `rogen-intake` with aliases NA24385 and GM24385,
aligns with minimap2, calls with untuned bcftools, and benchmarks with
`rtg vcfeval` against GIAB NISTv4.2.1. Metrics:
`analysis/dryrun/hg002_shortread_metrics.json`. Narrative: `docs/HG002_DRYRUN.md`
(table values rounded to three decimals; raw values stay in the JSON). Mapping
rate is near 100% by construction of the window realignment, not a whole-genome
rate.

Status: DONE-PUBLIC-DATA.

### ONT versus bisulfite at chr20 clock CpGs

`scripts/clock/ont_clock_site_coverage.py` compared ONT CRAM methylation
(modkit 0.6.4) to public bisulfite coverage on `chr20:10000000-20000000` at
three levels: all CpGs, HM450, Hannum non-zero. Metrics:
`analysis/dryrun/hg002_ont_methylation_metrics.json`. Docs:
`docs/METHYLATION_CLOCK_VALIDATION.md`. Hannum-in-window has n = 4 CpGs — too
few for a correlation (r kept only in the JSON). The ONT and bisulfite
materials were cultured and sequenced years apart.

Status: DONE-PUBLIC-DATA for the chr20 window.

### Published-clock benchmark

`scripts/clock/benchmark_published_clocks.py` scored pyaging 0.5.7 clocks
(`horvath2013`, `hannum`, `dnamphenoage`, `dunedinpace`, `grimage2`) plus
`rogen_hannum_retrain` on GSE87571 and the GSE40279 held-out split. Commit
`d1b45e8` added the benchmark; `98ebcf0` repaired GrimAge2 features and the GEO
matrix. Outputs under `analysis/clocks_benchmark/`; doc
`docs/CLOCK_BENCHMARK.md`. Four aggregate tables passed `rogen-release-check`.

Status: DONE.

### Clock transfer check

No `docs/CLOCK_TRANSFER_CHECK.md` and no transfer/dropout outputs under
`analysis/clock/`.

Status: not run.

### Genome-wide ONT clock sites

No genome-wide Hannum-site metrics file is in the tree.

Status: not run.

### Clock sex check (chrX / chrY)

Commit `441ba92` exported coefficients, retrained excluding chrX/chrY before the
age filter, and wrote `docs/CLOCK_SEX_CHECK.md`. Main model has 0 non-zero
chrX/chrY coefficients; neither main nor noXY shows a sex difference on
GSE40279 held-out or GSE87571. Later cleanup deleted the byte-identical noXY
joblib and CSVs; provenance JSON kept.

Status: DONE.

### Fix C — allele decisions, GWAS 47, integrated evidence

Allele rule (DT, 2026-10-28): alt = allele reported at consistent frequency by
gnomAD v4 NFE and 1000 Genomes CEU/GBR; LongevityMap supporting only. All 47
rows in `data/processed/variants_47_input_v2.csv` are `alt_status=resolved`.
Decision notes in `docs/PENDING_ALLELES.md`. GWAS lookup: `n_resolved=47`,
`n_not_tested=0`. Integrated evidence refreshed after counting significant
palindromic eQTL hits for rs6435326 as eQTL support with direction
`ambiguous`. Panel sanity: `flag_count=0`.

Status: DONE for allele resolution and GWAS/integrated refresh. Longevity-allele
curation from papers still open.

### Panel rebuild and TSI audit

Live Ensembl REST returned intermittent HTTP 500/503 during panel work.
Release-116 variation cache was prefilled. Panel provenance
`generated_at_utc` 2026-10-28T18:10:11Z: Ensembl `cache_hits` 47,
`new_queries` 0, `request_failures` 0, release 116.

TSI versus NFE allele audit (`analysis/panel/tsi_nfe_allele_audit.tsv`): all 47
rows `alt_match=true`, action `keep`. Five q < 0.05 rsIDs unchanged. TSI
ref-complement only at biallelic TSI sites (`624a559`).

Status: panel rebuild DONE. TSI audit DONE.

## 3. Results

Every number below is taken from the cited file by
`scripts/reports/build_results_table.py`.

| Item | Value | Source |
|------|-------|--------|
| HG002 region | chr20:10000000-20000000 | `analysis/dryrun/hg002_shortread_metrics.json` `region` |
| Mean depth | 28.761419123858087 | `analysis/dryrun/hg002_shortread_metrics.json` `mean_depth` |
| Percent mapped | 99.39 | `analysis/dryrun/hg002_shortread_metrics.json` `percent_mapped` |
| Percent duplicates | 0.25508340435894317 | `analysis/dryrun/hg002_shortread_metrics.json` `percent_duplicates` |
| SNV TP | 13214 | `analysis/dryrun/hg002_shortread_metrics.json` `snv.TP` |
| SNV FP | 71 | `analysis/dryrun/hg002_shortread_metrics.json` `snv.FP` |
| SNV FN | 93 | `analysis/dryrun/hg002_shortread_metrics.json` `snv.FN` |
| SNV precision | 0.9947 | `analysis/dryrun/hg002_shortread_metrics.json` `snv.precision` |
| SNV recall | 0.993 | `analysis/dryrun/hg002_shortread_metrics.json` `snv.recall` |
| SNV F1 | 0.9938 | `analysis/dryrun/hg002_shortread_metrics.json` `snv.f1` |
| Indel TP | 1898 | `analysis/dryrun/hg002_shortread_metrics.json` `indel.TP` |
| Indel FP | 182 | `analysis/dryrun/hg002_shortread_metrics.json` `indel.FP` |
| Indel FN | 172 | `analysis/dryrun/hg002_shortread_metrics.json` `indel.FN` |
| Indel precision | 0.9159 | `analysis/dryrun/hg002_shortread_metrics.json` `indel.precision` |
| Indel recall | 0.9169 | `analysis/dryrun/hg002_shortread_metrics.json` `indel.recall` |
| Indel F1 | 0.9164 | `analysis/dryrun/hg002_shortread_metrics.json` `indel.f1` |
| Short-read metrics time | 2026-10-28T10:02:01Z | `analysis/dryrun/hg002_shortread_metrics.json` `generated_at_utc` |
| Short-read metrics git commit | 339c81435e352ad01423dc5387b24fe70bd4444e | `analysis/dryrun/hg002_shortread_metrics.json` `git_commit` |
| Benchmark VCF | HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz | `analysis/dryrun/hg002_shortread_metrics.json` `benchmark_vcf` |
| Benchmark BED | HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed | `analysis/dryrun/hg002_shortread_metrics.json` `benchmark_bed` |
| Reference sha256 | 9cce8b926416dd96b152deea85188495b75f7ac8d634cc723a017067be8702b7 | `analysis/dryrun/hg002_shortread_metrics.json` `sha256.reference` |
| Benchmark VCF sha256 | adb4d4a50048aa13353a06b84fcfcbca09a5d17525efaa4cea44f8822e81175c | `analysis/dryrun/hg002_shortread_metrics.json` `sha256.benchmark_vcf` |
| Benchmark BED sha256 | fba9a57c36ec88e5d14ea3e259c8866c7935f4998d3ec0fa2d6c3962da5b5575 | `analysis/dryrun/hg002_shortread_metrics.json` `sha256.benchmark_bed` |
| Wall time fetch (s) | 3474.107 | `analysis/dryrun/hg002_shortread_metrics.json` `wall_time_seconds.fetch` |
| Wall time intake (s) | 3.169 | `analysis/dryrun/hg002_shortread_metrics.json` `wall_time_seconds.intake` |
| Wall time align (s) | 25.96 | `analysis/dryrun/hg002_shortread_metrics.json` `wall_time_seconds.align` |
| Wall time call (s) | 30.933 | `analysis/dryrun/hg002_shortread_metrics.json` `wall_time_seconds.call` |
| Wall time benchmark (s) | 5.522 | `analysis/dryrun/hg002_shortread_metrics.json` `wall_time_seconds.benchmark` |
| ONT all_cpgs n | 80896 | `analysis/dryrun/hg002_ont_methylation_metrics.json` `levels[level=all_cpgs].n_cpgs` |
| ONT all_cpgs depth≥10 | 80333 | `analysis/dryrun/hg002_ont_methylation_metrics.json` `levels[level=all_cpgs].n_ont_depth_ge_min` |
| ONT all_cpgs Pearson r | 0.9755608294915357 | `analysis/dryrun/hg002_ont_methylation_metrics.json` `levels[level=all_cpgs].pearson_r` |
| ONT all_cpgs mean ONT−bisulfite | 0.013932325591353035 | `analysis/dryrun/hg002_ont_methylation_metrics.json` `levels[level=all_cpgs].mean_ont_minus_bisulfite` |
| ONT hm450 n | 653 | `analysis/dryrun/hg002_ont_methylation_metrics.json` `levels[level=hm450].n_cpgs` |
| ONT hm450 depth≥10 | 634 | `analysis/dryrun/hg002_ont_methylation_metrics.json` `levels[level=hm450].n_ont_depth_ge_min` |
| ONT hm450 Pearson r | 0.9883002412213961 | `analysis/dryrun/hg002_ont_methylation_metrics.json` `levels[level=hm450].pearson_r` |
| ONT hm450 mean diff | 0.013610999427427344 | `analysis/dryrun/hg002_ont_methylation_metrics.json` `levels[level=hm450].mean_ont_minus_bisulfite` |
| ONT hannum_nonzero n | 4 | `analysis/dryrun/hg002_ont_methylation_metrics.json` `levels[level=hannum_nonzero].n_cpgs` |
| ONT hannum_nonzero depth≥10 | 4 | `analysis/dryrun/hg002_ont_methylation_metrics.json` `levels[level=hannum_nonzero].n_ont_depth_ge_min` |
| ONT hannum_nonzero Pearson r | n = 4, too few for a correlation | `analysis/dryrun/hg002_ont_methylation_metrics.json` `levels[level=hannum_nonzero].pearson_r` |
| ONT hannum_nonzero mean diff | 0.03086132666268903 | `analysis/dryrun/hg002_ont_methylation_metrics.json` `levels[level=hannum_nonzero].mean_ont_minus_bisulfite` |
| modkit | 0.6.4 --cpg --combine-strands --modified-bases 5mC 5hmC | `analysis/dryrun/hg002_ont_methylation_metrics.json` `modkit` |
| Panel variants | 47 | `analysis/panel/ro_longevity_panel_provenance.json` `input_variants.n_rows` |
| Panel provenance time | 2026-10-28T18:10:11Z | `analysis/panel/ro_longevity_panel_provenance.json` `generated_at_utc` |
| Ensembl release | 116 | `analysis/panel/ro_longevity_panel_provenance.json` `ensembl.release` |
| Ensembl host | https://rest.ensembl.org | `analysis/panel/ro_longevity_panel_provenance.json` `ensembl.host` |
| Ensembl cache hits | 47 | `analysis/panel/ro_longevity_panel_provenance.json` `ensembl.cache_hits` |
| Ensembl new queries | 0 | `analysis/panel/ro_longevity_panel_provenance.json` `ensembl.new_queries` |
| Ensembl request failures | 0 | `analysis/panel/ro_longevity_panel_provenance.json` `ensembl.request_failures` |
| gnomAD cache hits | 47 | `analysis/panel/ro_longevity_panel_provenance.json` `gnomad.cache_hits` |
| gnomAD new queries | 0 | `analysis/panel/ro_longevity_panel_provenance.json` `gnomad.new_queries` |
| TSI vs NFE n tested | 47 | `analysis/panel/ro_longevity_panel_provenance.json` `tsi_vs_nfe.n_tests` |
| TSI vs NFE n_chi_square | 46 | `analysis/panel/ro_longevity_panel_provenance.json` `tsi_vs_nfe.n_chi_square` |
| TSI vs NFE n_fisher | 1 | `analysis/panel/ro_longevity_panel_provenance.json` `tsi_vs_nfe.n_fisher` |
| TSI vs NFE q_below_fdr05 count | 5 | `analysis/panel/ro_longevity_panel_provenance.json` `len:tsi_vs_nfe.q_below_fdr05` |
| Panel release check status (first) | passed | `analysis/panel/ro_longevity_panel_provenance.json` `release_check[0].status` |
| Sanity flag_count | 0 | `analysis/panel/ro_longevity_panel_provenance.json` `sanity_check.flag_count` |
| GWAS n_resolved | 47 | `analysis/longevity/longevity_gwas_lookup_provenance.json` `panel.n_resolved` |
| GWAS n_not_tested | 0 | `analysis/longevity/longevity_gwas_lookup_provenance.json` `panel.n_not_tested` |
| GWAS n_tests | 138 | `analysis/longevity/longevity_gwas_lookup_provenance.json` `multiple_testing.n_tests` |
| GWAS n_p_lt_0_05 | 11 | `analysis/longevity/longevity_gwas_lookup_provenance.json` `multiple_testing.n_p_lt_0_05` |
| GWAS n_expected_p_lt_0_05_under_null | 6.9 | `analysis/longevity/longevity_gwas_lookup_provenance.json` `multiple_testing.n_expected_p_lt_0_05_under_null` |
| GWAS n_q_bh_lt_0_05 | 3 | `analysis/longevity/longevity_gwas_lookup_provenance.json` `multiple_testing.n_q_bh_lt_0_05` |
| GWAS n_bonferroni_sig | 3 | `analysis/longevity/longevity_gwas_lookup_provenance.json` `multiple_testing.n_bonferroni_sig` |
| GWAS lookup time | 2026-10-29T10:51:09Z | `analysis/longevity/longevity_gwas_lookup_provenance.json` `generated_at_utc` |
| GWAS release check | passed | `analysis/longevity/longevity_gwas_lookup_provenance.json` `release_check.status` |
| Integrated n_variants | 47 | `analysis/longevity/integrated_evidence_provenance.json` `n_variants` |
| Integrated GWAS support count | 1 | `analysis/longevity/integrated_evidence_provenance.json` `label_counts.GWAS support` |
| Integrated eQTL support count | 16 | `analysis/longevity/integrated_evidence_provenance.json` `label_counts.eQTL support` |
| Integrated population difference count | 5 | `analysis/longevity/integrated_evidence_provenance.json` `label_counts.population difference` |
| Integrated eQTL input status | ok | `analysis/longevity/integrated_evidence_provenance.json` `inputs.eqtl_input_status` |
| eQTL datasets locked | 14 | `analysis/longevity/integrated_evidence_provenance.json` `inputs.eqtl_datasets_used.n_datasets` |
| Integrated eqtl_n_ambiguous (sum) | 14 | `analysis/longevity/integrated_evidence.csv` `sum:eqtl_n_ambiguous` |
| Integrated eqtl_n_flipped (sum) | 0 | `analysis/longevity/integrated_evidence.csv` `sum:eqtl_n_flipped` |
| Integrated release check | passed | `analysis/longevity/integrated_evidence_provenance.json` `release_check.status` |
| Integrated time | 2026-10-29T13:17:04Z | `analysis/longevity/integrated_evidence_provenance.json` `generated_at_utc` |
| Clocks benchmark generated_at_utc | 2026-10-28T21:43:52Z | `analysis/clocks_benchmark/provenance.json` `generated_at_utc` |
| Clocks benchmark finished_at_utc | 2026-10-28T21:47:30Z | `analysis/clocks_benchmark/provenance.json` `finished_at_utc` |
| pyaging version | 0.5.7 | `analysis/clocks_benchmark/provenance.json` `pyaging_version` |
| CpG columns used | 28724 | `analysis/clocks_benchmark/provenance.json` `cpg_restriction.n_probes` |
| GSE87571 horvath2013 MAE | 3.89887351011231 | `analysis/clocks_benchmark/metrics_by_clock_cohort.csv` `clock=horvath2013,cohort=GSE87571:mae` |
| GSE87571 horvath2013 r | 0.971268432811853 | `analysis/clocks_benchmark/metrics_by_clock_cohort.csv` `clock=horvath2013,cohort=GSE87571:pearson_r` |
| GSE87571 horvath2013 mean signed error | 0.3329914019919036 | `analysis/clocks_benchmark/metrics_by_clock_cohort.csv` `clock=horvath2013,cohort=GSE87571:mean_signed_error` |
| GSE87571 hannum MAE | 5.132376501451235 | `analysis/clocks_benchmark/metrics_by_clock_cohort.csv` `clock=hannum,cohort=GSE87571:mae` |
| GSE87571 hannum r | 0.9749768268038559 | `analysis/clocks_benchmark/metrics_by_clock_cohort.csv` `clock=hannum,cohort=GSE87571:pearson_r` |
| GSE87571 hannum mean signed error | 4.011251363514884 | `analysis/clocks_benchmark/metrics_by_clock_cohort.csv` `clock=hannum,cohort=GSE87571:mean_signed_error` |
| GSE87571 dnamphenoage MAE | 9.337906601562581 | `analysis/clocks_benchmark/metrics_by_clock_cohort.csv` `clock=dnamphenoage,cohort=GSE87571:mae` |
| GSE87571 dnamphenoage r | 0.9681415391343388 | `analysis/clocks_benchmark/metrics_by_clock_cohort.csv` `clock=dnamphenoage,cohort=GSE87571:pearson_r` |
| GSE87571 dnamphenoage mean signed error | -8.898725968061047 | `analysis/clocks_benchmark/metrics_by_clock_cohort.csv` `clock=dnamphenoage,cohort=GSE87571:mean_signed_error` |
| GSE87571 dunedinpace r | 0.5926512212446541 | `analysis/clocks_benchmark/metrics_by_clock_cohort.csv` `clock=dunedinpace,cohort=GSE87571:pearson_r` |
| GSE87571 grimage2 MAE | 6.579675233207164 | `analysis/clocks_benchmark/metrics_by_clock_cohort.csv` `clock=grimage2,cohort=GSE87571:mae` |
| GSE87571 grimage2 r | 0.9715219859184909 | `analysis/clocks_benchmark/metrics_by_clock_cohort.csv` `clock=grimage2,cohort=GSE87571:pearson_r` |
| GSE87571 grimage2 mean signed error | 5.408548010663117 | `analysis/clocks_benchmark/metrics_by_clock_cohort.csv` `clock=grimage2,cohort=GSE87571:mean_signed_error` |
| GSE87571 rogen_hannum_retrain MAE | 3.760192970835657 | `analysis/clocks_benchmark/metrics_by_clock_cohort.csv` `clock=rogen_hannum_retrain,cohort=GSE87571:mae` |
| GSE87571 rogen_hannum_retrain r | 0.9867447892597446 | `analysis/clocks_benchmark/metrics_by_clock_cohort.csv` `clock=rogen_hannum_retrain,cohort=GSE87571:pearson_r` |
| GSE87571 rogen_hannum_retrain mean signed error | 3.1032341482232693 | `analysis/clocks_benchmark/metrics_by_clock_cohort.csv` `clock=rogen_hannum_retrain,cohort=GSE87571:mean_signed_error` |
| GSE40279 held-out rogen_hannum_retrain MAE | 3.524863951133959 | `analysis/clocks_benchmark/metrics_by_clock_cohort.csv` `clock=rogen_hannum_retrain,cohort=GSE40279_heldout:mae` |
| GSE40279 held-out rogen_hannum_retrain r | 0.9566776321483395 | `analysis/clocks_benchmark/metrics_by_clock_cohort.csv` `clock=rogen_hannum_retrain,cohort=GSE40279_heldout:pearson_r` |
| GSE40279 held-out rogen_hannum_retrain mean signed error | -0.054869463949492485 | `analysis/clocks_benchmark/metrics_by_clock_cohort.csv` `clock=rogen_hannum_retrain,cohort=GSE40279_heldout:mean_signed_error` |
| GSE87571 horvath2013 sex coef | -0.9985615117003562 | `analysis/clocks_benchmark/sex_acceleration_models.csv` `clock=horvath2013,cohort=GSE87571:sex_coef_female_minus_male` |
| GSE87571 horvath2013 sex Bonf p | 0.04443506161549815 | `analysis/clocks_benchmark/sex_acceleration_models.csv` `clock=horvath2013,cohort=GSE87571:sex_p_bonferroni` |
| GSE87571 hannum sex coef | -2.089089525560717 | `analysis/clocks_benchmark/sex_acceleration_models.csv` `clock=hannum,cohort=GSE87571:sex_coef_female_minus_male` |
| GSE87571 hannum sex Bonf p | 1.6008054287839974e-09 | `analysis/clocks_benchmark/sex_acceleration_models.csv` `clock=hannum,cohort=GSE87571:sex_p_bonferroni` |
| GSE87571 dnamphenoage sex coef | 0.21940909788234944 | `analysis/clocks_benchmark/sex_acceleration_models.csv` `clock=dnamphenoage,cohort=GSE87571:sex_coef_female_minus_male` |
| GSE87571 dnamphenoage sex Bonf p | 1.0 | `analysis/clocks_benchmark/sex_acceleration_models.csv` `clock=dnamphenoage,cohort=GSE87571:sex_p_bonferroni` |
| GSE87571 grimage2 sex coef | -0.8227415397694006 | `analysis/clocks_benchmark/sex_acceleration_models.csv` `clock=grimage2,cohort=GSE87571:sex_coef_female_minus_male` |
| GSE87571 grimage2 sex Bonf p | 0.04688537609513253 | `analysis/clocks_benchmark/sex_acceleration_models.csv` `clock=grimage2,cohort=GSE87571:sex_p_bonferroni` |
| GSE87571 rogen_hannum_retrain sex coef | -0.23503770604855034 | `analysis/clocks_benchmark/sex_acceleration_models.csv` `clock=rogen_hannum_retrain,cohort=GSE87571:sex_coef_female_minus_male` |
| GSE87571 rogen_hannum_retrain sex Bonf p | 1.0 | `analysis/clocks_benchmark/sex_acceleration_models.csv` `clock=rogen_hannum_retrain,cohort=GSE87571:sex_p_bonferroni` |
| v2 alt_status=resolved count | 47 | `data/processed/variants_47_input_v2.csv` `count_where:alt_status=resolved` |
| Allele audit keep count | 47 | `analysis/panel/tsi_nfe_allele_audit.tsv` `count_where:action=keep` |
| Sanity TSV row count | 47 | `analysis/panel/ro_longevity_panel_sanity.tsv` `row_count` |
| Sanity complement_skipped count | 10 | `analysis/panel/ro_longevity_panel_sanity.tsv` `count_where:complement_skipped=true` |

## 4. Status

| Item | Status | Evidence | Next step |
|------|--------|----------|-----------|
| September leftovers | DONE | `974539b`, sanity TSV, `PENDING_ALLELES.md` | — |
| Merge feature tree on follow-up branch | DONE | `8f13043` | Open PR to `main` when ready |
| HG002 short-read dry run | DONE-PUBLIC-DATA | `hg002_shortread_metrics.json` | Do not label GIAB as Romanian/ROGEN |
| ONT vs bisulfite (chr20) | DONE-PUBLIC-DATA | `hg002_ont_methylation_metrics.json` | Genome-wide ONT clock sites still open |
| Genome-wide ONT clock sites | not run | no metrics artefact | Run or defer in writing |
| Published-clock benchmark | DONE | `d1b45e8`, `98ebcf0` | — |
| Clock transfer check | not run | no transfer outputs | Run if sample-size planning needs it |
| Clock sex check / noXY cleanup | DONE | `441ba92`, `CLOCK_SEX_CHECK.md` | — |
| Fix C allele decisions | DONE | v2 all resolved | — |
| Panel rebuild | DONE | provenance 2026-10-28T18:10:11Z | Prefer cache when REST unstable |
| TSI vs NFE allele audit | DONE | 47 keep; q list unchanged | Keep biallelic-only complement |
| Longevity GWAS lookup | DONE-PUBLIC-DATA | 47 tested | Fill curated longevity alleles |
| Integrated evidence | DONE-PUBLIC-DATA | eQTL support includes rs6435326 ambiguous | Rebuild eQTL tables vs v2 if needed |
| Intake on Romanian delivery | DONE-SYNTHETIC | `docs/DATA_FLOW.md` | Data-protection review pending |
| PR to `main` | Not opened | — | Human opener |

## 5. Decisions and problems

Dated September report text stays at `664881d`.

HG002 short-read reads came from a pre-aligned GIAB BAM window and were
realigned to chr20 only. Percent mapped near 100% is by construction.
Precision and recall are window-scoped. Metrics JSON now records
`git_commit`, `benchmark_vcf`, `benchmark_bed`, and `sha256` (verified against
files under `data/dryrun/hg002_shortread/`).

Ensembl live REST returned HTTP 500/503 during the panel rebuild. Release-116
cache prefilled; panel rebuild completed with `cache_hits` 47, `new_queries` 0.

TSI audit: panel alt matches gnomAD counts for all 47 rows. “5 of 47” q < 0.05
unchanged. After Fix C those five no longer rest on unconfirmed alts.

GWAS direction uses only curated longevity alleles; empty curation ⇒
`needs_curation`. Only rs4420638 reaches q < 0.05 (three traits).

rs6435326: 14 eQTL association rows are palindromic (MAF > 0.42); 4 of them
have `fdr_bh` < 0.05. They count as eQTL support with
`eqtl_alt_effect_direction=ambiguous` (4 significant datasets).

Integrated eQTL inputs are local gitignored `results/` files.

POD5: with pod5 0.3.49, `Reader` has no `run_infos`; intake still reports
metadata not scanned.

## 6. Open items

1. Open PR for `chore/sept-followup-to-main` → `main`.
2. Fill `data/curation/longevity_allele_direction.csv` from papers (not
   LongevityMap), especially rs28366003 and rs882696.
3. Genome-wide ONT clock sites — run or defer explicitly.
4. Clock transfer check — run if sample-size planning is in scope.
5. Rebuild or document local eQTL Catalogue tables against
   `variants_47_input_v2.csv`.
6. Data-protection review of `docs/DATA_FLOW.md`; Romanian delivery intake
   still blocked.
7. POD5 metadata scan versus pod5 Reader API (`run_infos` missing on 0.3.49).
8. Four August branches not merged:
   `origin/fix/ci-omitted-tests-mypy` (`bb6249c`),
   `origin/feature/july-clock-console-entries` (`689d301`),
   `origin/fix/ensembl-clock-eval-dedup` (`cbf355c`),
   `origin/chore/ruff-unify-tooling` (`704df36`).

## 7. Provenance

Selected commits on `chore/sept-followup-to-main` (`git log --date=iso`):

| Hash | Time | Message |
|------|------|---------|
| 82f8f98 | 2026-10-29 12:53:28 +0200 | Apply Fix C: resolve 47 panel alleles and rerun GWAS. |
| 441ba92 | 2026-10-29 00:24:53 +0200 | Make the Hannum clock reproducible and test chrX/chrY sex sensitivity. |
| 98ebcf0 | 2026-10-28 23:49:45 +0200 | Repair published-clock benchmark: GrimAge2 features and GEO matrix. |
| 31b2886 | 2026-10-28 21:17:34 +0200 | Orient integrated eQTL evidence to the panel alt. |
| 624a559 | 2026-10-28 20:42:28 +0200 | Audit TSI vs NFE alleles and restrict TSI complement to biallelic sites. |
| c419c9e | 2026-10-28 20:25:33 +0200 | Fix GWAS direction to use curated longevity alleles only. |
| d1b45e8 | 2026-10-28 20:08:15 +0200 | Add published-clock benchmark and integrated longevity evidence. |
| 8f13043 | 2026-10-28 19:40:25 +0200 | Record main as merged; main's tree equals e0a2233, already in this branch |
| cafe70b | 2026-10-28 15:28:45 +0200 | Document the longevity GWAS lookup inputs and script paths. |
| d732b60 | 2026-10-28 15:28:14 +0200 | Add candidate-gene lookup against public longevity GWAS. |
| 11cc8b8 | 2026-10-28 14:52:39 +0200 | Record ONT vs bisulfite HG002 chr20 comparison results. |
| 339c814 | 2026-10-28 12:29:33 +0200 | Note the completed HG002 short-read dry run in INPUT_MANIFEST. |
| e5e8033 | 2026-10-28 12:29:18 +0200 | Add HG002 short-read dry run and preserve paired FASTQ names. |
| b6936fa | 2026-10-28 10:46:36 +0200 | Document September follow-up: seq extra, sanity TSV, pending alleles. |
| 974539b | 2026-10-28 07:59:17 +0000 | Close September leftovers: panel sanity, alleles, heatmap, seq extra. |
