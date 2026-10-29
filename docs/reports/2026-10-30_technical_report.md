# Technical report, 2026-10-30 (updated 2026-10-29)

Branch: `chore/sept-followup-to-main` tip `82f8f98`. Numbers below are taken from the cited files and `git log`; analyses for this note were not re-run beyond what Fix C already wrote on 2026-10-29.

The dated report `docs/reports/2026-09-27_technical_report.md` stays at its text from commit `664881d`. October status belongs here.

Repo-improvement Fixes 1–6 (August review) and September Prompts 5–6 appear only where they actually ran. Fix 4 (shim deletion), Fix 6 (optional-extras cleanup), and Prompt 6 (no discrete labeled run found) are not claimed as done.

## 1. Summary

September leftovers missing from `main` after the squash merge are on this branch. Public HG002 short-read and ONT-versus-bisulfite dry runs finished on the chr20 window. A pyaging published-clock benchmark sits next to the retrained Hannum model (`d1b45e8`, repaired in `98ebcf0`). Commit `8f13043` recorded that main’s tree already equals `e0a2233` on this branch. The 47-variant panel was rebuilt after Ensembl HTTP 500/503 (cache prefilled). Fix C (`82f8f98`) recorded DT allele decisions: all 47 `alt_status=resolved`; GWAS lookup now tests 47 variants (`n_not_tested=0`); integrated evidence and panel sanity provenance were refreshed. Longevity-allele curation is still empty → direction `needs_curation`. ChrX/chrY sensitivity produced a byte-identical clock; noXY joblib/CSVs were deleted, provenance kept (`docs/CLOCK_SEX_CHECK.md`). No Romanian sequencing delivery is in the tree. PR to `main` is not opened.

## 2. Sub-activities

### September leftovers on `chore/sept-followup-to-main`

An earlier merge attempt of `origin/feature/sept-secure-pipeline` into a branch from `main` hit conflicts. Commit `8f13043` recorded that main’s tree already equals `e0a2233` and is already on this branch (merge resolved). Commit `974539b` closed the leftover checklist: sanity TSV + flag count in provenance, power-heatmap caption, `analysis/eqtl/eqtl_datasets_used.tsv` (14 `dataset_id` rows), `docs/PENDING_ALLELES.md`, and optional seq extra. `b6936fa` synced docs. INPUT_MANIFEST later recorded the completed HG002 dry runs (`339c814`, `11cc8b8`).

Status: DONE for the leftover doc/script items listed above.

### A.2.1.8.1 HG002 short-read dry run (public GIAB)

`scripts/dryrun/hg002_shortread_dryrun.py` downloads a chr20 window from a pre-aligned Illumina BAM, runs `rogen-intake` with aliases NA24385 and GM24385, aligns with minimap2, calls with untuned bcftools, and benchmarks with `rtg vcfeval` against GIAB NISTv4.2.1. Metrics: `analysis/dryrun/hg002_shortread_metrics.json`. Narrative: `docs/HG002_DRYRUN.md` (table values rounded to three decimals; raw values stay in the JSON). Mapping rate is near 100% by construction of the window realignment, not a whole-genome rate.

Status: DONE-PUBLIC-DATA.

### A.2.1.10.1 ONT versus bisulfite at chr20 clock CpGs

`scripts/clock/ont_clock_site_coverage.py` compared ONT CRAM methylation (modkit 0.6.4) to public bisulfite coverage on `chr20:10000000-20000000` at three levels: all CpGs, HM450, Hannum non-zero. Metrics: `analysis/dryrun/hg002_ont_methylation_metrics.json`. Docs: `docs/METHYLATION_CLOCK_VALIDATION.md`. The ONT and bisulfite materials were cultured and sequenced years apart. A later task to extend the comparison genome-wide over 2314 Hannum sites was opened; no genome-wide metrics file from that extension is in the tree for this report.

Status: DONE-PUBLIC-DATA for the chr20 window. BUILT-NOT-RUN for the genome-wide extension.

### Published-clock benchmark

`scripts/clock/benchmark_published_clocks.py` scored pyaging 0.5.7 clocks (`horvath2013`, `hannum`, `dnamphenoage`, `dunedinpace`, `grimage2`) plus `rogen_hannum_retrain` on GSE87571 (n=729) and the GSE40279 held-out split (n=132). Commit `d1b45e8` added the benchmark and integrated evidence; `98ebcf0` repaired GrimAge2 features and the GEO matrix. Current provenance: `generated_at_utc` 2026-10-28T21:43:52Z, `finished_at_utc` 2026-10-28T21:47:30Z. Outputs under `analysis/clocks_benchmark/`; doc `docs/CLOCK_BENCHMARK.md`. Four aggregate tables passed `rogen-release-check`. Clock transfer / missing-CpG dropout: no `docs/CLOCK_TRANSFER_CHECK.md` / transfer outputs in tree.

Status: DONE for the published-clock benchmark. BUILT-NOT-RUN for the transfer/dropout check.

### Clock sex check (chrX / chrY)

Commit `441ba92` exported coefficients, retrained excluding chrX/chrY before the age filter, and wrote `docs/CLOCK_SEX_CHECK.md`. Main model has 0 non-zero chrX/chrY coefficients; neither main nor noXY shows a sex difference on GSE40279 held-out or GSE87571. Fix C (`82f8f98`) deleted the byte-identical noXY joblib and CSVs; provenance JSON kept. Sentence in the doc: removing chrX/chrY changed two candidate probes, neither in the top 8,000, so the retrained model was byte-identical and is not stored separately.

Status: DONE.

### Fix C — allele decisions, GWAS 47, integrated evidence

Commit `82f8f98`. Rule (DT, 2026-10-28): alt = allele reported at consistent frequency by gnomAD v4 NFE and 1000 Genomes CEU/GBR; LongevityMap supporting only. Twelve alts confirmed without changing alleles; `alt_status` set to `resolved` for all 47 rows in `data/processed/variants_47_input_v2.csv`. Decision notes in `docs/PENDING_ALLELES.md` (strand flips for rs28366003 / rs882696; palindrome note for rs155979). GWAS lookup re-run: `n_resolved=47`, `n_not_tested=0`. Integrated evidence re-run. Panel sanity re-run; provenance `sanity_check.flag_count=0`. Curation CSV: notes only for the two LongevityMap strand flips; longevity alleles not filled from papers.

Status: DONE for Fix C steps that edit tree / rerun lookup. Longevity-allele curation from papers still open.

### A.2.1.11.1 Panel rebuild (Prompt 5), TSI audit

Live Ensembl REST returned intermittent HTTP 500/503 during panel work. Release-116 variation cache was prefilled. Prompt 5 provenance `generated_at_utc` 2026-10-28T18:10:11Z: Ensembl `cache_hits` 47, `new_queries` 0, `request_failures` 0, release 116.

TSI versus NFE allele audit (`analysis/panel/tsi_nfe_allele_audit.tsv`): all 47 rows `alt_match=true`, action `keep`. Five q < 0.05 rsIDs unchanged (`tsi_nfe_q_before_after.tsv`). TSI ref-complement only at biallelic TSI sites (`624a559`).

Prompt 6: no discrete labeled execution found; not marked done.

Status: panel Prompt 5 DONE. TSI audit DONE.

### Repo-improvement Fixes 1–6 (August; status only)

| Fix | Intent | Ran? | Evidence |
|-----|--------|------|----------|
| 1 | CI: config + methylation preflight tests; mypy on `src/rogen_aging` | yes | `bb6249c` |
| 2 | Promote July annotate / GSE87571 eval to package CLIs | yes | `689d301` |
| 3 | Route VEP and GSE87571 loading through shared package APIs | yes | `cbf355c` |
| 4 | Delete root/scripts runpy shims | no | skipped (session went to Fix 5) |
| 5 | Docs lint/format on Ruff; drop black/isort/flake8 | yes | `704df36` |
| 6 | Optional extras / lower-priority cleanup from the same review | no | no matching commit in that sequence |

## 3. Results

Every number below is taken from the cited file.

| Item | Value | Source |
|------|-------|--------|
| HG002 region | chr20:10000000-20000000 | `analysis/dryrun/hg002_shortread_metrics.json` `region` |
| Mean depth | 28.761419123858087 | same file, `mean_depth` |
| Percent mapped | 99.39 | same file, `percent_mapped` |
| Percent duplicates | 0.25508340435894317 | same file, `percent_duplicates` |
| SNV TP / FP / FN | 13214 / 71 / 93 | same file, `snv` |
| SNV precision / recall / F1 | 0.9947 / 0.993 / 0.9938 | same file, `snv` |
| Indel TP / FP / FN | 1898 / 182 / 172 | same file, `indel` |
| Indel precision / recall / F1 | 0.9159 / 0.9169 / 0.9164 | same file, `indel` |
| Short-read metrics time | 2026-10-28T10:02:01Z | same file, `generated_at_utc` |
| Short-read metrics git commit | 339c81435e352ad01423dc5387b24fe70bd4444e | same file, `git_commit` |
| Benchmark VCF | HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz | same file, `benchmark_vcf` |
| Benchmark BED | HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed | same file, `benchmark_bed` |
| Reference sha256 | 9cce8b926416dd96b152deea85188495b75f7ac8d634cc723a017067be8702b7 | same file, `sha256.reference` |
| Benchmark VCF sha256 | adb4d4a50048aa13353a06b84fcfcbca09a5d17525efaa4cea44f8822e81175c | same file, `sha256.benchmark_vcf` |
| Benchmark BED sha256 | fba9a57c36ec88e5d14ea3e259c8866c7935f4998d3ec0fa2d6c3962da5b5575 | same file, `sha256.benchmark_bed` |
| Wall time fetch / intake / align / call / benchmark (s) | 3474.107 / 3.169 / 25.96 / 30.933 / 5.522 | same file, `wall_time_seconds` |
| ONT all_cpgs n / depth≥10 / Pearson r / mean ONT−bisulfite | 80896 / 80333 / 0.9755608294915357 / 0.013932325591353035 | `analysis/dryrun/hg002_ont_methylation_metrics.json` |
| ONT hm450 n / depth≥10 / r / mean diff | 653 / 634 / 0.9883002412213961 / 0.013610999427427344 | same file |
| ONT hannum_nonzero n / depth≥10 / r / mean diff | 4 / 4 / 0.9979509346301134 / 0.03086132666268903 | same file |
| modkit | 0.6.4 --cpg --combine-strands --modified-bases 5mC 5hmC | same file, `modkit` |
| Panel variants | 47 | `analysis/panel/ro_longevity_panel_provenance.json` `input_variants.n_rows` |
| Panel provenance time | 2026-10-28T18:10:11Z | same file, `generated_at_utc` |
| Ensembl release / host | 116 / https://rest.ensembl.org | same file, `ensembl` |
| Ensembl cache hits / new queries / failures | 47 / 0 / 0 | same file, `ensembl` |
| gnomAD cache hits / new queries | 47 / 0 | same file, `gnomad` |
| TSI vs NFE tests | 47 (46 chi-square, 1 Fisher) | same file, `tsi_vs_nfe` |
| Rows with `tsi_differs_fdr05` | 5 of 47 | same file, `tsi_vs_nfe.q_below_fdr05` |
| Allele audit alt_match / keep | 47 / 47 | row tally of `analysis/panel/tsi_nfe_allele_audit.tsv` |
| v2 `alt_status=resolved` | 47 of 47 | row tally of `data/processed/variants_47_input_v2.csv` |
| Sanity flagged rows | 0 of 47 | `analysis/panel/ro_longevity_panel_sanity.tsv`; provenance `sanity_check.flag_count` |
| Sanity complement_skipped | 10 | same TSV |
| Panel / power-grid release check | passed, return code 0 | provenance `release_check` |
| GWAS resolved / not tested | 47 / 0 | `analysis/longevity/longevity_gwas_lookup_provenance.json` `panel` |
| GWAS n_tests / n_p&lt;0.05 / expected under null | 138 / 11 / 6.9 | same file, `multiple_testing` |
| GWAS n_q_bh&lt;0.05 / Bonferroni sig | 3 / 3 (all rs4420638) | same file, `multiple_testing` |
| GWAS lookup time | 2026-10-29T10:51:09Z | same file, `generated_at_utc` |
| GWAS release check | passed | same file, `release_check.status` |
| Integrated rows | 47 | row tally of `analysis/longevity/integrated_evidence.csv` |
| Integrated label counts | GWAS support 1; eQTL support 15; population difference 5 | `analysis/longevity/integrated_evidence_provenance.json` `label_counts` |
| Integrated eQTL input status | ok | same file, `inputs.eqtl_input_status` |
| eQTL datasets locked | 14 | same file; also `analysis/eqtl/eqtl_datasets_used.tsv` |
| Integrated eqtl_n_ambiguous (sum) | 14 | sum of column in `integrated_evidence.csv` |
| Integrated eqtl_n_flipped (sum) | 0 | same CSV |
| Integrated release check | passed | provenance `release_check` |
| Integrated time | 2026-10-29T10:51:51Z | provenance `generated_at_utc` |
| Clocks benchmark generated / finished | 2026-10-28T21:43:52Z / 2026-10-28T21:47:30Z | `analysis/clocks_benchmark/provenance.json` |
| pyaging version | 0.5.7 | same file |
| CpG columns used | 28724 | same file, `cpg_restriction.n_probes` |
| GSE87571 horvath2013 MAE / r / mean signed error | 3.8986034530126403 / 0.9712699385193299 / 0.33326145909157345 | `analysis/clocks_benchmark/metrics_by_clock_cohort.csv` |
| GSE87571 hannum MAE / r / mean signed error | 5.131737360643968 / 0.9749669350179596 / 4.010110995632912 | same CSV |
| GSE87571 dnamphenoage MAE / r / mean signed error | 9.33871368354254 / 0.9681562373823652 / −8.899621138285726 | same CSV |
| GSE87571 dunedinpace r / mean signed error | 0.5926332174273056 / −46.41822103394301 (MAE omitted: rate) | same CSV |
| GSE87571 grimage2 MAE / r / mean signed error | 18.8237875950521 / 0.7540762407262465 / 15.200331096666627 | same CSV |
| GSE87571 rogen_hannum_retrain MAE / r / mean signed error | 3.760192970835657 / 0.9867447892597446 / 3.1032341482232693 | same CSV |
| GSE40279 held-out rogen_hannum_retrain MAE / r / mean signed error | 3.524863951133959 / 0.9566776321483395 / −0.054869463949492485 | same CSV |
| GSE87571 hannum sex coef (female−male) / Bonf p | −2.08670560467607 / 1.6904162092147052e-09 | `analysis/clocks_benchmark/sex_acceleration_models.csv` |
| GSE87571 horvath2013 sex coef / Bonf p | −0.9989383079492378 / 0.04429368968127986 | same CSV |
| GSE87571 rogen_hannum_retrain sex Bonf p | 1.0 | same CSV |
| Clock aggregate release checks | 4 passed | provenance `release_check` |
| noXY joblib stored | no (deleted; byte-identical to main) | `82f8f98`; provenance JSON kept |

## 4. Status

| Item | Status | Evidence | Next step |
|------|--------|----------|-----------|
| September leftovers | DONE | `974539b`, sanity TSV, `PENDING_ALLELES.md` | — |
| Merge feature→main on follow-up branch | DONE (tree) | `8f13043` | Open PR to `main` when ready |
| HG002 short-read dry run | DONE-PUBLIC-DATA | `hg002_shortread_metrics.json` | Do not label GIAB as Romanian/ROGEN |
| ONT vs bisulfite (chr20) | DONE-PUBLIC-DATA | `hg002_ont_methylation_metrics.json` | Genome-wide extension still open |
| Genome-wide ONT clock extension | BUILT-NOT-RUN | no metrics artefact | Run or defer in writing |
| Published-clock benchmark | DONE | `d1b45e8`, `98ebcf0`, provenance 2026-10-28T21:43:52Z | Transfer/dropout still open |
| Clock sex check / noXY cleanup | DONE | `441ba92`, `82f8f98`, `CLOCK_SEX_CHECK.md` | — |
| Fix C allele decisions | DONE | `82f8f98`, v2 all resolved | — |
| Panel Prompt 5 rerun | DONE | provenance 2026-10-28T18:10:11Z | Prefer cache when REST unstable |
| Prompt 6 | not run | — | Identify intent before claiming |
| TSI vs NFE allele audit | DONE | 47 keep; q list unchanged | Keep biallelic-only complement |
| Longevity GWAS lookup | DONE-PUBLIC-DATA | provenance 2026-10-29T10:51:09Z; **47 tested** | Fill curated longevity alleles |
| Integrated evidence | DONE-PUBLIC-DATA | provenance 2026-10-29T10:51:51Z; eQTL `ok` | Rebuild eQTL tables vs v2 if needed |
| Intake on Romanian delivery | DONE-SYNTHETIC | `docs/DATA_FLOW.md` | Data-protection review pending |
| Fix 1 / 2 / 3 / 5 (August) | DONE | `bb6249c`, `689d301`, `cbf355c`, `704df36` | — |
| Fix 4 / Fix 6 (August) | not run | — | Only if still wanted |
| PR to `main` | Not opened | — | Human opener |

## 5. Decisions and problems

Dated September report text stays at `664881d`.

HG002 short-read reads came from a pre-aligned GIAB BAM window and were realigned to chr20 only. Percent mapped near 100% is by construction. Precision and recall are window-scoped.

Ensembl live REST returned HTTP 500/503 during the October panel rebuild. Release-116 cache prefilled; Prompt 5 completed with `cache_hits` 47, `new_queries` 0.

TSI audit: panel alt matches gnomAD counts for all 47 rows. “5 of 47” q < 0.05 unchanged. After Fix C those five no longer rest on unconfirmed alts (TSI already used these alts).

GWAS direction uses only curated longevity alleles; empty curation ⇒ `needs_curation`. Only rs4420638 reaches q < 0.05 (three traits). LongevityMap must not fill curation for rs28366003 / rs882696 (opposite strand).

Integrated eQTL inputs are local gitignored `results/` files. Orientation recorded 14 ambiguous palindromic rows and 0 flips.

POD5: with pod5 0.3.49, `Reader` has no `run_infos`; intake still reports metadata not scanned.

Earlier chat reports that said the merge failed or the clock benchmark was not run were wrong: `8f13043` resolved the merge; `d1b45e8` and `98ebcf0` ran and repaired the benchmark. This file uses git and output files only.

## 6. Open items

1. Open PR for `chore/sept-followup-to-main` → `main`.
2. Fill `data/curation/longevity_allele_direction.csv` from papers (not LongevityMap), especially rs28366003 and rs882696.
3. Genome-wide ONT versus bisulfite over Hannum non-zero CpGs, or explicit deferral.
4. Clock transfer / missing-CpG dropout outputs if sample-size planning is in scope.
5. Rebuild or document local eQTL Catalogue tables against `variants_47_input_v2.csv`.
6. Data-protection review of `docs/DATA_FLOW.md`; Romanian delivery intake still blocked.
7. Decide whether August Fix 4 / Fix 6 are still wanted.
8. Identify Prompt 6 if it was meant to be a specific October run.
9. POD5 metadata scan versus pod5 Reader API (`run_infos` missing on 0.3.49).

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

August Fix commits cited above (not re-listed as October work): `bb6249c`, `689d301`, `cbf355c`, `704df36`.
