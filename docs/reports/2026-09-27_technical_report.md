# Technical report, 2026-09-27

Branch: `feature/sept-secure-pipeline`. Read from git history and the files named below. No analysis was rerun for this note.

`git status` on this branch is clean. Nothing is uncommitted.

`outputs/clock_metrics.json` is not in the tree. It is not used as a result from today. The copy under `outputs/archive/fixture_2026-07/` is the July fixture. See Decisions.

## 1. Summary

This morning the branch did not have a sequencing intake command, a release guard, a GSE40279 clock with Illumina probe names, or a public-frequency table for the 47 prioritized variants. Those pieces are now committed. `rogen-intake` checks checksums, writes pseudonyms, and treats manifest aliases as identifiers. It has been exercised on synthetic files only. A Hannum ElasticNet pipeline was trained on public GSE40279 data, with a recorded train and test split, and a later commit wrote GSE87571 evaluation metrics for that model. The old `cg_test_*` clock was separated as a fixture and its July metrics were archived. A panel script stored gnomAD v4 NFE and 1000 Genomes phase 3 EUR frequencies for 47 variants, an illustrative power grid, and two release-check reports that passed. Comparison of ONT and bisulfite methylation at chr20 clock CpGs is coded and not run. No HG002 or GM24385 files were downloaded.

## 2. Sub-activities

### A.2.1.8.1 Secure intake and dry-run preparation

Built `rogen-intake` (`src/rogen_aging/cli/intake.py`, package `src/rogen_aging/intake/`) and `docs/DATA_FLOW.md`. Intake verifies checksums before it writes a processing directory, derives an `RO-` pseudonym from a key that must sit outside the checkout, rewrites sample identifiers in FASTQ, BAM, CRAM, and VCF outputs, and fails if an original id or a manifest alias remains. Aliases were added in a later commit so a header id that is not the primary sample id still fails the leak check. POD5 files are renamed after the checksum check. If the `pod5` package is missing, the report says `pod5 metadata not scanned`.

`docs/DATA_FLOW.md` states that no Romanian sequencing data has arrived and that intake is exercised on small synthetic files the tests create. The design is not yet reviewed by the institutional data protection officer.

Dry-run preparation is documentation, not a completed download. `docs/METHYLATION_CLOCK_VALIDATION.md` lists public URLs for an HG002 CRAM, a bisulfite coverage file, an HM450 hg38 manifest, and an Illumina chr20 BAM. `docs/TOOLCHAIN_SEPT2026.md` records samtools, bcftools, minimap2, modkit, tabix, and bgzip in a micromamba environment named `rogen-tools`, not on the default `PATH`. `pysam` and `pod5` are absent from the project environment and from system Python. `INPUT_MANIFEST.md` section "September 2026: secure intake and dry runs" still says "Filled in by later steps."

Test file: `tests/test_intake.py`. A recorded pytest pass count is not available.

Status: DONE-SYNTHETIC for intake and alias handling. BUILT-NOT-RUN for HG002 short-read and ONT dry runs.

### A.2.1.10.1 Clock audit, Hannum training, fixture split, ONT site script, doc sync

`scripts/audit_clock_artifacts.py` opens clock serializations read-only and writes `docs/CLOCK_ARTIFACT_AUDIT.md`. A MAE row is TRACEABLE only when that file records the same feature count as an artifact and either names that artifact or uses its probe-id pattern. The current report lists one CANDIDATE REAL model, `models/gse40279_hannum450k_elasticnet.joblib` (8000 Illumina probe names), and three FIXTURE models with `cg_test_*` names.

`rogen-clock prepare-gse40279`, `train`, and `evaluate` live in `src/rogen_aging/cli/clock.py`, `src/rogen_aging/clock/gse40279.py`, and `src/rogen_aging/clock/hannum.py`. Preparation writes the GEO wide table. Training filters probes on the training split only and writes `models/gse40279_hannum450k_elasticnet.provenance.json`. Evaluate refuses a model whose features are all `cg_test_*` unless `--demo` is passed. The demo clock path is `models/fixtures/fixture_clock_cg_test.pkl`.

Held-out MAE, RMSE, and Pearson r in the provenance file are from the GSE40279 train/test split. They are not external validation.

Commit `56740c1` archived the July fixture outputs and wrote `outputs/validation_metrics.json` (8000 features, 729 samples). That file is the GSE87571 evaluation of the Hannum model. `docs/GSE40279_CLOCK_TRAINING.md` now quotes `n_samples`, `mae_overall`, `median_ae`, and `pearson_r` from that file.

The age-correlation filter and the 8,000-probe cap are in `select_probes_from_parquet` (`src/rogen_aging/clock/hannum.py`). That function is what `train_published_gse40279` calls. `train_test_split` builds `train_idx` and `test_idx`. Missingness and Pearson r with age are computed on `batch[train_idx]` and `train_age` only (`_batch_train_stats`). Probes with absolute correlation at least `AGE_ABS_CORRELATION_MIN` (0.2) are kept, then `trained = correlated[:MAX_PROBES_FOR_ELASTICNET]` with `MAX_PROBES_FOR_ELASTICNET` set to 8,000. `test_idx` is counted and is not passed into the filter. The in-memory twin `select_training_probes` does the same on `frame.iloc[train_idx]`. Held-out GSE40279 metrics are the train/test split from `train_clock`. They are not optimistic from probe selection on test rows.

`outputs/validation_metrics.json` stores `mae_by_decade`. It does not store mean signed error (predicted minus chronological). `evaluate_external` computes that residual in memory and writes the JSON and three PNG figures. No predictions file is under `outputs/`.

`scripts/clock/ont_clock_site_coverage.py` defines an ONT versus bisulfite comparison at chr20 clock CpGs. `docs/METHYLATION_CLOCK_VALIDATION.md` says the method is defined, the ONT run is still pending, and nothing has been downloaded. The CLI tools are installed in micromamba `rogen-tools` (`docs/TOOLCHAIN_SEPT2026.md`).

Clock docs updated today: `docs/CLOCK_LIBRARY.md`, `docs/GSE40279_CLOCK_TRAINING.md`, `docs/METHYLATION_CLOCK_VALIDATION.md`, `docs/WORKFLOWS.md`, `INPUT_MANIFEST.md`, `README.md`.

Test files: `tests/test_audit_clock_artifacts.py`, `tests/test_gse40279_clock.py`, `tests/test_ont_clock_site_coverage.py`. Full-suite `uv run pytest -q` on 2026-09-27: 179 passed, 0 failed, 2 skipped.

Status: Hannum training DONE. GSE87571 evaluation DONE. ONT comparison BUILT-NOT-RUN.

### A.2.1.11.1 Release guard, longevity-variant panel, sanity check

`rogen-release-check` (`src/rogen_aging/cli/release_check.py`, `src/rogen_aging/privacy/release_guard.py`) blocks individual-level columns, direct identifiers, small allele-count cells, and recoverable strata cells. Summary tables skip the cell-count checks. The project minimum cell count in the release-check JSON is 5.

`scripts/panel/ro_longevity_panel_power.py` reads the 47 prioritized variants and writes `analysis/panel/ro_longevity_panel.csv`, `analysis/panel/ro_longevity_power_grid.csv`, a heatmap, and `analysis/panel/ro_longevity_panel_provenance.json`. Frequencies are gnomAD v4 NFE from a local cache and 1000 Genomes phase 3 TSI, IBS, CEU, GBR, and FIN from Ensembl REST. They are not averaged and are not labeled as Romanian frequencies. Arm sizes and odds ratios in the power grid are illustrative. The provenance disclaimer says ROGEN sample size and design are not fixed.

`tsi_nfe_maf_abs_diff_gt_0_05_descriptive` is an absolute MAF gap above 0.05. The provenance text says it is not a test. The separate count test uses a two-sided Fisher exact test, or Pearson chi-square when every expected count is at least 5, with Benjamini-Hochberg q-values. `tsi_differs_fdr05` is that test at q below 0.05. It compares alternate-allele frequency between public reference panels. It is not a Romanian or Southeast European frequency. Three rows are true. Eleven rows have an empty p-value because the gnomAD cache entry has no `ac_nfe` or `an_nfe`. The lists are in the results table.

`scripts/panel/ro_longevity_panel_sanity.py` prints NFE, CEU, and GBR alternate-allele frequencies. It flags an absolute NFE minus CEU gap above 0.15, and a complement mismatch, except when NFE AF is between 0.4 and 0.6. In that band both the frequency and its complement are near 0.5, so a smaller distance to 1 minus CEU is not treated as a swapped allele. No saved sanity log is in the tree, so a count of flags from a real run is not available.

Both panel tables were passed to `rogen-release-check` with `--table-kind summary`. Provenance records return code 0 and status `passed` for each. The JSON reports set `release_blocked` to false.

Ensembl host: `config/default.yaml` keys `apis.ensembl_rest` and `apis.ensembl_rest_dated` are `https://rest.ensembl.org`. Provenance records release 116 from `GET /info/data`. Commit `49cc8e6` says the June 2026 dated mirror returned HTTP 503 and that variation records omit population frequencies unless `pops=1`.

`.mailmap` maps `4mitya@gmail.com` and Cursor agent addresses to Dmitri Toren `52277457+mitya-toren@users.noreply.github.com`.

Test files: `tests/test_release_guard.py`, `tests/test_ro_longevity_panel_power.py` (includes the complement-skip test). Full-suite counts are in the results table.

Status: release guard DONE on the two panel tables. Panel DONE on public reference data. TSI versus NFE count test DONE where counts were present (`n_tests` 36 in provenance). Sanity check BUILT-NOT-RUN as a saved execution (the script and unit test exist; no output file).

### Other

`87fadcf` and `9ae9094` also added `scripts/build_eqtl_coverage.py` and a section in `docs/GENOMICS_ANALYSIS.md`. The script builds an allele by eQTL Catalogue `dataset_id` matrix from git-ignored `results/` tables and refuses a 13-dataset collapse. Test file: `tests/test_eqtl_coverage.py`. This does not match the three sub-activities above.

`10e5c01` keeps the UK Biobank compliance audit and the overlap enrichment test from failing when ignored local files are present. Test files: `tests/test_ukbb_ci_compliance_audit.py`, `tests/test_overlap_enrichment.py`.

`e0a2233` merges `origin/main`. The only content conflict was one README index row. Both the artifact-audit link and the Romanian-clock link were kept. Main also added `.cursor/` to `.gitignore` and removed duplicate root docs. That merge is housekeeping, not a sub-activity.

## 3. Results

Every number below is taken from the cited file, from a row tally of a cited CSV, or from the pytest command recorded in the last row. Counts of `nfe_maf_below_0_01` and `tsi_differs_fdr05` are row tallies of `analysis/panel/ro_longevity_panel.csv`.

| Item | Value | Source |
|------|-------|--------|
| GSE40279 train samples | 524 | `models/gse40279_hannum450k_elasticnet.provenance.json` key `n_train_samples` |
| GSE40279 test samples (held-out split) | 132 | same file, key `n_test_samples` |
| Input CpGs | 470043 | same file, `probe_counts.input_cg` |
| CpGs overlapping GSE87571 probe set during training filters | 470043 | same file, `probe_counts.overlap_gse87571` |
| CpGs after missingness at most 5 percent on the train split | 470043 | same file, `probe_counts.after_missingness_le_5pct_on_train` |
| CpGs after absolute age correlation on the train split | 51165 | same file, `probe_counts.after_abs_age_correlation_on_train` |
| CpGs used to train | 8000 | same file, `probe_counts.trained` and `n_cpgs_features` |
| Missing-fraction maximum | 0.05 | same file, `filters.missing_fraction_max` |
| Absolute age-correlation minimum | 0.2 | same file, `filters.age_abs_correlation_min` |
| Max probes | 8000 | same file, `filters.max_probes` |
| Non-zero coefficients | 2314 | same file, `n_nonzero_coefficients`; also `docs/CLOCK_ARTIFACT_AUDIT.md` for the Hannum joblib |
| Held-out MAE (GSE40279 split, not external) | 3.5248637199401855 | same file, `holdout_metrics_from_train_clock.test_mae` |
| Held-out RMSE | 4.727581024169922 | same file, `holdout_metrics_from_train_clock.test_rmse` |
| Held-out Pearson r | 0.9566776752471924 | same file, `holdout_metrics_from_train_clock.test_pearson_r` |
| Training CV folds | 5 | same file, `hyperparameters.cv` |
| Chosen alpha | 0.0034306552100727577 | same file, `hyperparameters.alpha` |
| Chosen l1_ratio | 0.5 | same file, `hyperparameters.l1_ratio` |
| sklearn version recorded for that fit | 1.8.0 | same file, `sklearn_version` |
| Provenance git commit | d509f0ebe3fd1eb6a41659ec0db276de4870446d | same file, `git_commit` |
| Provenance time | 2026-09-27T08:42:19Z | same file, `timestamp_utc` |
| GSE87571 evaluation samples | 729 | `outputs/validation_metrics.json` key `n_samples` |
| GSE87571 evaluation features used | 8000 | same file, `n_features_used` |
| GSE87571 MAE | 3.760192551894116 | same file, `mae_overall` |
| GSE87571 median absolute error | 3.2548561096191406 | same file, `median_ae` |
| GSE87571 Pearson r | 0.9867447913584797 | same file, `pearson_r` |
| GSE87571 missing CpGs imputed | 0 | same file, `n_missing_imputed` |
| GSE87571 MAE, age &lt;20 | 2.0683278465270996 | same file, `mae_by_decade["<20"]` |
| GSE87571 MAE, 20-29 | 2.7305740977442543 | same file, `mae_by_decade["20-29"]` |
| GSE87571 MAE, 30-39 | 4.614984880010765 | same file, `mae_by_decade["30-39"]` |
| GSE87571 MAE, 40-49 | 4.905645311795748 | same file, `mae_by_decade["40-49"]` |
| GSE87571 MAE, 50-59 | 4.463715764095909 | same file, `mae_by_decade["50-59"]` |
| GSE87571 MAE, 60-69 | 3.806238769203104 | same file, `mae_by_decade["60-69"]` |
| GSE87571 MAE, 70-79 | 3.4909284037928425 | same file, `mae_by_decade["70-79"]` |
| GSE87571 MAE, 80-89 | 3.6742114191469937 | same file, `mae_by_decade["80-89"]` |
| GSE87571 MAE, 90+ | 2.345769246419271 | same file, `mae_by_decade["90+"]` |
| GSE87571 mean signed error (predicted minus chronological), overall and by age stratum | not available | `outputs/validation_metrics.json` has no signed-error key. No predictions file under `outputs/` |
| Audit artifacts inspected | 4 | `docs/CLOCK_ARTIFACT_AUDIT.md` Counts |
| Audit FIXTURE | 3 | same section |
| Audit CANDIDATE REAL | 1 | same section |
| Audit MAE rows | 9 | same section |
| Fixture non-zero coefficients | 3 | same report, FIXTURE rows |
| Fixture feature count | 12 | same report, FIXTURE rows |
| Panel variant rows | 47 | `analysis/panel/ro_longevity_panel_provenance.json` `input_variants.n_rows` |
| Ensembl release | 116 | same file, `ensembl.release` |
| Ensembl host | https://rest.ensembl.org | same file, `ensembl.host` |
| Ensembl cache hits | 47 | same file, `ensembl.cache_hits` |
| Ensembl new queries | 0 | same file, `ensembl.new_queries` |
| Ensembl request failures | 0 | same file, `ensembl.request_failures` |
| gnomAD cache hits | 47 | same file, `gnomad.cache_hits` |
| gnomAD new queries | 0 | same file, `gnomad.new_queries` |
| gnomAD cache entries with AC and AN | 36 | same file, `gnomad.allele_counts.entries_with_counts` |
| TSI vs NFE tests | 36 | same file, `tsi_vs_nfe.n_tests` |
| Of those, chi-square | 35 | same file, `tsi_vs_nfe.n_chi_square` |
| Of those, Fisher exact | 1 | same file, `tsi_vs_nfe.n_fisher` |
| FDR threshold for that test | 0.05 | same file, `tsi_vs_nfe.fdr` |
| Rows with `tsi_differs_fdr05` true | 3 | row tally of `analysis/panel/ro_longevity_panel.csv` |
| rs3804474 | af_tsi 0.593457943925234, af_nfe 0.4935718278367803, af_diff 0.09988611608845366, q_tsi_vs_nfe 0.04226800304311284 | same CSV |
| rs2116538 | af_tsi 0.878504672897196, af_nfe 0.7795259381249264, af_diff 0.09897873477226959, q_tsi_vs_nfe 0.00872929668450113 | same CSV |
| rs9530108 | af_tsi 0.308411214953271, af_nfe 0.1631007140881239, af_diff 0.1453105008651471, q_tsi_vs_nfe 3.169302692610224e-07 | same CSV |
| Variants not tested (empty `p_tsi_vs_nfe`) | 11: rs139170, rs1981429, rs9916344, rs41383, rs1800774, rs1801318, rs524533, rs7207422, rs155979, rs882696, rs28366003 | same CSV, rows with null `p_tsi_vs_nfe`. Each note is "TSI versus NFE allele-count test left empty; counts were not filled from another source". `data/geo/gnomad_r4_nfe_cache.json` has no `ac_nfe` or `an_nfe` on these 11 entries (`ac_nfe_matches_af` false). Provenance `gnomad.allele_counts.entries_with_counts` is 36 of 47 |
| Variants with NFE MAF below 0.01 | 0 true, 47 false | row tally of column `nfe_maf_below_0_01` in `analysis/panel/ro_longevity_panel.csv` |
| Bonferroni alpha 0.05/47 | 0.0010638297872340426 | provenance `power.alpha_bonferroni_0.05_over_47` |
| Genome-wide alpha | 5e-08 | provenance `power.alpha_genomewide` |
| Panel release check | passed, return code 0, `release_blocked` false | provenance `release_check` for `ro_longevity_panel.csv`; `analysis/panel/release_checks/ro_longevity_panel.json` |
| Power-grid release check | passed, return code 0, `release_blocked` false | provenance `release_check` for `ro_longevity_power_grid.csv`; `analysis/panel/release_checks/ro_longevity_power_grid.json` |
| Release minimum cell count | 5 | both release-check JSON files, key `min_cell_count` |
| Panel provenance time | 2026-09-27T09:21:56Z | provenance `generated_at_utc` |
| samtools version in micromamba | 1.24 | `docs/TOOLCHAIN_SEPT2026.md` |
| bcftools version | 1.24 | same file |
| minimap2 version | 2.31-r1302 | same file |
| modkit version | 0.6.4 | same file |
| tabix and bgzip | htslib 1.24 | same file |
| pysam and pod5 | not installed | same file, Found = no |
| Orientation flags from a saved sanity run | not available | `scripts/panel/ro_longevity_panel_sanity.py` prints flags and does not write a file |
| Pytest | 179 passed, 0 failed, 2 skipped | `uv run pytest -q` on 2026-09-27 (14.51s, 1 warning) |

## 4. Status

| Item | Status | Evidence | Next step |
|------|--------|----------|-----------|
| Intake | DONE-SYNTHETIC | `docs/DATA_FLOW.md`: no Romanian sequencing data; tests create the files. Code in `src/rogen_aging/intake/`. | Run on a real delivery only after the data protection review. |
| Release guard | DONE | `rogen-release-check` return code 0 on both panel CSVs in `analysis/panel/ro_longevity_panel_provenance.json`. | Keep using `--table-kind summary` for frequency tables. |
| Hannum clock training | DONE | `models/gse40279_hannum450k_elasticnet.provenance.json` timestamp 2026-09-27T08:42:19Z, GSE40279, 524 train and 132 test samples. | None for the fit itself. Refresh the audit headline if the "may exist" wording is still wanted. |
| GSE87571 external validation | DONE | `outputs/validation_metrics.json` written 2026-09-27, `n_features_used` 8000, `n_samples` 729, commit `56740c1`. `docs/GSE40279_CLOCK_TRAINING.md` quotes n, MAE, median AE, and r from that file. | None for the recorded metrics. Mean signed error was not saved. |
| Short-read HG002 dry run | BUILT-NOT-RUN | Illumina chr20 BAM URL is listed in `docs/METHYLATION_CLOCK_VALIDATION.md`. That page says nothing has been downloaded. | Download only when the dry run is started. Do not label the GIAB files as Romanian or ROGEN data. |
| ONT HG002 methylation dry run | BUILT-NOT-RUN | HG002 CRAM URL is listed in the same doc. Run pending. | Same as the short-read item. `pod5` is still not installed if POD5 metadata must be scanned. |
| ONT vs bisulfite comparison | BUILT-NOT-RUN | `scripts/clock/ont_clock_site_coverage.py` exists. Validation doc: ONT run still pending, nothing downloaded. Toolchain is in micromamba `rogen-tools`. | Run the comparison. The install is recorded; the download has not started. |
| Longevity panel | DONE | `analysis/panel/ro_longevity_panel.csv` and provenance `generated_at_utc` 2026-09-27T09:21:56Z. Public gnomAD and Ensembl frequencies. | Do not describe the power ranks as a longevity result. |
| TSI vs NFE statistical test | DONE | Provenance `tsi_vs_nfe.n_tests` 36, `n_chi_square` 35, `n_fisher` 1, `fdr` 0.05. The column `tsi_nfe_maf_abs_diff_gt_0_05_descriptive` is not this test. | Keep the descriptive gap column separate from `tsi_differs_fdr05` in any write-up. |
| Alias-ID support in intake | DONE-SYNTHETIC | Commit `ddd9ac7`. `docs/DATA_FLOW.md` names HG002 with aliases NA24385 and GM24385 as the leak-check example. Tests are synthetic. | Confirm on a real header only when a delivery exists. |

## 5. Decisions and problems

The files named `gse40279` that existed before this training were the demo clock. `docs/CLOCK_ARTIFACT_AUDIT.md` flags `analysis/gse40279_elasticnet_clock.pkl`, `models/fixtures/fixture_clock_cg_test.pkl`, and `models/methylation_clock_v1.joblib` as FIXTURE, each with 12 `cg_test_*` features and 3 non-zero coefficients. Evaluate now refuses that probe set unless `--demo` is passed. The Hannum joblib is the only CANDIDATE REAL artifact in the audit (8000 `cg` plus eight digits, 2314 non-zero coefficients).

The June 2026 Ensembl mirror returned HTTP 503 on 2026-09-27 (commit `49cc8e6`). Requests go to `https://rest.ensembl.org`. Provenance records release 116. Variation calls use `pops=1` because population frequencies are otherwise omitted. The script stops if `GET /info/data` does not report release 116.

`outputs/clock_metrics.json` is not in the tree today, so it is not a result from this date. `outputs/archive/fixture_2026-07/README.md` says that archived metrics file was produced 2026-07-27/28 by `models/ro_clock_elasticnet_gse40279.pkl` and is not a GSE40279 clock result. The archived JSON `model_path` is that fixture pickle, `n_features_used` is 12, and `n_samples` is 40. The audit report still lists a row for `outputs/clock_metrics.json` with MAE 3.535416311865679 marked UNTRACEABLE. Treat that row as the stale fixture file, not as the Hannum evaluation. The Hannum GSE87571 numbers are in `outputs/validation_metrics.json`.

The allele-orientation check is `scripts/panel/ro_longevity_panel_sanity.py`. A gap above 0.15 between NFE and CEU alternate-allele frequency is flagged. The complement check is skipped when NFE AF is between 0.4 and 0.6, because the frequency and its complement are both near 0.5 and a smaller distance to the complement would be a false swap flag. A count of such flags from a saved run is not available. The unit test checks that AF 0.52 against CEU 0.40 returns no flag.

This morning samtools, bcftools, minimap2, modkit, tabix, and bgzip were absent on the default `PATH`. `docs/TOOLCHAIN_SEPT2026.md` records them in micromamba `rogen-tools` (versions in the results table). `pysam` and `pod5` are still missing. `docs/METHYLATION_CLOCK_VALIDATION.md` now says the same: the tools are in `rogen-tools`, and the ONT run itself is still pending.

The age-correlation filter and the 8,000-probe cap use the training split only. In `src/rogen_aging/clock/hannum.py`, `AGE_ABS_CORRELATION_MIN` is 0.2 and `MAX_PROBES_FOR_ELASTICNET` is 8,000. `select_probes_from_parquet` (called by `train_published_gse40279`) passes `batch[train_idx]` into `_batch_train_stats` and then takes `correlated[:MAX_PROBES_FOR_ELASTICNET]`. `select_training_probes` uses `frame.iloc[train_idx]` the same way. Test rows are split off and are not used to rank or cap probes. Held-out GSE40279 MAE, RMSE, and Pearson r stay the recorded split metrics. They are not marked optimistic.

`.mailmap` (commit `6ff6132`) maps `4mitya@gmail.com` and the Cursor agent addresses onto Dmitri Toren `52277457+mitya-toren@users.noreply.github.com`.

## 6. Open items for tomorrow

1. Run `scripts/clock/ont_clock_site_coverage.py` on the listed public HG002 and bisulfite URLs, or leave a one-line note if the download is still blocked. Do not call those files Romanian, ROGEN, or cohort data. The CLI tools are already in micromamba `rogen-tools`.
2. Fill the `INPUT_MANIFEST.md` section "September 2026: secure intake and dry runs", which still says "Filled in by later steps."
3. Decide whether to add `pysam` and `pod5`. Until `pod5` is installed, intake will not scan POD5 metadata.
4. Send `docs/DATA_FLOW.md` for data-protection review. The file says it has not been reviewed.
5. If a sanity log is required, run `scripts/panel/ro_longevity_panel_sanity.py` and keep the output. There is no saved flag count today.
6. The eQTL coverage script is committed and unmentioned in the September manifest section. Say whether that matrix is in scope for the next write-up.
7. Mean signed error on GSE87571 is not in `outputs/validation_metrics.json` and no predictions file was saved. A later evaluation can write it. Do not recompute it from a new model fit if the point is to keep the 2026-09-27 numbers.

## 7. Provenance

Commits on `feature/sept-secure-pipeline` with author dates on 2026-09-27, from `git log --date=iso`:

| Hash | Time | Message |
|------|------|---------|
| 56740c1 | 2026-09-27 12:58:59 +0300 | Archive fixture clock outputs; evaluate Hannum clock on GSE87571 |
| 6ff6132 | 2026-09-27 12:46:27 +0300 | Add longevity-variant priority panel with reference frequencies, power grid, and sanity check |
| 211352c | 2026-09-27 12:41:50 +0300 | Point the clock docs at the trained Hannum model and record that the ONT comparison has not run. |
| e0a2233 | 2026-09-27 12:25:04 +0300 | Merge remote-tracking branch 'origin/main' into feature/sept-secure-pipeline |
| 3921f3b | 2026-09-27 12:17:42 +0300 | Record the micromamba environment that supplies the September sequencing tools. |
| b4ce416 | 2026-09-27 12:17:40 +0300 | Compare ONT and bisulfite methylation at chr20 clock CpGs. |
| 4b6ce99 | 2026-09-27 12:17:37 +0300 | Rank the 47 prioritized variants by illustrative European-frequency power. |
| 49cc8e6 | 2026-09-27 12:17:33 +0300 | Request Ensembl population frequencies from the live REST host. |
| 55da420 | 2026-09-27 12:07:04 +0300 | Let summary tables pass the release guard without small-cell suppression. |
| ddd9ac7 | 2026-09-27 12:06:59 +0300 | Scrub sample aliases during intake so a header id that is not the primary sample id still fails the leak check. |
| b9e3e69 | 2026-09-27 11:47:23 +0300 | Add a traceable GSE40279 clock and keep the demo fixture out of evaluation. |
| d509f0e | 2026-09-27 11:04:32 +0300 | Add output release guard with small-cell and complementary suppression. |
| 10e5c01 | 2026-09-27 10:42:25 +0300 | Keep compliance and overlap checks stable when ignored local data is present. |
| 9ae9094 | 2026-09-27 10:41:39 +0300 | Document the clock artifact audit and the eQTL coverage matrix. |
| 87fadcf | 2026-09-27 10:33:02 +0300 | Make the clock audit follow inspected artifacts and lock eQTL coverage checks. |
| f70cd6b | 2026-09-27 10:28:37 +0300 | Add secure intake: checksum verification, keyed pseudonymization, leak check, data-flow document. |
| 2565cba | 2026-09-27 10:07:33 +0300 | Prepare September branch: toolchain inventory and ignore rules for intake data. |
