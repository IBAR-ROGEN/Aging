# SYNTHETIC — clock harness synthetic preflight

## 1. Purpose

This file is a pass/fail record of how the existing evaluator treats
failure modes that the real external methylation matrix is expected to
have. It is not a clock validation result.

## 2. Preflight artefact

The run used `models/ro_clock_elasticnet_gse40279.pkl`. Despite that
filename, the artefact is a **fixture**, not a GSE40279-trained 450K
clock: its features are `cg_test_0000` … `cg_test_0011` (12 sites), of
which 3 coefficients are non-zero. Pass/fail findings below apply to
this pipeline object. They do not describe a production 450K/EPIC clock.

## 3. Failure modes

| Failure mode | Outcome | Missing-CpG behaviour | What the harness did |
|--------------|---------|----------------------|----------------------|
| (a) drop 5% of clock CpGs entirely | handled | mean-imputed | Scored 200 samples; mean-imputed 1 absent clock CpG(s) ['cg_test_0004'] (training SimpleImputer column means; not dropped, not errored). |
| (b) set 2% of values to NaN | handled | n/a (none missing) | Scored 200 samples. No clock CpG was dropped. NaN cells were mean-imputed at predict time by the training SimpleImputer (not dropped, not errored). |
| (c) 20 sample IDs absent from metadata | handled | n/a (none missing) | Inner join on sample_id dropped matrix rows absent from metadata; scored 180 of 200 samples. No clock CpG was missing. |
| combined (a)+(b)+(c) | handled | mean-imputed | Scored 180 samples; mean-imputed 1 absent clock CpG(s) ['cg_test_0000'] (training SimpleImputer column means; not dropped, not errored). |

## 4. Unmatched sample IDs

Sample IDs present in the matrix but absent from metadata are dropped
by an inner join on `sample_id` (`load_validation_cohort`). They are
not scored and they do not crash the run when any IDs still overlap.

## 5. Missing CpG imputation

The evaluator does **not drop** clock CpGs that are absent from the
beta matrix, and it does **not error**. It **mean-imputes** them.

Named behaviour: **training-column-mean imputation** via the fitted
`sklearn.impute.SimpleImputer(strategy='mean')` step on the clock
pipeline. `rogen_aging.clock.evaluate.build_feature_matrix` inserts a
column of NaN for each missing probe and leaves those NaNs for the
pipeline imputer, which fills them from training `statistics_`
(the per-CpG mean on the training matrix). This is **not** a test-set
mean and **not** a row mean.

That filling rule is a property of the pipeline object, so the
**mechanism** is established for any later input the harness scores:
a missing clock probe is replaced by that site's training-cohort mean.
This run did not measure the **magnitude** of the resulting bias, and
a 12-feature fixture cannot estimate that magnitude for a full clock.
Whether predicted age would shift, and by how much, remains unmeasured.

## 6. Cell-level NaNs

Cell-level NaNs in CpGs that *are* present are handled the same way:
left as NaN and filled by the same training-mean SimpleImputer. The
mechanism is the same as in section 5; the size of any effect on
predicted age was not measured here.

## 7. Alignment diagnostics (synthetic data, not a performance estimate)

- Rows after ID join: 180
- Clock CpGs mean-imputed because they were absent: ['cg_test_0000']
- NaN cells left in the aligned matrix for the training imputer: 221

## 8. Generating parameters

Generating parameters: `/Users/mityatoren/code/rogen_aging/synthetic/GENERATION_PARAMS.json`.

## 9. Implications for a later matrix

The **mechanism** in section 5 will apply if a real 450K/EPIC matrix is
scored with this same pipeline object: absent clock CpGs are
mean-imputed from fitted training `statistics_`, not dropped and not
errored. That is a property of the estimator, not of this synthetic
table.

The **magnitude** of any bias that fill would induce on a full clock
is unmeasured. Nothing in this run estimated it. The fixture has 12
features, three of them non-zero; it cannot stand in for a 450K-trained
model, and no figure from this session should be read as the size of
the shift on a later external cohort.
