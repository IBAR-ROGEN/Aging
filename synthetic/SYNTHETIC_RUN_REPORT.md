# SYNTHETIC — clock harness synthetic preflight

This file is a pass/fail record of how the existing evaluator treats
failure modes that the real external methylation matrix is expected to
have. It is not a clock validation result.

## Failure modes

| Failure mode | Outcome | Missing-CpG behaviour | What the harness did |
|--------------|---------|----------------------|----------------------|
| (a) drop 5% of clock CpGs entirely | handled | mean-imputed | Scored 200 samples; mean-imputed 1 absent clock CpG(s) ['cg_test_0004'] (training SimpleImputer column means; not dropped, not errored). |
| (b) set 2% of values to NaN | handled | n/a (none missing) | Scored 200 samples. No clock CpG was dropped. NaN cells were mean-imputed at predict time by the training SimpleImputer (not dropped, not errored). |
| (c) 20 sample IDs absent from metadata | handled | n/a (none missing) | Inner join on sample_id dropped matrix rows absent from metadata; scored 180 of 200 samples. No clock CpG was missing. |
| combined (a)+(b)+(c) | handled | mean-imputed | Scored 180 samples; mean-imputed 1 absent clock CpG(s) ['cg_test_0000'] (training SimpleImputer column means; not dropped, not errored). |

## Missing CpG imputation (will bias the real run)

The evaluator does **not drop** clock CpGs that are absent from the
beta matrix, and it does **not error**. It **mean-imputes** them.

Named behaviour: **training-column-mean imputation** via the fitted
`sklearn.impute.SimpleImputer(strategy='mean')` step on the clock
pipeline. `rogen_aging.clock.evaluate.build_feature_matrix` inserts a
column of NaN for each missing probe and leaves those NaNs for the
pipeline imputer, which fills them from training `statistics_`
(the per-CpG mean on the training matrix). This is **not** a test-set
mean and **not** a row mean. On a real external dataset, every missing
clock probe is replaced by the training-cohort mean methylation for
that site, which shrinks predicted age toward the training intercept.

Cell-level NaNs in CpGs that *are* present are handled the same way:
left as NaN and filled by the same training-mean SimpleImputer.

Sample IDs present in the matrix but absent from metadata are dropped
by an inner join on `sample_id` (`load_validation_cohort`). They are
not scored and they do not crash the run when any IDs still overlap.

## Alignment diagnostics (synthetic data, not a performance estimate)

- Rows after ID join: 180
- Clock CpGs mean-imputed because they were absent: ['cg_test_0000']
- NaN cells left in the aligned matrix for the training imputer: 221

Generating parameters: `/Users/mityatoren/code/rogen_aging/synthetic/GENERATION_PARAMS.json`.
