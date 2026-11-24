# Clock transfer check (ONT vs bisulfite, same sample)

Agreement between **ONT** and **bisulfite** methylation inputs on the same
public sample (HG002 / NA24385 / GM24385). This is **not** accuracy against
chronological age. A single sample cannot support an age-accuracy claim.
Do not report MAE, Pearson *r* vs age, or similar from this check.

## Status

**Not run.** Genome-wide ONT bedMethyl for this sample is not in the tree.

| Input | Path | Present |
|-------|------|---------|
| Chr20 metrics (thresholds only) | `analysis/dryrun/hg002_ont_methylation_metrics.json` | yes |
| Bisulfite (genome-wide, same cell line) | `data/dryrun/ont_vs_bisulfite/CpG.gz.bismark.cov.gz` | yes |
| ONT bedMethyl (genome-wide) | *(required CLI `--ont-bedmethyl`)* | **no** |
| Chr20-only ONT pileup (not valid here) | `data/dryrun/ont_vs_bisulfite/modkit_cpg_chr20.bed` | yes — **refused** as a substitute |

The October report already marked genome-wide ONT clock sites and this transfer
check as not run. Romanian sequencing data are not used.

## Thresholds (no new tuning)

Read from `analysis/dryrun/hg002_ont_methylation_metrics.json` before any
comparison. Every `levels[*].min_depth` in that file is `10`. The script
asserts they agree and uses that value. Do not change depth cutoffs here.

Pearson *r* is reported only when **n ≥ 30**. If **n < 30**, the script writes
`pearson_r: "not estimated"` and still reports `n`.

## What is compared

1. **CpG fractions** (ONT `(5mC + 5hmC) / valid` vs bisulfite strand-merged
   fraction), at:
   - all shared CpGs
   - HM450 sites (Zhou HM450 hg38 manifest)
   - sites of each clock: `horvath2013`, `hannum`, `rogen_hannum_retrain`
2. **Clock outputs** on the same sample: age from the ONT beta vector versus
   age from the bisulfite beta vector for those three clocks (pyaging
   **0.5.7** for Horvath and Hannum; retrained ElasticNet joblib for
   `rogen_hannum_retrain`).

pyaging fills missing CpGs (`df_to_adata(..., imputer_strategy="knn")`). The
metrics file must report, for each clock and each input, the fraction of
sites missing before imputation. Do not omit that.

## How to run

```bash
uv run python scripts/clock/ont_vs_bisulfite_genomewide.py \
  --ont-bedmethyl /path/to/genomewide_HG002.bedMethyl \
  --bisulfite-cov data/dryrun/ont_vs_bisulfite/CpG.gz.bismark.cov.gz
```

`--ont-bedmethyl` is required. Passing the chr20 window bed fails with an
assertion. Outputs (when the run succeeds):

| File | Role |
|------|------|
| `analysis/clock/ont_vs_bisulfite_genomewide_metrics.json` | Levels, Pearson *n*/r, transfer ages, missing fractions |
| `analysis/clock/ont_vs_bisulfite_genomewide_provenance.json` | Tool versions, input sha256, git commit, timestamp |

## Results

No genome-wide metrics yet. Fill this section only from
`analysis/clock/ont_vs_bisulfite_genomewide_metrics.json` after a successful
run. Until then, do not invent numbers.

## Related

- Chr20 dry run: [METHYLATION_CLOCK_VALIDATION.md](METHYLATION_CLOCK_VALIDATION.md)
- Helpers / URLs: `scripts/clock/ont_clock_site_coverage.py`
- Published-clock benchmark (age accuracy on cohorts): [CLOCK_BENCHMARK.md](CLOCK_BENCHMARK.md)
