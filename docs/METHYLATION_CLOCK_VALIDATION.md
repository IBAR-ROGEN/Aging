# Methylation clock validation (`evaluate_methylation_clock.py`)

**Project:** IBAR-ROGEN Aging  
**Activity:** 2.1.10.1 — methylation aging clock (GSE40279 train, GSE87571 validate)  
**Script:** [`scripts/clock/evaluate_methylation_clock.py`](../scripts/clock/evaluate_methylation_clock.py) (deprecated shim: [`evaluate_methylation_clock.py`](../evaluate_methylation_clock.py))  
**Library:** `rogen_aging.clock` (`load_model`, `build_feature_matrix`)  
**Related:** [CLOCK_LIBRARY.md](CLOCK_LIBRARY.md) · [CLOCK_EVAL_FIGURES.md](CLOCK_EVAL_FIGURES.md) · [GSE40279_CLOCK_TRAINING.md](GSE40279_CLOCK_TRAINING.md) · [CLOCK_ARTIFACT_AUDIT.md](CLOCK_ARTIFACT_AUDIT.md)

## Purpose

Publication-oriented **external validation** of a trained ElasticNet DNAm-age clock on an independent GSE87571 cohort. The script:

1. Joins a processed methylation beta matrix with phenotype ages.
2. Aligns CpGs to the model’s training features (mean-imputing missing probes).
3. Computes MAE, median absolute error, Pearson *r*, and age-stratified MAE (`<30`, `30–60`, `>60`).
4. Writes a three-panel figure: predicted-vs-chronological scatter (with 95% CI), residual plot, and top-25 CpG weights labeled by nearest gene.

This complements `rogen-clock evaluate` (decade MAE / residual PNGs) and `plot_clock_eval.py` (two-panel figure under `figures/`).

## Required inputs

| Input | Default path | Description |
|-------|--------------|-------------|
| Methylation matrix | `data/methylation/GSE87571_processed.parquet` | Wide beta table (`sample_id` + `cg*` columns; probes×samples also accepted) |
| Phenotype metadata | `data/methylation/GSE87571_meta.csv` | `sample_id` (or GEO accession) + `chronological_age` / `age` |
| Trained model | `models/gse40279_hannum450k_elasticnet.joblib` | GSE40279 pipeline (`SimpleImputer` → `ElasticNetCV`). [CLOCK_ARTIFACT_AUDIT.md](CLOCK_ARTIFACT_AUDIT.md) flags it CANDIDATE REAL |
| Probe→gene annotation (optional) | `data/methylation/HM450_probe_annotation.csv` | `IlmnID` + `UCSC_RefGene_Name` (falls back to Horvath S3 table) |

Preflight: the script reads [`INPUT_MANIFEST.md`](../INPUT_MANIFEST.md) and aborts if any required file is missing.

## Outputs

| Output | Path |
|--------|------|
| Metrics JSON | `outputs/clock_metrics.json` |
| Figure (raster) | `outputs/figures/Figure_Epigenetic_Clock_Panels.png` (300 dpi) |
| Figure (vector) | `outputs/figures/Figure_Epigenetic_Clock_Panels.pdf` |
| Markdown summary | Printed to stdout |

## CLI usage

```bash
uv sync

uv run python scripts/clock/evaluate_methylation_clock.py

# Explicit paths / options
uv run python scripts/clock/evaluate_methylation_clock.py \
  --methylation data/methylation/GSE87571_processed.parquet \
  --meta data/methylation/GSE87571_meta.csv \
  --model models/gse40279_hannum450k_elasticnet.joblib \
  --metrics-out outputs/clock_metrics.json \
  --figure-stem outputs/figures/Figure_Epigenetic_Clock_Panels \
  --annotation data/methylation/HM450_probe_annotation.csv \
  --top-n 25

# Skip INPUT_MANIFEST.md preflight when overriding paths in tests/CI
uv run python scripts/clock/evaluate_methylation_clock.py \
  --model /tmp/ro_clock.pkl \
  --methylation /tmp/meth.parquet \
  --meta /tmp/meta.csv \
  --skip-manifest-check
```

## Technical notes

- **Estimator contract:** a bare `sklearn.linear_model.ElasticNet`, or a `Pipeline` whose last step is `ElasticNet` or `ElasticNetCV`. `rogen-clock evaluate` refuses a model whose features are all `cg_test_*` unless `--demo` is passed.
- **Feature alignment:** uses `rogen_aging.clock.evaluate.build_feature_matrix` to reorder CpGs to `feature_names_in_` and mean-impute probes missing from GSE87571.
- **Strata:** middle bin is closed `[30, 60]`; empty strata report `MAE = null` in JSON / `NA` in the markdown summary.
- **Tests:** `uv run pytest tests/test_evaluate_methylation_clock.py -q`

## Metrics definition

| Metric | Definition |
|--------|------------|
| MAE | Mean \|chronological − predicted\| (years) |
| Median Absolute Error | Median \|chronological − predicted\| |
| Pearson *r* | Correlation of chronological vs predicted age |
| Age-stratified MAE | MAE within `<30`, `30–60` (inclusive), and `>60` years |
| Residual (panel B) | chronological − predicted |

## Model provenance

The production clock is `models/gse40279_hannum450k_elasticnet.joblib`. The demo clock with `cg_test_*` probe names is `models/fixtures/fixture_clock_cg_test.pkl`. Before treating any other serialization as the GSE40279 clock, regenerate the audit:

```bash
uv run python scripts/audit_clock_artifacts.py
# → docs/CLOCK_ARTIFACT_AUDIT.md
```

[CLOCK_ARTIFACT_AUDIT.md](CLOCK_ARTIFACT_AUDIT.md) is the list of what is loaded in this checkout. A MAE is TRACEABLE only when that metrics file records the same feature count as an artifact and either names that artifact or uses its probe-id pattern. Tests: `uv run pytest tests/test_audit_clock_artifacts.py -q`.

## ONT and bisulfite on HG002 chr20

Ran on region `chr20:10000000-20000000` (same window as the short-read dry run). ONT reads from `SAMPLE.haplotagged.cram` (+ `.crai`) over HTTPS against `GCA_000001405.15_GRCh38_no_alt_analysis_set`; the CRAM subset went through `rogen-intake` (aliases NA24385, GM24385). Own `modkit` 0.6.4 pileup used `--cpg --combine-strands --modified-bases 5mC 5hmC -r REF --region REGION`. Do not use `gm24385_mod_2021.09/extra_analysis/all.bam`.

Metrics: [analysis/dryrun/hg002_ont_methylation_metrics.json](../analysis/dryrun/hg002_ont_methylation_metrics.json). `MIN_DEPTH=10`. The ONT run is dated 2024-06-19 and the bisulfite file was published in the 2021.09 bucket. The samples were cultured and sequenced years apart, so part of any difference is biological.

| Level | n CpGs (both) | n with ONT depth ≥ 10 | Pearson r | mean (ONT − bisulfite) | Scatter |
|-------|---------------|------------------------|-----------|-------------------------|---------|
| All CpGs in region | 80896 | 80333 | 0.976 | 0.014 | [all_cpgs](../analysis/dryrun/hg002_ont_vs_bisulfite_all_cpgs.png) |
| HM450 probes in region | 653 | 634 | 0.988 | 0.014 | [hm450](../analysis/dryrun/hg002_ont_vs_bisulfite_hm450.png) |
| Hannum nonzero coefficients in region | 4 | 4 | 0.998 | 0.031 | [hannum_nonzero](../analysis/dryrun/hg002_ont_vs_bisulfite_hannum_nonzero.png) |

The CLI tools are installed in the micromamba environment `rogen-tools`, not on the default `PATH` ([TOOLCHAIN_SEPT2026.md](TOOLCHAIN_SEPT2026.md)): samtools 1.24, bcftools 1.24, minimap2 2.31-r1302, modkit 0.6.4, tabix and bgzip from htslib 1.24. `pysam` and `pod5` are the optional project extra `seq` (`uv sync --extra seq --extra dev`).

`scripts/clock/ont_clock_site_coverage.py` holds the URL decisions, strand merge, and comparison helpers.

Confirmed input URLs, copied from that script:

- `https://ont-open-data.s3.amazonaws.com/giab_2025.01/analysis/wf-human-variation/sup/HG002/PAW70337/output/SAMPLE.haplotagged.cram`
- `https://ont-open-data.s3.amazonaws.com/gm24385_mod_2021.09/bisulphite/cpg/CpG.gz.bismark.cov.gz`
- `https://zhouserver.research.chop.edu/InfiniumAnnotation/current/HM450/HM450.hg38.manifest.tsv.gz`
- `https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.300x_chr20.bam`

## See also

- [CLOCK_LIBRARY.md](CLOCK_LIBRARY.md) — package API and `rogen-clock` CLI
- [CLOCK_EVAL_FIGURES.md](CLOCK_EVAL_FIGURES.md) — two-panel `plot_clock_eval.py` figure
- [CLOCK_ARTIFACT_AUDIT.md](CLOCK_ARTIFACT_AUDIT.md) — which serializations are fixtures
- [ACTIVITIES.md](ACTIVITIES.md#21101--methylation-aging-clock) — activity 2.1.10.1 map

---

**Last updated:** October 28, 2026
