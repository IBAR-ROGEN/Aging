# ROGEN Aging — Activity map

Maps IBAR-ROGEN activity IDs to console entry points, source modules, outputs, and documentation.

## Summary

| Activity | Title | Primary CLI | Package / module | Documentation |
|----------|-------|-------------|------------------|---------------|
| **2.1.7.1** | AlphaGenome LA-SNP regulatory comparison | `scripts/alphagenome/alphagenome_sequence_comparer.py` | — (repo scripts) | [ALPHAGENOME_ANALYSIS_EXPLANATION.md](ALPHAGENOME_ANALYSIS_EXPLANATION.md) |
| **2.1.7.1** | LA-SNP VEP functional consequences (manuscript) | `scripts/ukb/annotate_la_snps_vep.py` | — (repo script) | [LA_SNP_VEP_ANNOTATION.md](LA_SNP_VEP_ANNOTATION.md) |
| **2.1.7.1** | LA-SNP pathway network figure | `scripts/figures/generate_network_fig.py` | — | [FIGURES.md](FIGURES.md) |
| **2.1.8.1** | Methylation calling pipeline (ONT) | `scripts/dev/pipeline_validation.sh`, `scripts/dev/downstream_analysis.R` | `rogen_aging.methylation_visualizations` | [METHYLATION_PIPELINE_README.md](METHYLATION_PIPELINE_README.md) |
| **2.1.8.1** | LA-SNP manifest + public AF validation | `rogen-ukb-manifest`, `rogen-compare-af-gnomad` | `rogen_aging.ukb.manifest`, `rogen_aging.ukb.gnomad` | [LA_SNP_PUBLIC_FREQUENCY_PIPELINE.md](LA_SNP_PUBLIC_FREQUENCY_PIPELINE.md) |
| **2.1.8.1** | Synthetic UKB-RAP mock folder | `rogen-ukb-mock-rap` | `rogen_aging.ukb.mock_rap`, `rogen_aging.vcf` | [SYNTHETIC_UKB_RAP_GENERATOR.md](SYNTHETIC_UKB_RAP_GENERATOR.md) |
| **2.1.10.1** | Methylation aging clock (GSE40279 train, GSE87571 validate) | `rogen-clock`, `python -m rogen_aging.clock.external_data` | `rogen_aging.clock` | [CLOCK_LIBRARY.md](CLOCK_LIBRARY.md), [GSE40279_CLOCK_TRAINING.md](GSE40279_CLOCK_TRAINING.md) |
| **2.1.11.1** | Synthetic UKB integrative validation | `rogen-ukb-integrate` | `rogen_aging.integration` | [UKB_INTEGRATION_PIPELINE.md](UKB_INTEGRATION_PIPELINE.md) |
| — | Multi-omics EDA dashboard | `streamlit run src/rogen_aging/eda_dashboard/app.py` | `rogen_aging.eda_dashboard` | [EDA_DASHBOARD.md](EDA_DASHBOARD.md) |

---

## 2.1.7.1 — AlphaGenome LA-SNP regulatory comparison

**Title:** AlphaGenome LA-SNP regulatory comparison (sequence variant effect on RNA-seq tracks).

| | |
|---|---|
| **Console entry points** | `uv run python scripts/alphagenome/alphagenome_sequence_comparer.py` (primary) · `uv run python scripts/alphagenome/analyze_alphagenome_results.py` · `uv run python scripts/alphagenome/visualize_alphagenome_results.py` |
| **Deprecated shims** | `scripts/alphagenome_sequence_comparer.py` → forwards to `scripts/alphagenome/` |
| **Source modules** | Standalone scripts under `scripts/alphagenome/` (uses `alphagenome` package + Ensembl REST) |
| **Output files** | `analysis/alphagenome/alphagenome_comparison_results.csv` · `analysis/alphagenome/alphagenome_impact_analysis.csv` · `figures/alphagenome/alphagenome_impact_bar_plot.png` · `figures/alphagenome/alphagenome_ref_vs_alt_scatter.png` |
| **Documentation** | [ALPHAGENOME_ANALYSIS_EXPLANATION.md](ALPHAGENOME_ANALYSIS_EXPLANATION.md) |

### 2.1.7.1 — LA-SNP VEP functional consequences (manuscript)

**Title:** Inferred functional consequences table for LA-SNPs via Ensembl VEP REST (GRCh38).

| | |
|---|---|
| **Console entry points** | `uv run python scripts/ukb/annotate_la_snps_vep.py` (deprecated shim: `annotate_la_snps_vep.py` at repo root) |
| **Source modules** | `scripts/ukb/annotate_la_snps_vep.py` (Ensembl VEP REST + local JSON cache) |
| **Input files** | `analysis/alphagenome/alphagenome_impact_analysis.csv` (default `snp` column) · or plain-text / CSV rsID list |
| **Output files** | `analysis/vep_annotation/la_snp_vep_annotations.csv` · `.xlsx` · `la_snp_vep_not_found.txt` · cache under `analysis/vep_cache/` |
| **Documentation** | [LA_SNP_VEP_ANNOTATION.md](LA_SNP_VEP_ANNOTATION.md) |

### 2.1.7.1 — LA-SNP pathway network figure

**Title:** LA-SNP pathway network figure (manuscript supplementary).

| | |
|---|---|
| **Console entry points** | `uv run python scripts/figures/generate_network_fig.py` · deprecated shim: `scripts/generate_network_fig.py` |
| **Source modules** | `scripts/figures/generate_network_fig.py` (networkx + matplotlib) |
| **Output files** | `analysis/Fig_LA_SNP_network.png`, `analysis/Fig_LA_SNP_network.pdf` (committed snapshots; regenerate under `figures/` — see [FIGURES.md](FIGURES.md)) |
| **Documentation** | [FIGURES.md](FIGURES.md) |

---

## 2.1.8.1 — Methylation calling pipeline (ONT)

**Title:** Oxford Nanopore methylation calling pipeline validation.

| | |
|---|---|
| **Console entry points** | `./pipeline_validation.sh` · `Rscript downstream_analysis.R` (or source from R) |
| **Source modules** | `src/rogen_aging/methylation_visualizations.py` (downstream plots) |
| **Output files** | `basecalled_methylation.bam`, `methylation_calls.bedMethyl` (from `pipeline_validation.sh`; paths configurable in script) |
| **Documentation** | [METHYLATION_PIPELINE_README.md](METHYLATION_PIPELINE_README.md) |

### 2.1.8.1 — LA-SNP manifest + public AF validation

**Title:** LA-SNP manifest builder, 1KG extraction, gnomAD v4 NFE comparison.

| | |
|---|---|
| **Console entry points** | `uv run rogen-ukb-manifest build …` · `uv run rogen-ukb-manifest extract …` · `uv run rogen-compare-af-gnomad …` |
| **Source modules** | `src/rogen_aging/ukb/manifest.py` · `src/rogen_aging/ukb/gnomad.py` · CLI: `src/rogen_aging/cli/ukb_manifest.py`, `src/rogen_aging/cli/compare_af_gnomad.py` |
| **Output files** | `analysis/ukb_snp_manifest_v0.1.csv` (build) · `analysis/la_snp_1kg_frequencies.csv` (extract) · `analysis/la_snp_af_1kg_vs_gnomad.csv`, `figures/af_1kg_vs_gnomad_scatter.png`, `data/geo/gnomad_r4_nfe_cache.json` (compare) |
| **Documentation** | [LA_SNP_PUBLIC_FREQUENCY_PIPELINE.md](LA_SNP_PUBLIC_FREQUENCY_PIPELINE.md) |

Example:

```bash
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

### 2.1.8.1 — Synthetic UKB-RAP mock folder

**Title:** Synthetic UKB-RAP phenotype table + LA-SNP VCF for pipeline QA.

| | |
|---|---|
| **Console entry points** | `uv run rogen-ukb-mock-rap` (preferred) · `uv run python scripts/ukb/mock_rap_folder.py` (thin wrapper) |
| **Source modules** | `src/rogen_aging/ukb/mock_rap.py` · `src/rogen_aging/vcf/synthetic.py` · CLI: `src/rogen_aging/cli/ukb_mock_rap.py` |
| **Output files** | `test_data/mock_ukb_rap/phenotypes/ukb_phenotypes.csv` · `test_data/mock_ukb_rap/genotypes/ukb_la_snps.vcf` (default `--output-dir`) |
| **Documentation** | [SYNTHETIC_UKB_RAP_GENERATOR.md](SYNTHETIC_UKB_RAP_GENERATOR.md) |

---

## 2.1.10.1 — Methylation aging clock

**Title:** Epigenetic clock training on GSE40279-style wide tables and external validation on GSE87571.

| | |
|---|---|
| **Console entry points** | `uv run rogen-clock train …` · `uv run rogen-clock evaluate …` · `uv run python -m rogen_aging.clock.external_data --output …` |
| **Script equivalents** | `scripts/clock/run_clock.py` (same as `rogen-clock`) · deprecated: `scripts/clock/train_clock_on_gse40279.py`, `scripts/clock/validate_clock.py` |
| **Source modules** | `src/rogen_aging/clock/data.py` (wide-table I/O) · `model.py` (`make_clock_pipeline`) · `train.py` (`train_clock`) · `evaluate.py` (`evaluate_clock`) · `external_data.py` (`load_gse87571`, GSE87571 Parquet builder) · CLI: `src/rogen_aging/cli/clock.py` |
| **Figure script** | [`scripts/figures/plot_clock_eval.py`](../scripts/figures/plot_clock_eval.py) — two-panel external-validation figure (scatter + top CpG weights); see [CLOCK_EVAL_FIGURES.md](CLOCK_EVAL_FIGURES.md) |
| **Output files** | **Train:** `--output_model` (`.pkl`/`.joblib` pipeline), `--output_metrics` (training JSON with `test_mae`, `selected_cpgs`, …) · **External data:** `--output` Parquet (e.g. `data/gse87571.parquet`) · **Evaluate:** `{output_dir}/validation_metrics.json`, `Fig_Clock_Residuals.png`, `Fig_Clock_MAE_by_decade.png` · **Eval figure:** `figures/validation_gse87571/clock_eval_gse87571.png` + `.pdf` from `plot_clock_eval.py` |
| **Documentation** | [CLOCK_LIBRARY.md](CLOCK_LIBRARY.md) · [GSE40279_CLOCK_TRAINING.md](GSE40279_CLOCK_TRAINING.md) · [CLOCK_EVAL_FIGURES.md](CLOCK_EVAL_FIGURES.md) · [ROMANIAN_EPIGENETIC_CLOCK.md](ROMANIAN_EPIGENETIC_CLOCK.md) |

Example (GSE40279 train → GSE87571 external validate):

```bash
# Train on a GSE40279-style wide Parquet (cg* + chronological_age)
uv run rogen-clock train \
  --input_data data/gse40279_beta_age.parquet \
  --output_model analysis/gse40279_elasticnet_clock.pkl \
  --output_metrics analysis/gse40279_train_metrics.json

# Build GSE87571 evaluation cohort
uv run python -m rogen_aging.clock.external_data \
  --output data/gse87571.parquet \
  --geo-cache-dir data/geo

# External validation
uv run rogen-clock evaluate \
  --model_path analysis/gse40279_elasticnet_clock.pkl \
  --test_data data/gse87571.parquet \
  --output_dir analysis/validation_gse87571

# Publication figure (predicted vs chronological + top CpG weights)
uv run python scripts/figures/plot_clock_eval.py
```

---

## 2.1.11.1 — Synthetic UKB integrative validation

**Title:** Join mock phenotype CSV to LA-SNP VCF and run LA-SNP association scans (synthetic QA only).

| | |
|---|---|
| **Console entry points** | `uv run rogen-ukb-integrate` (preferred) · `uv run python scripts/ukb/run_integration.py` |
| **Deprecated shims** | `scripts/run_integration.py` → forwards to `scripts/ukb/run_integration.py` |
| **Source modules** | `src/rogen_aging/integration/ukb_joiner.py` (join + Fisher scan) · `src/rogen_aging/integration/run_cli.py` · CLI: `src/rogen_aging/cli/ukb_integrate.py` |
| **Output files** | `{output_dir}/assoc_la_snp_parental_longevity.csv` · `{output_dir}/assoc_la_snp_ad.csv` (default `analysis/`) |
| **Documentation** | [UKB_INTEGRATION_PIPELINE.md](UKB_INTEGRATION_PIPELINE.md) |

Example:

```bash
uv run rogen-ukb-mock-rap   # if mock inputs missing
uv run rogen-ukb-integrate \
  --pheno test_data/mock_ukb_rap/phenotypes/ukb_phenotypes.csv \
  --vcf test_data/mock_ukb_rap/genotypes/ukb_la_snps.vcf \
  --output-dir analysis/
```

---

## Script layout (by workflow)

```
scripts/
├── clock/           run_clock.py; deprecated train_clock_on_gse40279.py, validate_clock.py
├── ukb/             mock_rap_folder.py, run_integration.py (thin wrappers)
├── vcf/             synthetic Romanian cohort VCF
├── figures/         matplotlib / networkx manuscript renders
├── alphagenome/     sequence comparer, analyze, visualize
├── eda/             mock epigenetic EDA
└── dev/             security hook, CI audit, R bootstrap, utilities
```

Deprecated one-line shims at the old flat `scripts/*.py` paths forward to `scripts/<workflow>/` (including `scripts/figures/`).
