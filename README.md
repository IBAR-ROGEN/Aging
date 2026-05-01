# rogen_aging

Project scaffold for genomic notebooks and analysis, managed with `uv`.

## Quickstart

1) Install uv (one-time):

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

2) Create the environment, install dependencies, and install this repo as a package (editable):

```bash
uv sync
```

For running tests, include dev dependencies:

```bash
uv sync --extra dev
uv run pytest
```

3) Register a Jupyter kernel for this project (optional but recommended):

```bash
uv run python -m ipykernel install --user --name rogen-aging --display-name "Python (rogen-aging)"
```

4) Launch JupyterLab and start working in `notebooks/`:

```bash
uv run jupyter lab
```

## Common tasks

- Add a runtime dependency:

```bash
uv add <package>
```

- Add a dev-only dependency (linters, tests, etc.):

```bash
uv add --dev <package>
```

- Install optional genomics extras:

```bash
uv add .[genomics]
```

## Multi-omics EDA dashboard (Streamlit)

Interactive dashboard for the **merged integration Parquet** (clinical, methylation-derived ages, LA-SNP columns). Run from the repo root:

```bash
uv run streamlit run src/rogen_aging/eda_dashboard/app.py
```

By default the app expects `data/merged_cohort.parquet`, or set `ROGEN_MERGED_COHORT_PARQUET`. Full schema, tabs, and module map: **[docs/EDA_DASHBOARD.md](docs/EDA_DASHBOARD.md)**.

Static **manuscript figure** (matplotlib mockup of the React component in `components/DashboardFigureMockup.tsx`):

```bash
uv run python scripts/render_dashboard_figure_mockup.py
```

Writes `analysis/dashboard_figure_mockup.png`. With Node.js, `components/dashboard-figure-render` can render the TSX via Vite + Playwright (`npm install` then `npm run capture`).

## Manuscript figure renders (matplotlib / networkx)

Static figures aligned with the React/Vite mockups (no participant data):

```bash
# Longevity conceptual network (default: figures/longevity_network_diagram.png; dir is git-ignored unless you pass --output elsewhere)
uv run python scripts/render_longevity_network_diagram.py

# Figure 1C — LA-SNP mechanisms hub-and-spoke (PNG + PDF under analysis/)
uv run python scripts/render_figure1c_mechanisms_network.py

# EDA dashboard mockup (matplotlib twin of DashboardFigureMockup.tsx)
uv run python scripts/render_dashboard_figure_mockup.py

# Supplementary figure — LA-SNPs per gene (defaults: Gene / SNP_rsID columns)
uv run python scripts/generate_la_snp_per_gene_plot.py
```

Writes `analysis/Fig_LA_SNPs_per_gene.png`. If your Excel uses different headers (for example `Gene Symbol` / `SNP Identifier`), pass `--gene-column` and `--snp-column`.

See **[docs/CODE_MODULES_REFERENCE.md](docs/CODE_MODULES_REFERENCE.md)** (section 8 and manuscript figure script sections, including 3.12–3.18) for parameters and dependencies.

## UK Biobank SNP manifest (offline, Ensembl)

Build a CSV manifest from an Excel overlap table (`Gene`, `SNP_rsID`): resolve rs IDs to **GRCh38** coordinates via the Ensembl Variation REST API, add an imputed-bulk chunk label for downstream UKB genotype extraction planning (no DNAnexus or dx-toolkit calls).

```bash
uv run python scripts/ukb_la_snp_lookup.py \
  --input overlapping_genes_with_snps.xlsx \
  --output analysis/ukb_snp_manifest_v0.1.csv
```

The Git pre-commit hook exempts this script from the generic `UKB_` content scan because it intentionally emits manifest column names for extraction metadata; it must never contain real participant IDs. See **[docs/UKB_PRE_COMMIT_HOOK.md](docs/UKB_PRE_COMMIT_HOOK.md)**.

After building the CSV, run the exploratory notebook **`notebooks/05_ukb_exploration/UKB_LA_SNP_FirstContact.ipynb`** (from the repo root with `uv run jupyter lab`) to sanity-check coverage, chromosome distribution, gene-level counts, and per-chromosome GRCh38 spans before spending UK Biobank extraction credits.

## Layout

- `src/rogen_aging/` — Installable Python package (`import rogen_aging`, `from rogen_aging import …`)
- `tests/` — Pytest tests (`uv run pytest` after `uv sync --extra dev`), including `test_mock_clinical_csv.py` and `test_synthetic_vcf.py` for the synthetic generators (`pyproject.toml` adds `scripts/` to pytest’s `pythonpath`)
- `scripts/` — CLI scripts (AlphaGenome, mock tabular/VCF generators, figure renders, LA-SNP per-gene supplementary plot, UKB manifest builder, security hook)
- `notebooks/` — Genomic analysis notebooks (including `05_ukb_exploration/` for UKB manifest QA)
- `docs/` — Project documentation
- `components/` — React/TypeScript manuscript figure mockups (e.g. dashboard) and a small Vite capture app under `components/dashboard-figure-render/`
- `frontend/` — Vite + React app for the longevity network diagram (`LongevityNetworkDiagram.tsx`) and headless capture script
- `test_data/` — Synthetic test data (versioned)
- `analysis/` — Committed or local figure exports (see figure scripts below); large ad hoc outputs may remain git-ignored
- `data/` — Large/local data (git-ignored)
- `repo_structure.txt` — Optional snapshot of tracked paths (regenerate with `git ls-files` when refreshing the tree)
- `setup.py` / `pyproject.toml` — Packaging (setuptools + uv); `requirements.txt` is an optional pip-oriented pin list

## Python version

This project targets Python 3.12 (configured in `pyproject.toml`).

## Running the AlphaGenome Notebook

The `notebooks/01_genomics_analysis/AlphaGenome.ipynb` notebook performs a comprehensive analysis of gene lists for Alzheimer's and Parkinson's diseases. To run it, you'll need to set up your environment with the necessary API keys and data files.

### 1. API Keys

The notebook requires API keys for NCBI and AlphaGenome. Follow these steps to set them up:

1. Copy the `.env.example` file to `.env`:

```bash
cp .env.example .env
```

2. Edit the `.env` file and add your actual API keys:

```bash
NCBI_API_KEY=your_actual_ncbi_api_key_here
ALPHA_GENOME_API_KEY=your_actual_alpha_genome_api_key_here
```

**Important:** The `.env` file is git-ignored and will not be committed to version control, keeping your API keys secure.

### 2. Data Files

The notebook needs two data files:

- `Supplementary Table 3.xlsx`: An Excel file with gene lists.
- `longevitymap.sqlite`: A SQLite database from the Longevity Map.

You should place these files in the `data/` directory. This directory is included in `.gitignore` to prevent large data files from being committed to your repository.

### 3. Launch JupyterLab

Once you've set up your API keys and data files, you can launch JupyterLab:

```bash
uv run jupyter lab
```

Now, you can open `notebooks/01_genomics_analysis/AlphaGenome.ipynb` and run the cells. The notebook is configured to read the data files from the `data/` directory and the API keys from your environment.

## AlphaGenome Sequence Comparer (Activity 2.1.7.1)

Automated tool for batch submission to AlphaGenome API to compare regulatory scores of Reference vs. Alternate alleles for longevity-associated SNPs.

### Documentation & Analysis
- **[docs/ALPHAGENOME_ANALYSIS_EXPLANATION.md](docs/ALPHAGENOME_ANALYSIS_EXPLANATION.md)** - **Detailed methodology and logic explanation**
- **Results:**
  - `alphagenome_comparison_results.csv` - Raw API outputs
  - `alphagenome_impact_analysis.csv` - Processed regulatory impact scores
  - `alphagenome_impact_bar_plot.png` - Visualization of top impacts
  - `alphagenome_ref_vs_alt_scatter.png` - Scatter plot of Ref vs Alt scores

### Usage
```bash
# 1. Ensure API keys are set in .env
# 2. Run the sequence comparer
uv run scripts/alphagenome_sequence_comparer.py

# 3. Analyze and visualize results
uv run scripts/analyze_alphagenome_results.py
uv run scripts/visualize_alphagenome_results.py
```

## Methylation Calling Pipeline (ROGEN Activity 2.1.8.1)

This repository includes a complete methylation calling pipeline for Oxford Nanopore sequencing data.

### Quick Links

- **[METHYLATION_PIPELINE_README.md](METHYLATION_PIPELINE_README.md)** - Comprehensive user guide
- **[docs/METHYLATION_PIPELINE_USAGE.md](docs/METHYLATION_PIPELINE_USAGE.md)** - Detailed step-by-step usage guide
- **[docs/UKB_COMPLIANCE_AUDITOR.md](docs/UKB_COMPLIANCE_AUDITOR.md)** - Compliance and safety documentation
- **[docs/UKBB_CI_COMPLIANCE_AUDIT.md](docs/UKBB_CI_COMPLIANCE_AUDIT.md)** - CI/CD repository audit (`scripts/ukbb_ci_compliance_audit.sh`): genomic files, large tabular data, hardcoded paths
- **[docs/UKB_PRE_COMMIT_HOOK.md](docs/UKB_PRE_COMMIT_HOOK.md)** — Git pre-commit hook blocks commits containing `patient_id`, `UKB_`, or `.vcf`/`.bed`. Install: `./scripts/install_pre_commit_hook.sh`
- **Pipeline Scripts:**
  - `pipeline_validation.sh` - Basecalling and methylation extraction
  - `downstream_analysis.R` - DMR calling and analysis
  - `notebooks/02_methylation_pipeline/DownstreamMethylationAnalysis.ipynb` - Interactive R notebook
  - `notebooks/03_validation_and_compliance/UKB_Compliance_Auditor.ipynb` - Compliance auditing tool

### Pipeline Overview

The methylation pipeline integrates three main tools:
1. **Dorado** - Basecalling with methylation-aware models
2. **Modkit** - BAM to bedMethyl conversion
3. **DMRcaller** - Differential methylation analysis

### Quick Start

```bash
# 1. Run basecalling and methylation extraction
./pipeline_validation.sh

# 2. Run downstream analysis
Rscript downstream_analysis.R

# Or use the interactive notebook
uv run jupyter lab
# Open: notebooks/02_methylation_pipeline/DownstreamMethylationAnalysis.ipynb
```

For detailed instructions, see [METHYLATION_PIPELINE_README.md](METHYLATION_PIPELINE_README.md).

## Synthetic UK Biobank Data Generator

Generate fake UK Biobank-style tabular data for pipeline development before real data arrives. Safe for GitHub (no real participant IDs).

```bash
uv run scripts/mock_ukb_generator.py
# Output: test_data/mock_clinical_data.csv (1000 samples by default)
```

See **[docs/SYNTHETIC_UKB_GENERATOR.md](docs/SYNTHETIC_UKB_GENERATOR.md)** for columns, CLI options, and pytest coverage in **`tests/test_mock_clinical_csv.py`**.

## Synthetic Romanian cohort VCF

Generate a **streaming VCF v4.2** with Hardy–Weinberg diploid genotypes, **EUR-like** allele frequencies, **`GT:AD:DP:GQ`** per sample, and GRCh38 **chr1–chr22** contigs (bcftools-friendly sort order). Output path is your choice; large VCFs should live under **`data/`** (git-ignored).

```bash
uv run scripts/generate_synthetic_romanian_vcf.py \
  --samples 100 --variants 5000 --seed 42 \
  --output data/mock_romanian_eur.vcf
```

See **[docs/SYNTHETIC_ROMANIAN_VCF_GENERATOR.md](docs/SYNTHETIC_ROMANIAN_VCF_GENERATOR.md)** for options, headers, compliance notes, and **`tests/test_synthetic_vcf.py`** for lightweight pytest checks.

## EDA Mock Integration

Exploratory analysis for mock epigenetic aging data (chronological vs epigenetic age, EAA residuals).

```bash
uv run python scripts/eda_mock_integration.py
# Input: test_data/mock_epigenetic_clinical.csv
# Output: results/mock_eaa_plot.png
```

See **[docs/EDA_MOCK_INTEGRATION.md](docs/EDA_MOCK_INTEGRATION.md)** for details.

## Romanian cohort mock epigenetic clock (Elastic Net)

Train a **custom epigenetic aging clock** with **`ElasticNetCV`** (5-fold CV over `alpha` and `l1_ratio`), **Polars**-loaded methylation matrix + metadata (`chronological_age`), test-set **MAE** and **Pearson r**, and a **scatter plot** (chronological vs predicted age with best-fit and identity lines).

```bash
uv run python scripts/train_romanian_epigenetic_clock.py
# Mock CSVs: data/mock_romanian_cohort/ (created if missing; git-ignored)
# Plot: figures/romanian_mock_epigenetic_clock_scatter.png (git-ignored)
```

See **[docs/ROMANIAN_EPIGENETIC_CLOCK.md](docs/ROMANIAN_EPIGENETIC_CLOCK.md)** for file formats, options, and validation notes.

## Epigenetic clock validation (held-out test set)

Evaluate a **pre-trained** scikit-learn model (`.joblib` / `.pkl`) on a single table with **`chronological_age`** and CpG columns whose names start with **`cg`**. Reports overall **MAE** and **Pearson r**, **MAE by age decade** (`pandas.cut`), and saves residual / decade bar figures plus **`validation_metrics.json`**.

```bash
uv run python scripts/validate_clock.py \
  --model_path path/to/clock.joblib \
  --test_data path/to/test.parquet \
  --output_dir path/to/validation_out
```

Supported test formats: **`.parquet`**, **`.csv`** (also **`.tsv`**). Missing model CpGs (when the estimator stores `feature_names_in_`) are mean-imputed from available test `cg*` values. Details: **[docs/ROMANIAN_EPIGENETIC_CLOCK.md](docs/ROMANIAN_EPIGENETIC_CLOCK.md#held-out-validation-validate_clockpy)**.

## R environment (optional)

Bootstrap a local R toolchain under `.r-env/` (micromamba, conda-forge `r-base`):

```bash
./scripts/bootstrap_r_env.sh
```

See comments in the script for IDE integration; `.r-env/` is git-ignored.

## Documentation

| Document | Description |
|----------|-------------|
| [docs/PROJECT_STRUCTURE.md](docs/PROJECT_STRUCTURE.md) | Bioinformatics project directory layout |
| [docs/UKB_PRE_COMMIT_HOOK.md](docs/UKB_PRE_COMMIT_HOOK.md) | Git pre-commit security hook |
| [docs/SYNTHETIC_UKB_GENERATOR.md](docs/SYNTHETIC_UKB_GENERATOR.md) | Mock UK Biobank data generator |
| [docs/SYNTHETIC_ROMANIAN_VCF_GENERATOR.md](docs/SYNTHETIC_ROMANIAN_VCF_GENERATOR.md) | Synthetic Romanian cohort VCF (VCF 4.2) |
| [docs/EDA_MOCK_INTEGRATION.md](docs/EDA_MOCK_INTEGRATION.md) | EDA mock epigenetic aging script |
| [docs/ROMANIAN_EPIGENETIC_CLOCK.md](docs/ROMANIAN_EPIGENETIC_CLOCK.md) | Romanian cohort Elastic Net clock (`train_romanian_epigenetic_clock.py`) and held-out validation (`validate_clock.py`) |
| [docs/UKB_COMPLIANCE_AUDITOR.md](docs/UKB_COMPLIANCE_AUDITOR.md) | UK Biobank compliance tool |
| [docs/UKBB_CI_COMPLIANCE_AUDIT.md](docs/UKBB_CI_COMPLIANCE_AUDIT.md) | CI/CD UKB-oriented repo audit script and usage |
| [docs/CODE_MODULES_REFERENCE.md](docs/CODE_MODULES_REFERENCE.md) | Code modules reference (`ukb_la_snp_lookup.py`, figure render scripts, `components/` / `frontend/`) |
