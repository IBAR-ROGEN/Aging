# rogen_aging

Project scaffold for genomic notebooks and analysis, managed with `uv`.

## Quickstart

1) Install uv (one-time):

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

2) Create the environment and install deps:

```bash
uv sync
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

## Layout

- `src/rogen_aging/` — Python package for shared code
- `scripts/` — CLI scripts (AlphaGenome, mock data generator, security hook)
- `notebooks/` — Genomic analysis notebooks
- `docs/` — Project documentation
- `test_data/` — Synthetic test data (versioned)
- `data/` — Large/local data (git-ignored)

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

See **[docs/SYNTHETIC_UKB_GENERATOR.md](docs/SYNTHETIC_UKB_GENERATOR.md)** for usage and options.

## EDA Mock Integration

Exploratory analysis for mock epigenetic aging data (chronological vs epigenetic age, EAA residuals).

```bash
uv run python scripts/eda_mock_integration.py
# Input: test_data/mock_epigenetic_clinical.csv
# Output: results/mock_eaa_plot.png
```

See **[docs/EDA_MOCK_INTEGRATION.md](docs/EDA_MOCK_INTEGRATION.md)** for details.

## Documentation

| Document | Description |
|----------|-------------|
| [docs/PROJECT_STRUCTURE.md](docs/PROJECT_STRUCTURE.md) | Bioinformatics project directory layout |
| [docs/UKB_PRE_COMMIT_HOOK.md](docs/UKB_PRE_COMMIT_HOOK.md) | Git pre-commit security hook |
| [docs/SYNTHETIC_UKB_GENERATOR.md](docs/SYNTHETIC_UKB_GENERATOR.md) | Mock UK Biobank data generator |
| [docs/EDA_MOCK_INTEGRATION.md](docs/EDA_MOCK_INTEGRATION.md) | EDA mock epigenetic aging script |
| [docs/UKB_COMPLIANCE_AUDITOR.md](docs/UKB_COMPLIANCE_AUDITOR.md) | UK Biobank compliance tool |
| [docs/CODE_MODULES_REFERENCE.md](docs/CODE_MODULES_REFERENCE.md) | Code modules reference |
