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

- `src/rogen_aging/`: Python package for shared code
- `notebooks/`: genomic analysis notebooks
- `data/`: put large/local data here (git-ignored)

## Python version

This project targets Python 3.12 (configured in `pyproject.toml`).

## Running the AlphaGenome Notebook

The `notebooks/AlphaGenome.ipynb` notebook performs a comprehensive analysis of gene lists for Alzheimer's and Parkinson's diseases. To run it, you'll need to set up your environment with the necessary API keys and data files.

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

Now, you can open `notebooks/AlphaGenome.ipynb` and run the cells. The notebook is configured to read the data files from the `data/` directory and the API keys from your environment.

## Methylation Calling Pipeline (ROGEN Activity 2.1.8.1)

This repository includes a complete methylation calling pipeline for Oxford Nanopore sequencing data.

### Quick Links

- **[METHYLATION_PIPELINE_README.md](METHYLATION_PIPELINE_README.md)** - Comprehensive user guide
- **[docs/METHYLATION_PIPELINE_USAGE.md](docs/METHYLATION_PIPELINE_USAGE.md)** - Detailed step-by-step usage guide
- **Pipeline Scripts:**
  - `pipeline_validation.sh` - Basecalling and methylation extraction
  - `downstream_analysis.R` - DMR calling and analysis
  - `notebooks/DownstreamMethylationAnalysis.ipynb` - Interactive R notebook

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
# Open: notebooks/DownstreamMethylationAnalysis.ipynb
```

For detailed instructions, see [METHYLATION_PIPELINE_README.md](METHYLATION_PIPELINE_README.md).
