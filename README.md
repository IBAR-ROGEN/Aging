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

The notebook requires API keys for NCBI and AlphaGenome. You should set these as environment variables:

```bash
export NCBI_API_KEY="your_ncbi_api_key"
export ALPHA_GENOME_API_KEY="your_alpha_genome_api_key"
```

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
