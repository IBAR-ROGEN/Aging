# Bioinformatics Project Directory Structure

**Project:** IBAR-ROGEN Aging  
**Related:** [Code Modules Reference](CODE_MODULES_REFERENCE.md)

## Overview

This document describes the directory layout of the rogen_aging bioinformatics project. The structure separates source code, analysis, data, and documentation to support reproducible workflows while keeping sensitive data out of version control. Small **Vite + React** apps under `components/` and `frontend/` support manuscript figures alongside Python render scripts.

## Top-Level Structure

```
rogen_aging/
├── src/rogen_aging/          # Installable Python package (shared code)
├── tests/                    # Pytest package smoke / unit tests
├── scripts/                  # CLI and utility scripts
├── notebooks/                # Jupyter notebooks by functional area
├── docs/                     # Documentation
├── components/               # React/TSX manuscript mockups + Vite capture app
├── frontend/                 # Vite + React longevity network diagram + capture
├── analysis/                 # Generated figures and reports (some committed)
├── test_data/                # Small synthetic/test datasets (versioned)
├── data/                     # Large/local data (git-ignored)
├── results/                  # Pipeline outputs (git-ignored)
├── outputs/                  # Other outputs (git-ignored)
├── repo_structure.txt        # Optional snapshot of tracked paths (git ls-files)
├── pyproject.toml            # Project metadata, deps, build (uv / setuptools)
├── setup.py                  # Setuptools shim (reads pyproject.toml)
├── requirements.txt          # Optional pip-style pin list for non-uv workflows
├── .env.example              # Example environment variables
└── README.md
```

## Directory Descriptions

### `src/rogen_aging/`

Installable package (`rogen-aging` on the environment path after `uv sync` or `uv pip install -e .`). Shared analysis logic and visualizations.

| Path | Purpose |
|------|---------|
| `__init__.py` | Public API: re-exports visualization helpers and submodules (`methylation_visualizations`, `network_visualizer`, `pipeline`) |
| `pipeline/` | Placeholder subpackage for shared pipeline steps (grow as scripts move here) |
| `methylation_visualizations.py` | Pipeline diagrams, DMR plots, clock validation figures |
| `network_visualizer.py` | Protein interaction network visualization |
| `eda_dashboard/` | Streamlit EDA app for merged multi-omics Parquet ([EDA_DASHBOARD.md](EDA_DASHBOARD.md)) |

### `tests/`

Pytest tests. Run with `uv run pytest` (install dev extras first: `uv sync --extra dev`). `pyproject.toml` sets `pythonpath = ["src"]` so imports resolve without a manual install in many setups.

### `scripts/`

Executable scripts and shell utilities. Run with `uv run scripts/<script>.py` or `./scripts/<script>.sh`.

| Script | Purpose |
|--------|---------|
| `mock_ukb_generator.py` | Synthetic UK Biobank-style tabular data |
| `ukb_la_snp_lookup.py` | Offline UKB genotype manifest: Excel overlap → Ensembl GRCh38 → CSV (no dx-toolkit) |
| `render_longevity_network_diagram.py` | Matplotlib twin of `frontend` longevity network TSX |
| `render_figure1c_mechanisms_network.py` | Figure 1C mechanisms network (networkx + matplotlib; PNG/PDF) |
| `render_dashboard_figure_mockup.py` | Matplotlib twin of `components/DashboardFigureMockup.tsx` |
| `bootstrap_r_env.sh` | Optional micromamba R base under `.r-env/` |
| `security_check.sh` | UK Biobank pre-commit security scan |
| `install_pre_commit_hook.sh` | Install Git pre-commit hook |
| `alphagenome_sequence_comparer.py` | AlphaGenome API batch submission |
| `analyze_alphagenome_results.py` | Process AlphaGenome outputs |
| `visualize_alphagenome_results.py` | AlphaGenome visualizations |
| `generate_*.py` | Figure generation (pipeline, heatmaps, agent schema, etc.) |
| `validate_clock.py` | Held-out epigenetic clock validation (see ROMANIAN_EPIGENETIC_CLOCK.md) |

### `notebooks/`

Jupyter notebooks grouped by analysis type. Run with `uv run jupyter lab`.

| Folder | Purpose |
|--------|---------|
| `01_genomics_analysis/` | AlphaGenome, gene lists, network analysis |
| `02_methylation_pipeline/` | Methylation downstream analysis, clocks |
| `03_validation_and_compliance/` | UKB compliance auditor, validations |
| `04_exploratory_visualizations/` | Publication figures, exploratory plots |

### `docs/`

Project documentation (Markdown).

| Document | Purpose |
|----------|---------|
| `PROJECT_STRUCTURE.md` | This file — directory layout |
| `CODE_MODULES_REFERENCE.md` | Code-level reference |
| `UKB_PRE_COMMIT_HOOK.md` | Pre-commit security hook |
| `UKB_COMPLIANCE_AUDITOR.md` | UKB compliance tool |
| `SYNTHETIC_UKB_GENERATOR.md` | Mock data generator |
| `METHYLATION_PIPELINE_*.md` | Methylation pipeline usage |
| `ALPHAGENOME_ANALYSIS_EXPLANATION.md` | AlphaGenome methodology |

### `components/` and `frontend/`

TypeScript/React **manuscript figure** utilities (not the Streamlit EDA app). Each folder is a small **Vite** project: install with `npm install`, develop with `npm run dev`, and (where provided) export PNG via Playwright (`npm run capture` or `node scripts/capture*.mjs`). Python scripts under `scripts/render_*.py` produce matplotlib equivalents for CI and paper workflows without Node.

| Path | Purpose |
|------|---------|
| `components/DashboardFigureMockup.tsx` | Wide multi-panel dashboard mockup for figures |
| `components/dashboard-figure-render/` | Vite wrapper + capture script for the dashboard mockup |
| `frontend/src/components/LongevityNetworkDiagram.tsx` | Longevity conceptual network diagram |
| `frontend/scripts/capture-diagram.mjs` | Headless PNG export for the longevity diagram |

### `analysis/`

Generated figures, reports, and downstream outputs. Some publication-ready PNG/PDF assets are versioned; large or sensitive outputs should stay git-ignored.

| Subfolder / file (examples) | Purpose |
|-----------------------------|---------|
| `Figure1C_Mechanisms.*` | LA-SNP mechanisms panel from `render_figure1c_mechanisms_network.py` |
| `dashboard_figure_mockup.png` | Dashboard mockup from `render_dashboard_figure_mockup.py` |
| `methylation/` | Methylation pipeline outputs |
| `aging_related_datasets/` | Aging dataset search results |

### `test_data/`

Small synthetic datasets suitable for version control and CI.

| File | Purpose |
|------|---------|
| `mock_clinical_data.csv` | Synthetic UK Biobank-style clinical data (Sample_ID, Age, EAA, SNPs) |
| `gb-2013-14-10-r115-S3.csv` | Example/public dataset |

### `data/` (git-ignored)

Large or sensitive data. Not committed.

- `raw/` — Raw sequencing data
- `processed/` — Processed intermediate files
- `fasta/` — Reference genomes, etc.

Place Supplementary tables, longevitymap.sqlite, and similar files here.

### `results/`, `outputs/` (git-ignored)

Pipeline outputs and exported results. Kept local only.

## Data Flow Principles

1. **Sensitive data** → `data/` (git-ignored)
2. **Synthetic/mock data** → `test_data/` (versioned; safe for GitHub)
3. **Generated outputs** → `analysis/`, `results/`, `outputs/` (typically git-ignored)
4. **No real UK Biobank data** — Use `mock_ukb_generator.py` for pipeline development.

## Configuration

| File | Purpose |
|------|---------|
| `pyproject.toml` | Dependencies, Python version (≥3.12), `[build-system]` for setuptools, optional `dev` extra (pytest) |
| `setup.py` | Thin `setuptools.setup()` entry point; metadata lives in `pyproject.toml` |
| `requirements.txt` | Loose lower bounds for core bioinformatics tools; prefer `uv sync` for the full graph |
| `.env` | API keys, local paths (git-ignored; copy from `.env.example`) |
| `.gitignore` | Excludes `data/`, `results/`, `outputs/`, most `.csv`, `.vcf`, `.bed` |

---

**Last Updated:** May 1, 2026
