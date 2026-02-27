# Bioinformatics Project Directory Structure

**Project:** IBAR-ROGEN Aging  
**Related:** [Code Modules Reference](CODE_MODULES_REFERENCE.md)

## Overview

This document describes the directory layout of the rogen_aging bioinformatics project. The structure separates source code, analysis, data, and documentation to support reproducible workflows while keeping sensitive data out of version control.

## Top-Level Structure

```
rogen_aging/
├── src/rogen_aging/          # Python package (shared code)
├── scripts/                  # CLI and utility scripts
├── notebooks/                # Jupyter notebooks by functional area
├── docs/                     # Documentation
├── analysis/                 # Generated figures and reports
├── test_data/                # Small synthetic/test datasets (versioned)
├── data/                     # Large/local data (git-ignored)
├── results/                  # Pipeline outputs (git-ignored)
├── outputs/                  # Other outputs (git-ignored)
├── pyproject.toml            # Project config (uv)
├── .env.example              # Example environment variables
└── README.md
```

## Directory Descriptions

### `src/rogen_aging/`

Core Python package for shared analysis logic and visualizations.

| Module | Purpose |
|--------|---------|
| `methylation_visualizations.py` | Pipeline diagrams, DMR plots, clock validation figures |
| `network_visualizer.py` | Protein interaction network visualization |

### `scripts/`

Executable scripts and shell utilities. Run with `uv run scripts/<script>.py` or `./scripts/<script>.sh`.

| Script | Purpose |
|--------|---------|
| `mock_ukb_generator.py` | Synthetic UK Biobank-style tabular data |
| `security_check.sh` | UK Biobank pre-commit security scan |
| `install_pre_commit_hook.sh` | Install Git pre-commit hook |
| `alphagenome_sequence_comparer.py` | AlphaGenome API batch submission |
| `analyze_alphagenome_results.py` | Process AlphaGenome outputs |
| `visualize_alphagenome_results.py` | AlphaGenome visualizations |
| `generate_*.py` | Figure generation (pipeline, heatmaps, etc.) |

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

### `analysis/`

Generated figures, reports, and downstream outputs. May be git-ignored for large artifacts.

| Subfolder | Purpose |
|-----------|---------|
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
| `pyproject.toml` | Dependencies, Python version (≥3.12) |
| `.env` | API keys, local paths (git-ignored; copy from `.env.example`) |
| `.gitignore` | Excludes `data/`, `results/`, `outputs/`, most `.csv`, `.vcf`, `.bed` |

---

**Last Updated:** February 27, 2026
