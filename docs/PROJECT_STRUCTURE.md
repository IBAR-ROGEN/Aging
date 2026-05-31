# Bioinformatics Project Directory Structure

**Project:** IBAR-ROGEN Aging  
**Navigation:** [WORKFLOWS.md](WORKFLOWS.md) · [ACTIVITIES.md](ACTIVITIES.md) · [FIGURES.md](FIGURES.md)

## Overview

The structure separates installable Python code (`src/`), grouped CLI scripts (`scripts/`), notebooks, documentation, and git-ignored data. Synthetic fixtures live in `test_data/`; real or large data in `data/`.

## Top-Level Structure

```
rogen_aging/
├── src/rogen_aging/          # Installable package
├── scripts/                  # Grouped CLIs (+ flat deprecation shims)
├── tests/
├── notebooks/
├── docs/
├── components/               # React dashboard figure mockup + Vite capture
├── frontend/                 # Longevity network diagram + Vite capture
├── analysis/                 # Committed figure exports (PNG/PDF)
├── test_data/                # Versioned synthetic fixtures
├── data/                     # Large/local data (git-ignored)
├── .github/workflows/        # CI (pytest + UKB audit)
├── pyproject.toml
└── README.md
```

## `src/rogen_aging/`

| Path | Purpose |
|------|---------|
| `clock/` | Epigenetic clock train/eval/data ([CLOCK_LIBRARY.md](CLOCK_LIBRARY.md)) |
| `ukb/` | LA-SNP manifest, gnomAD compare, mock clinical CSV ([LA_SNP_PUBLIC_FREQUENCY_PIPELINE.md](LA_SNP_PUBLIC_FREQUENCY_PIPELINE.md)) |
| `vcf/` | Synthetic VCF utilities ([SYNTHETIC_ROMANIAN_VCF_GENERATOR.md](SYNTHETIC_ROMANIAN_VCF_GENERATOR.md)) |
| `integration/` | Synthetic UKB join + associations ([UKB_INTEGRATION_PIPELINE.md](UKB_INTEGRATION_PIPELINE.md)) |
| `eda_dashboard/` | Streamlit merged-cohort EDA ([EDA_DASHBOARD.md](EDA_DASHBOARD.md)) |
| `cli/` | Console entry points (`rogen-clock`, `rogen-ukb-manifest`, `rogen-ukb-integrate`, …) |
| `methylation_visualizations.py`, `network_visualizer.py` | Shared visualization helpers |
| `pipeline/` | Placeholder for shared pipeline steps |

## `scripts/`

Grouped by workflow. Flat `scripts/*.py` paths are **deprecation shims** that forward to these folders.

| Folder | Contents |
|--------|----------|
| `clock/` | `run_clock.py` (canonical), deprecated `validate_clock.py` / `train_clock_on_gse40279.py`, Romanian demo |
| `ukb/` | `la_snp_lookup.py`, `compare_af_gnomad.py`, `mock_clinical_csv.py`, `mock_rap_folder.py`, `run_integration.py` |
| `vcf/` | `generate_synthetic_romanian_vcf.py` |
| `figures/` | `render_*`, `generate_*` manuscript figure scripts |
| `alphagenome/` | Sequence comparer + analysis + visualize |
| `eda/` | Mock epigenetic EDA |
| `dev/` | `security_check.sh`, `install_pre_commit_hook.sh`, `summarize_af_comparison.py`, R bootstrap, utilities |

Shell wrappers at `scripts/security_check.sh` (etc.) delegate to `scripts/dev/`.

## `tests/`

Run with `uv run pytest` after `uv sync --extra dev`. Imports use `rogen_aging.*` only (no `scripts/` on `pythonpath`).

## `docs/`

| Document | Purpose |
|----------|---------|
| [WORKFLOWS.md](WORKFLOWS.md) | Workflow index and console commands |
| [ACTIVITIES.md](ACTIVITIES.md) | Activity ID → code map |
| [FIGURES.md](FIGURES.md) | Manuscript figure assets |
| [CODE_MODULES_REFERENCE.md](CODE_MODULES_REFERENCE.md) | Package and CLI reference |
| [LA_SNP_PUBLIC_FREQUENCY_PIPELINE.md](LA_SNP_PUBLIC_FREQUENCY_PIPELINE.md) | LA-SNP manifest, 1KG extract, gnomAD compare + summarize |
| Per-topic guides | Clock, UKB, methylation, compliance, … |

## Data flow

1. **Sensitive data** → `data/` (git-ignored)
2. **Synthetic fixtures** → `test_data/` (versioned)
3. **Figure exports** → `analysis/` (selected PNG/PDF committed)
4. **No real UKB data** — use `rogen-ukb-mock-clinical`, `rogen-ukb-mock-rap`, `rogen-ukb-integrate` ([UKB_INTEGRATION_PIPELINE.md](UKB_INTEGRATION_PIPELINE.md))

## Configuration

| File | Purpose |
|------|---------|
| `pyproject.toml` | Dependencies, `[project.scripts]` console entry points, pytest |
| `.env` | API keys (git-ignored; copy from `.env.example`) |

---

**Last updated:** May 31, 2026
