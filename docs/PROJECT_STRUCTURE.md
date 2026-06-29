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
├── analysis/                 # Committed figure snapshots, alphagenome tables, pipeline CSVs
├── figures/                  # Local regenerated plots (git-ignored; .gitkeep tracked)
├── outputs/                  # Optional scratch for ad-hoc runs (git-ignored)
├── data/                     # Large/local data (git-ignored)
├── test_data/                # Versioned synthetic fixtures
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
| `ukb/` | `la_snp_lookup.py`, `compare_af_gnomad.py`, `mock_clinical_csv.py`, `mock_rap_folder.py`, `run_integration.py`, `annotate_la_snps_vep.py` |
| `vcf/` | `generate_synthetic_romanian_vcf.py` |
| `figures/` | `render_*`, `generate_*` manuscript figure scripts |
| `alphagenome/` | Sequence comparer + analysis + visualize |
| `eda/` | Mock epigenetic EDA |
| `dev/` | `security_check.sh`, `install_pre_commit_hook.sh`, `summarize_af_comparison.py`, R bootstrap, utilities |

Shell wrappers at `scripts/security_check.sh` (etc.) delegate to `scripts/dev/`.

## `tests/`

Run with `uv run pytest` after `uv sync --extra dev`. Imports use `rogen_aging.*` only (no `scripts/` on `pythonpath`).

| File | Role |
|------|------|
| `test_clock_regression.py` | Refactored clock vs legacy ElasticNet metrics on `test_data/mock_clock_wide.csv` |
| `test_package_imports.py` | Smoke imports for installable subpackages |
| `test_ukb_integration.py` / `test_ukb_mock_gen.py` | Synthetic UKB join and mock RAP layout |
| `test_synthetic_vcf.py` / `test_mock_clinical_csv.py` | VCF generator and mock clinical CSV |
| `test_af_comparison_summary.py` | gnomAD comparison summarizer |

**CI:** `.github/workflows/ci.yml` runs `uv sync --extra dev`, `uv run pytest -q`, then `./scripts/dev/ukbb_ci_compliance_audit.sh` on every push/PR to `main`. There is no committed `uv.lock`; CI resolves the latest compatible dependency versions (including scikit-learn), so clock tests must track current sklearn APIs (for example `ElasticNetCV(alphas=20)` since 1.9).

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

1. **Sensitive / large inputs** → `data/` (git-ignored; `.gitkeep` tracked)
2. **Synthetic fixtures** → `test_data/` (versioned)
3. **Regenerated plots** → `figures/` (git-ignored; `.gitkeep` tracked) — default for all matplotlib/networkx figure scripts
4. **Pipeline tables & models** → `analysis/` (manifests, comparison CSVs, `.pkl` models, VEP cache — most regenerable artifacts git-ignored)
5. **Committed manuscript snapshots** → selected PNG/PDF under `analysis/` (historical exports; regenerate fresh copies under `figures/`)
6. **Ad-hoc scratch** → `outputs/` (git-ignored; `.gitkeep` tracked)
7. **No real UKB data** — use `rogen-ukb-mock-clinical`, `rogen-ukb-mock-rap`, `rogen-ukb-integrate` ([UKB_INTEGRATION_PIPELINE.md](UKB_INTEGRATION_PIPELINE.md))

### Output path quick reference

| Artifact type | Default location | Example |
|---------------|------------------|---------|
| Matplotlib / networkx figures | `figures/` | `figures/Fig_LA_SNP_network.png`, `figures/validation_gse87571/clock_eval_gse87571.png` |
| gnomAD scatter plot | `figures/` | `figures/af_1kg_vs_gnomad_scatter.png` |
| AlphaGenome tables | `analysis/alphagenome/` | `analysis/alphagenome/alphagenome_impact_analysis.csv` |
| AlphaGenome plots (new runs) | `figures/alphagenome/` | `figures/alphagenome/alphagenome_impact_bar_plot.png` |
| LA-SNP manifest / AF CSVs | `analysis/` | `analysis/ukb_snp_manifest_v0.1.csv` |
| VEP annotation + cache | `analysis/` | `analysis/vep_annotation/`, `analysis/vep_cache/` |
| GTEx eQTL annotation + cache | `analysis/` | `analysis/gtex_annotation/`, `analysis/gtex_cache/` |
| Clock model + metrics | `analysis/` | `analysis/gse40279_elasticnet_clock.pkl` |

### Root script shims (deprecated paths)

These one-line forwards remain at the repo root for backward compatibility; canonical paths live under `scripts/`:

| Shim | Canonical |
|------|-----------|
| `annotate_la_snps_vep.py` | `scripts/ukb/annotate_la_snps_vep.py` |
| `annotate_la_snps_gtex.py` | `scripts/ukb/annotate_la_snps_gtex.py` |
| `plot_clock_eval.py` | `scripts/figures/plot_clock_eval.py` |
| `pipeline_validation.sh` | `scripts/dev/pipeline_validation.sh` |
| `find_r.sh` | `scripts/dev/find_r.sh` |

`downstream_analysis.R` moved to `scripts/dev/downstream_analysis.R` (update notebook references if any).

## Configuration

| File | Purpose |
|------|---------|
| `pyproject.toml` | Dependencies, `[project.scripts]` console entry points, pytest |
| `.env` | API keys (git-ignored; copy from `.env.example`) |

---

**Last updated:** June 29, 2026
