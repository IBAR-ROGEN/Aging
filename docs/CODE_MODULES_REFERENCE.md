# ROGEN Aging — Code Modules Reference

**Navigation:** [WORKFLOWS.md](WORKFLOWS.md) · [PROJECT_STRUCTURE.md](PROJECT_STRUCTURE.md) · [ACTIVITIES.md](ACTIVITIES.md)

This document summarizes installable modules and CLI entry points. Per-file inventories live in source docstrings; avoid duplicating README content here.

## Installable package (`src/rogen_aging/`)

| Module | Responsibility |
|--------|----------------|
| `clock/` | `train_clock`, `evaluate_clock`, `load_wide_table`, GSE87571 loader — [CLOCK_LIBRARY.md](CLOCK_LIBRARY.md) |
| `ukb/manifest.py` | Ensembl manifest build + 1KG VCF extract |
| `ukb/gnomad.py` | 1KG vs gnomAD v4 NFE comparison |
| `ukb/mock_clinical.py` | Synthetic clinical CSV generator |
| `integration/ukb_joiner.py` | Mock phenotype–genotype join + LA-SNP associations |
| `vcf/synthetic.py` | Streaming synthetic VCF (Hardy–Weinberg, GRCh38) |
| `eda_dashboard/` | Streamlit multi-omics EDA — [EDA_DASHBOARD.md](EDA_DASHBOARD.md) |
| `cli/` | Console entry points wired in `pyproject.toml` |
| `methylation_visualizations.py` | Pipeline / DMR / clock validation plots |
| `network_visualizer.py` | Protein interaction networks |

## Console entry points

| Command | Module |
|---------|--------|
| `rogen-clock` | `rogen_aging.cli.clock` |
| `rogen-ukb-manifest` | `rogen_aging.cli.ukb_manifest` → `ukb.manifest` |
| `rogen-compare-af-gnomad` | `rogen_aging.cli.compare_af_gnomad` → `ukb.gnomad` (`compare` + `summarize` subcommands) |
| `rogen-ukb-mock-clinical` | `rogen_aging.cli.ukb_mock_clinical` |
| `rogen-ukb-mock-rap` | `rogen_aging.cli.ukb_mock_rap` |
| `rogen-ukb-integrate` | `rogen_aging.cli.ukb_integrate` → `integration.run_cli` |
| `rogen-vcf-synthetic` | `rogen_aging.cli.vcf_synthetic` |

**Canonical clock CLI:** `uv run rogen-clock train|evaluate` or `scripts/clock/run_clock.py`.  
`scripts/clock/validate_clock.py` and `scripts/clock/train_clock_on_gse40279.py` are deprecated Typer wrappers.

## Scripts by folder

See [ACTIVITIES.md](ACTIVITIES.md) for the full tree. Highlights:

- **`scripts/ukb/`** — manifest, gnomAD, mock generators, integration
- **`scripts/figures/`** — matplotlib/networkx renders ([FIGURES.md](FIGURES.md))
- **`scripts/alphagenome/`** — AlphaGenome batch + analysis
- **`scripts/dev/`** — `security_check.sh`, CI audit, R bootstrap

## Root-level non-Python

| File | Purpose |
|------|---------|
| `pipeline_validation.sh` | ONT methylation pipeline validation |
| `downstream_analysis.R` | DMRcaller downstream workflow |
| `find_r.sh` | Locate R for IDE setup |

## Notebooks

| Folder | Focus |
|--------|-------|
| `01_genomics_analysis/` | AlphaGenome, gene lists |
| `02_methylation_pipeline/` | Methylation downstream |
| `03_validation_and_compliance/` | UKB compliance |
| `04_exploratory_visualizations/` | Publication figures |
| `05_ukb_exploration/` | LA-SNP manifest QA |

## TypeScript figure apps

| Path | Role |
|------|------|
| `components/DashboardFigureMockup.tsx` | Dashboard manuscript mockup |
| `components/dashboard-figure-render/` | Vite + Playwright capture |
| `frontend/` | Longevity network diagram + capture |

## Tests

| File | Covers |
|------|--------|
| `test_clock_regression.py` | Clock refactor vs legacy metrics |
| `test_ukb_integration.py` | Synthetic join + association scan |
| `test_ukb_mock_gen.py` | Mock RAP folder layout |
| `test_synthetic_vcf.py` | `rogen_aging.vcf` |
| `test_mock_clinical_csv.py` | `rogen_aging.ukb.mock_clinical` |
| `test_package_imports.py` | Smoke imports |

---

**Last updated:** May 31, 2026
