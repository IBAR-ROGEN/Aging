# ROGEN Aging — Code Modules Reference

**Navigation:** [WORKFLOWS.md](WORKFLOWS.md) · [PROJECT_STRUCTURE.md](PROJECT_STRUCTURE.md) · [ACTIVITIES.md](ACTIVITIES.md)

This document summarizes installable modules and CLI entry points. Per-file inventories live in source docstrings; avoid duplicating README content here.

## Installable package (`src/rogen_aging/`)

| Module | Responsibility |
|--------|----------------|
| `clock/` | `train_clock`, `evaluate_clock`, `load_wide_table`, `external_data.load_gse87571` (Activity **2.1.10.1**) — [CLOCK_LIBRARY.md](CLOCK_LIBRARY.md) |
| `ukb/manifest.py` | Ensembl manifest build + 1KG VCF extract |
| `ukb/gnomad.py` | 1KG vs gnomAD v4 NFE comparison |
| `ukb/mock_clinical.py` | Synthetic clinical CSV generator |
| `ukb/mock_rap.py` | Synthetic UKB-RAP phenotype + LA-SNP VCF folder |
| `integration/ukb_joiner.py` | Mock phenotype–genotype join + LA-SNP associations |
| `vcf/synthetic.py` | Streaming synthetic VCF (Hardy–Weinberg, GRCh38) |
| `eda_dashboard/` | Streamlit multi-omics EDA — [EDA_DASHBOARD.md](EDA_DASHBOARD.md) |
| `cli/` | Console entry points wired in `pyproject.toml` |
| `methylation_visualizations.py` | Pipeline / DMR / clock validation plots (default output: `figures/`) |
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

- **`scripts/ukb/`** — manifest, gnomAD, VEP + GTEx annotation, mock generators, integration
- **`scripts/figures/`** — matplotlib/networkx renders + `plot_clock_eval.py` ([FIGURES.md](FIGURES.md)); flat `scripts/generate_*.py` shims forward here
- **`scripts/alphagenome/`** — AlphaGenome batch + analysis (tables → `analysis/alphagenome/`, plots → `figures/alphagenome/`)
- **`scripts/dev/`** — `security_check.sh`, CI audit, ONT pipeline validation, R bootstrap, `find_r.sh`

## Root-level shims (deprecated)

| Shim | Forwards to |
|------|-------------|
| `annotate_la_snps_vep.py` | `scripts/ukb/annotate_la_snps_vep.py` |
| `annotate_la_snps_gtex.py` | `scripts/ukb/annotate_la_snps_gtex.py` |
| `plot_clock_eval.py` | `scripts/figures/plot_clock_eval.py` |
| `plot_af_comparison.py` | `scripts/figures/plot_af_comparison.py` |
| `pipeline_validation.sh` | `scripts/dev/pipeline_validation.sh` |
| `find_r.sh` | `scripts/dev/find_r.sh` |

## Root-level non-Python (canonical under `scripts/dev/`)

| File | Purpose |
|------|---------|
| `scripts/dev/pipeline_validation.sh` | ONT methylation pipeline validation |
| `scripts/dev/downstream_analysis.R` | DMRcaller downstream workflow |
| `scripts/dev/find_r.sh` | Locate R for IDE setup |

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
| `components/dashboard-figure-render/` | Vite + Playwright capture → `figures/dashboard_figure_mockup.png` |
| `frontend/` | Longevity network diagram + capture → `figures/longevity_network_diagram.png` |

## Tests

| File | Covers |
|------|--------|
| `test_clock_regression.py` | Clock refactor vs legacy metrics (`ElasticNetCV(alphas=20)`; sklearn 1.9+ compatible) |
| `test_ukb_integration.py` | Synthetic join + association scan |
| `test_ukb_mock_gen.py` | Mock RAP folder layout |
| `test_synthetic_vcf.py` | `rogen_aging.vcf` |
| `test_mock_clinical_csv.py` | `rogen_aging.ukb.mock_clinical` |
| `test_package_imports.py` | Smoke imports |

---

**Last updated:** June 29, 2026
