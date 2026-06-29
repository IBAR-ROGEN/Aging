# rogen_aging

IBAR-ROGEN aging analysis: epigenetic clocks, UK Biobank LA-SNP public-frequency tooling, synthetic cohort integration, and manuscript figure pipelines. Python ≥3.12, managed with [uv](https://docs.astral.sh/uv/).

## Install

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
uv sync --extra dev
uv run pytest
```

Optional Jupyter kernel:

```bash
uv run python -m ipykernel install --user --name rogen-aging --display-name "Python (rogen-aging)"
uv run jupyter lab
```

Install the UK Biobank pre-commit security hook:

```bash
./scripts/dev/install_pre_commit_hook.sh
```

## Continuous integration

GitHub Actions (`.github/workflows/ci.yml`) on every push/PR to `main`:

1. `uv sync --extra dev`
2. `uv run pytest -q`
3. `./scripts/dev/ukbb_ci_compliance_audit.sh`

See [docs/UKBB_CI_COMPLIANCE_AUDIT.md](docs/UKBB_CI_COMPLIANCE_AUDIT.md) for the audit rules. Clock training uses `ElasticNetCV(alphas=20)` for scikit-learn 1.9+ compatibility — details in [docs/GSE40279_CLOCK_TRAINING.md](docs/GSE40279_CLOCK_TRAINING.md#scikit-learn-compatibility).

## Package map

Installable code lives under `src/rogen_aging/`. Console entry points are registered in `pyproject.toml`.

| Module | Purpose |
|--------|---------|
| **`rogen_aging.clock`** | Train and evaluate ElasticNet epigenetic clocks; load GSE87571 and Romanian mock cohorts |
| **`rogen_aging.ukb`** | Build LA-SNP manifest CSVs, extract 1KG allele frequencies, compare to gnomAD, generate synthetic UKB-RAP mocks |
| `rogen_aging.vcf` | Synthetic Romanian VCF generation for pipeline testing |
| `rogen_aging.integration` | Join synthetic UKB phenotypes/genotypes and run LA-SNP association summaries |
| `rogen_aging.eda_dashboard` | Streamlit EDA on merged mock clinical / epigenetic-age tables |
| `rogen_aging.cli` | Typer wrappers: `rogen-clock`, `rogen-ukb-manifest`, `rogen-ukb-integrate`, … |

## Directory layout

| Path | Role |
|------|------|
| `src/rogen_aging/` | Installable package |
| `scripts/` | Grouped CLIs (`clock/`, `ukb/`, `figures/`, `alphagenome/`, `dev/`); flat paths are deprecation shims |
| `tests/` | pytest (`uv run pytest`) |
| `notebooks/` | Jupyter by analysis area |
| `docs/` | Workflow guides — start at [docs/WORKFLOWS.md](docs/WORKFLOWS.md) |
| `analysis/` | Committed manuscript figure snapshots, alphagenome tables, pipeline CSVs/manifests |
| `figures/` | **Local regenerated plots** (git-ignored; `.gitkeep` only in git) |
| `data/` | Large/local inputs and caches (git-ignored) |
| `outputs/` | Optional scratch for ad-hoc pipeline runs (git-ignored) |
| `results/` | EDA and ad-hoc plot output (git-ignored; see [docs/EDA_MOCK_INTEGRATION.md](docs/EDA_MOCK_INTEGRATION.md)) |
| `frontend/` | React/Vite longevity network diagram (layout + PNG capture) |
| `components/` | React/Vite dashboard figure mockup (layout + PNG capture) |
| `test_data/` | Versioned synthetic fixtures |

## Common commands

### Epigenetic clock

```bash
uv run rogen-clock train \
  --input_data data/gse40279.parquet \
  --output_model analysis/gse40279_elasticnet_clock.pkl \
  --output_metrics analysis/gse40279_train_metrics.json

uv run python -m rogen_aging.clock.external_data --output data/gse87571.parquet

uv run rogen-clock evaluate \
  --model_path analysis/gse40279_elasticnet_clock.pkl \
  --test_data data/gse87571.parquet \
  --output_dir figures/validation_gse87571

uv run python scripts/figures/plot_clock_eval.py
# → figures/validation_gse87571/clock_eval_gse87571.png + .pdf
```

### UK Biobank LA-SNP pipeline (public data only)

```bash
uv run rogen-ukb-manifest build --input overlapping_genes_with_snps.xlsx \
  --output analysis/ukb_snp_manifest_v0.1.csv
uv run rogen-compare-af-gnomad \
  --input analysis/la_snp_1kg_frequencies.csv \
  --output analysis/la_snp_af_1kg_vs_gnomad.csv \
  --scatter figures/af_1kg_vs_gnomad_scatter.png
uv run python scripts/figures/plot_af_comparison.py
# → figures/af_1kg_vs_gnomad_comparison.png + .pdf
uv run rogen-ukb-mock-rap --n-samples 1000 --output-dir test_data/mock_ukb_rap/
uv run rogen-ukb-integrate --output-dir analysis/
```

### Figures and annotation scripts

```bash
uv run python scripts/figures/generate_network_fig.py      # → figures/Fig_LA_SNP_network.*
uv run python scripts/ukb/annotate_la_snps_vep.py          # VEP table + cache under analysis/
uv run python scripts/ukb/annotate_la_snps_gtex.py         # GTEx eQTL table + cache under analysis/
uv run python scripts/alphagenome/alphagenome_sequence_comparer.py
uv run streamlit run src/rogen_aging/eda_dashboard/app.py
```

Flat paths such as `plot_clock_eval.py`, `plot_af_comparison.py`, `annotate_la_snps_vep.py`, and `annotate_la_snps_gtex.py` at the repo root forward to `scripts/` with a deprecation warning.

## Documentation index

| Topic | Link |
|-------|------|
| Workflow index | [docs/WORKFLOWS.md](docs/WORKFLOWS.md) |
| Activity map | [docs/ACTIVITIES.md](docs/ACTIVITIES.md) |
| Code reference | [docs/CODE_MODULES_REFERENCE.md](docs/CODE_MODULES_REFERENCE.md) |
| Directory layout | [docs/PROJECT_STRUCTURE.md](docs/PROJECT_STRUCTURE.md) |
| Manuscript figures | [docs/FIGURES.md](docs/FIGURES.md) |
| LA-SNP VEP / GTEx annotation | [docs/LA_SNP_VEP_ANNOTATION.md](docs/LA_SNP_VEP_ANNOTATION.md) · [docs/LA_SNP_GTEX_ANNOTATION.md](docs/LA_SNP_GTEX_ANNOTATION.md) |
| Epigenetic clock | [docs/CLOCK_LIBRARY.md](docs/CLOCK_LIBRARY.md) · [eval figure](docs/CLOCK_EVAL_FIGURES.md) |
| LA-SNP public AF validation | [docs/LA_SNP_PUBLIC_FREQUENCY_PIPELINE.md](docs/LA_SNP_PUBLIC_FREQUENCY_PIPELINE.md) · [comparison figure](docs/AF_COMPARISON_FIGURES.md) |
| Synthetic UKB integration | [docs/UKB_INTEGRATION_PIPELINE.md](docs/UKB_INTEGRATION_PIPELINE.md) |
| Methylation pipeline | [docs/METHYLATION_PIPELINE_README.md](docs/METHYLATION_PIPELINE_README.md) |

## License

MIT — see [LICENSE](LICENSE).
