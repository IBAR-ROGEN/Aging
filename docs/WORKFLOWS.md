# ROGEN Aging — Workflow index

Navigation hub for installable packages, CLI entry points, and detailed guides.  
**Project layout:** [PROJECT_STRUCTURE.md](PROJECT_STRUCTURE.md) · **Activity map:** [ACTIVITIES.md](ACTIVITIES.md)

## Quickstart

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
uv sync --extra dev
uv run pytest
uv run jupyter lab   # optional
```

Install the Git pre-commit hook (UKB data protection):

```bash
./scripts/dev/install_pre_commit_hook.sh
```

## Console entry points (`uv run …`)

| Command | Purpose |
|---------|---------|
| `rogen-clock train …` / `rogen-clock evaluate …` | Epigenetic clock train & validation ([CLOCK_LIBRARY.md](CLOCK_LIBRARY.md)) |
| `rogen-ukb-manifest build …` / `… extract …` | LA-SNP manifest + 1KG AF ([LA_SNP_PUBLIC_FREQUENCY_PIPELINE.md](LA_SNP_PUBLIC_FREQUENCY_PIPELINE.md)) |
| `rogen-compare-af-gnomad …` | 1KG vs gnomAD v4 NFE ([LA_SNP_PUBLIC_FREQUENCY_PIPELINE.md](LA_SNP_PUBLIC_FREQUENCY_PIPELINE.md)) |
| `rogen-compare-af-gnomad summarize …` | Markdown summary + top-|ΔAF| table from comparison CSV |
| `rogen-ukb-mock-clinical …` | Synthetic clinical CSV ([SYNTHETIC_UKB_GENERATOR.md](SYNTHETIC_UKB_GENERATOR.md)) |
| `rogen-ukb-mock-rap …` | Synthetic UKB-RAP folder (phenotypes + LA-SNP VCF) ([SYNTHETIC_UKB_RAP_GENERATOR.md](SYNTHETIC_UKB_RAP_GENERATOR.md)) |
| `rogen-ukb-integrate …` | Mock phenotype–genotype join + LA-SNP associations ([UKB_INTEGRATION_PIPELINE.md](UKB_INTEGRATION_PIPELINE.md)) |
| `rogen-vcf-synthetic …` | Streaming synthetic VCF ([SYNTHETIC_ROMANIAN_VCF_GENERATOR.md](SYNTHETIC_ROMANIAN_VCF_GENERATOR.md)) |

Legacy script paths under `scripts/*.py` remain as **deprecation shims** forwarding to `scripts/<workflow>/`.

## Workflows

### Epigenetic clock

- **Package:** `src/rogen_aging/clock/`
- **Canonical CLI:** `uv run rogen-clock train|evaluate` or `scripts/clock/run_clock.py`
- **Romanian mock demo** (separate StandardScaler path): `scripts/clock/train_romanian_epigenetic_clock.py`
- **Docs:** [CLOCK_LIBRARY.md](CLOCK_LIBRARY.md), [GSE40279_CLOCK_TRAINING.md](GSE40279_CLOCK_TRAINING.md), [ROMANIAN_EPIGENETIC_CLOCK.md](ROMANIAN_EPIGENETIC_CLOCK.md)

### UK Biobank (synthetic + LA-SNP)

- **Packages:** `src/rogen_aging/ukb/`, `src/rogen_aging/integration/`
- **Mock clinical CSV:** `scripts/ukb/mock_clinical_csv.py` · **Mock RAP folder:** `scripts/ukb/mock_rap_folder.py`
- **Integration scan:** `uv run rogen-ukb-integrate` or `scripts/ukb/run_integration.py`
- **Docs:** [SYNTHETIC_UKB_GENERATOR.md](SYNTHETIC_UKB_GENERATOR.md), [SYNTHETIC_UKB_RAP_GENERATOR.md](SYNTHETIC_UKB_RAP_GENERATOR.md), [UKB_INTEGRATION_PIPELINE.md](UKB_INTEGRATION_PIPELINE.md), [LA_SNP_PUBLIC_FREQUENCY_PIPELINE.md](LA_SNP_PUBLIC_FREQUENCY_PIPELINE.md)

### Methylation (Oxford Nanopore)

- **Root scripts:** `pipeline_validation.sh`, `downstream_analysis.R`
- **Docs:** [METHYLATION_PIPELINE_README.md](METHYLATION_PIPELINE_README.md), [METHYLATION_PIPELINE_USAGE.md](METHYLATION_PIPELINE_USAGE.md)

### Multi-omics EDA dashboard

```bash
uv run streamlit run src/rogen_aging/eda_dashboard/app.py
```

See [EDA_DASHBOARD.md](EDA_DASHBOARD.md).

### AlphaGenome (Activity 2.1.7.1)

```bash
uv run python scripts/alphagenome/alphagenome_sequence_comparer.py
uv run python scripts/alphagenome/analyze_alphagenome_results.py
uv run python scripts/alphagenome/visualize_alphagenome_results.py
```

See [ALPHAGENOME_ANALYSIS_EXPLANATION.md](ALPHAGENOME_ANALYSIS_EXPLANATION.md).

### Manuscript figures

See [FIGURES.md](FIGURES.md) for React vs matplotlib assets and render commands.

### Compliance & CI

- Pre-commit: [UKB_PRE_COMMIT_HOOK.md](UKB_PRE_COMMIT_HOOK.md)
- CI audit: [UKBB_CI_COMPLIANCE_AUDIT.md](UKBB_CI_COMPLIANCE_AUDIT.md) · `./scripts/dev/ukbb_ci_compliance_audit.sh`
- GitHub Actions: `.github/workflows/ci.yml`

## Notebooks

Grouped under `notebooks/01_*` … `05_*`. See [notebooks/README.md](../notebooks/README.md).
