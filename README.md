# rogen_aging

Project scaffold for genomic notebooks and analysis, managed with `uv`.

## Quickstart

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

## Where to go next

| Topic | Link |
|-------|------|
| **Workflow index** | [docs/WORKFLOWS.md](docs/WORKFLOWS.md) |
| **Activity map** | [docs/ACTIVITIES.md](docs/ACTIVITIES.md) |
| **Code reference** | [docs/CODE_MODULES_REFERENCE.md](docs/CODE_MODULES_REFERENCE.md) |
| **Directory layout** | [docs/PROJECT_STRUCTURE.md](docs/PROJECT_STRUCTURE.md) |
| **Manuscript figures** | [docs/FIGURES.md](docs/FIGURES.md) |
| **Epigenetic clock** | [docs/CLOCK_LIBRARY.md](docs/CLOCK_LIBRARY.md) |
| **Methylation pipeline** | [docs/METHYLATION_PIPELINE_README.md](docs/METHYLATION_PIPELINE_README.md) |

## Layout (summary)

- `src/rogen_aging/` — installable package (`clock`, `ukb`, `vcf`, `integration`, `eda_dashboard`, …)
- `scripts/` — grouped CLIs (`clock/`, `ukb/`, `figures/`, …); flat paths are deprecation shims
- `tests/` — pytest (`uv run pytest`)
- `notebooks/` — Jupyter by analysis area
- `docs/` — guides and workflow index
- `components/`, `frontend/` — React/Vite manuscript figure mockups
- `test_data/` — versioned synthetic fixtures
- `data/`, `results/`, `outputs/` — local/large data (git-ignored)

## Common console commands

```bash
uv run rogen-clock train --input_data … --output_model … --output_metrics …
uv run rogen-clock evaluate --model_path … --test_data … --output_dir …
uv run rogen-ukb-manifest build --input overlap.xlsx --output analysis/ukb_snp_manifest_v0.1.csv
uv run streamlit run src/rogen_aging/eda_dashboard/app.py
```

See [docs/WORKFLOWS.md](docs/WORKFLOWS.md) for full command tables and legacy script paths.

## Python version

Python ≥3.12 (`pyproject.toml`).

## License

MIT — see [LICENSE](LICENSE).
