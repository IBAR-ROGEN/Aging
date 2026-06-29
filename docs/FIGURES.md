# Manuscript figures — React vs matplotlib

Static figures for the IBAR-ROGEN Aging manuscript. Each concept has a **React/Vite mockup** (for layout) and/or a **Python render script** (for CI and publication exports without Node).

**Navigation:** [ACTIVITIES.md](ACTIVITIES.md) · [WORKFLOWS.md](WORKFLOWS.md) · [PROJECT_STRUCTURE.md](PROJECT_STRUCTURE.md)

All matplotlib/networkx renders live under **`scripts/figures/`**. Legacy paths at `scripts/render_*.py` and `scripts/generate_*.py` forward there with a `DeprecationWarning`.

## Dashboard mockup

| Asset | Tool | Output |
|-------|------|--------|
| React component | `components/DashboardFigureMockup.tsx` | — |
| Vite capture | `components/dashboard-figure-render/` (`npm run capture`) | `figures/dashboard_figure_mockup.png` |
| Matplotlib twin | `scripts/figures/render_dashboard_figure_mockup.py` | `figures/dashboard_figure_mockup.png` |

```bash
uv run python scripts/figures/render_dashboard_figure_mockup.py
```

## Longevity conceptual network

| Asset | Tool | Output |
|-------|------|--------|
| React | `frontend/src/components/LongevityNetworkDiagram.tsx` | — |
| Vite capture | `frontend/scripts/capture-diagram.mjs` | `figures/longevity_network_diagram.png` (git-ignored by default) |
| Matplotlib twin | `scripts/figures/render_longevity_network_diagram.py` | same default path |

```bash
uv run python scripts/figures/render_longevity_network_diagram.py
```

## LA-SNP mechanisms (Figure 1C)

```bash
uv run python scripts/figures/render_figure1c_mechanisms_network.py
# → figures/Figure1C_Mechanisms.png + .pdf
```

## LA-SNPs per gene (supplementary bar chart)

```bash
uv run python scripts/figures/generate_la_snp_per_gene_plot.py
# → figures/Fig_LA_SNPs_per_gene.png
```

Pass `--gene-column` / `--snp-column` if your Excel headers differ.

## LA-SNP pathway network (Activity 2.1.7.1)

```bash
uv run python scripts/figures/generate_network_fig.py
# → figures/Fig_LA_SNP_network.png + .pdf
```

Optional `--pathway-map` CSV with `Gene`, `Pathway` columns. Requires `overlapping_genes_with_snps.xlsx` at repo root (or pass `--input`).

## Methylation pipeline & clock validation plots

| Script | Output |
|--------|--------|
| `scripts/figures/generate_methylation_visualizations.py` | Multiple pipeline QA PNGs under `figures/` |
| `scripts/figures/generate_bimodal_heatmap.py` | `figures/Fig2_Risk_Heatmap.png` |
| `scripts/figures/generate_clock_validation.py` | `figures/Fig3_Clock_Validation.png` |

```bash
uv run python scripts/figures/generate_methylation_visualizations.py
uv run python scripts/figures/generate_bimodal_heatmap.py
uv run python scripts/figures/generate_clock_validation.py
```

Clock **train/evaluate** figures (`Fig_Clock_Residuals.png`, etc.) come from **`rogen-clock evaluate`** — see [CLOCK_LIBRARY.md](CLOCK_LIBRARY.md) (Activity **2.1.10.1**).

### Real-data external validation (GSE87571)

Two-panel figure from a trained model and held-out cohort (predicted vs chronological scatter + top CpG weights):

```bash
uv run python scripts/figures/plot_clock_eval.py
# → figures/validation_gse87571/clock_eval_gse87571.png + .pdf
```

Full configuration and input options: [CLOCK_EVAL_FIGURES.md](CLOCK_EVAL_FIGURES.md).

## Architecture diagrams

| Script | Output | Notes |
|--------|--------|-------|
| `scripts/figures/generate_agent_system_schema.py` | `figures/Fig4_Agent_System_Schema.png` | Uses Graphviz when `dot` is on `PATH`; matplotlib fallback otherwise |
| `scripts/figures/generate_pipeline_diagram.py` | `figures/Bioinformatics_Pipeline_Diagram.png` | Requires Graphviz (`brew install graphviz`) |

```bash
uv run python scripts/figures/generate_agent_system_schema.py
uv run python scripts/figures/generate_pipeline_diagram.py
```

## `scripts/figures/` inventory

| Script | Activity / role |
|--------|-----------------|
| `generate_network_fig.py` | 2.1.7.1 LA-SNP pathway network |
| `render_figure1c_mechanisms_network.py` | Figure 1C mechanisms |
| `generate_la_snp_per_gene_plot.py` | LA-SNPs per gene bar chart |
| `render_dashboard_figure_mockup.py` | Dashboard layout mockup |
| `render_longevity_network_diagram.py` | Longevity conceptual network |
| `generate_methylation_visualizations.py` | ONT methylation pipeline QA |
| `generate_bimodal_heatmap.py` | Risk heatmap |
| `generate_clock_validation.py` | Clock validation scatter |
| `generate_agent_system_schema.py` | Agent system schema |
| `generate_agent_system_schema_fallback.py` | Matplotlib fallback (called by schema script) |
| `plot_clock_eval.py` | GSE87571 external-validation two-panel figure |
| `generate_pipeline_diagram.py` | Bioinformatics pipeline architecture |

## Which asset is canonical?

- **Publication PNG/PDF:** regenerate under **`figures/`** via **`scripts/figures/`** (reproducible, no Node in CI). Committed snapshots under **`analysis/`** are historical exports for the manuscript.
- **Layout exploration:** use React apps under `components/` and `frontend/`.
- **Streamlit EDA** (live data, not a static figure): [EDA_DASHBOARD.md](EDA_DASHBOARD.md).
