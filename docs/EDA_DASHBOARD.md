# Multi-Omics Aging EDA Dashboard (Streamlit)

Interactive exploratory data analysis (EDA) for the **merged cohort Parquet** produced by the ROGEN multi-omics integration pipeline. The app is implemented as the `rogen_aging.eda_dashboard` package and launched with Streamlit.

## Prerequisites

- Python 3.12+ and a synced uv environment (`uv sync`).
- Dependencies include **Streamlit**, **Plotly**, **Seaborn**, **Polars**, **Pandas**, **SciPy**, and **PyArrow** (see `pyproject.toml`).

## Launch

From the repository root:

```bash
uv run streamlit run src/rogen_aging/eda_dashboard/app.py
```

Streamlit opens a local URL (default `http://localhost:8501`).

## Data input

### Default path

The app looks for:

`<repo_root>/data/merged_cohort.parquet`

The `data/` directory is typically git-ignored; place your pipeline output there or override the path.

### Environment variable

Set **`ROGEN_MERGED_COHORT_PARQUET`** to an absolute or user-expanded path to use that file instead of the default (sidebar path still initializes from this when set in the environment for the default display; the sidebar text field can override interactively).

### Missing file

If the Parquet file is absent, the UI shows an error and a **Generate in-memory mock cohort** button. That loads a cached synthetic dataset (same schema as the reference design) so layouts and filters can be exercised without real data.

## Expected schema (and aliases)

After load, the dashboard **normalizes** common aliases to canonical names where possible:

| Concept | Canonical column | Accepted aliases (examples) |
|--------|-------------------|-----------------------------|
| Sample identifier | `Sample_ID` | `sample_id`, `IID`, `participant_id` |
| Chronological age | `Chronological_Age` | `chronological_age`, `Age`, `age` |
| Sex | `Sex` | `sex`, `gender`, … |
| Disease / group | `Disease_Status` | `disease_status`, `Case_Control`, … |
| Epigenetic age | `Epigenetic_Age` | `epigenetic_age`, `DNAmAge`, … |

**Epigenetic age acceleration** (`Epigenetic_Age_Acceleration`) is added automatically when chronological and epigenetic ages are both present: `Epigenetic_Age - Chronological_Age`.

**SNP / LA-SNP columns** are detected by a case-insensitive `rs\d+` substring in the column name (for example `rs5882_CETP`). Genotypes may be encoded as **dosage** `0`, `1`, `2` (shown as 0/0, 0/1, 1/1) or as strings such as `0/0`.

## Sidebar (global filters)

All tabs use the same filtered cohort:

- **Age range** — slider on chronological age (numeric).
- **Sex** — multiselect (all unique values in the column).
- **Disease status** — multiselect.

A caption under the title reports filtered **N** versus total **N**.

## Tabs

### 1. Clinical & phenotypic overview

- **Stacked histogram** of chronological age by disease status (Plotly Express).
- **Correlation heatmap** of continuous clinical numeric columns (Seaborn / Pearson). SNP-like and ID-like columns are excluded from the heatmap set.

### 2. Epigenetic clock validation

- Scatter of **chronological age** vs **epigenetic age** (Plotly Express) with an **OLS** line (`numpy.polyfit`, degree 1).
- **Metrics:** mean absolute error (MAE) of OLS fit, Pearson *r*, and sample count.

### 3. Genotype–phenotype (LA-SNPs)

- Choose a **SNP column** and a **numeric trait** (for example HDL cholesterol or epigenetic age acceleration).
- **Box plot** by genotype group; **Kruskal–Wallis** omnibus *p*-value (and *H*) in the caption when there are at least two groups with sufficient data.

## Module layout (`src/rogen_aging/eda_dashboard/`)

| Module | Role |
|--------|------|
| `app.py` | Page config, data source selection, tabs, session state. |
| `data.py` | `@st.cache_data` Parquet load (Polars → pandas) and synthetic cohort builder. |
| `schema.py` | Column resolution, SNP discovery, heatmap column selection, EAA helper. |
| `sidebar.py` | Global filters dataclass, sidebar UI, `apply_global_filters()`. |
| `tabs.py` | `render_tab_clinical_overview`, `render_tab_epigenetic_clock_validation`, `render_tab_la_snp_genotype_phenotype`. |

## Related documentation

- **[EDA_MOCK_INTEGRATION.md](EDA_MOCK_INTEGRATION.md)** — Earlier CSV-based mock EDA script (complementary, not the Streamlit app).
- **[PROJECT_STRUCTURE.md](PROJECT_STRUCTURE.md)** — Repository layout.
