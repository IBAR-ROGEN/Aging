# ROGEN Aging Project — Code Modules Reference

This document describes each code file and module in the IBAR-ROGEN Aging project and what it is responsible for.

---

## Directory structure (overview)

See **[docs/PROJECT_STRUCTURE.md](PROJECT_STRUCTURE.md)** for the full directory layout.

```
Aging/
├── downstream_analysis.R          # Root R script
├── find_r.sh
├── pipeline_validation.sh
├── pyproject.toml
├── setup.py                       # Setuptools entry (metadata in pyproject.toml)
├── requirements.txt               # Optional pip-style dependencies
├── src/rogen_aging/               # Installable Python package
│   ├── __init__.py                # Public API + submodule exports
│   ├── pipeline/                  # Shared pipeline code (placeholder)
│   │   └── __init__.py
│   ├── methylation_visualizations.py
│   ├── network_visualizer.py
│   └── eda_dashboard/            # Streamlit merged-cohort EDA (see docs/EDA_DASHBOARD.md)
├── tests/                         # Pytest (package smoke + synthetic UKB/VCF generators)
├── scripts/                       # Entry-point scripts
│   ├── mock_ukb_generator.py      # Synthetic UK Biobank data
│   ├── generate_synthetic_romanian_vcf.py  # Synthetic EUR-style cohort VCF v4.2
│   ├── ukb_la_snp_lookup.py       # Offline UKB SNP manifest (Ensembl GRCh38)
│   ├── render_*_network*.py       # Manuscript networks (matplotlib / networkx)
│   ├── render_dashboard_figure_mockup.py
│   ├── bootstrap_r_env.sh         # Optional micromamba R under .r-env/
│   ├── security_check.sh          # UK Biobank pre-commit hook
│   ├── install_pre_commit_hook.sh
│   ├── alphagenome_sequence_comparer.py
│   ├── analyze_alphagenome_results.py
│   ├── visualize_alphagenome_results.py
│   ├── generate_*.py
│   ├── generate_la_snp_per_gene_plot.py  # LA-SNPs per gene (supplementary bar chart)
│   ├── train_romanian_epigenetic_clock.py  # Elastic Net mock clock
│   ├── validate_clock.py          # Held-out validation for saved clock models
│   ├── install_graphviz.sh
│   └── README_GRAPHVIZ.md
├── components/                    # React/Vite dashboard figure mockup + capture
├── frontend/                      # React/Vite longevity network + capture
├── notebooks/                     # Notebooks by functional area
│   ├── 01_genomics_analysis/
│   ├── 02_methylation_pipeline/
│   ├── 03_validation_and_compliance/
│   ├── 04_exploratory_visualizations/
│   ├── 05_ukb_exploration/        # UKB manifest first-contact QA
│   └── README.md
├── docs/                          # Documentation
├── test_data/                     # Synthetic test data (versioned)
├── analysis/                      # Generated figures and reports
└── .vscode/, .env.example, ...
```

---

## 1. Project root (Aging/)

### 1.1 downstream_analysis.R

**Purpose:** Downstream methylation analysis with DMRcaller (Activity 2.1.8.1).

**Responsibilities:**
- Imports bedMethyl files produced by Modkit (from the methylation pipeline).
- Prepares methylation data for the DMRcaller Bioconductor package.
- Identifies Differentially Methylated Regions (DMRs) between groups (e.g., Old vs Young).
- Provides a full workflow: data import → preparation → DMR calling → (optional) export.

**Main functions:**
- `import_bedmethyl()` — Reads bedMethyl (BED 9+3) via rtracklayer into GRanges.
- `prepare_dmrcaller_data()` — Converts bedMethyl data into coverage and methylation matrices for DMRcaller.
- `calculate_dmrs()` — Calls DMRs between two groups with configurable coverage, p-value, and min CpG options.
- `run_dmrcaller_workflow()` — Runs the end-to-end analysis workflow.

**Dependencies:** DMRcaller, GenomicRanges, rtracklayer (R/Bioconductor).

**Input:** bedMethyl files from Modkit.  
**Output:** DMRs (GRanges) and optional BED/CSV exports.

---

### 1.2 find_r.sh

**Purpose:** Locate the R installation on macOS for IDE/VS Code configuration.

**Responsibilities:**
- Searches common paths for the R executable (e.g. `/Library/Frameworks/R.framework`, `/opt/homebrew/bin/R`, `/usr/local/bin/R`).
- Prints the path and version of the first R binary found.
- Suggests the exact `r.rpath.mac` setting for `.vscode/settings.json`.
- If R is not found, suggests how to install R (CRAN or Homebrew).

**Usage:** Run from the project root; use the printed path in VS Code/Cursor R extension settings.

---

### 1.3 pipeline_validation.sh

**Purpose:** Validate the Oxford Nanopore methylation calling pipeline (Activity 2.1.8.1).

**Responsibilities:**
- **Step 1:** Checks that Dorado and Modkit are installed and reports versions.
- **Step 2:** Downloads the official ONT/Epi2Me Labs test dataset (wf-basecalling-demo).
- **Step 3:** Runs (or documents) Dorado basecalling with a methylation-aware model (5mC/5hmC), producing BAM with MM/ML tags. GPU (CUDA) or CPU.
- **Step 4:** Documents/runs Modkit to convert BAM to bedMethyl; notes that a reference FASTA is required.
- **Step 5:** Prints a validation summary and next steps (uncomment Dorado/Modkit when ready, then run `downstream_analysis.R`).

**Output:** Test data in `wf-basecalling-demo/`, optional `basecalled_methylation.bam`, optional `methylation_calls.bedMethyl`, or placeholder files when steps are commented out.

---

## 2. Python package: src/rogen_aging/

### 2.1 __init__.py

**Purpose:** Package root for `rogen_aging`; defines the stable import surface for notebooks and scripts.

**Responsibilities:**
- Re-exports visualization entrypoints so callers can use `from rogen_aging import create_pipeline_workflow_diagram`, `create_network_visualization`, `generate_all_visualizations`, etc.
- Binds submodules for explicit use: `from rogen_aging import methylation_visualizations`, `network_visualizer`, `pipeline`.
- Declares `__all__` for `from rogen_aging import *` (discouraged in application code, but supported).

**After install:** With `uv sync` or `uv pip install -e .`, use absolute imports only (`from rogen_aging…`), not `sys.path` hacks.

---

### 2.2 pipeline/ (subpackage)

**Purpose:** Home for shared pipeline logic as it is refactored out of `scripts/` and notebooks.

**Responsibilities:** Currently a placeholder (`__all__ = []`). Add modules here (e.g. QC, I/O, model steps) and export names from `pipeline/__init__.py`.

---

### 2.3 methylation_visualizations.py

**Purpose:** Central module for all methylation-related figures and pipeline diagrams.

**Responsibilities:**
- **Pipeline workflow diagram** — Draws the flow: POD5 → Dorado → BAM → Modkit → bedMethyl → DMRcaller → DMR results (with arrows and legend).
- **Example DMR visualizations** — Simulated DMR data: Manhattan plot, DMR width distribution, methylation difference distribution, p-value distribution, CpG sites per DMR.
- **Pipeline summary diagram** — High-level overview of tools (Dorado, Modkit, DMRcaller) and their roles.
- **Bimodal risk heatmap (Figure 2)** — Protective vs risk effects of longevity genes across conditions (Activity 2.1.7).
- **Clock validation plot (Figure 3)** — Scatter of chronological age vs DNAm predicted age with MAE ~2.1 years (Activity 2.1.10).
- **generate_all_visualizations()** — Runs pipeline workflow, example DMR, and pipeline summary in one call.

**Main functions:**
- `create_pipeline_workflow_diagram()`
- `create_example_dmr_visualizations()`
- `create_pipeline_summary_diagram()`
- `create_bimodal_risk_heatmap()`
- `create_clock_validation_plot()`
- `generate_all_visualizations()`

**Dependencies:** matplotlib, seaborn, numpy, pandas.  
**Output:** PNGs under `analysis/` (or paths passed in).

---

### 2.4 network_visualizer.py

**Purpose:** Protein interaction network visualization (“Resilience Core”).

**Responsibilities:**
- Builds a NetworkX graph of longevity/neuro genes (Hubs, Longevity, Neuro categories).
- Adds edges representing interactions (e.g. APOE–MAPT, TP53–SIRT1).
- Sizes nodes by degree centrality.
- Draws the graph with category-based colors and saves as PNG.

**Main function:** `create_network_visualization(output_path="Network_Analysis_Nov.png")`

**Dependencies:** networkx, matplotlib, numpy.  
**Output:** Single PNG (default `Network_Analysis_Nov.png`).

---

### 2.5 eda_dashboard/ (package subfolder)

**Purpose:** Official **Streamlit** exploratory dashboard for the multi-omics aging **merged Parquet** cohort.

**Responsibilities:**
- Load merged table with `@st.cache_data` (Polars read, pandas for Plotly/Seaborn).
- Sidebar **global filters** (age range, sex, disease status) applied across all tabs.
- **Tab 1** — Stacked age histogram by disease; clinical correlation heatmap.
- **Tab 2** — Chronological vs epigenetic age scatter (Plotly Express) with OLS line; MAE and Pearson *r* metrics.
- **Tab 3** — LA-SNP genotype vs continuous trait boxplots; Kruskal–Wallis caption.
- Synthetic in-memory cohort when Parquet is missing.

**Entry point:** `uv run streamlit run src/rogen_aging/eda_dashboard/app.py`

**Documentation:** [docs/EDA_DASHBOARD.md](EDA_DASHBOARD.md)

**Dependencies:** streamlit, plotly, seaborn, polars, pandas, pyarrow, scipy, numpy.

---

## 3. Scripts (scripts/)

Scripts are entry points that call into `src/rogen_aging` or external tools. Run them from the `Aging/` directory (e.g. `python scripts/generate_*.py` or `uv run python scripts/...`).

### 3.1 generate_agent_system_schema.py

**Purpose:** Generate Figure 4 — LongevityForest agent system architecture diagram.

**Responsibilities:**
- Uses the `diagrams` library (Graphviz) to draw: Researcher → Cursor IDE → MCP → LongevityForest cluster (BioMART, AlphaFold, STRING).
- Locates the Graphviz `dot` executable (PATH or common install paths).
- If Graphviz is missing, delegates to `generate_agent_system_schema_fallback.py` (matplotlib version).

**Output:** `analysis/Fig4_Agent_System_Schema.png` (or Graphviz source in `analysis/`).

---

### 3.2 generate_agent_system_schema_fallback.py

**Purpose:** Fallback for Figure 4 when Graphviz is not installed.

**Responsibilities:**
- Draws the same agent system schema using matplotlib (boxes, arrows, labels): Researcher, Cursor IDE, MCP, BioMART, AlphaFold, STRING, cluster.

**Output:** `analysis/Fig4_Agent_System_Schema.png`.

---

### 3.3 generate_bimodal_heatmap.py

**Purpose:** Generate the bimodal risk heatmap (Figure 2).

**Responsibilities:** Adds project root to `sys.path`, then calls `create_bimodal_risk_heatmap()` from `methylation_visualizations.py`.

**Output:** Figure 2 PNG in `analysis/` (e.g. `Fig2_Risk_Heatmap.png`).

---

### 3.4 generate_clock_validation.py

**Purpose:** Generate the methylation clock validation plot (Figure 3).

**Responsibilities:** Adds project root to `sys.path`, then calls `create_clock_validation_plot()` from `methylation_visualizations.py`.

**Output:** Figure 3 PNG in `analysis/` (e.g. `Fig3_Clock_Validation.png`).

---

### 3.5 generate_methylation_visualizations.py

**Purpose:** Generate the main methylation pipeline visualizations.

**Responsibilities:** Calls `generate_all_visualizations()` from `methylation_visualizations.py`, which produces:
1. Pipeline workflow diagram  
2. Example DMR visualizations  
3. Pipeline summary diagram  

**Output:** Multiple PNGs in `analysis/` (workflow, example DMR, summary).

---

### 3.6 generate_pipeline_diagram.py

**Purpose:** Generate the bioinformatics pipeline architecture diagram (diagrams library).

**Responsibilities:**
- Uses the `diagrams` library (Graphviz) to draw: Nanopore pod5 → Dorado/Modkit-style flow, Parquet storage, Polars/DuckDB, Dagster orchestration.
- Saves diagram in `analysis/` (e.g. `Bioinformatics_Pipeline_Diagram.png`).

**Output:** Pipeline diagram PNG (and optional Graphviz source) in `analysis/`.

---

### 3.7 install_graphviz.sh

**Purpose:** Install Graphviz on macOS for diagram generation.

**Responsibilities:**
- Detects Homebrew (`brew`, or `/opt/homebrew/bin/brew`, `/usr/local/bin/brew`).
- Runs `brew install graphviz`.
- If Homebrew is missing, prints install instructions and exits with an error.

**Usage:** Run once on macOS when Graphviz is needed for `generate_agent_system_schema.py` or `generate_pipeline_diagram.py`.

---

### 3.8 README_GRAPHVIZ.md

**Purpose:** Short documentation for Graphviz setup and diagram scripts (location and usage).

---

### 3.9 mock_ukb_generator.py

**Purpose:** Generate synthetic UK Biobank-style tabular data for pipeline development and testing.

**Responsibilities:**
- Produces fake data with `Sample_ID`, `Age`, `Sex`, `BMI`, `AD_diagnosis`, `EAA` (Epigenetic Age Acceleration), and 5 dummy SNP columns (0, 1, 2).
- Uses `MOCK_` prefix for Sample_IDs (whitelisted by UK Biobank pre-commit security hook).
- Typer CLI: `--n-samples`, `--output`, `--seed`. Python API `generate_synthetic_ukb_data()` accepts `min_age` / `max_age` and other sampling parameters.

**Dependencies:** pandas, numpy, typer.  
**Output:** `test_data/mock_clinical_data.csv` (default).  
**Tests:** `tests/test_mock_clinical_csv.py`.  
**Related doc:** [docs/SYNTHETIC_UKB_GENERATOR.md](SYNTHETIC_UKB_GENERATOR.md).

---

### 3.10 train_romanian_epigenetic_clock.py

**Purpose:** Train a custom epigenetic aging clock with Elastic Net (`ElasticNetCV`) on a methylation feature matrix and chronological age metadata; Romanian-style mock cohort IDs when generating synthetic data.

**Responsibilities:**
- Load `methylation_matrix.csv` and `metadata.csv` (inner join on `sample_id`) with Polars, or write deterministic mock CSVs under `data/mock_romanian_cohort/` when files are missing.
- Fit `Pipeline(StandardScaler, ElasticNetCV)` with 5-fold CV over `l1_ratio` and `alpha`.
- Report test-set MAE and Pearson r; save scatter plot (chronological vs predicted age) under `figures/` (git-ignored).

**Dependencies:** polars, numpy, scipy, scikit-learn, matplotlib, typer.  
**Related doc:** [docs/ROMANIAN_EPIGENETIC_CLOCK.md](ROMANIAN_EPIGENETIC_CLOCK.md).

---

### 3.11 validate_clock.py

**Purpose:** Evaluate a **pre-trained** epigenetic clock (`joblib` / `pickle`) on a held-out table with `chronological_age` and `cg*` feature columns.

**Responsibilities:**
- Load model and test Parquet/CSV; align features to `feature_names_in_` when present; mean-impute missing expected CpGs from test-set statistics.
- Compute overall MAE and Pearson r; stratify MAE by age decade (`pandas.cut`).
- Write `Fig_Clock_Residuals.png`, `Fig_Clock_MAE_by_decade.png`, and `validation_metrics.json` under `--output_dir`.

**Dependencies:** pandas, numpy, scipy, scikit-learn, matplotlib, seaborn, joblib (via sklearn), pyarrow for Parquet.  
**Related doc:** [docs/ROMANIAN_EPIGENETIC_CLOCK.md](ROMANIAN_EPIGENETIC_CLOCK.md#held-out-validation-validate_clockpy).

---

### 3.12 ukb_la_snp_lookup.py

**Purpose:** Offline manifest builder for UK Biobank-style genotype extraction planning: read Excel (`Gene`, `SNP_rsID`), query the Ensembl Variation REST API for **GRCh38** coordinates, emit CSV with chromosome, position, and an imputed-bulk chunk label (no DNAnexus / dx-toolkit).

**Responsibilities:** argparse CLI, rate-limited HTTP with retries on 429/5xx, pandas + openpyxl input, CSV output under `analysis/` by default.

**Dependencies:** pandas, requests, openpyxl.  
**Related:** [README.md](../README.md), [UKB_PRE_COMMIT_HOOK.md](UKB_PRE_COMMIT_HOOK.md), exploratory QA notebook `notebooks/05_ukb_exploration/UKB_LA_SNP_FirstContact.ipynb` (§4.5).

---

### 3.13 render_longevity_network_diagram.py

**Purpose:** Render the longevity conceptual network diagram with matplotlib to mirror `frontend/src/components/LongevityNetworkDiagram.tsx`.

**Responsibilities:** Typer CLI; default output `figures/longevity_network_diagram.png` (typically git-ignored; override path as needed).

**Dependencies:** matplotlib, typer.

---

### 3.14 render_figure1c_mechanisms_network.py

**Purpose:** Figure 1, panel C — static hub-and-spoke network of LA-SNP mechanisms (networkx + matplotlib + seaborn + adjustText).

**Responsibilities:** Writes `Figure1C_Mechanisms.png` and `.pdf` under `analysis/` (default `--out-dir`).

**Dependencies:** networkx, matplotlib, seaborn, adjustText, typer.

---

### 3.15 render_dashboard_figure_mockup.py

**Purpose:** Matplotlib static export mirroring `components/DashboardFigureMockup.tsx` when Node/Playwright is unavailable.

**Responsibilities:** Typer CLI; default `analysis/dashboard_figure_mockup.png`.

**Dependencies:** matplotlib, numpy, typer.

---

### 3.16 bootstrap_r_env.sh

**Purpose:** Create a local micromamba environment with conda-forge `r-base` under `.r-env/` for reproducible R workflows without modifying system R.

**Usage:** `./scripts/bootstrap_r_env.sh` from repo root. `.r-env/` is git-ignored.

---

### 3.17 security_check.sh & install_pre_commit_hook.sh

**Purpose:** UK Biobank pre-commit security hook — blocks commits containing `patient_id`, `UKB_`, or `.vcf`/`.bed` files.

**Responsibilities (exemptions):** Content scanning skips `docs/*`, root `README.md`, `notebooks/README.md`, `notebooks/05_ukb_exploration/*`, the hook scripts themselves, `scripts/mock_ukb_generator.py`, `test_data/mock_clinical_data.csv`, **`scripts/ukb_la_snp_lookup.py`** (offline manifest only; must never hold participant IDs), and **`repo_structure.txt`** (generated tree listing). Reinstall the hook after editing `security_check.sh`.

**Usage:** Run `./scripts/install_pre_commit_hook.sh` to install.  
**Related doc:** [docs/UKB_PRE_COMMIT_HOOK.md](UKB_PRE_COMMIT_HOOK.md).

---

### 3.18 generate_la_snp_per_gene_plot.py

**Purpose:** Supplementary manuscript figure — horizontal bar chart of **unique LA-SNPs per gene** from an Excel overlap table.

**Responsibilities:** Read `.xlsx` with pandas; group by gene column and count unique SNP identifiers; sort by count (descending); highlight genes with at least three SNPs; annotate bar ends; journal-style matplotlib (sans-serif, no top/right spines); write `analysis/Fig_LA_SNPs_per_gene.png` at 300 DPI.

**CLI:** `argparse` defaults: `--input` repo-root `overlapping_genes_with_snps.xlsx`, `--output` `analysis/Fig_LA_SNPs_per_gene.png`, `--gene-column` `Gene`, `--snp-column` `SNP_rsID` (override when the spreadsheet uses other headers, for example `Gene Symbol` / `SNP Identifier`).

**Dependencies:** pandas, openpyxl, matplotlib.

---

## 4. Jupyter notebooks (notebooks/)

Notebooks are grouped by functional area in numbered subfolders. Run with `uv run jupyter lab` from the project root. Large data should live in the root `data/` directory (git-ignored).

### 4.1 notebooks/01_genomics_analysis/

Notebooks for genomic data analysis, gene-list exploration, and network analysis.

| File | Purpose |
|------|--------|
| **AlphaGenome.ipynb** | Comprehensive analysis of AD/PD gene lists using the AlphaGenome API; longevity genes, network analysis, and visualization. |
| **AlphaGenome_updated.ipynb** | Updated version with enhanced network visualizations and AlphaGenome API usage. |

---

### 4.2 notebooks/02_methylation_pipeline/

Notebooks for processing and analyzing DNA methylation data from Oxford Nanopore sequencing.

| File | Purpose |
|------|--------|
| **DownstreamMethylationAnalysis.ipynb** | Interactive downstream analysis and DMR calling; complements `downstream_analysis.R`. |
| **MethylationClocks.ipynb** | Epigenetic clock exploration and validation: chronological vs DNAm age, MAE/RMSE/R²; overview of Horvath, Hannum, PhenoAge, GrimAge, DunedinPACE; foundation for 450K/EPIC array data. |

---

### 4.3 notebooks/03_validation_and_compliance/

Tools for data quality, code correctness, and regulatory compliance.

| File | Purpose |
|------|--------|
| **UKB_Compliance_Auditor.ipynb** | Scanner for restricted UK Biobank identifiers (EIDs) before public sharing; run before pushing analysis to public portals. |
| **Validations.ipynb** | General pipeline validation and quality-control checks. |

**Related doc:** `docs/UKB_COMPLIANCE_AUDITOR.md` — Documentation for the UKB compliance auditor.

---

### 4.4 notebooks/04_exploratory_visualizations/

Notebooks for project-wide visualizations and heatmaps.

| File | Purpose |
|------|--------|
| **Visualizations.ipynb** | Centralized notebook for generating publication-ready figures and exploratory plots. |

---

### 4.5 notebooks/05_ukb_exploration/

Exploratory QA for the offline UK Biobank longevity-associated SNP manifest (`scripts/ukb_la_snp_lookup.py` → `analysis/ukb_snp_manifest_v0.1.csv`). No UKB participant data.

| File | Purpose |
|------|--------|
| **UKB_LA_SNP_FirstContact.ipynb** | Manifest shape and sample rows; GRCh38 resolution coverage vs failures; SNPs per chromosome (bar plot); gene-level SNP counts and chromosomes; min–max GRCh38 position per chromosome for extraction planning; short Markdown summary. |

**Related:** `scripts/ukb_la_snp_lookup.py` (§3.12), [README.md](../README.md) (UK Biobank SNP manifest).

---

### 4.6 notebooks/README.md

**Purpose:** Describes the notebook directory structure, the role of each subfolder, guidelines (data locality, environment, compliance), and a short summary of each notebook.

---

## 5. Documentation (docs/)

| File | Purpose |
|------|--------|
| **CODE_MODULES_REFERENCE.md** | This document — reference for all code files and modules. |
| **PROJECT_STRUCTURE.md** | Bioinformatics project directory layout. |
| **UKB_PRE_COMMIT_HOOK.md** | Git pre-commit security hook for UK Biobank compliance. |
| **SYNTHETIC_UKB_GENERATOR.md** | Synthetic UK Biobank data generator (`mock_ukb_generator.py`). |
| **ROMANIAN_EPIGENETIC_CLOCK.md** | Romanian cohort Elastic Net clock (`train_romanian_epigenetic_clock.py`) and held-out validation (`validate_clock.py`). |
| **UKB_COMPLIANCE_AUDITOR.md** | UK Biobank compliance auditor (used with `03_validation_and_compliance/UKB_Compliance_Auditor.ipynb`). |
| **METHYLATION_PIPELINE_QUICK_REFERENCE.md** | Quick reference for the methylation pipeline. |
| **METHYLATION_PIPELINE_USAGE.md** | Detailed usage and workflow for the methylation pipeline. |

---

## 6. Configuration and environment

| File | Purpose |
|------|--------|
| **pyproject.toml** | Python project config (uv): package name `rogen-aging`, `[build-system]` (setuptools), dependencies including pysam and streamlit, optional `dev` extra (pytest), test paths / `pythonpath` (`src` and `scripts`) for pytest, Python ≥3.12. |
| **setup.py** | Setuptools shim; metadata and dependencies are defined in `pyproject.toml`. |
| **requirements.txt** | Optional pip-oriented list with loose pins for a subset of scientific stack. |
| **.env.example** | Example environment variables (e.g. API keys or paths) for local runs. |
| **.vscode/settings.json** | Editor settings (e.g. R path, Python interpreter) for the workspace. |

---

## 7. Summary: module vs responsibility

| Module / File | Main responsibility |
|---------------|---------------------|
| `downstream_analysis.R` | DMR calling from bedMethyl (DMRcaller workflow). |
| `find_r.sh` | Find R on macOS for IDE. |
| `pipeline_validation.sh` | Validate ONT pipeline (Dorado + Modkit + test data). |
| `src/rogen_aging/__init__.py` | Public API and submodule exports for `import rogen_aging`. |
| `src/rogen_aging/pipeline/` | Placeholder subpackage for shared pipeline code. |
| `src/rogen_aging/methylation_visualizations.py` | All methylation pipeline and clock figures. |
| `src/rogen_aging/network_visualizer.py` | Protein interaction network figure. |
| `src/rogen_aging/eda_dashboard/` | Streamlit EDA dashboard for merged multi-omics Parquet. |
| `tests/` | Pytest tests (e.g. package import smoke tests). |
| `scripts/generate_agent_system_schema*.py` | Figure 4 — agent system architecture. |
| `scripts/generate_bimodal_heatmap.py` | Figure 2 — bimodal risk heatmap. |
| `scripts/generate_clock_validation.py` | Figure 3 — clock validation. |
| `scripts/generate_methylation_visualizations.py` | Pipeline workflow + example DMR + summary diagrams. |
| `scripts/generate_pipeline_diagram.py` | Bioinformatics pipeline (diagrams/Graphviz). |
| `scripts/install_graphviz.sh` | Install Graphviz on macOS. |
| `scripts/mock_ukb_generator.py` | Synthetic UK Biobank-style mock data. |
| `scripts/train_romanian_epigenetic_clock.py` | Elastic Net epigenetic clock (mock Romanian cohort / custom CSVs). |
| `scripts/validate_clock.py` | Held-out validation for a saved clock model (MAE, r, decade MAE, figures). |
| `scripts/ukb_la_snp_lookup.py` | Offline UKB SNP manifest via Ensembl GRCh38 (CSV). |
| `scripts/generate_la_snp_per_gene_plot.py` | Supplementary LA-SNPs-per-gene bar chart (Excel → PNG). |
| `scripts/render_longevity_network_diagram.py` | Matplotlib longevity network (twin of `frontend/` TSX). |
| `scripts/render_figure1c_mechanisms_network.py` | Figure 1C mechanisms network (PNG/PDF). |
| `scripts/render_dashboard_figure_mockup.py` | Matplotlib dashboard mockup (twin of `components/` TSX). |
| `scripts/bootstrap_r_env.sh` | Optional micromamba R under `.r-env/`. |
| `scripts/security_check.sh` | UK Biobank pre-commit security hook. |
| `components/`, `frontend/` | Vite + React manuscript figure apps (see §8). |
| `notebooks/01_genomics_analysis/` | AlphaGenome AD/PD gene analysis and networks. |
| `notebooks/02_methylation_pipeline/` | Methylation downstream analysis and epigenetic clocks. |
| `notebooks/03_validation_and_compliance/` | UKB compliance auditor and pipeline validations. |
| `notebooks/04_exploratory_visualizations/` | Publication-ready and exploratory figures. |
| `notebooks/05_ukb_exploration/` | UKB LA-SNP CSV manifest sanity checks before extraction. |

For pipeline usage and quick commands, see **docs/METHYLATION_PIPELINE_QUICK_REFERENCE.md** and **docs/METHYLATION_PIPELINE_USAGE.md**.

---

## 8. TypeScript figure apps (`components/`, `frontend/`)

Small **Vite + React** projects for publication-quality diagrams. Use `npm install`, `npm run dev` for interactive editing, and package-specific capture scripts (often `npm run capture` or `node scripts/capture-diagram.mjs`) where Playwright is configured. Python `scripts/render_*.py` counterparts exist for matplotlib-only pipelines.

| Path | Role |
|------|------|
| `components/DashboardFigureMockup.tsx` | Multi-panel dashboard mockup component |
| `components/dashboard-figure-render/` | Vite app + Playwright capture for the dashboard mockup |
| `frontend/src/components/LongevityNetworkDiagram.tsx` | Longevity network layout |
| `frontend/scripts/capture-diagram.mjs` | Headless PNG export for the longevity diagram |
