# ROGEN Aging Project — Code Modules Reference

This document describes each code file and module in the IBAR-ROGEN Aging project and what it is responsible for.

---

## 1. Project Root (Aging/)

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

## 2. Python Package: src/rogen_aging/

### 2.1 __init__.py

**Purpose:** Package initializer for `rogen_aging`.

**Responsibilities:** Declares the package and its public API (currently `__all__ = []`). Visualization and network code live in the other modules and are run via scripts or notebooks.

---

### 2.2 methylation_visualizations.py

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

### 2.3 network_visualizer.py

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

## 3. Scripts (scripts/)

Scripts are entry points that call into `src/rogen_aging` or external tools. They are run from the `Aging/` directory (e.g. `python scripts/generate_*.py` or `uv run python scripts/...`).

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

## 4. Jupyter Notebooks (notebooks/)

Notebooks contain interactive analyses and may generate figures or tables; they are not listed as “modules” in the same way as scripts but are part of the codebase.

| Notebook | Purpose |
|----------|--------|
| **AlphaGenome.ipynb** | Longevity gene analysis using the AlphaGenome API; genes linked to Alzheimer’s and Parkinson’s; network analysis and visualization. |
| **AlphaGenome_updated.ipynb** | Updated version of the AlphaGenome longevity/AD/PD analysis and network exploration. |
| **DownstreamMethylationAnalysis.ipynb** | Downstream methylation analysis (complements `downstream_analysis.R` in an interactive, notebook form). |
| **MethylationClocks.ipynb** | Epigenetic clock analysis: chronological vs DNAm age, MAE/RMSE/R², overview of Horvath, Hannum, PhenoAge, GrimAge, DunedinPACE; basis for array data (450K/EPIC). |
| **Validations.ipynb** | Validation analyses (e.g. method or result checks). |
| **Visualizations.ipynb** | Ad-hoc visualizations and plots. |

**notebooks/README.md** — Describes suggested naming, data location (`data/`), and the available notebooks (AlphaGenome, MethylationClocks).

---

## 5. Configuration and Environment

| File | Purpose |
|------|--------|
| **pyproject.toml** | Python project config (uv): package name `rogen-aging`, dependencies (alphagenome, pandas, matplotlib, requests, polars, python-dotenv, networkx, numpy, scipy, scikit-learn, seaborn, diagrams), dev deps (ipykernel, jupyterlab), Python ≥3.12. |
| **.env.example** | Example environment variables (e.g. API keys or paths) for local runs. |
| **.vscode/settings.json** | Editor settings (e.g. R path, Python interpreter) for the workspace. |

---

## 6. Summary: Module vs Responsibility

| Module / File | Main responsibility |
|---------------|---------------------|
| `downstream_analysis.R` | DMR calling from bedMethyl (DMRcaller workflow). |
| `find_r.sh` | Find R on macOS for IDE. |
| `pipeline_validation.sh` | Validate ONT pipeline (Dorado + Modkit + test data). |
| `src/rogen_aging/methylation_visualizations.py` | All methylation pipeline and clock figures. |
| `src/rogen_aging/network_visualizer.py` | Protein interaction network figure. |
| `scripts/generate_agent_system_schema*.py` | Figure 4 — agent system architecture. |
| `scripts/generate_bimodal_heatmap.py` | Figure 2 — bimodal risk heatmap. |
| `scripts/generate_clock_validation.py` | Figure 3 — clock validation. |
| `scripts/generate_methylation_visualizations.py` | Pipeline workflow + example DMR + summary diagrams. |
| `scripts/generate_pipeline_diagram.py` | Bioinformatics pipeline (diagrams/Graphviz). |
| `scripts/install_graphviz.sh` | Install Graphviz on macOS. |
| **Notebooks** | Interactive analyses (AlphaGenome, methylation clocks, downstream methylation, validations, visualizations). |

For pipeline usage and quick commands, see **docs/METHYLATION_PIPELINE_QUICK_REFERENCE.md** and **docs/METHYLATION_PIPELINE_USAGE.md**.
