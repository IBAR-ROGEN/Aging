# Analysis Notebooks

This directory contains Jupyter notebooks for the ROGEN Aging Research project, organized by functional area.

## Directory Structure

### `01_genomics_analysis/`
Notebooks focused on genomic data analysis, gene list exploration, and network analysis.
- **AlphaGenome.ipynb**: Comprehensive analysis of AD/PD gene lists using the AlphaGenome API.
- **AlphaGenome_updated.ipynb**: Updated version with enhanced network visualizations.

### `02_methylation_pipeline/`
Notebooks for processing and analyzing DNA methylation data from Oxford Nanopore sequencing.
- **DownstreamMethylationAnalysis.ipynb**: Interactive R notebook for DMR calling and downstream analysis.
- **MethylationClocks.ipynb**: Exploration and validation of epigenetic clocks for biological age prediction.

### `03_validation_and_compliance/`
Tools for ensuring data quality, code correctness, and regulatory compliance.
- **UKB_Compliance_Auditor.ipynb**: Scanner for identifying restricted UK Biobank identifiers (EIDs) before public sharing.
- **Validations.ipynb**: General pipeline validation and quality control checks.

### `04_exploratory_visualizations/`
Notebooks dedicated to generating project-wide visualizations and heatmaps.
- **Visualizations.ipynb**: Centralized notebook for generating publication-ready figures.

### `05_ukb_exploration/`
Exploratory checks on the offline UK Biobank SNP manifest produced by `scripts/ukb_la_snp_lookup.py` (no participant data).
- **UKB_LA_SNP_FirstContact.ipynb**: Loads `analysis/ukb_snp_manifest_v0.1.csv`, reports GRCh38 resolution failures, chromosome and gene summaries, and per-chromosome position ranges to support chunk-based extraction planning.

## Guidelines
- **Data Locality**: Keep large datasets in the root `data/` directory (git-ignored).
- **Environment**: Use `uv run jupyter lab` to ensure all dependencies are available.
- **Compliance**: Always run the `UKB_Compliance_Auditor` before pushing new analysis to public portals.
