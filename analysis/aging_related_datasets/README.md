# Aging-related NGS datasets

This folder contains tools to discover NGS sequencing datasets relevant to aging research.

## search_aging_ngs_datasets.py

Searches **NCBI SRA** (Sequence Read Archive) for aging-related sequencing studies and writes a CSV table with:

- **accession** – SRA run/experiment/study ID (e.g. SRR, SRX, SRP)
- **title** – Record title
- **description** – Dataset description (from SRA Summary, Abstract, or Design when available)
- **link** – URL to the dataset on NCBI SRA
- **number_of_samples** – Number of samples/runs (from SRA when available, otherwise 1 for run-level records)
- **platform** – Normalized platform name (e.g. Oxford Nanopore, Illumina, PacBio)
- **platform_raw** – Raw platform string from SRA
- **is_oxford_nanopore** – `True` if the dataset was generated with Oxford Nanopore; `False` otherwise
- **organism** – Organism (if available)
- **study_accession** – Parent study accession (if available)

### Usage

From the **Aging** project root:

```bash
uv run python analysis/aging_related_datasets/search_aging_ngs_datasets.py
```

The script runs **two SRA searches** and writes **two CSV files**:

1. **All datasets** – Search with aging + sequencing terms (no platform filter) → `aging_ngs_datasets.csv`
2. **Oxford Nanopore only** – Same aging terms **plus** `AND (nanopore OR "Oxford Nanopore")` in the query, so the database returns only Nanopore datasets → `aging_ngs_datasets_nanopore_only.csv`

The Nanopore-only file is filled from the second search, not by filtering the first result set in memory.

Options:

- `--max-results N` – Maximum number of records to fetch (default: 100)
- `--output FILE` – Output CSV path for all datasets (default: `aging_ngs_datasets.csv`)
- `--nanopore-output FILE` – Output CSV path for Oxford Nanopore–only datasets (default: same dir as `--output`, suffix `_nanopore_only.csv`)
- `--query "..."` – Override the search query (default: aging + sequencing terms)

### Requirements

- **requests** (in project dependencies)
- Optional: **NCBI_API_KEY** in `.env` for higher rate limits ([get a key](https://www.ncbi.nlm.nih.gov/account/settings/))

### Search terms

The default query combines aging-related terms (aging, ageing, longevity, senescence, biological age, epigenetic age, DNA methylation age) with sequencing context (sequencing, RNA-Seq, whole genome, WGS, methylation). Use `--query` to customize.
