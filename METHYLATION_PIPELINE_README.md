# ROGEN Methylation Calling Pipeline - User Guide

**Project:** ROGEN Aging Research  
**Activity:** 2.1.8.1 - Methylation Calling Pipeline  
**Last Updated:** January 8, 2025

## Table of Contents

1. [Overview](#overview)
2. [Pipeline Architecture](#pipeline-architecture)
3. [Prerequisites](#prerequisites)
4. [Installation](#installation)
5. [Quick Start](#quick-start)
6. [Detailed Usage](#detailed-usage)
7. [Output Files](#output-files)
8. [Troubleshooting](#troubleshooting)
9. [Advanced Usage](#advanced-usage)
10. [References](#references)

---

## Overview

This pipeline provides a complete workflow for analyzing DNA methylation patterns from Oxford Nanopore sequencing data. It integrates three main tools:

1. **Dorado** - Basecalling with methylation-aware models
2. **Modkit** - Conversion of BAM files to bedMethyl format
3. **DMRcaller** - Identification of Differentially Methylated Regions (DMRs)

### Pipeline Workflow

```
POD5 files
    ↓
[Dorado] Basecalling with methylation model (5mC/5hmC)
    ↓
BAM files with MM/ML tags
    ↓
[Modkit] Extract methylation calls
    ↓
bedMethyl files
    ↓
[DMRcaller] Differential methylation analysis
    ↓
DMRs (Differentially Methylated Regions)
```

---

## Pipeline Architecture

### Components

1. **`pipeline_validation.sh`** - Bash script for basecalling and methylation extraction
2. **`downstream_analysis.R`** - R script for DMR calling and analysis
3. **`notebooks/DownstreamMethylationAnalysis.ipynb`** - Interactive R notebook version

### Data Flow

```
Input: POD5 files (raw Nanopore data)
  ↓
Step 1: Dorado basecalling → BAM with methylation tags
  ↓
Step 2: Modkit extraction → bedMethyl files
  ↓
Step 3: DMRcaller analysis → DMR results
  ↓
Output: BED files, CSV summaries, visualizations
```

---

## Prerequisites

### Required Software

1. **Dorado** (v5.0.0 or later)
   - Download from: https://github.com/nanoporetech/dorado
   - Required for basecalling with methylation models

2. **Modkit**
   - Download from: https://github.com/nanoporetech/modkit
   - Required for BAM to bedMethyl conversion

3. **R** (v4.0 or later)
   - Required for downstream analysis
   - Install from: https://www.r-project.org/

4. **Bioconductor Packages**
   - DMRcaller
   - GenomicRanges
   - rtracklayer

5. **Optional Tools**
   - `samtools` - For BAM file manipulation
   - `wget` or `curl` - For downloading test datasets

### System Requirements

- **GPU**: NVIDIA GPU recommended for Dorado basecalling (CUDA support)
- **RAM**: Minimum 16GB recommended for large datasets
- **Storage**: Sufficient space for POD5 files, BAM files, and bedMethyl outputs

---

## Installation

### 1. Install Dorado

```bash
# Download and install Dorado
# Follow instructions at: https://github.com/nanoporetech/dorado

# Verify installation
dorado --version
```

### 2. Install Modkit

```bash
# Download and install Modkit
# Follow instructions at: https://github.com/nanoporetech/modkit

# Verify installation
modkit --version
```

### 3. Install R and Bioconductor Packages

```bash
# Install R (if not already installed)
# macOS: brew install r
# Ubuntu: sudo apt-get install r-base

# Launch R and install Bioconductor
R

# In R console:
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("DMRcaller", "GenomicRanges", "rtracklayer"))
```

### 4. Download Methylation Model (for Dorado)

```bash
# Download the R10.4.1 methylation model
dorado download --model dna_r10.4.1_e8.2_400bps_fast@v5.0.0
```

---

## Quick Start

### Step 1: Run Pipeline Validation Script

```bash
# Make script executable (if not already)
chmod +x pipeline_validation.sh

# Run validation script
./pipeline_validation.sh
```

This script will:
- Check if Dorado and Modkit are installed
- Download a test dataset
- Run basecalling (if GPU available)
- Convert BAM to bedMethyl format

### Step 2: Run Downstream Analysis

**Option A: Using R Script**

```bash
Rscript downstream_analysis.R
```

**Option B: Using R Notebook**

```bash
# Launch JupyterLab
uv run jupyter lab

# Open notebooks/DownstreamMethylationAnalysis.ipynb
# Run cells sequentially
```

---

## Detailed Usage

### Part 1: Basecalling and Methylation Extraction (`pipeline_validation.sh`)

#### Basic Usage

```bash
./pipeline_validation.sh
```

#### What the Script Does

1. **Prerequisites Check**
   - Verifies Dorado and Modkit installation
   - Displays tool versions

2. **Test Data Download**
   - Downloads ONT test dataset from Epi2Me Labs
   - Extracts POD5 files

3. **Basecalling** (commented out by default)
   - Uses model: `dna_r10.4.1_e8.2_400bps_fast@v5.0.0`
   - Enables 5mC/5hmC detection with `--modified-bases 5mC_5hmC`
   - Outputs BAM file with MM/ML tags

4. **Methylation Extraction** (commented out by default)
   - Converts BAM to bedMethyl format
   - Requires reference genome FASTA

#### Customizing the Script

**To run actual basecalling**, uncomment lines ~100-110:

```bash
dorado basecaller \
    "$DORADO_MODEL" \
    --modified-bases 5mC_5hmC \
    --device cuda:0 \
    "$POD5_DIR" \
| samtools view -bS - \
> "$OUTPUT_BAM"
```

**To run Modkit**, uncomment lines ~140-150 and provide reference:

```bash
REFERENCE_FASTA="path/to/reference.fasta"

modkit extract \
    "$OUTPUT_BAM" \
    "$OUTPUT_BEDMETHYL" \
    --ref "$REFERENCE_FASTA" \
    --bedgraph \
    --combine-strands \
    --filter-threshold 0.0
```

#### Script Parameters

| Variable | Default | Description |
|----------|---------|-------------|
| `DORADO_MODEL` | `dna_r10.4.1_e8.2_400bps_fast@v5.0.0` | Dorado methylation model |
| `OUTPUT_BAM` | `basecalled_methylation.bam` | Output BAM filename |
| `OUTPUT_BEDMETHYL` | `methylation_calls.bedMethyl` | Output bedMethyl filename |
| `DEVICE` | `cuda:0` (or `cpu`) | Computing device for basecalling |

#### Flag Explanations

**Dorado Flags:**
- `--modified-bases 5mC_5hmC`: Enables detection of 5-methylcytosine and 5-hydroxymethylcytosine modifications
- `--device cuda:0`: Uses NVIDIA GPU (use `cpu` if no GPU available)

**Modkit Flags:**
- `--ref`: Reference genome FASTA file (required)
- `--bedgraph`: Also generate bedGraph format output
- `--combine-strands`: Combine forward and reverse strand calls
- `--filter-threshold 0.0`: Minimum modification probability (0.0 = include all)

---

### Part 2: Downstream Analysis (`downstream_analysis.R`)

#### Basic Usage

```bash
# Run the complete workflow
Rscript downstream_analysis.R

# Or launch R interactively
R
source("downstream_analysis.R")
```

#### Step-by-Step Workflow

**Step 1: Load Libraries**

```r
library(DMRcaller)
library(GenomicRanges)
library(rtracklayer)
```

**Step 2: Import bedMethyl Files**

```r
# Import bedMethyl files for each sample
young_sample_1 <- import_bedmethyl("data/young_sample1.bedMethyl")
young_sample_2 <- import_bedmethyl("data/young_sample2.bedMethyl")
old_sample_1 <- import_bedmethyl("data/old_sample1.bedMethyl")
old_sample_2 <- import_bedmethyl("data/old_sample2.bedMethyl")
```

**Step 3: Prepare Data for DMRcaller**

```r
# Prepare data matrices
young_samples <- list(
    prepare_dmrcaller_data(young_sample_1, "Young_1"),
    prepare_dmrcaller_data(young_sample_2, "Young_2")
)

old_samples <- list(
    prepare_dmrcaller_data(old_sample_1, "Old_1"),
    prepare_dmrcaller_data(old_sample_2, "Old_2")
)
```

**Step 4: Calculate DMRs**

```r
dmrs <- calculate_dmrs(
    group1_data = young_samples,
    group2_data = old_samples,
    group1_name = "Young",
    group2_name = "Old",
    min_coverage = 5,
    p_value_threshold = 0.01,
    min_cpg = 3
)
```

**Step 5: Export Results**

```r
# Export to BED file
export(dmrs, "dmrs_young_vs_old.bed", format = "bed")

# Export summary statistics
dmr_summary <- data.frame(
    chr = seqnames(dmrs),
    start = start(dmrs),
    end = end(dmrs),
    width = width(dmrs),
    n_cpg = mcols(dmrs)$n_cpg,
    p_value = mcols(dmrs)$p_value,
    mean_meth_diff = mcols(dmrs)$mean_meth_diff
)
write.csv(dmr_summary, "dmr_summary.csv", row.names = FALSE)
```

#### Function Reference

**`import_bedmethyl(bedmethyl_path)`**
- Imports bedMethyl file into GRanges object
- Parameters:
  - `bedmethyl_path`: Path to bedMethyl file
- Returns: GRanges object with methylation data

**`prepare_dmrcaller_data(bedmethyl_granges, sample_name)`**
- Prepares data for DMRcaller analysis
- Parameters:
  - `bedmethyl_granges`: GRanges object from `import_bedmethyl()`
  - `sample_name`: Name identifier for the sample
- Returns: List with coverage and methylation matrices

**`calculate_dmrs(group1_data, group2_data, ...)`**
- Identifies DMRs between two groups
- Parameters:
  - `group1_data`: List of prepared data for group 1
  - `group2_data`: List of prepared data for group 2
  - `group1_name`: Name of group 1 (default: "Group1")
  - `group2_name`: Name of group 2 (default: "Group2")
  - `min_coverage`: Minimum coverage per position (default: 5)
  - `p_value_threshold`: P-value threshold for significance (default: 0.01)
  - `min_cpg`: Minimum CpG sites per DMR (default: 3)
- Returns: GRanges object with identified DMRs

---

### Part 3: Using the R Notebook

#### Launching the Notebook

```bash
# Start JupyterLab
uv run jupyter lab

# Navigate to notebooks/DownstreamMethylationAnalysis.ipynb
```

#### Notebook Structure

1. **Introduction** - Overview and workflow steps
2. **Section 1** - Package installation and loading
3. **Section 2** - Data import functions
4. **Section 3** - Data preparation functions
5. **Section 4** - DMR calling function
6. **Section 5** - Complete workflow example
7. **Section 6** - Helper functions

#### Running the Notebook

1. Execute cells sequentially from top to bottom
2. Uncomment code blocks when you have data available
3. Modify file paths to match your data structure
4. Adjust parameters based on your experimental design

---

## Output Files

### From `pipeline_validation.sh`

- **`basecalled_methylation.bam`** - BAM file with MM/ML methylation tags
- **`methylation_calls.bedMethyl`** - bedMethyl format file with methylation calls
- **`wf-basecalling-demo/`** - Test dataset directory (if downloaded)

### From `downstream_analysis.R`

- **`dmrs_young_vs_old.bed`** - BED file with DMR coordinates
- **`dmr_summary.csv`** - CSV file with DMR statistics
- **`dmr_visualization.pdf`** - Visualization plots (if generated)

### Pipeline Visualizations

The repository includes visualization scripts that generate diagrams and example plots:

**Generated Visualizations:**
- **`analysis/Methylation_Pipeline_Workflow.png`** - Workflow diagram showing data flow through the pipeline
- **`analysis/Example_DMR_Visualizations.png`** - Example DMR analysis plots (Manhattan plot, distributions)
- **`analysis/Methylation_Pipeline_Summary.png`** - Component overview diagram

**To Generate Visualizations:**

```bash
# Option 1: Using the script
uv run python scripts/generate_methylation_visualizations.py

# Option 2: Direct Python module
uv run python -m src.rogen_aging.methylation_visualizations
```

These visualizations are useful for:
- Understanding the pipeline architecture
- Demonstrating expected output formats
- Documentation and presentations
- Quality control and result interpretation

### File Formats

**bedMethyl Format (BED9+3):**
```
chrom  start  end  name  score  strand  coverage  percent_methylated  count_methylated  count_unmethylated  ...
```

**DMR Summary CSV:**
```
chr, start, end, width, n_cpg, p_value, mean_meth_diff
chr1, 1000, 2000, 1000, 5, 0.001, 0.25
...
```

---

## Troubleshooting

### Common Issues

#### 1. Dorado Not Found

**Error:** `dorado: command not found`

**Solution:**
```bash
# Add Dorado to PATH or use full path
export PATH=$PATH:/path/to/dorado/bin
```

#### 2. GPU Not Available

**Error:** CUDA errors or slow basecalling

**Solution:**
```bash
# Edit pipeline_validation.sh, change device to CPU:
DEVICE="cpu"
```

#### 3. Modkit Requires Reference

**Error:** `modkit: error: the following arguments are required: --ref`

**Solution:**
```bash
# Download reference genome and specify path
wget https://example.com/reference.fasta
# Update REFERENCE_FASTA variable in script
```

#### 4. R Package Installation Fails

**Error:** Bioconductor packages fail to install

**Solution:**
```r
# Update BiocManager
BiocManager::install(version = "3.18")

# Install packages individually
BiocManager::install("DMRcaller")
BiocManager::install("GenomicRanges")
BiocManager::install("rtracklayer")
```

#### 5. bedMethyl Import Fails

**Error:** `Error in import(): invalid format`

**Solution:**
```r
# Check bedMethyl file format
head -5 your_file.bedMethyl

# Ensure it's BED9+3 format
# Adjust import parameters if needed
```

#### 6. Memory Issues

**Error:** Out of memory errors with large datasets

**Solution:**
- Process samples in batches
- Filter by chromosome or region
- Increase system RAM or use a computing cluster

---

## Advanced Usage

### Batch Processing Multiple Samples

```bash
# Process multiple POD5 directories
for pod5_dir in sample1 sample2 sample3; do
    dorado basecaller \
        "$DORADO_MODEL" \
        --modified-bases 5mC_5hmC \
        --device cuda:0 \
        "$pod5_dir" \
    | samtools view -bS - \
    > "${pod5_dir}_methylation.bam"
done
```

### Custom DMR Parameters

```r
# More stringent DMR calling
dmrs_strict <- calculate_dmrs(
    group1_data = young_samples,
    group2_data = old_samples,
    min_coverage = 10,        # Higher coverage requirement
    p_value_threshold = 0.001, # More stringent p-value
    min_cpg = 5               # More CpG sites per DMR
)

# Less stringent DMR calling
dmrs_relaxed <- calculate_dmrs(
    group1_data = young_samples,
    group2_data = old_samples,
    min_coverage = 3,         # Lower coverage requirement
    p_value_threshold = 0.05, # Less stringent p-value
    min_cpg = 2               # Fewer CpG sites per DMR
)
```

### Parallel Processing

```r
# Use parallel processing for multiple samples
library(parallel)

# Prepare data in parallel
cl <- makeCluster(4)
young_prepared <- parLapply(cl, young_samples, prepare_dmrcaller_data)
stopCluster(cl)
```

### Custom Visualizations

```r
# Create custom DMR visualizations
library(ggplot2)

# DMR size distribution
ggplot(dmr_summary, aes(x = width)) +
    geom_histogram(bins = 50) +
    labs(title = "DMR Size Distribution", x = "DMR Width (bp)")

# P-value distribution
ggplot(dmr_summary, aes(x = p_value)) +
    geom_histogram(bins = 50) +
    labs(title = "DMR P-value Distribution")
```

---

## References

### Tools and Packages

- **Dorado**: https://github.com/nanoporetech/dorado
- **Modkit**: https://github.com/nanoporetech/modkit
- **DMRcaller**: https://bioconductor.org/packages/DMRcaller/
- **GenomicRanges**: https://bioconductor.org/packages/GenomicRanges/
- **rtracklayer**: https://bioconductor.org/packages/rtracklayer/

### Documentation

- Dorado Documentation: https://github.com/nanoporetech/dorado#readme
- Modkit Documentation: https://github.com/nanoporetech/modkit#readme
- DMRcaller Vignette: Available in R via `vignette("DMRcaller")`

### Test Data

- ONT Test Dataset: https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/wf-basecalling/wf-basecalling-demo.tar.gz

### Citation

If you use this pipeline, please cite:

- Dorado: Nanopore Technologies
- Modkit: Nanopore Technologies  
- DMRcaller: Catoni et al., 2018 (see package citation in R)

---

## Support

For issues or questions:

1. Check the troubleshooting section above
2. Review tool-specific documentation
3. Check GitHub issues for known problems
4. Contact the ROGEN project team

---

**Last Updated:** January 8, 2025  
**Version:** 1.0  
**Maintained by:** ROGEN Project Team
