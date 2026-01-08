# Methylation Pipeline - Detailed Usage Guide

This document provides step-by-step instructions for using the ROGEN methylation calling pipeline.

## Table of Contents

1. [Before You Begin](#before-you-begin)
2. [Step-by-Step: Complete Workflow](#step-by-step-complete-workflow)
3. [Understanding the Output](#understanding-the-output)
4. [Customization Guide](#customization-guide)
5. [Best Practices](#best-practices)

---

## Before You Begin

### Check Your Environment

```bash
# Verify all tools are installed
which dorado
which modkit
which R
which samtools

# Check versions
dorado --version
modkit --version
R --version
```

### Prepare Your Data

1. **POD5 Files**: Ensure your POD5 files are in a directory
2. **Reference Genome**: Download reference FASTA file for your organism
3. **Sample Groups**: Organize samples into groups (e.g., Young vs Old)

### Directory Structure

Recommended structure:

```
project/
├── pod5_data/
│   ├── sample1/
│   ├── sample2/
│   └── sample3/
├── reference/
│   └── genome.fasta
├── output/
│   ├── bam/
│   ├── bedmethyl/
│   └── dmrs/
└── scripts/
    ├── pipeline_validation.sh
    └── downstream_analysis.R
```

---

## Step-by-Step: Complete Workflow

### Step 1: Basecalling with Dorado

#### 1.1 Download Methylation Model

```bash
# Download the R10.4.1 methylation model
dorado download --model dna_r10.4.1_e8.2_400bps_fast@v5.0.0
```

#### 1.2 Run Basecalling

```bash
# Single sample
dorado basecaller \
    dna_r10.4.1_e8.2_400bps_fast@v5.0.0 \
    --modified-bases 5mC_5hmC \
    --device cuda:0 \
    /path/to/pod5_files \
| samtools view -bS - \
> output/sample1_methylation.bam

# Multiple samples (using a loop)
for sample in sample1 sample2 sample3; do
    dorado basecaller \
        dna_r10.4.1_e8.2_400bps_fast@v5.0.0 \
        --modified-bases 5mC_5hmC \
        --device cuda:0 \
        pod5_data/${sample} \
    | samtools view -bS - \
    > output/bam/${sample}_methylation.bam
done
```

#### 1.3 Verify BAM Files

```bash
# Check BAM file has methylation tags
samtools view output/sample1_methylation.bam | head -1 | grep -o "MM\|ML"

# Count reads
samtools view -c output/sample1_methylation.bam
```

### Step 2: Extract Methylation Calls with Modkit

#### 2.1 Prepare Reference Genome

```bash
# Index reference genome (if not already indexed)
samtools faidx reference/genome.fasta
```

#### 2.2 Run Modkit

```bash
# Single sample
modkit extract \
    output/bam/sample1_methylation.bam \
    output/bedmethyl/sample1.bedMethyl \
    --ref reference/genome.fasta \
    --bedgraph \
    --combine-strands \
    --filter-threshold 0.0

# Multiple samples
for sample in sample1 sample2 sample3; do
    modkit extract \
        output/bam/${sample}_methylation.bam \
        output/bedmethyl/${sample}.bedMethyl \
        --ref reference/genome.fasta \
        --bedgraph \
        --combine-strands \
        --filter-threshold 0.0
done
```

#### 2.3 Verify bedMethyl Files

```bash
# Check file format
head -5 output/bedmethyl/sample1.bedMethyl

# Count methylation calls
wc -l output/bedmethyl/sample1.bedMethyl
```

### Step 3: Downstream Analysis with R

#### 3.1 Launch R Environment

```bash
# Option A: R console
R

# Option B: RStudio
# Open RStudio and set working directory

# Option C: R Notebook
uv run jupyter lab
```

#### 3.2 Load Required Libraries

```r
# Install if needed
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("DMRcaller", "GenomicRanges", "rtracklayer"))

# Load libraries
library(DMRcaller)
library(GenomicRanges)
library(rtracklayer)
library(ggplot2)  # For visualization
```

#### 3.3 Import bedMethyl Files

```r
# Set working directory
setwd("/path/to/project")

# Import files for each group
# Group 1: Young samples
young_samples <- list(
    import_bedmethyl("output/bedmethyl/young_sample1.bedMethyl"),
    import_bedmethyl("output/bedmethyl/young_sample2.bedMethyl"),
    import_bedmethyl("output/bedmethyl/young_sample3.bedMethyl")
)

# Group 2: Old samples
old_samples <- list(
    import_bedmethyl("output/bedmethyl/old_sample1.bedMethyl"),
    import_bedmethyl("output/bedmethyl/old_sample2.bedMethyl"),
    import_bedmethyl("output/bedmethyl/old_sample3.bedMethyl")
)
```

#### 3.4 Prepare Data for DMRcaller

```r
# Prepare data for each sample
young_prepared <- lapply(seq_along(young_samples), function(i) {
    prepare_dmrcaller_data(young_samples[[i]], paste0("Young_", i))
})

old_prepared <- lapply(seq_along(old_samples), function(i) {
    prepare_dmrcaller_data(old_samples[[i]], paste0("Old_", i))
})
```

#### 3.5 Calculate DMRs

```r
# Run DMR calling
dmrs <- calculate_dmrs(
    group1_data = young_prepared,
    group2_data = old_prepared,
    group1_name = "Young",
    group2_name = "Old",
    min_coverage = 5,
    p_value_threshold = 0.01,
    min_cpg = 3
)

# Check results
print(paste("Number of DMRs found:", length(dmrs)))
```

#### 3.6 Export Results

```r
# Export DMRs to BED file
export(dmrs, "output/dmrs/dmrs_young_vs_old.bed", format = "bed")

# Create summary table
dmr_summary <- data.frame(
    chr = seqnames(dmrs),
    start = start(dmrs),
    end = end(dmrs),
    width = width(dmrs),
    n_cpg = mcols(dmrs)$n_cpg,
    p_value = mcols(dmrs)$p_value,
    mean_meth_diff = mcols(dmrs)$mean_meth_diff
)

# Save summary
write.csv(dmr_summary, "output/dmrs/dmr_summary.csv", row.names = FALSE)
```

---

## Understanding the Output

### BAM File Structure

BAM files contain:
- **MM tag**: Modification probabilities for each base
- **ML tag**: Modification calls (methylated/unmethylated)

### bedMethyl File Format

Standard BED9+3 format:

```
chrom  start  end  name  score  strand  coverage  percent_methylated  count_methylated  count_unmethylated  ...
```

**Example:**
```
chr1  1000  1001  .  0  +  25  80.0  20  5  ...
```

- **chrom**: Chromosome name
- **start/end**: Genomic coordinates (0-based)
- **coverage**: Number of reads covering this position
- **percent_methylated**: Percentage of reads showing methylation
- **count_methylated**: Number of methylated reads
- **count_unmethylated**: Number of unmethylated reads

### DMR Results

DMR summary contains:
- **chr**: Chromosome
- **start/end**: DMR boundaries
- **width**: DMR size in base pairs
- **n_cpg**: Number of CpG sites in DMR
- **p_value**: Statistical significance
- **mean_meth_diff**: Mean methylation difference between groups

---

## Customization Guide

### Adjusting Basecalling Parameters

```bash
# Use different model
DORADO_MODEL="dna_r10.4.1_e8.2_400bps_sup@v5.0.0"  # Super accuracy model

# Use CPU instead of GPU
--device cpu

# Limit number of reads (for testing)
dorado basecaller ... | head -1000 | samtools view -bS - > test.bam
```

### Adjusting Modkit Parameters

```bash
# Filter by minimum coverage
--filter-threshold 0.5  # Only sites with >50% methylation probability

# Don't combine strands (keep separate)
# Remove --combine-strands flag

# Generate additional formats
--bedgraph  # BedGraph format
--wiggle    # Wiggle format
```

### Adjusting DMR Parameters

```r
# More stringent (fewer, higher confidence DMRs)
dmrs_strict <- calculate_dmrs(
    group1_data = young_prepared,
    group2_data = old_prepared,
    min_coverage = 10,        # Higher coverage
    p_value_threshold = 0.001, # More stringent
    min_cpg = 5               # More CpG sites
)

# Less stringent (more DMRs, lower confidence)
dmrs_relaxed <- calculate_dmrs(
    group1_data = young_prepared,
    group2_data = old_prepared,
    min_coverage = 3,         # Lower coverage
    p_value_threshold = 0.05, # Less stringent
    min_cpg = 2               # Fewer CpG sites
)
```

---

## Best Practices

### 1. Quality Control

- **Check BAM file quality**: Use `samtools stats` to verify read counts
- **Verify methylation tags**: Ensure MM/ML tags are present
- **Check coverage**: Ensure sufficient coverage per sample (>5x recommended)

### 2. Data Organization

- Use consistent naming conventions
- Keep raw data separate from processed data
- Document sample metadata (age, condition, etc.)

### 3. Reproducibility

- Record all parameters used
- Save R session info: `sessionInfo()`
- Version control your scripts
- Document any manual interventions

### 4. Performance Optimization

- Use GPU for basecalling when available
- Process samples in parallel when possible
- Filter low-coverage sites early in pipeline

### 5. Validation

- Compare results with known DMRs (if available)
- Visualize DMRs in genome browser
- Check for batch effects between samples

---

## Example: Complete Analysis Script

Save this as `run_complete_analysis.sh`:

```bash
#!/bin/bash
set -euo pipefail

# Configuration
POD5_DIR="pod5_data"
REFERENCE="reference/genome.fasta"
OUTPUT_DIR="output"
DORADO_MODEL="dna_r10.4.1_e8.2_400bps_fast@v5.0.0"

# Create output directories
mkdir -p ${OUTPUT_DIR}/{bam,bedmethyl,dmrs}

# Step 1: Basecalling
echo "Step 1: Basecalling..."
for sample in sample1 sample2 sample3; do
    echo "  Processing ${sample}..."
    dorado basecaller \
        ${DORADO_MODEL} \
        --modified-bases 5mC_5hmC \
        --device cuda:0 \
        ${POD5_DIR}/${sample} \
    | samtools view -bS - \
    > ${OUTPUT_DIR}/bam/${sample}_methylation.bam
done

# Step 2: Methylation extraction
echo "Step 2: Extracting methylation calls..."
for sample in sample1 sample2 sample3; do
    echo "  Processing ${sample}..."
    modkit extract \
        ${OUTPUT_DIR}/bam/${sample}_methylation.bam \
        ${OUTPUT_DIR}/bedmethyl/${sample}.bedMethyl \
        --ref ${REFERENCE} \
        --bedgraph \
        --combine-strands \
        --filter-threshold 0.0
done

# Step 3: DMR analysis (run R script)
echo "Step 3: Running DMR analysis..."
Rscript downstream_analysis.R

echo "Analysis complete!"
```

---

## Additional Resources

- See `METHYLATION_PIPELINE_README.md` for overview
- Check tool-specific documentation for advanced options
- Review example outputs in `analysis/` directory

---

**Last Updated:** January 8, 2025
