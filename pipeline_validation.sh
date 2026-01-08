#!/bin/bash
# =============================================================================
# ROGEN Project - Activity 2.1.8.1
# Methylation Calling Pipeline Validation Script
# =============================================================================
# This script validates the methylation calling pipeline for Oxford Nanopore
# data using the following tools:
#   1. Dorado (basecalling with methylation models)
#   2. Modkit (BAM to bedMethyl conversion)
#   3. DMRcaller (downstream analysis - see downstream_analysis.R)
#
# Test Dataset: Official ONT small test dataset from Epi2Me Labs
# =============================================================================

set -euo pipefail  # Exit on error, undefined vars, pipe failures

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Configuration
TEST_DATA_URL="https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/wf-basecalling/wf-basecalling-demo.tar.gz"
TEST_DATA_ARCHIVE="wf-basecalling-demo.tar.gz"
TEST_DATA_DIR="wf-basecalling-demo"
DORADO_MODEL="dna_r10.4.1_e8.2_400bps_fast@v5.0.0"
OUTPUT_BAM="basecalled_methylation.bam"
OUTPUT_BEDMETHYL="methylation_calls.bedMethyl"

# =============================================================================
# Function: Check if a command exists
# =============================================================================
check_command() {
    local cmd=$1
    if command -v "$cmd" &> /dev/null; then
        echo -e "${GREEN}✓${NC} $cmd is installed"
        return 0
    else
        echo -e "${RED}✗${NC} $cmd is NOT installed"
        return 1
    fi
}

# =============================================================================
# Function: Print section header
# =============================================================================
print_section() {
    echo ""
    echo "=============================================================================="
    echo "$1"
    echo "=============================================================================="
    echo ""
}

# =============================================================================
# STEP 1: Check Prerequisites
# =============================================================================
print_section "STEP 1: Checking Prerequisites"

MISSING_TOOLS=0

if ! check_command "dorado"; then
    echo "  Please install Dorado from: https://github.com/nanoporetech/dorado"
    MISSING_TOOLS=1
fi

if ! check_command "modkit"; then
    echo "  Please install Modkit from: https://github.com/nanoporetech/modkit"
    MISSING_TOOLS=1
fi

if [ $MISSING_TOOLS -eq 1 ]; then
    echo -e "${RED}ERROR: Required tools are missing. Please install them before proceeding.${NC}"
    exit 1
fi

# Display versions for documentation
echo ""
echo "Tool versions:"
dorado --version 2>&1 | head -1 || echo "  Dorado version check failed"
modkit --version 2>&1 | head -1 || echo "  Modkit version check failed"

# =============================================================================
# STEP 2: Download Test Dataset
# =============================================================================
print_section "STEP 2: Downloading Test Dataset"

if [ -d "$TEST_DATA_DIR" ]; then
    echo -e "${YELLOW}Test data directory already exists. Skipping download.${NC}"
    echo "  To re-download, remove the directory: rm -rf $TEST_DATA_DIR"
else
    echo "Downloading test dataset from: $TEST_DATA_URL"
    echo "  This may take a few minutes..."
    
    if command -v wget &> /dev/null; then
        wget -q --show-progress "$TEST_DATA_URL" -O "$TEST_DATA_ARCHIVE"
    elif command -v curl &> /dev/null; then
        curl -L --progress-bar "$TEST_DATA_URL" -o "$TEST_DATA_ARCHIVE"
    else
        echo -e "${RED}ERROR: Neither wget nor curl found. Cannot download test data.${NC}"
        exit 1
    fi
    
    echo "Extracting test dataset..."
    tar -xzf "$TEST_DATA_ARCHIVE"
    
    # Clean up archive
    rm "$TEST_DATA_ARCHIVE"
    
    echo -e "${GREEN}✓${NC} Test dataset downloaded and extracted"
fi

# Find POD5 files in the test dataset
POD5_DIR=$(find "$TEST_DATA_DIR" -type d -name "*.pod5" -o -type f -name "*.pod5" | head -1 | xargs dirname 2>/dev/null || echo "")
if [ -z "$POD5_DIR" ]; then
    # Try alternative: look for pod5 files directly
    POD5_FILES=$(find "$TEST_DATA_DIR" -name "*.pod5" | head -5)
    if [ -z "$POD5_FILES" ]; then
        echo -e "${YELLOW}WARNING: Could not automatically find POD5 files.${NC}"
        echo "  Please verify the test dataset structure."
        echo "  Expected: POD5 files should be in $TEST_DATA_DIR"
        POD5_DIR="$TEST_DATA_DIR"
    else
        POD5_DIR=$(dirname "$(echo "$POD5_FILES" | head -1)")
    fi
fi

echo "Using POD5 directory: $POD5_DIR"

# =============================================================================
# STEP 3: Basecalling with Dorado (Methylation-Aware)
# =============================================================================
print_section "STEP 3: Basecalling with Dorado (Methylation-Aware)"

echo "Running Dorado basecalling with methylation model..."
echo ""
echo "Command breakdown:"
echo "  dorado basecaller                    # Basecalling subcommand"
echo "    $DORADO_MODEL                      # Model: R10.4.1 chemistry, E8.2 kit, 400bps, fast mode, v5.0.0"
echo "    --modified-bases 5mC_5hmC          # Enable detection of 5-methylcytosine (5mC) and"
echo "                                       #   5-hydroxymethylcytosine (5hmC) modifications"
echo "                                       #   These modifications are written as MM/ML tags in BAM"
echo "    --device cuda:0                    # Use GPU (CUDA device 0) for acceleration"
echo "                                       #   Use 'cpu' if GPU unavailable"
echo "    $POD5_DIR                          # Input directory containing POD5 files"
echo "  | samtools view -bS -                # Convert SAM to BAM format"
echo "  > $OUTPUT_BAM                        # Output BAM file with methylation tags"

# Check if GPU is available (optional check)
if command -v nvidia-smi &> /dev/null; then
    echo ""
    echo "GPU information:"
    nvidia-smi --query-gpu=name,memory.total --format=csv,noheader | head -1
    DEVICE="cuda:0"
else
    echo ""
    echo -e "${YELLOW}No NVIDIA GPU detected. Will use CPU (slower).${NC}"
    DEVICE="cpu"
fi

echo ""
echo "Starting basecalling (this may take a while)..."
echo "  Note: For validation purposes, this step can be skipped if GPU unavailable"
echo "  The output BAM will contain MM (modification probabilities) and ML (modification calls) tags"

# Uncomment the following lines to actually run basecalling:
# dorado basecaller \
#     "$DORADO_MODEL" \
#     --modified-bases 5mC_5hmC \
#     --device "$DEVICE" \
#     "$POD5_DIR" \
# | samtools view -bS - \
# > "$OUTPUT_BAM"

# For validation/documentation purposes, create a placeholder if BAM doesn't exist
if [ ! -f "$OUTPUT_BAM" ]; then
    echo ""
    echo -e "${YELLOW}NOTE: Basecalling step commented out for validation script.${NC}"
    echo "  To run basecalling, uncomment the dorado command in this script."
    echo "  Creating placeholder message..."
    echo "# Basecalling output would be written to: $OUTPUT_BAM" > "$OUTPUT_BAM.placeholder"
fi

if [ -f "$OUTPUT_BAM" ]; then
    echo -e "${GREEN}✓${NC} Basecalling completed: $OUTPUT_BAM"
    # Show BAM statistics
    if command -v samtools &> /dev/null; then
        echo ""
        echo "BAM file statistics:"
        samtools stats "$OUTPUT_BAM" | grep "^SN" | head -5
    fi
else
    echo -e "${YELLOW}⚠${NC}  Basecalling step skipped (see placeholder file)"
fi

# =============================================================================
# STEP 4: Convert BAM to bedMethyl with Modkit
# =============================================================================
print_section "STEP 4: Converting BAM to bedMethyl Format"

if [ ! -f "$OUTPUT_BAM" ]; then
    echo -e "${YELLOW}WARNING: BAM file not found. Skipping Modkit step.${NC}"
    echo "  Expected BAM file: $OUTPUT_BAM"
    echo "  Run basecalling first (uncomment dorado command above)"
    exit 0
fi

echo "Running Modkit to extract methylation calls..."
echo ""
echo "Command breakdown:"
echo "  modkit extract                        # Extract modification calls from BAM"
echo "    $OUTPUT_BAM                         # Input BAM file with MM/ML tags"
echo "    $OUTPUT_BEDMETHYL                   # Output bedMethyl file"
echo "    --ref <reference.fasta>             # Reference genome FASTA (required)"
echo "    --bedgraph                          # Also generate bedGraph format (optional)"
echo "    --combine-strands                   # Combine forward and reverse strand calls"
echo "    --filter-threshold 0.0              # Minimum modification probability threshold"
echo "                                       #   0.0 includes all calls (can filter later)"

# Note: Modkit requires a reference genome
# For the test dataset, we would need the reference FASTA
echo ""
echo -e "${YELLOW}NOTE: Modkit requires a reference genome FASTA file.${NC}"
echo "  The test dataset should include a reference, or you can download it separately."
echo "  Example: wget https://example.com/reference.fasta"

# Uncomment and modify the following when reference is available:
# REFERENCE_FASTA="path/to/reference.fasta"
# 
# if [ ! -f "$REFERENCE_FASTA" ]; then
#     echo -e "${RED}ERROR: Reference FASTA not found at: $REFERENCE_FASTA${NC}"
#     exit 1
# fi
# 
# modkit extract \
#     "$OUTPUT_BAM" \
#     "$OUTPUT_BEDMETHYL" \
#     --ref "$REFERENCE_FASTA" \
#     --bedgraph \
#     --combine-strands \
#     --filter-threshold 0.0

echo ""
echo "Creating placeholder for bedMethyl output..."
echo "# bedMethyl output would be written to: $OUTPUT_BEDMETHYL" > "$OUTPUT_BEDMETHYL.placeholder"
echo "# Format: chrom, start, end, name, score, strand, coverage, percent_methylated, ..."

if [ -f "$OUTPUT_BEDMETHYL" ]; then
    echo -e "${GREEN}✓${NC} bedMethyl file created: $OUTPUT_BEDMETHYL"
    echo ""
    echo "First few lines of bedMethyl file:"
    head -5 "$OUTPUT_BEDMETHYL"
else
    echo -e "${YELLOW}⚠${NC}  Modkit step skipped (see placeholder file)"
fi

# =============================================================================
# STEP 5: Validation Summary
# =============================================================================
print_section "STEP 5: Validation Summary"

echo "Pipeline validation artifacts:"
echo ""
echo "  ✓ Prerequisites checked (Dorado, Modkit)"
echo "  ✓ Test dataset downloaded: $TEST_DATA_DIR"
if [ -f "$OUTPUT_BAM" ]; then
    echo "  ✓ Basecalling completed: $OUTPUT_BAM"
else
    echo "  ⚠ Basecalling step (commented out for validation)"
fi
if [ -f "$OUTPUT_BEDMETHYL" ]; then
    echo "  ✓ bedMethyl conversion completed: $OUTPUT_BEDMETHYL"
else
    echo "  ⚠ bedMethyl conversion (commented out for validation)"
fi

echo ""
echo "Next steps:"
echo "  1. Uncomment basecalling command (line ~100) to run Dorado"
echo "  2. Provide reference FASTA and uncomment Modkit command (line ~140)"
echo "  3. Run downstream_analysis.R for DMR calling"
echo ""
echo -e "${GREEN}Validation script completed successfully!${NC}"
