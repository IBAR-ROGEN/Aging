# Methylation Pipeline - Quick Reference

Quick reference guide for common tasks and commands.

## Installation Checklist

- [ ] Install Dorado
- [ ] Install Modkit
- [ ] Install R and Bioconductor packages
- [ ] Download methylation model: `dorado download --model dna_r10.4.1_e8.2_400bps_fast@v5.0.0`
- [ ] Download reference genome FASTA

## Common Commands

### Basecalling

```bash
# Single sample
dorado basecaller dna_r10.4.1_e8.2_400bps_fast@v5.0.0 \
    --modified-bases 5mC_5hmC --device cuda:0 pod5_dir \
| samtools view -bS - > output.bam

# CPU mode
dorado basecaller dna_r10.4.1_e8.2_400bps_fast@v5.0.0 \
    --modified-bases 5mC_5hmC --device cpu pod5_dir \
| samtools view -bS - > output.bam
```

### Methylation Extraction

```bash
modkit extract input.bam output.bedMethyl \
    --ref reference.fasta \
    --bedgraph \
    --combine-strands \
    --filter-threshold 0.0
```

### R Analysis

```r
# Load libraries
library(DMRcaller)
library(GenomicRanges)
library(rtracklayer)

# Import data
sample <- import_bedmethyl("sample.bedMethyl")

# Prepare data
prepared <- prepare_dmrcaller_data(sample, "Sample1")

# Calculate DMRs
dmrs <- calculate_dmrs(group1, group2, 
                       min_coverage=5, 
                       p_value_threshold=0.01, 
                       min_cpg=3)

# Export results
export(dmrs, "dmrs.bed", format="bed")
```

## File Formats

### bedMethyl (BED9+3)
```
chrom  start  end  name  score  strand  coverage  percent_methylated  ...
```

### DMR Summary CSV
```
chr, start, end, width, n_cpg, p_value, mean_meth_diff
```

## Parameter Reference

### Dorado
- `--modified-bases 5mC_5hmC`: Detect 5mC and 5hmC
- `--device cuda:0`: Use GPU (or `cpu`)
- Model: `dna_r10.4.1_e8.2_400bps_fast@v5.0.0`

### Modkit
- `--ref`: Reference genome (required)
- `--bedgraph`: Generate bedGraph output
- `--combine-strands`: Combine forward/reverse strands
- `--filter-threshold 0.0`: Minimum methylation probability

### DMRcaller
- `min_coverage`: Minimum reads per position (default: 5)
- `p_value_threshold`: Significance threshold (default: 0.01)
- `min_cpg`: Minimum CpG sites per DMR (default: 3)

## Troubleshooting Quick Fixes

| Issue | Solution |
|-------|----------|
| Dorado not found | Add to PATH or use full path |
| GPU errors | Use `--device cpu` |
| Modkit needs reference | Download reference FASTA |
| R packages fail | Update BiocManager, install individually |
| Out of memory | Process in batches, filter by chromosome |

## Output Files

- `*.bam` - BAM files with MM/ML tags
- `*.bedMethyl` - Methylation calls in bedMethyl format
- `dmrs_*.bed` - DMR coordinates
- `dmr_summary.csv` - DMR statistics

## Directory Structure

```
project/
├── pod5_data/          # Raw POD5 files
├── reference/          # Reference genome
├── output/
│   ├── bam/           # Basecalled BAM files
│   ├── bedmethyl/     # bedMethyl files
│   └── dmrs/          # DMR results
└── scripts/            # Pipeline scripts
```

## Useful Commands

```bash
# Check BAM has methylation tags
samtools view file.bam | head -1 | grep -o "MM\|ML"

# Count reads in BAM
samtools view -c file.bam

# Check bedMethyl format
head -5 file.bedMethyl

# Count methylation calls
wc -l file.bedMethyl

# Check R session info
Rscript -e "sessionInfo()"
```

## Links

- Full Documentation: [METHYLATION_PIPELINE_README.md](../METHYLATION_PIPELINE_README.md)
- Detailed Usage: [METHYLATION_PIPELINE_USAGE.md](METHYLATION_PIPELINE_USAGE.md)
- Dorado: https://github.com/nanoporetech/dorado
- Modkit: https://github.com/nanoporetech/modkit
- DMRcaller: https://bioconductor.org/packages/DMRcaller/

---

**Last Updated:** January 8, 2025
