# Mitochondrial (chrM) extract, Mutect2, HaploGrep 3

Validation pipeline for **public** GRCh38 alignments with a published
haplogroup. The script extracts mitochondrial reads, calls variants with GATK
Mutect2 `--mitochondria-mode`, writes depth/VAF tables, assigns a haplogroup
with HaploGrep 3, and reports major-clade concordance as `n/N`.

This is not a population study. Do not use Romanian sequencing data. Do not
claim population heteroplasmy rates or aging interpretations from these
outputs.

## Inputs

| Path | Role |
|------|------|
| `data/external/mtdna_validation_samples.csv` | You create this. Columns: `sample_id`, `path_or_url`, `published_haplogroup`, `source_citation`. Public samples only. |
| GRCh38 FASTA (`--reference`) | Same reference used for the alignments. Must be `samtools faidx`-indexed and contain contig `chrM`. |
| `config/mtdna_calling.yaml` | Fixed `min_depth` / `min_vaf`, contig name, HaploGrep tree, optional NUMT mask path. |

### Samples CSV

```text
sample_id,path_or_url,published_haplogroup,source_citation
HG002,/path/to/HG002.bam,H1a,citation or DOI
```

`path_or_url` is a local BAM/CRAM (preferred) or an `http(s)` URL downloaded
into the output `_downloads/` cache. Alignments must be GRCh38 with contig
`chrM` (asserted against the reference `.fai`).

### Thresholds (do not tune for concordance)

Defaults in [`config/mtdna_calling.yaml`](../../config/mtdna_calling.yaml):

| Key | Default | Meaning |
|-----|---------|---------|
| `thresholds.min_depth` | `10` | Minimum DP for a site/allele to enter the filtered TSV and HaploGrep input |
| `thresholds.min_vaf` | `0.05` | Minimum alt allele fraction (`AD_alt / DP`) for the same |

Raw Mutect2 VCFs keep all called sites. Thresholds only gate the filtered
variant TSV / HaploGrep input. **Do not change these defaults to improve
haplogroup concordance.**

### NUMT mask [UNVERIFIED]

Config key `numt_mask` defaults to `data/external/mtdna/numt_mask.bed`. That
path is **unverified** and may be absent. Behaviour:

- **File present:** variants overlapping the BED (`chrM` intervals) get
  `numt_overlap=true` in the TSV.
- **File absent:** every output header records `NUMT filtering not applied`.
  No substitute mask is invented.

## Tools

Required on `PATH` (resolved with `which`; missing tools stop the run):

- `samtools`
- `gatk` (Mutect2)
- `haplogrep3` (with tree `phylotree-fu-rcrs@1.2` installed)

## Run

```bash
uv run python scripts/mtdna/extract_and_call_mtdna.py \
  --samples data/external/mtdna_validation_samples.csv \
  --reference /path/to/GRCh38.fa \
  --config config/mtdna_calling.yaml \
  --out-dir analysis/mtdna
```

### Per sample

1. Assert BAM/CRAM header contig `chrM` matches the reference `.fai`.
2. Extract `chrM` reads → `{sample}.chrM.bam`; report read count and mean depth.
3. **Zero chrM reads:** status `skipped_zero_chrm_reads`, recorded in the
   summary, **no** filled VCF/haplogroup.
4. Mutect2 `--mitochondria-mode` → `{sample}.mutect2.vcf`.
5. TSV with depth and VAF → `{sample}.variants.tsv`.
6. HaploGrep 3 on threshold-passing alleles → `{sample}.haplogrep.tsv`.

### Concordance

Major clade = alphabetic prefix before the first digit (`H1a` → `H`,
`HV1` → `HV`). Among samples with status `ok`:

```text
major_clade_concordance = n/N
```

`N` excludes `skipped_zero_chrm_reads`. The string is printed and written in
`analysis/mtdna/haplogroup_concordance.tsv` and `run_provenance.json`.

## Outputs

| Path | Description |
|------|-------------|
| `analysis/mtdna/{sample_id}/` | chrM BAM, Mutect2 VCF, variants TSV, HaploGrep files, `metrics.json` |
| `analysis/mtdna/haplogroup_concordance.tsv` | Per-sample status + major-clade comparison; header has `n/N` |
| `analysis/mtdna/run_provenance.json` | Tool versions, input sha256, git commit, timestamp, NUMT status |

Every TSV/VCF header (or leading `#` comment) carries the NUMT note and a
provenance line (git commit, timestamp, input hashes, tool versions).

## Prohibited

- Romanian WGS/ONT or any non-public ROGEN delivery as input
- Population heteroplasmy claims
- Aging interpretation of mtDNA calls or concordance
