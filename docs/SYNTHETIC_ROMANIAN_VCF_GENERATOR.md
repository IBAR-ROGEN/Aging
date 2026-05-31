# Synthetic Romanian Cohort VCF Generator

**Project:** [IBAR-ROGEN/Aging](https://github.com/IBAR-ROGEN/Aging)  
**CLI:** `uv run rogen-vcf-synthetic`  
**Script:** `scripts/vcf/generate_synthetic_romanian_vcf.py`  
**Package:** `rogen_aging.vcf.synthetic`  
**Output:** User-specified path (uncompressed VCF v4.2 text)

## Overview

This tool writes a **synthetic, bcftools-friendly VCF v4.2** representing a **mock Romanian population cohort** with **European (EUR)-style** allele-frequency structure. It is intended for **pipeline testing, teaching, and tooling validation** when real cohort VCFs are unavailable or must not be committed.

Important properties:

- **Streaming output** — Variant lines are generated and written one at a time (no full VCF held in memory).
- **Hardy–Weinberg genotypes** — For each site, an alternate-allele frequency `q` is drawn uniformly in **[0.01, 0.5]**; diploid genotypes are sampled from `(p², 2pq, q²)` with `p = 1 - q`.
- **Standard FORMAT** — Per-sample fields are **`GT:AD:DP:GQ`** (unphased `GT`, integer `AD`/`DP`/`GQ`).
- **Sorted records** — Variants are emitted in **genome order** (chr1 → chr22, increasing `POS` within each chromosome) so **`bcftools index`** works without a separate sort step.
- **Reference metadata** — `##contig` lines use **GRCh38** lengths for autosomes **chr1–chr22**.

All sample IDs and variant IDs are **synthetic** (for example `RO_EUR_000001`, `RO_MOCK_00000001`). **No real participant data** is used.

## Requirements

- Python **3.12+** (see `pyproject.toml`).
- Dependencies installed via **`uv sync`** (uses **NumPy** for random sampling).
- Optional: **[HTSlib](https://www.htslib.org/)** / **`bcftools`** to validate or compress the output.

## Usage

### Minimal example

```bash
uv run rogen-vcf-synthetic \
  --samples 100 \
  --variants 1000 \
  --output data/mock_romanian_eur.vcf
```

### Reproducibility

```bash
uv run rogen-vcf-synthetic \
  --samples 50 --variants 500 --seed 42 \
  --output data/mock_romanian_eur.vcf
```

### Logging

Progress logs go to **stderr** at INFO by default. Use **`-v` / `--verbose`** for DEBUG.

```bash
uv run rogen-vcf-synthetic \
  --samples 10 --variants 100 --output data/mock.vcf -v
```

### Validate with bcftools (optional)

```bash
bcftools view -H data/mock_romanian_eur.vcf | head
bcftools index data/mock_romanian_eur.vcf
```

### All options

```bash
uv run rogen-vcf-synthetic --help
```

| Option | Default | Description |
|--------|---------|-------------|
| `--samples` | (required) | Number of diploid samples (columns after `FORMAT`). |
| `--variants` | (required) | Number of variant data rows to write. |
| `--output` | (required) | Output file path (uncompressed VCF). |
| `--seed` | unset | Random seed for reproducible output. |
| `--mean-depth` | `32.0` | Mean total read depth per sample per site (Poisson mean for `DP`; `AD` is simulated consistently with `GT`). |
| `--cohort-label` | `mock_RO_EUR_cohort` | Value for the `##synthetic_cohort` meta line. |
| `--sample-prefix` | `RO_EUR` | Prefix for sample column names (`RO_EUR_000001`, …). |
| `-v`, `--verbose` | off | Enable DEBUG logging. |

## VCF content summary

### Headers

- `##fileformat=VCFv4.2`
- `##FILTER=<ID=PASS,...>`
- `##INFO` — `AC`, `AN`, `AF`, `END` (biallelic SNPs; `END` equals `POS` for each record)
- `##FORMAT` — `GT`, `AD`, `DP`, `GQ`
- `##contig=<ID=chrN,length=...>` for chr1–chr22 (GRCh38 lengths)
- `#CHROM` … `FORMAT` plus one column per sample

### Variant rows

- **SNPs only** — Random `REF`/`ALT` single bases (`A`, `C`, `G`, `T`), always biallelic.
- **ID** — `RO_MOCK_` + zero-padded index.
- **QUAL** — Fixed `60`; **FILTER** — `PASS`.
- **INFO** — `AC`, `AN`, and `AF` derived from the **simulated genotypes** for that row (`AF = AC/AN`).

## Automated tests

**`tests/test_synthetic_vcf.py`** exercises `rogen_aging.vcf.synthetic.main()` with a tiny cohort written under pytest’s **`tmp_path`**.

```bash
uv run pytest tests/test_synthetic_vcf.py
```

## Security and compliance

- **Synthetic-only** — Safe for public repositories **as code**; generated `.vcf` files can be large and often should stay under **`data/`** (git-ignored).
- **Pre-commit hook** — This repository’s hook blocks committing raw **`.vcf`** / **`.bed`** in many cases.
- **Naming** — Sample names use the configurable prefix (default `RO_EUR`), not real study IDs.

## Related documentation

- [WORKFLOWS.md](WORKFLOWS.md)
- [Synthetic UKB-RAP Folder Generator](SYNTHETIC_UKB_RAP_GENERATOR.md) — reuses VCF helpers for LA-SNP mock VCF
- [UK Biobank Pre-Commit Hook](UKB_PRE_COMMIT_HOOK.md)

---

**Last updated:** May 31, 2026
