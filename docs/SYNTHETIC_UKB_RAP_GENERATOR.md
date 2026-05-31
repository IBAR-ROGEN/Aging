# Synthetic UKB-RAP Folder Generator

**Project:** IBAR-ROGEN Aging  
**Activity:** 2.1.8.1  
**Script:** `scripts/ukb_mock_gen.py`  
**Output:** `test_data/mock_ukb_rap/` (default)

## Overview

This tool writes a **UK Biobank RAP-style directory** with two joinable artefacts:

1. **Phenotype table** — CSV keyed by synthetic `eid`
2. **Genotype VCF** — LA-SNP VCF restricted to the ~70 SNPs in the offline manifest, with sample columns matching `eid`

All participant identifiers and genotypes are **strictly synthetic**. No real UKB EIDs or participant data are used. Safe for GitHub when only the **code** is committed; generated `.vcf` and phenotype CSV outputs are git-ignored by default (see [Security and compliance](#security-and-compliance)).

Genotype simulation reuses helpers from `scripts/generate_synthetic_romanian_vcf.py` (Hardy–Weinberg sampling, EUR-like allele frequencies, `GT:AD:DP:GQ` FORMAT fields).

## Output layout

```
test_data/mock_ukb_rap/
├── phenotypes/
│   └── ukb_phenotypes.csv
└── genotypes/
    └── ukb_la_snps.vcf
```

Both files include an **Activity 2.1.8.1** safety header noting the cohort is synthetic.

## Phenotype columns (v2 dictionary)

| Column | Description |
|--------|-------------|
| `eid` | Synthetic participant ID (`SYN_EID_00000001`, … — not real 7-digit UKB EIDs) |
| `age` | Age in years (uniform 40–80 by default) |
| `sex` | Binary sex code (`0` / `1`) |
| `parental_longevity` | Parental longevity flag (`0` / `1`) |
| `ad_diagnosis_code` | Synthetic ICD-10-style AD code, or empty string |
| `pd_diagnosis_code` | Synthetic ICD-10-style PD code, or empty string |
| `frailty_weight_loss` | Frailty component (`0` / `1`) |
| `frailty_exhaustion` | Frailty component (`0` / `1`) |
| `frailty_weakness` | Frailty component (`0` / `1`) |
| `frailty_slowness` | Frailty component (`0` / `1`) |
| `frailty_low_activity` | Frailty component (`0` / `1`) |

## Prerequisites

Build the LA-SNP manifest once (offline Ensembl lookup; no participant data):

```bash
uv run python scripts/ukb_la_snp_lookup.py build \
  --input overlapping_genes_with_snps.xlsx \
  --output analysis/ukb_snp_manifest_v0.1.csv
```

Optionally validate public allele frequencies (1KG extract + gnomAD comparison) before generating synthetic data — see **[LA_SNP_PUBLIC_FREQUENCY_PIPELINE.md](LA_SNP_PUBLIC_FREQUENCY_PIPELINE.md)**.

See `notebooks/05_ukb_exploration/UKB_LA_SNP_FirstContact.ipynb` for manifest sanity checks.

## Usage

### Basic (1000 samples, default output)

```bash
uv run scripts/ukb_mock_gen.py
```

Writes under `test_data/mock_ukb_rap/` using `analysis/ukb_snp_manifest_v0.1.csv`.

### Custom cohort size and seed

```bash
uv run scripts/ukb_mock_gen.py \
  --n-samples 500 \
  --snp-manifest analysis/ukb_snp_manifest_v0.1.csv \
  --output-dir data/mock_ukb_rap/ \
  --seed 42
```

### All options

```bash
uv run scripts/ukb_mock_gen.py --help
```

| Option | Short | Default | Description |
|--------|-------|---------|-------------|
| `--n-samples` | `-n` | 1000 | Number of synthetic participants |
| `--snp-manifest` | | `analysis/ukb_snp_manifest_v0.1.csv` | LA-SNP manifest CSV (rsID + GRCh38 coords) |
| `--output-dir` | `-o` | `test_data/mock_ukb_rap/` | Root directory for RAP-style layout |
| `--seed` | `-s` | unset | Random seed for reproducibility |
| `--mean-depth` | | 32.0 | Mean simulated read depth per site |
| `--verbose` | `-v` | off | Enable debug logging |

## VCF content

- **VCF v4.2** with GRCh38 `##contig` lines (chr1–chr22)
- One row per manifest SNP (`ID` = rsID), genome-sorted
- **FORMAT:** `GT:AD:DP:GQ` per sample; sample column names = `eid` values
- Hardy–Weinberg genotypes with EUR-like alternate-allele frequencies in [0.01, 0.5]

## Automated tests

**`tests/test_ukb_mock_gen.py`** builds a 70-row manifest fixture under pytest’s `tmp_path`, runs `generate_ukb_rap_mock()`, and checks:

- Phenotype and VCF share the same `eid` set
- VCF contains exactly **70** variant rows
- All v2 phenotype dictionary columns are present

```bash
uv run pytest tests/test_ukb_mock_gen.py
```

## Security and compliance

- **Synthetic EIDs:** `SYN_EID_*` prefix — not real UKB 7-digit participant IDs.
- **Whitelisted script:** `scripts/ukb_mock_gen.py` is excluded from the UK Biobank pre-commit content scan (see [UKB_PRE_COMMIT_HOOK.md](UKB_PRE_COMMIT_HOOK.md)).
- **Generated outputs:** `.vcf` files are blocked by the pre-commit hook; phenotype CSVs under `test_data/` are git-ignored except for explicitly whitelisted fixtures. Keep generated RAP folders under `data/` or local paths.

## Related documentation

- [LA-SNP Public Frequency Pipeline](LA_SNP_PUBLIC_FREQUENCY_PIPELINE.md) — manifest, 1KG extract, gnomAD comparison
- [Synthetic UK Biobank Data Generator](SYNTHETIC_UKB_GENERATOR.md) — tabular mock clinical CSV (`mock_ukb_generator.py`)
- [Synthetic Romanian Cohort VCF Generator](SYNTHETIC_ROMANIAN_VCF_GENERATOR.md) — general-purpose streaming VCF generator
- [UK Biobank Pre-Commit Hook](UKB_PRE_COMMIT_HOOK.md) — security checks and whitelisting
- [Code Modules Reference](CODE_MODULES_REFERENCE.md) — `ukb_la_snp_lookup.py` manifest builder

---

**Last updated:** May 31, 2026
