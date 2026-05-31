# Synthetic UKB Integrative Validation Pipeline

**Project:** IBAR-ROGEN Aging  
**Activity:** 2.1.11.1  
**Script:** `scripts/run_integration.py`  
**Library:** `src/rogen_aging/integration/ukb_joiner.py`  
**Input:** Mock UKB-RAP output from `scripts/ukb_mock_gen.py`

## Overview

Joins a **synthetic** UK Biobank RAP-style phenotype CSV and LA-SNP VCF on `eid`, then runs dominant-model association scans across all manifest SNPs for two binary outcomes:

1. **Parental longevity** (`parental_longevity`)
2. **AD diagnosis** (derived from non-empty `ad_diagnosis_code`)

This validates the integrative architecture end to end on **mock data only**. Outputs are for pipeline QA — **not** biological conclusions.

## Prerequisites

Generate mock RAP artefacts first:

```bash
uv run python scripts/ukb_la_snp_lookup.py build \
  --input overlapping_genes_with_snps.xlsx \
  --output analysis/ukb_snp_manifest_v0.1.csv

uv run scripts/ukb_mock_gen.py \
  --n-samples 1000 \
  --snp-manifest analysis/ukb_snp_manifest_v0.1.csv \
  --output-dir test_data/mock_ukb_rap/ \
  --seed 42
```

See [Synthetic UKB-RAP Generator](SYNTHETIC_UKB_RAP_GENERATOR.md) and [LA-SNP Public Frequency Pipeline](LA_SNP_PUBLIC_FREQUENCY_PIPELINE.md) for upstream steps.

## Usage

```bash
uv run python scripts/run_integration.py \
  --pheno test_data/mock_ukb_rap/phenotypes/ukb_phenotypes.csv \
  --vcf test_data/mock_ukb_rap/genotypes/ukb_la_snps.vcf \
  --output-dir analysis/
```

Defaults point at `test_data/mock_ukb_rap/` paths. Use `-v` for debug logging.

## Outputs

Written under `--output-dir` (default `analysis/`):

| File | Description |
|------|-------------|
| `assoc_la_snp_parental_longevity.csv` | Per-SNP dominant OR, 95% CI, Fisher p, n |
| `assoc_la_snp_ad.csv` | Same for AD diagnosis flag |

Each CSV includes a synthetic-data disclaimer header. Columns: `rsID`, `OR`, `CI_low`, `CI_high`, `p_value`, `n`.

## Implementation notes

- Genotypes loaded with **cyvcf2**; dosages 0/1/2 from `gt_type`.
- Inner join on `eid`; row count must match phenotype table (sample IDs in VCF columns must equal `eid` values).
- Association: 2×3 contingency (outcome × dosage) collapsed to dominant 2×2; **Fisher exact** two-sided p; Woolf log-OR CI with Haldane correction when needed.

## Python API

```python
from pathlib import Path
from rogen_aging.integration.ukb_joiner import run_integration_pipeline

joined, parental, ad = run_integration_pipeline(
    Path("test_data/mock_ukb_rap/phenotypes/ukb_phenotypes.csv"),
    Path("test_data/mock_ukb_rap/genotypes/ukb_la_snps.vcf"),
    Path("analysis"),
)
```

## Tests

```bash
uv run pytest tests/test_ukb_integration.py -q
```

Uses a 70-SNP mock manifest fixture and `ukb_mock_gen.generate_ukb_rap_mock`.

## Security and compliance

- **Synthetic data only** — no real UKB participant IDs or genotypes.
- Input VCF paths are read locally; do not commit `.vcf` files (git-ignored).
- `scripts/run_integration.py` and `src/rogen_aging/integration/` are whitelisted in the pre-commit hook because they reference mock `ukb_*` paths by design.

## Related documentation

- [Synthetic UKB-RAP Generator](SYNTHETIC_UKB_RAP_GENERATOR.md)
- [Code Modules Reference](CODE_MODULES_REFERENCE.md) — §2.4 `rogen_aging.integration`, §3.21 `run_integration.py`
- [UKB Pre-Commit Hook](UKB_PRE_COMMIT_HOOK.md)

---

**Last updated:** May 31, 2026
