# Git Pre-Commit Hook — UK Biobank Security

**Project:** IBAR-ROGEN Aging  
**Related:** [UKB Compliance Auditor](UKB_COMPLIANCE_AUDITOR.md), [Methylation Pipeline](METHYLATION_PIPELINE_USAGE.md)

## Overview

The UK Biobank pre-commit hook automatically blocks Git commits that could expose restricted data. It runs before every commit and enforces compliance with the UK Biobank Data Usage Agreement (DUA) and ROGEN data safety policies.

## Installation

Install the hook from the project root:

```bash
./scripts/install_pre_commit_hook.sh
```

This copies `scripts/security_check.sh` to `.git/hooks/pre-commit` and makes it executable. The hook runs automatically on every `git commit`.

To reinstall (e.g. after updating the security script):

```bash
./scripts/install_pre_commit_hook.sh
```

## What It Blocks

| Category | Patterns / Extensions | Action |
|----------|----------------------|--------|
| **Content patterns** | `patient_id`, `UKB_` (case-insensitive) | Blocks commit |
| **File extensions** | `.vcf`, `.vcf.gz`, `.vcf.bgz`, `.bed`, `.bed.gz` | Blocks commit |

These formats may contain participant identifiers or individual-level genomic data that must not be committed to public repositories.

## Exemptions (Whitelisted)

The following are excluded from content scanning:

| Path | Reason |
|------|--------|
| `docs/UKB_COMPLIANCE*` | Compliance documentation |
| `README*` | Project documentation |
| `*security_check*` | The security script itself |
| `*install_pre_commit*` | Hook installer |
| `scripts/mock_ukb_generator.py` | Synthetic data generator (no real UKB data) |
| `test_data/mock_clinical_data.csv` | Synthetic mock output (MOCK_ IDs) |

## Manual Run

Run the check without committing:

```bash
./scripts/security_check.sh
```

The script inspects **staged** files (`git diff --cached`). Stage files first if you want to test:

```bash
git add <files>
./scripts/security_check.sh
```

## Blocked Commit Output

When a violation is detected:

```
Running UK Biobank security check on staged files...

BLOCKED: Staged file contains sensitive pattern: path/to/file
  -> Found 'patient_id', 'UKB_', or similar. Remove before committing.

========================================
  UK BIOBANK DATA SECURITY WARNING
========================================

Your commit was BLOCKED to prevent accidental exposure of restricted data.
...
```

## Resolving Violations

1. **Content patterns**  
   Remove or replace `patient_id`, `UKB_`, or similar strings. For synthetic/mock data, use `MOCK_` IDs and ensure the mock generator and its output are whitelisted.

2. **Restricted file extensions**  
   Do not commit `.vcf` or `.bed` files. Place them in the git-ignored `data/` directory and reference them by path in code or config.

3. **Legitimate mock/test data**  
   Use `scripts/mock_ukb_generator.py` to produce synthetic data with `MOCK_` Sample_IDs. The generator and `test_data/mock_clinical_data.csv` are whitelisted.

## Related Documentation

- [UKB Compliance Auditor](UKB_COMPLIANCE_AUDITOR.md) — Full compliance tool and audit reporting
- [Synthetic UK Biobank Data Generator](SYNTHETIC_UKB_GENERATOR.md) — Generating safe mock data for pipeline development

---

**Last Updated:** February 27, 2026
