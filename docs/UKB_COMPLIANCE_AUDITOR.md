# UKB Compliance Auditor - Documentation

**Project:** ROGEN Aging Research  
**Activity:** 2.1.8.1 - UKB Compliance Auditing  
**Tool:** `notebooks/UKB_Compliance_Auditor.ipynb`

## Overview

The **UKB Compliance Auditor** is a security and compliance tool designed to prevent the accidental disclosure of restricted UK Biobank (UKB) data. It scans the `rogen_aging` repository for sensitive identifiers—specifically Participant IDs (EIDs)—and restricted file formats before any code or data is pushed to public repositories or shared portals.

This tool ensures adherence to the UK Biobank Data Usage Agreement (DUA) and the ROGEN project's data safety protocols.

## Key Features

- **Pattern-Based Scanning**: Detects 7-digit integers (the standard format for UKB EIDs) and sensitive keywords (`participant_id`, `eid`, `patient_name`).
- **File Type Enforcement**: Identifies restricted genomic and data formats (`.bam`, `.pod5`, `.sqlite`, `.xlsx`, etc.).
- **Automated Exclusion**: Skips non-relevant files like `uv.lock`, `.git/`, and dependency directories to reduce false positives.
- **Audit Trail**: Records the researcher's GitHub username and provides a structured report of all findings.
- **UKB API Integration**: Includes a wrapper for submitting audit reports to the UK Biobank Git Audit Tool API (mock implementation).

## Usage Instructions

### 1. Prerequisites

Ensure you have the project environment set up:

```bash
uv sync
```

### 2. Launching the Auditor

The auditor is provided as an interactive Jupyter notebook for easy review of findings.

```bash
uv run jupyter lab notebooks/UKB_Compliance_Auditor.ipynb
```

### 3. Running the Audit

1. Open the notebook.
2. Update the `GITHUB_USERNAME` variable in the configuration cell with your GitHub handle.
3. Run all cells (`Cell > Run All`).
4. Review the "Potential Compliance Issues Found" section.

## Configuration

You can customize the scan behavior in the first code cell of the notebook:

| Parameter | Description |
|-----------|-------------|
| `RESTRICTED_PATTERNS` | Regex patterns to search for in file contents. |
| `RESTRICTED_EXTENSIONS` | File extensions that are strictly prohibited. |
| `IGNORE_DIRS` | Directories to skip (e.g., `.venv`, `.git`). |
| `IGNORE_FILES` | Specific files to skip (e.g., `uv.lock`). |

## Compliance Policy

According to ROGEN Activity 2.1.8.1:
1. **Zero EID Policy**: No file containing 7-digit participant IDs should be committed to version control.
2. **Data Locality**: Raw Nanopore data (`.pod5`, `.bam`) must remain in the `data/` directory, which is git-ignored.
3. **Audit Requirement**: A compliance scan must be performed and reported before any major release or public push.

## Troubleshooting

### High False Positives
If the scan flags many 7-digit numbers that are not EIDs (e.g., library versions in `uv.lock` or hashes), add those files to the `IGNORE_FILES` list in the configuration cell.

### Binary File Scanning
The tool defaults to scanning files as UTF-8 text. Binary files (like `.png` or `.bam`) are checked for extension only to avoid performance issues and decoding errors.

---
**Last Updated:** January 30, 2026  
**Status:** Active  
**Maintained by:** ROGEN Compliance Team
