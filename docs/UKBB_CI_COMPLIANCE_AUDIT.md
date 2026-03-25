# UK Biobank–oriented CI repository audit

**Script:** `scripts/ukbb_ci_compliance_audit.sh`  
**Purpose:** Fail CI when the repository tree contains high-risk genomic/clinical files or disallowed hardcoded absolute paths in Python and shell code.

This check complements the [pre-commit hook](UKB_PRE_COMMIT_HOOK.md) (`scripts/security_check.sh`), which inspects **staged** files, and the [notebook compliance auditor](UKB_COMPLIANCE_AUDITOR.md), which focuses on patterns such as EIDs. The CI script scans the **checked-out tree** (what the pipeline actually builds from), which is the right place to catch files that are already committed.

## What it is not

- It does **not** replace legal review, your institutional policy, or [UK Biobank](https://www.ukbiobank.ac.uk/)’s own guidance and agreements.
- It does **not** prove absence of sensitive **content** inside allowed file types (for example, a `.json` or small `.csv` could still hold restricted fields). Use content-oriented checks (pre-commit, notebook auditor, manual review) alongside this tool.

## Requirements

- Bash 4+ (typical on GitHub Actions `ubuntu-latest`).
- `find` and `grep` available on `PATH` (standard on Linux runners).

No extra Python packages are required.

## How to run

From the repository root:

```bash
./scripts/ukbb_ci_compliance_audit.sh
```

With optional environment overrides:

```bash
MIN_SENSITIVE_BYTES=2097152 WORK=/data/ukb_isolated ./scripts/ukbb_ci_compliance_audit.sh
```

## Exit codes

| Code | Meaning |
|------|---------|
| `0` | No blocking violations. |
| `1` | At least one violation; the job should fail. |

Warnings (for example, no `.py`/`.sh` files under `REPO_ROOT`) do not by themselves change the exit code unless you extend the script to treat them as fatal.

## Environment variables

| Variable | Default | Description |
|----------|---------|-------------|
| `REPO_ROOT` | `.` | Root directory to scan. In CI, run from the repo root or set this to the clone path. |
| `MIN_SENSITIVE_BYTES` | `1048576` (1 MiB) | Minimum file size in **bytes** for `.csv`, `.tsv`, and `.ped` to be reported. Smaller files are ignored for these extensions only. |
| `WORK` | *(unset)* | If set, quoted absolute paths under this directory are **allowed** in `.py`/`.sh` (isolated work area, e.g. HPC `$WORK`). Trailing slashes are normalised away. |

The script reads **`WORK`**, not `WORK_ROOT`; `WORK_ROOT` is an internal normalisation of `WORK`.

## Check 1: Genomic and clinical-style files

**Rationale:** Individual-level or bulk genomic/clinical artefacts must not live in version control under typical data-use agreements; committing them breaks separation between public code and restricted data.

**Method:** `find` lists files under `REPO_ROOT` (with common junk paths pruned), then `grep` filters paths by extension (null-terminated stream).

**Always reported (any size):** paths whose names match restricted genomic-style suffixes, including for example:

- `.vcf`, `.vcf.gz`, `.vcf.bgz`, `.bcf`, `.bcf.csi`
- `.bam`, `.cram`, `.bai`, `.crai`
- `.bed`, `.bed.gz`, `.ped`

**Reported only if size ≥ `MIN_SENSITIVE_BYTES`:** `.csv`, `.tsv`, `.ped` (again matched on the path). Small tracked fixtures may be acceptable under local policy; large tabular or pedigree files are treated as high risk.

**Pruned directories (not scanned as part of the tree walk logic):** `.git`, `node_modules`, `.venv`, `venv`, `__pycache__`.

## Check 2: Hardcoded absolute paths in `.py` and `.sh`

**Rationale:** Literals such as `/home/...`, `/data/...`, or `/Users/...` encode machine-specific layout and can signal where restricted copies might live. CI should push contributors toward environment variables, config, or allowed roots.

**Method:** For each `*.py` and `*.sh` file, `grep` extracts **double-quoted or single-quoted** substrings that start with `/`. Each candidate is classified as allowed or not.

**Allowed path prefixes (after unquoting):**

- `/tmp`, `/scratch`, `/mnt/scratch`
- Any path under `$WORK` when `WORK` is set
- Any string containing a literal `$WORK` (runtime resolution)
- Typical OS and package roots: `/usr`, `/bin`, `/sbin`, `/lib`, `/lib64`, `/opt`, `/etc`, `/dev`, `/proc`, `/sys`, `/run`
- macOS application locations: `/Applications`, `/Library`

Anything else that looks like an absolute path in a quoted string is a **violation**.

**Limitations:**

- Only **quoted** paths are detected; unquoted absolutes in shell (for example) are not matched by this heuristic.
- URLs, documentation, or comments that contain quoted strings starting with `/` may occasionally produce noise; prefer rewording or using allowed roots.
- The audit script avoids certain patterns in its own source so it does not false-positive on itself.

## Output

The script prints a short banner, each check as a subsection, lines prefixed with `[ OK ]`, `[FAIL]`, or `[WARN]`, and a summary. ANSI colours are used so warnings and errors stand out in GitHub Actions logs.

## GitHub Actions example

```yaml
jobs:
  ukb-repo-audit:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4

      - name: UKB-oriented repository compliance audit
        run: ./scripts/ukbb_ci_compliance_audit.sh
        env:
          # Optional: allow literals under your isolated work filesystem
          # WORK: ${{ vars.UKB_ISOLATED_WORK_ROOT }}
          MIN_SENSITIVE_BYTES: '1048576'
```

Ensure the step runs from the repository root (default for `actions/checkout`).

## When the job fails

1. **Restricted or large sensitive files:** Remove them from the branch; if they were ever committed, purge them from Git history per your security process and rotate any exposed credentials if applicable.
2. **Disallowed absolute paths:** Replace literals with environment variables, configuration files (not committed if they contain secrets), or paths under `/tmp`, `/scratch`, `/mnt/scratch`, or `$WORK` as appropriate.

## Related tooling

| Tool | Scope |
|------|--------|
| `scripts/ukbb_ci_compliance_audit.sh` | Full tree on CI; files + quoted paths in `.py`/`.sh` |
| `scripts/security_check.sh` | Staged files; content patterns and some extensions |
| `notebooks/03_validation_and_compliance/UKB_Compliance_Auditor.ipynb` | Interactive audit (EIDs, patterns, extensions) |

---

**Maintained with:** ROGEN compliance and DevSecOps practices. Update this document when the script’s behaviour or variables change.
