#!/usr/bin/env bash
#
# UK Biobank / sensitive-data compliance audit for CI/CD (e.g. GitHub Actions).
#
# Rationale: Version control must not host individual-level genomic/clinical
# artefacts or machine-specific paths that reveal where restricted data live.
# This script automates checks aligned with typical DUA-style obligations;
# it does not replace legal review or UK Biobank's own guidance.
#
# Usage:
#   ./scripts/ukbb_ci_compliance_audit.sh
#   MIN_SENSITIVE_BYTES=2097152 ./scripts/ukbb_ci_compliance_audit.sh
#
# Exit codes: 0 = no violations, 1 = one or more violations (fail the pipeline).
#
# Documentation: docs/UKBB_CI_COMPLIANCE_AUDIT.md
#

set -euo pipefail

# --- Terminal formatting (visible in GitHub Actions logs; harmless when not a TTY) ---
RED=$'\033[0;31m'
YELLOW=$'\033[1;33m'
GREEN=$'\033[0;32m'
BLUE=$'\033[0;34m'
BOLD=$'\033[1m'
NC=$'\033[0m'

# Minimum size (bytes) for .csv / .ped to count as a violation.
# Small tracked examples (fixtures, schemas) may be allowed by policy; large files
# are high risk for containing participant-level rows. Genomic extensions are
# always flagged regardless of size — they must not live in git for UKB-style work.
MIN_SENSITIVE_BYTES="${MIN_SENSITIVE_BYTES:-1048576}"

# Optional: repository root (default: current directory). CI should run from repo root.
REPO_ROOT="${REPO_ROOT:-.}"

# Normalised optional isolated work root from environment (e.g. HPC $WORK).
WORK_ROOT="${WORK:-}"
WORK_ROOT="${WORK_ROOT%/}"

violations=0
warn_count=0

section() {
  printf '\n%s%s%s\n' "$BLUE" "$1" "$NC"
}

subsection() {
  printf '%s%s%s\n' "$BOLD" "$1" "$NC"
}

pass() {
  printf '  %s[ OK ]%s %s\n' "$GREEN" "$NC" "$1"
}

warn_line() {
  printf '  %s[WARN]%s %s\n' "$YELLOW" "$NC" "$1"
  warn_count=$((warn_count + 1))
}

fail_line() {
  printf '  %s[FAIL]%s %s\n' "$RED" "$NC" "$1"
  violations=$((violations + 1))
}

# Security rationale: do not traverse VCS metadata or ephemeral dependency trees;
# those are not "our" code/data but inflate scans and produce false paths.
should_prune_path() {
  case "$1" in
    */.git/* | */.git | */node_modules/* | */.venv/* | */venv/* | */__pycache__/*) return 0 ;;
    *) return 1 ;;
  esac
}

# --- Check 1: sensitive / large genomic and clinical filenames in the tree ---
# Security rationale: .vcf/.bam/.bed (and friends) can carry individual-level
# genetic or functional data; committing them breaks segregation from the public
# repo and typical data-use agreements. .csv/.ped are flagged when large because
# they often hold phenotypes, pedigrees, or joined tables that must stay off-git.
scan_sensitive_files() {
  subsection "1. Genomic / clinical file presence (find + grep on paths)"

  local -a hits=()
  local f rel

  while IFS= read -r -d '' f; do
    rel="${f#./}"
    if should_prune_path "$f"; then
      continue
    fi
    hits+=("$rel")
  done < <(
    # find: locate candidates; grep: used on the path stream for extension matching
    # (explicit requirement — avoids relying on shell glob edge cases across find versions).
    find "$REPO_ROOT" \
      \( -path "$REPO_ROOT/.git" -o -path "*/.git/*" -o -path "*/node_modules/*" \
      -o -path "*/.venv/*" -o -path "*/venv/*" -o -path "*/__pycache__/*" \) -prune -o \
      -type f -print0 |
      grep -zE \
        '\.(vcf|vcf\.gz|vcf\.bgz|bcf|bcf\.csi|bam|cram|bai|crai|bed|bed\.gz|ped)(\.[^/]*)?$' || true
  )

  local csv_ped_hits=()
  while IFS= read -r -d '' f; do
    rel="${f#./}"
    if should_prune_path "$f"; then
      continue
    fi
    csv_ped_hits+=("$rel")
  done < <(
    find "$REPO_ROOT" \
      \( -path "$REPO_ROOT/.git" -o -path "*/.git/*" -o -path "*/node_modules/*" \
      -o -path "*/.venv/*" -o -path "*/venv/*" -o -path "*/__pycache__/*" \) -prune -o \
      -type f -size "+${MIN_SENSITIVE_BYTES}c" -print0 |
      grep -zE '\.(csv|tsv|ped)(\.[^/]*)?$' || true
  )

  if ((${#hits[@]})); then
    for rel in "${hits[@]}"; do
      fail_line "Restricted genomic/clinical-style file in tree: ${rel}"
    done
  else
    pass "No restricted genomic extensions (.vcf, .bam, .bed, .ped, etc.) found under ${REPO_ROOT}."
  fi

  if ((${#csv_ped_hits[@]})); then
    for rel in "${csv_ped_hits[@]}"; do
      fail_line "Large tabular / pedigree file (>= ${MIN_SENSITIVE_BYTES} bytes): ${rel}"
    done
  else
    pass "No large .csv/.tsv/.ped files (>= ${MIN_SENSITIVE_BYTES} bytes) found."
  fi
}

# Return 0 if path $1 is an allowed hardcoded absolute reference, 1 otherwise.
# Security rationale: only ephemeral (/tmp), shared scratch, or the lab-declared
# isolated work area should appear as literals; other absolutes often encode
# home directories, local mount points, or project paths that leak data layout.
is_allowed_abs_path() {
  local p="$1"
  local work_prefix

  [[ "$p" == /tmp/* || "$p" == /tmp ]] && return 0
  [[ "$p" == /scratch/* || "$p" == /scratch ]] && return 0
  # Many clusters expose scratch under /mnt; treat /mnt/scratch like /scratch.
  [[ "$p" == /mnt/scratch/* || "$p" == /mnt/scratch ]] && return 0

  # Use a separate prefix variable so the script source never contains a closing
  # double-quote of WORK_ROOT immediately followed by a slash. Otherwise the path
  # scanner grep matches that token and raises a spurious self-violation.
  if [[ -n "$WORK_ROOT" ]]; then
    work_prefix="${WORK_ROOT}/"
    [[ "$p" == "${work_prefix}"* || "$p" == "$WORK_ROOT" ]] && return 0
  fi

  # Literal use of $WORK inside a string is acceptable (resolved at runtime).
  [[ "$p" == *\$WORK* ]] && return 0

  # Common OS / packaging roots — not UKB data locations; excluding reduces noise in CI.
  [[ "$p" == /usr/* || "$p" == /usr ]] && return 0
  [[ "$p" == /bin/* || "$p" == /bin ]] && return 0
  [[ "$p" == /sbin/* || "$p" == /sbin ]] && return 0
  [[ "$p" == /lib/* || "$p" == /lib ]] && return 0
  [[ "$p" == /lib64/* || "$p" == /lib64 ]] && return 0
  [[ "$p" == /opt/* || "$p" == /opt ]] && return 0
  [[ "$p" == /etc/* || "$p" == /etc ]] && return 0
  [[ "$p" == /dev/* || "$p" == /dev ]] && return 0
  [[ "$p" == /proc/* || "$p" == /proc ]] && return 0
  [[ "$p" == /sys/* || "$p" == /sys ]] && return 0
  [[ "$p" == /run/* || "$p" == /run ]] && return 0
  # macOS standard install locations (tooling, not project data directories).
  [[ "$p" == /Applications/* || "$p" == /Applications ]] && return 0
  [[ "$p" == /Library/* || "$p" == /Library ]] && return 0

  return 1
}

# Strip one layer of surrounding quotes from captured path token.
strip_quotes() {
  local s="$1"
  s="${s#\"}"
  s="${s%\"}"
  s="${s#\'}"
  s="${s%\'}"
  printf '%s' "$s"
}

# --- Check 2: hardcoded absolute paths in Python and shell scripts ---
# Security rationale: absolute paths in scripts tend to encode operator-specific
# filesystem layout; in regulated workflows that implies where restricted copies
# might exist and complicates reproducibility without exposing infrastructure.
scan_hardcoded_paths() {
  subsection "2. Hardcoded absolute paths in *.py and *.sh"

  local scanned=0
  local file path candidate

  while IFS= read -r -d '' file; do
    if should_prune_path "$file"; then
      continue
    fi
    scanned=$((scanned + 1))

    while IFS= read -r path; do
      [[ -z "$path" ]] && continue
      candidate="$(strip_quotes "$path")"
      if [[ "$candidate" != /* ]]; then
        continue
      fi
      if is_allowed_abs_path "$candidate"; then
        continue
      fi
      fail_line "${file#./}: disallowed absolute path literal: ${candidate}"
    done < <(
      # grep -oE: pull quoted substrings that start with / from each line.
      # shellcheck disable=SC2016
      # ERE-safe: avoid PCRE-only groups; tolerate simple escaped chars inside strings.
      grep -Eho "['\"]/([^'\"]|\\\\.)+['\"]" "$file" 2>/dev/null || true
    )
  done < <(find "$REPO_ROOT" \
    \( -path "$REPO_ROOT/.git" -o -path "*/.git/*" -o -path "*/node_modules/*" \
    -o -path "*/.venv/*" -o -path "*/venv/*" -o -path "*/__pycache__/*" \) -prune -o \
    \( -name '*.py' -o -name '*.sh' \) -type f -print0)

  if [[ "$scanned" -eq 0 ]]; then
    warn_line "No .py/.sh files found under ${REPO_ROOT} (check REPO_ROOT)."
  else
    pass "Scanned ${scanned} Python/shell files for quoted absolute paths."
  fi
}

# --- Banner / summary ---
printf '%s%s%s\n' "$BOLD" "UK Biobank-oriented repository compliance audit" "$NC"
printf 'Repository root: %s\n' "$REPO_ROOT"
if [[ -n "$WORK_ROOT" ]]; then
  printf 'WORK (allowed prefix): %s\n' "$WORK_ROOT"
else
  printf '%sWORK not set — only /tmp, /scratch, /mnt/scratch, and system roots are allowed as absolutes.%s\n' "$YELLOW" "$NC"
fi
printf 'Large tabular threshold: %s bytes\n' "$MIN_SENSITIVE_BYTES"

section "Checks"
scan_sensitive_files
scan_hardcoded_paths

section "Summary"
if [[ "$violations" -eq 0 ]]; then
  printf '%sNo blocking violations.%s\n' "$GREEN" "$NC"
  if [[ "$warn_count" -gt 0 ]]; then
    printf '%sWarnings emitted: %s%s\n' "$YELLOW" "$warn_count" "$NC"
  fi
  exit 0
fi

printf '%sAudit failed: %s violation(s).%s\n' "$RED" "$violations" "$NC"
printf '%sRemediate before merging: remove sensitive files from git history if needed, replace absolutes with config/env vars or allowed roots.%s\n' "$YELLOW" "$NC"
exit 1
