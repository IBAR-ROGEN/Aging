#!/usr/bin/env bash
#
# Git Pre-Commit Security Hook - UK Biobank Data Protection
#
# Scans staged files for potential patient identifiers, UKB-specific patterns,
# and restricted genomic file formats (.vcf, .bed). Blocks commits that violate
# UK Biobank Data Usage Agreement (DUA) and ROGEN compliance policies.
#
# Usage: Run automatically via pre-commit hook, or manually:
#   ./scripts/security_check.sh
#
# Install as pre-commit hook:
#   cp scripts/security_check.sh .git/hooks/pre-commit && chmod +x .git/hooks/pre-commit
#   # Or run: ./scripts/install_pre_commit_hook.sh
#

set -e

RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m'

RESTRICTED_EXTENSIONS='\.(vcf|vcf\.gz|vcf\.bgz|bed|bed\.gz)$'
RESTRICTED_PATTERNS='patient_id|UKB_'

violations=0
staged_files=$(git diff --cached --name-only 2>/dev/null || true)

if [ -z "$staged_files" ]; then
  exit 0
fi

echo "Running UK Biobank security check on staged files..."
echo ""

# Check file names for restricted extensions
while IFS= read -r file; do
  [ -z "$file" ] && continue

  # Check restricted file extensions (no exceptions)
  if echo "$file" | grep -qiE "$RESTRICTED_EXTENSIONS"; then
    echo -e "${RED}BLOCKED:${NC} Staged file uses restricted format: $file"
    echo "  -> .vcf and .bed files may contain UK Biobank participant data."
    violations=$((violations + 1))
    continue
  fi

  # Skip content scan for compliance docs and the security tool itself
  case "$file" in
    docs/UKB_COMPLIANCE*|README*|*security_check*|*install_pre_commit*) continue ;;
  esac

  # Scan file contents for restricted patterns (-a treats binary as text for grep)
  if git show ":$file" 2>/dev/null | grep -aqiE "$RESTRICTED_PATTERNS"; then
    echo -e "${RED}BLOCKED:${NC} Staged file contains sensitive pattern: $file"
    echo "  -> Found 'patient_id', 'UKB_', or similar. Remove before committing."
    violations=$((violations + 1))
  fi
done <<< "$staged_files"

if [ $violations -gt 0 ]; then
  echo ""
  echo -e "${YELLOW}========================================${NC}"
  echo -e "${YELLOW}  UK BIOBANK DATA SECURITY WARNING${NC}"
  echo -e "${YELLOW}========================================${NC}"
  echo ""
  echo "Your commit was BLOCKED to prevent accidental exposure of restricted data."
  echo ""
  echo "UK Biobank Data Usage Agreement prohibits committing:"
  echo "  - Participant IDs (EIDs) or patient identifiers"
  echo "  - UKB-specific column names or paths"
  echo "  - Raw genomic files (.vcf, .bed) that may contain individual-level data"
  echo ""
  echo "Keep sensitive data in the git-ignored data/ directory."
  echo "See docs/UKB_COMPLIANCE_AUDITOR.md for compliance details."
  echo ""
  exit 1
fi

exit 0
