#!/usr/bin/env bash
#
# Install the UK Biobank security check as a Git pre-commit hook.
# Run from the repository root.
#

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(git rev-parse --show-toplevel 2>/dev/null)" || { echo "Not in a git repo."; exit 1; }
HOOK_PATH="$REPO_ROOT/.git/hooks/pre-commit"

cp "$SCRIPT_DIR/security_check.sh" "$HOOK_PATH"
chmod +x "$HOOK_PATH"
echo "Pre-commit security hook installed at $HOOK_PATH"
