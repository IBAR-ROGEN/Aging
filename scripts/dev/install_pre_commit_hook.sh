#!/usr/bin/env bash
#
# Install the repository pre-commit framework (Black, isort, flake8, mypy,
# genomics schema check, UK Biobank security scan).
# Run from the repository root.
#
# This repo uses core.hooksPath=.githooks, so the hook is installed as
# .githooks/pre-commit rather than via `pre-commit install` into .git/hooks.
#

set -euo pipefail

REPO_ROOT="$(git rev-parse --show-toplevel 2>/dev/null)" || {
  echo "Not in a git repo."
  exit 1
}
cd "$REPO_ROOT"

if ! command -v uv >/dev/null 2>&1; then
  echo "uv is required. Install: https://docs.astral.sh/uv/"
  exit 1
fi

uv sync --extra dev

HOOK_PATH="$REPO_ROOT/.githooks/pre-commit"
chmod +x "$HOOK_PATH"
chmod +x "$REPO_ROOT/.githooks/prepare-commit-msg" 2>/dev/null || true

# Warm hook environments / confirm config parses.
uv run pre-commit validate-config
uv run pre-commit install-hooks

echo "Pre-commit hooks ready at $HOOK_PATH"
echo "See docs/CONTRIBUTING.md for usage and bypass options."
