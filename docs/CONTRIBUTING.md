# Contributing to rogen_aging

Python ≥3.12, managed with [uv](https://docs.astral.sh/uv/). Quality gates run
through [pre-commit](https://pre-commit.com/) (Black, isort, flake8, mypy, plus
a lightweight genomics schema check and the UK Biobank security scan).

## Install

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
uv sync --extra dev
./scripts/dev/install_pre_commit_hook.sh
```

This repository sets `core.hooksPath=.githooks`, so the installer enables
`.githooks/pre-commit` (it does not call `pre-commit install` into `.git/hooks`).
Re-run `./scripts/dev/install_pre_commit_hook.sh` after cloning if hooks are missing.

## Git workflow

### Branch naming

Use standard git-flow prefixes with a concise hyphen-separated description.
Do **not** use a `cursor/` prefix.

| Prefix | Use for |
|--------|---------|
| `feature/` | New functionality |
| `fix/` | Bug fixes |
| `refactor/` | Restructuring without behavior change |
| `docs/` | Documentation only |
| `chore/` | Tooling, CI, deps, housekeeping |

Examples: `feature/methylation-matrix-preflight`, `fix/ukb-eid-join-audit`,
`chore/pre-commit-hooks`.

### Commits and attribution

- Configure `user.name` / `user.email` to match a **verified** email on your GitHub
  account (or your GitHub noreply address) so contributions attribute correctly.
- Do not add AI / agent attribution trailers (`Co-authored-by: Cursor`,
  `Made-with: Cursor`, and similar). Subject and body only.
- `.githooks/prepare-commit-msg` strips common Cursor attribution trailers if they
  are injected by the environment.

## What runs on commit

| Hook | Purpose |
|------|---------|
| `black` | Format Python to the project style (`line-length = 100`) |
| `isort` | Sort imports (Black-compatible profile) |
| `flake8` | Lint for unused imports and common errors |
| `mypy` | Strict static typing on `src/rogen_aging/` (pandas schemas via `pandas_schemas.py`) |
| `validate-genomics-tables` | Offline schema check of the gene–LA-SNP overlap table on `test_data/genomics_overlap_minimal.csv` |
| `ukb-security-check` | Block staged UKB / patient-identifying content (see [UKB_PRE_COMMIT_HOOK.md](UKB_PRE_COMMIT_HOOK.md)) |

## Useful commands

```bash
# Run every hook against the whole tree
uv run pre-commit run --all-files

# Run a single hook
uv run pre-commit run mypy --all-files
uv run pre-commit run validate-genomics-tables --all-files

# Manual genomics schema check (same as the hook)
uv run python scripts/dev/validate_genomics_tables.py \
  --fixture test_data/genomics_overlap_minimal.csv

# Full (networked) genomics audit — not part of pre-commit
uv run python analysis/validate_genomics_tables/validate_genomics_tables.py \
  --input overlapping_genes_with_snps.xlsx \
  --output-dir results
```

## Bypass hooks (emergency only)

Hooks exist to catch formatting, type, schema, and UKB compliance issues before
they land on `main`. Skip them only when you understand the risk:

```bash
# Skip all pre-commit hooks for one commit
git commit --no-verify -m "emergency: reason"

# Equivalent short flag
git commit -n -m "emergency: reason"
```

Do **not** use `--no-verify` to push past UK Biobank security failures. Fix the
staged content instead (see [UKB_PRE_COMMIT_HOOK.md](UKB_PRE_COMMIT_HOOK.md)).

To re-enable hooks after a local disable, re-run:

```bash
./scripts/dev/install_pre_commit_hook.sh
```

## Typing notes

`mypy` runs in **strict** mode on `src/rogen_aging/`. Pandas DataFrame column
contracts are declared in `src/rogen_aging/pandas_schemas.py` (`TypedDict` +
runtime `assert_*_schema` helpers) because mypy cannot see DataFrame columns
directly. Prefer those helpers when loading tabular inputs.
