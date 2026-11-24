# August branches status vs `main`

Read-only comparison. No merge, push, or branch deletion was performed.

## Provenance

| Field | Value |
|-------|-------|
| timestamp_utc | 2026-11-24T09:20:15Z |
| origin/main | `f7c4c3d4f59b57cfbe648012735768ebadf678cb` |
| workspace HEAD | `430f0558a0ad899c0f7dab0ff63fc7128ee0bf71` |
| git | 2.50.1 (Apple Git-155) |
| method | `git rev-list --left-right --count origin/main...<tip>`; `git diff --name-status origin/main...<tip>`; `git merge-tree --write-tree origin/main <tip>` |

Ahead/behind: left count = commits on `origin/main` not in tip (**behind**); right count = commits on tip not in `origin/main` (**ahead**).

## Summary table

| Branch | Tip | Ahead | Behind | Files changed vs main (`...`) | `git merge-tree --write-tree` |
|--------|-----|------:|-------:|-------------------------------|-------------------------------|
| `fix/ci-omitted-tests-mypy` | `bb6249c` | 1 | 51 | 7 | **conflicts** (exit 1) |
| `feature/july-clock-console-entries` | `689d301` | 1 | 51 | 26 | **conflicts** (exit 1) |
| `fix/ensembl-clock-eval-dedup` | `cbf355c` | 2 | 51 | 34 | **conflicts** (exit 1) |
| `chore/ruff-unify-tooling` | `704df36` | 1 | 51 | 67 | **clean** (exit 0, tree `a1fb45c…`) |

## Per branch

### `fix/ci-omitted-tests-mypy` (`bb6249c`)

- Unique commit: `bb6249c` Run config and methylation preflight tests plus mypy in CI.
- Ahead 1 / behind 51.
- Files changed:
  - `M` `.github/workflows/ci.yml`
  - `M` `README.md`
  - `M` `docs/CLOCK_LIBRARY.md`
  - `M` `docs/PROJECT_STRUCTURE.md`
  - `M` `docs/UKB_PRE_COMMIT_HOOK.md`
  - `M` `docs/WORKFLOWS.md`
  - `M` `scripts/dev/security_check.sh`
- Merge-tree conflicts:
  - `docs/CLOCK_LIBRARY.md`

### `feature/july-clock-console-entries` (`689d301`)

- Unique commit: `689d301` Promote July annotation and GSE87571 clock eval to package CLIs.
- Ahead 1 / behind 51.
- Files changed (26):
  - `M` `INPUT_MANIFEST.md`
  - `M` `README.md`
  - `M` `docs/ACTIVITIES.md`
  - `M` `docs/CLOCK_LIBRARY.md`
  - `M` `docs/CODE_MODULES_REFERENCE.md`
  - `M` `docs/FIGURES.md`
  - `M` `docs/GENOMICS_ANALYSIS.md`
  - `M` `docs/JULY_ANNOTATION_PIPELINE.md`
  - `M` `docs/METHYLATION_CLOCK_VALIDATION.md`
  - `M` `docs/PROJECT_STRUCTURE.md`
  - `M` `docs/WORKFLOWS.md`
  - `M` `pyproject.toml`
  - `M` `scripts/clock/evaluate_methylation_clock.py`
  - `M` `scripts/dev/write_pipeline_fixtures.py`
  - `M` `scripts/ukb/run_july_annotation_pipeline.py`
  - `A` `src/rogen_aging/annotation/__init__.py`
  - `A` `src/rogen_aging/annotation/july.py`
  - `M` `src/rogen_aging/cli/clock.py`
  - `A` `src/rogen_aging/cli/july_annotate.py`
  - `A` `src/rogen_aging/clock/gse87571.py`
  - `M` `src/rogen_aging/integrative/io.py`
  - `M` `src/rogen_aging/integrative/variant_tissue_mapper.py`
  - `M` `src/rogen_aging/pipeline_fixtures.py`
  - `M` `tests/test_evaluate_methylation_clock.py`
  - `M` `tests/test_package_imports.py`
  - `M` `tests/test_root_shims.py`
- Merge-tree conflicts:
  - `INPUT_MANIFEST.md`
  - `docs/ACTIVITIES.md`
  - `docs/CODE_MODULES_REFERENCE.md`
  - `docs/WORKFLOWS.md`
  - `scripts/clock/evaluate_methylation_clock.py`
  - `src/rogen_aging/cli/clock.py`

### `fix/ensembl-clock-eval-dedup` (`cbf355c`)

- Unique commits: `cbf355c` Route VEP and GSE87571 clock loading through shared package APIs.; `689d301` Promote July annotation and GSE87571 clock eval to package CLIs.
- Ahead 2 / behind 51.
- Files changed (34): same set as `689d301` plus:
  - `M` `docs/LA_SNP_VEP_ANNOTATION.md`
  - `M` `scripts/ukb/annotate_la_snps_vep.py`
  - `M` `src/rogen_aging/clock/__init__.py`
  - `M` `src/rogen_aging/clock/evaluate.py`
  - `M` `src/rogen_aging/ensembl/__init__.py`
  - `M` `src/rogen_aging/ensembl/client.py`
  - `A` `src/rogen_aging/ensembl/vep.py`
  - `M` `tests/test_ensembl_client.py`
- Merge-tree conflicts (same six paths as `689d301`):
  - `INPUT_MANIFEST.md`
  - `docs/ACTIVITIES.md`
  - `docs/CODE_MODULES_REFERENCE.md`
  - `docs/WORKFLOWS.md`
  - `scripts/clock/evaluate_methylation_clock.py`
  - `src/rogen_aging/cli/clock.py`

### `chore/ruff-unify-tooling` (`704df36`)

- Unique commit: `704df36` Unify lint and format on Ruff; drop Black, isort, and flake8.
- Ahead 1 / behind 51.
- Files changed (67): includes deleting `.flake8`, updating `.github/workflows/ci.yml`, `.pre-commit-config.yaml`, `pyproject.toml`, `uv.lock`, many `scripts/**` reformats, and selected `src/` / `tests/` files (full list via `git diff --name-status origin/main...704df36`).
- Merge-tree: **no conflicts** (`git merge-tree --write-tree` exit 0).
