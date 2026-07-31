"""Shared I/O helpers for integrative CLI scripts (Activity A.2.1.11.1).

Production defaults consume the July annotation workbook (or its parquet
equivalents) and write results under ``analysis/integrative/results/``.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl

from rogen_aging.config import cfg_path, find_repo_root, get_config

REPO_ROOT = find_repo_root()

# July annotation artefacts (production). Seeded from config; refreshed by set_config.
_cfg = get_config()
DEFAULT_JULY_XLSX = cfg_path(_cfg, "paths", "integrative", "july_xlsx")
DEFAULT_VARIANTS_PARQUET = cfg_path(_cfg, "paths", "integrative", "variants_parquet")
DEFAULT_EQTLS_PARQUET = cfg_path(_cfg, "paths", "integrative", "eqtls_parquet")
DEFAULT_EQTLS_CSV = cfg_path(_cfg, "paths", "integrative", "eqtls_csv")
DEFAULT_ALPHAGENOME = cfg_path(_cfg, "paths", "integrative", "alphagenome")
DEFAULT_OUTPUT_DIR = cfg_path(_cfg, "paths", "integrative", "output_dir")
del _cfg

# Pre-summarised GTEx columns on Combined_Master; dropped so the mapper rebuilds them.
_GTEX_SUMMARY_COLS = (
    "gtex_variant_id",
    "gtex_n_eqtls",
    "gtex_best_tissue",
    "gtex_best_gene",
    "gtex_best_slope",
    "gtex_best_p_value",
    "gtex_tissues",
)

_COLUMN_ALIASES = {
    "rsID": "rsid",
    "RsID": "rsid",
    "SNP": "rsid",
    "snp": "rsid",
    "slope": "nes",
    "NES": "nes",
}


def resolve_default_variants() -> Path:
    """Return the preferred production annotated-variant table path.

    Prefers the Combined_Master parquet when present, otherwise the July Excel
    workbook.

    Returns:
        Path to the production annotated-variant table.

    Raises:
        FileNotFoundError: If neither parquet nor Excel artefact exists.
    """
    if DEFAULT_VARIANTS_PARQUET.is_file():
        return DEFAULT_VARIANTS_PARQUET
    if DEFAULT_JULY_XLSX.is_file():
        return DEFAULT_JULY_XLSX
    raise FileNotFoundError(
        "July production annotation missing. Expected "
        f"{DEFAULT_VARIANTS_PARQUET} or {DEFAULT_JULY_XLSX}. "
        "Run: uv run python scripts/ukb/run_july_annotation_pipeline.py"
    )


def resolve_default_eqtls() -> Path:
    """Return the preferred production long eQTL table path.

    Prefers the July GTEx parquet, then the analysis GTEx CSV.

    Returns:
        Path to the production eQTL table.

    Raises:
        FileNotFoundError: If no known eQTL artefact exists.
    """
    if DEFAULT_EQTLS_PARQUET.is_file():
        return DEFAULT_EQTLS_PARQUET
    if DEFAULT_EQTLS_CSV.is_file():
        return DEFAULT_EQTLS_CSV
    if DEFAULT_JULY_XLSX.is_file():
        return DEFAULT_JULY_XLSX
    raise FileNotFoundError(
        "Production eQTL table missing. Expected "
        f"{DEFAULT_EQTLS_PARQUET} or {DEFAULT_EQTLS_CSV}."
    )


def ensure_july_parquet_cache(*, force: bool = False) -> tuple[Path, Path]:
    """Materialise Combined_Master + GTEx parquet siblings from the July Excel.

    Args:
        force: Rewrite parquet even when files already exist.

    Returns:
        ``(variants_parquet, eqtls_parquet)`` paths.

    Raises:
        FileNotFoundError: If the July Excel workbook is absent.
    """
    if not force and DEFAULT_VARIANTS_PARQUET.is_file() and DEFAULT_EQTLS_PARQUET.is_file():
        return DEFAULT_VARIANTS_PARQUET, DEFAULT_EQTLS_PARQUET
    if not DEFAULT_JULY_XLSX.is_file():
        raise FileNotFoundError(
            f"July workbook not found: {DEFAULT_JULY_XLSX}. "
            "Run: uv run python scripts/ukb/run_july_annotation_pipeline.py"
        )
    import pandas as pd

    DEFAULT_VARIANTS_PARQUET.parent.mkdir(parents=True, exist_ok=True)
    master = pl.from_pandas(pd.read_excel(DEFAULT_JULY_XLSX, sheet_name="Combined_Master"))
    eqtl = pl.from_pandas(pd.read_excel(DEFAULT_JULY_XLSX, sheet_name="GTEx_eQTL_Summary"))
    master.write_parquet(DEFAULT_VARIANTS_PARQUET)
    normalize_integrative_columns(eqtl).write_parquet(DEFAULT_EQTLS_PARQUET)
    return DEFAULT_VARIANTS_PARQUET, DEFAULT_EQTLS_PARQUET


def normalize_integrative_columns(frame: pl.DataFrame) -> pl.DataFrame:
    """Rename common production aliases onto integrative column names.

    Args:
        frame: Raw input table.

    Returns:
        Frame with ``rsid`` / ``nes`` aliases applied when needed.
    """
    renames = {
        src: dst
        for src, dst in _COLUMN_ALIASES.items()
        if src in frame.columns and dst not in frame.columns
    }
    return frame.rename(renames) if renames else frame


def strip_gtex_summary_columns(frame: pl.DataFrame) -> pl.DataFrame:
    """Drop pre-joined GTEx summary columns so the mapper can rebuild them.

    Args:
        frame: Annotated variants (e.g. Combined_Master).

    Returns:
        Frame without ``gtex_*`` summary columns.
    """
    drop = [c for c in _GTEX_SUMMARY_COLS if c in frame.columns]
    return frame.drop(drop) if drop else frame


def read_table(path: Path, *, sheet_name: str | None = None) -> pl.DataFrame:
    """Load a CSV, TSV, Excel, or Parquet table and normalise column aliases.

    For multi-sheet July Excel workbooks, ``sheet_name`` selects the sheet.
    When ``path`` is the July workbook and ``sheet_name`` is omitted, defaults
    to ``Combined_Master``.

    Args:
        path: Input table path.
        sheet_name: Optional Excel sheet name.

    Returns:
        Polars DataFrame with integrative column aliases applied.
    """
    suffix = path.suffix.lower()
    if suffix == ".parquet":
        frame = pl.read_parquet(path)
    elif suffix in {".xlsx", ".xls"}:
        import pandas as pd

        resolved_sheet = sheet_name
        if resolved_sheet is None and path.resolve() == DEFAULT_JULY_XLSX.resolve():
            resolved_sheet = "Combined_Master"
        frame = pl.from_pandas(pd.read_excel(path, sheet_name=resolved_sheet or 0))
    elif suffix == ".tsv":
        frame = pl.read_csv(path, separator="\t")
    else:
        frame = pl.read_csv(path)
    return normalize_integrative_columns(frame)


def load_production_variants(path: Path | None = None) -> pl.DataFrame:
    """Load production annotated variants, stripping stale GTEx summaries.

    Args:
        path: Optional override path. Defaults to :func:`resolve_default_variants`.

    Returns:
        Annotated variant frame ready for :class:`VariantTissueMapper`.
    """
    resolved = path if path is not None else resolve_default_variants()
    sheet = "Combined_Master" if resolved.suffix.lower() in {".xlsx", ".xls"} else None
    return strip_gtex_summary_columns(read_table(resolved, sheet_name=sheet))


def load_production_eqtls(path: Path | None = None) -> pl.DataFrame:
    """Load the production long eQTL table.

    Args:
        path: Optional override path. Defaults to :func:`resolve_default_eqtls`.

    Returns:
        Long eQTL frame with ``rsid`` / ``nes`` / ``tissue`` / ``p_value``.
    """
    resolved = path if path is not None else resolve_default_eqtls()
    sheet = None
    if resolved.suffix.lower() in {".xlsx", ".xls"}:
        sheet = "GTEx_eQTL_Summary"
    return read_table(resolved, sheet_name=sheet)


__all__ = [
    "DEFAULT_ALPHAGENOME",
    "DEFAULT_EQTLS_CSV",
    "DEFAULT_EQTLS_PARQUET",
    "DEFAULT_JULY_XLSX",
    "DEFAULT_OUTPUT_DIR",
    "DEFAULT_VARIANTS_PARQUET",
    "REPO_ROOT",
    "ensure_july_parquet_cache",
    "load_production_eqtls",
    "load_production_variants",
    "normalize_integrative_columns",
    "read_table",
    "resolve_default_eqtls",
    "resolve_default_variants",
    "strip_gtex_summary_columns",
]
