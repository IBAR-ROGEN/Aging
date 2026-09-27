#!/usr/bin/env python3
"""Build the 73-allele × 14-dataset eQTL Catalogue coverage matrix.

Reads ``results/eqtl_associations_full.csv`` and ``results/annotation_layer.csv``
and writes a categorical coverage table plus heatmap. Cells are one of four
statements; absence is never encoded as 0, NA-as-no-effect, or an imputed value.

Example:
    uv run python scripts/build_eqtl_coverage.py
"""

from __future__ import annotations

import hashlib
import sys
from datetime import UTC, datetime
from importlib import metadata
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import polars as pl  # noqa: E402
import typer  # noqa: E402
from matplotlib.colors import BoundaryNorm, ListedColormap  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

from rogen_aging.config import find_repo_root  # noqa: E402

REPO_ROOT = find_repo_root()

DEFAULT_ASSOCIATIONS = REPO_ROOT / "results" / "eqtl_associations_full.csv"
DEFAULT_ANNOTATION = REPO_ROOT / "results" / "annotation_layer.csv"
DEFAULT_MATRIX_OUT = REPO_ROOT / "results" / "eqtl_coverage_matrix.csv"
DEFAULT_FIGURE_OUT = REPO_ROOT / "results" / "fig_eqtl_coverage.png"

# Ground-truth counts from the annotation run (do not "fix" by collapsing tissues).
EXPECTED_ASSOC_ROWS = 1650
EXPECTED_ASSOC_DATASETS = 14
EXPECTED_ASSOC_RSIDS = 43
EXPECTED_ALLELE_ROWS = 73
EXPECTED_ANNOTATION_RSIDS = 47
EXPECTED_STATUS_COUNTS: dict[str, int] = {
    "tested": 43,
    "allele_not_represented_in_eQTL_Catalogue": 24,
    "variant_absent_from_all_14_datasets": 6,
}
EXPECTED_ABSENT_RSIDS = frozenset({"rs17044026", "rs1805329", "rs1868778", "rs1898422"})
EXPECTED_ABSENT_ALLELE_ROWS = 6
EXPECTED_FDR_VARIANTS = 16
FDR_THRESHOLD = 0.05

CELL_ASSOC = "association at FDR<0.05"
CELL_TESTED_NS = "tested, no association at FDR<0.05"
CELL_ALLELE_ABSENT = "allele not represented in resource"
CELL_VARIANT_ABSENT = "variant absent from dataset"

# Discrete legend order. Tested-not-significant is pale sand so it cannot be
# read as missing; the two absence statements are saturated orange and purple.
CELL_ORDER: tuple[str, ...] = (
    CELL_ASSOC,
    CELL_TESTED_NS,
    CELL_ALLELE_ABSENT,
    CELL_VARIANT_ABSENT,
)
CELL_COLORS: dict[str, str] = {
    CELL_ASSOC: "#1b9e77",
    CELL_TESTED_NS: "#f6e8c3",
    CELL_ALLELE_ABSENT: "#e66101",
    CELL_VARIANT_ABSENT: "#5e3c99",
}

INDEX_COLUMNS: tuple[str, ...] = (
    "gene",
    "canonical_rsid",
    "alt_allele",
    "eqtl_query_status",
)

# Anatomical / catalogue order. Cerebellum appears twice on purpose.
DATASET_ORDER: tuple[str, ...] = (
    "QTD000146",
    "QTD000151",
    "QTD000156",
    "QTD000161",
    "QTD000166",
    "QTD000171",
    "QTD000176",
    "QTD000181",
    "QTD000186",
    "QTD000191",
    "QTD000196",
    "QTD000201",
    "QTD000206",
    "QTD000356",
)

STATUS_TESTED = "tested"
STATUS_ALLELE_ABSENT = "allele_not_represented_in_eQTL_Catalogue"
STATUS_VARIANT_ABSENT = "variant_absent_from_all_14_datasets"


class CoverageCheckError(Exception):
    """An eQTL coverage input failed a snapshot check."""


def _fail(message: str) -> None:
    """Abort when a coverage input disagrees with the catalogue snapshot."""
    raise CoverageCheckError(message)


def _sha256_hex(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 16), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _package_version(name: str) -> str:
    try:
        return metadata.version(name)
    except metadata.PackageNotFoundError:
        return "not-installed"


def _require_columns(frame: pl.DataFrame, required: tuple[str, ...], label: str) -> None:
    missing = [c for c in required if c not in frame.columns]
    if missing:
        _fail(f"{label} is missing required columns: {missing}. Found: {frame.columns}")


def load_associations(path: Path) -> pl.DataFrame:
    """Load and assert the eQTL association table against ground-truth counts."""
    if not path.is_file():
        raise FileNotFoundError(f"Associations file not found: {path}")
    frame = pl.read_csv(path)
    _require_columns(
        frame,
        (
            "canonical_rsid",
            "eqtl_alt_allele",
            "tissue",
            "dataset_id",
            "n_samples",
            "fdr_bh",
        ),
        "eqtl_associations_full.csv",
    )
    n_datasets = frame["dataset_id"].n_unique()
    n_rsids = frame["canonical_rsid"].n_unique()
    n_tissues = frame["tissue"].n_unique()
    if frame.height != EXPECTED_ASSOC_ROWS:
        _fail(
            f"eqtl_associations_full.csv has {frame.height} rows; "
            f"expected {EXPECTED_ASSOC_ROWS}."
        )
    if n_datasets != EXPECTED_ASSOC_DATASETS:
        _fail(
            f"eqtl_associations_full.csv has {n_datasets} dataset_id values; "
            f"expected {EXPECTED_ASSOC_DATASETS}. Group by dataset_id, not tissue "
            f"(tissue n={n_tissues}; a 13-column tissue matrix is wrong because "
            "cerebellum maps to QTD000161 and QTD000166)."
        )
    if n_rsids != EXPECTED_ASSOC_RSIDS:
        _fail(
            f"eqtl_associations_full.csv has {n_rsids} unique rsIDs; "
            f"expected {EXPECTED_ASSOC_RSIDS}."
        )
    if n_tissues == 13 and n_datasets != 14:
        _fail("Refusing a tissue-collapsed matrix: cerebellum is two datasets.")
    observed_ids = set(frame["dataset_id"].unique().to_list())
    expected_ids = set(DATASET_ORDER)
    if observed_ids != expected_ids:
        _fail(
            "dataset_id set does not match the 14 GTEx Catalogue datasets: "
            f"extra={sorted(observed_ids - expected_ids)} "
            f"missing={sorted(expected_ids - observed_ids)}"
        )
    return frame


def load_annotation(path: Path) -> pl.DataFrame:
    """Load and assert the 73-row allele annotation table."""
    if not path.is_file():
        raise FileNotFoundError(f"Annotation file not found: {path}")
    frame = pl.read_csv(path)
    _require_columns(
        frame,
        INDEX_COLUMNS,
        "annotation_layer.csv",
    )
    if frame.height != EXPECTED_ALLELE_ROWS:
        _fail(
            f"annotation_layer.csv has {frame.height} rows; "
            f"expected {EXPECTED_ALLELE_ROWS} allele rows."
        )
    n_rsids = frame["canonical_rsid"].n_unique()
    if n_rsids != EXPECTED_ANNOTATION_RSIDS:
        _fail(
            f"annotation_layer.csv has {n_rsids} unique rsIDs; "
            f"expected {EXPECTED_ANNOTATION_RSIDS}."
        )
    status_counts = {
        str(row["eqtl_query_status"]): int(row["len"])
        for row in frame.group_by("eqtl_query_status").len().iter_rows(named=True)
    }
    if set(status_counts) != set(EXPECTED_STATUS_COUNTS):
        _fail(
            "eqtl_query_status must be exactly the three Catalogue statements "
            f"{sorted(EXPECTED_STATUS_COUNTS)}; got {sorted(status_counts)}."
        )
    for status, expected in EXPECTED_STATUS_COUNTS.items():
        got = status_counts.get(status, 0)
        if got != expected:
            _fail(f"eqtl_query_status={status!r} has {got} rows; expected {expected}.")

    absent = frame.filter(pl.col("eqtl_query_status") == STATUS_VARIANT_ABSENT)
    absent_rsids = set(absent["canonical_rsid"].unique().to_list())
    if absent.height != EXPECTED_ABSENT_ALLELE_ROWS:
        _fail(
            f"variant_absent_from_all_14_datasets has {absent.height} allele rows; "
            f"expected {EXPECTED_ABSENT_ALLELE_ROWS}."
        )
    if absent_rsids != set(EXPECTED_ABSENT_RSIDS):
        _fail(
            "Variants absent from all 14 datasets must be "
            f"{sorted(EXPECTED_ABSENT_RSIDS)}; got {sorted(absent_rsids)}. "
            "Both the 6-row and 4-variant counts are correct at their own level."
        )
    return frame


def dataset_columns(
    associations: pl.DataFrame,
    *,
    enforce_catalogue: bool = True,
) -> pl.DataFrame:
    """One row per dataset_id with tissue, sample size, and display label."""
    if enforce_catalogue:
        order = list(DATASET_ORDER)
    else:
        order = associations.get_column("dataset_id").unique(maintain_order=True).to_list()
    meta = (
        associations.select(["dataset_id", "tissue", "n_samples"])
        .unique()
        .with_columns(pl.col("dataset_id").cast(pl.Enum(order)).alias("_ord"))
        .sort("_ord")
        .drop("_ord")
    )
    if enforce_catalogue and meta.height != EXPECTED_ASSOC_DATASETS:
        _fail(
            f"Grouped dataset metadata has {meta.height} rows; "
            f"expected {EXPECTED_ASSOC_DATASETS} (cerebellum is two columns)."
        )
    labels = [
        f"{tissue} | {dataset_id} | n={n_samples}"
        for dataset_id, tissue, n_samples in meta.select(
            ["dataset_id", "tissue", "n_samples"]
        ).iter_rows()
    ]
    if len(set(labels)) != len(labels):
        _fail("Dataset column labels are not unique; studies that share a tissue name collided.")
    return meta.with_columns(pl.Series("column_label", labels, dtype=pl.String))


def classify_cells(
    annotation: pl.DataFrame,
    associations: pl.DataFrame,
    *,
    enforce_catalogue: bool = True,
) -> pl.DataFrame:
    """Return the long (allele × dataset) table with a categorical cell statement."""
    datasets = dataset_columns(associations, enforce_catalogue=enforce_catalogue)
    grid = annotation.select(list(INDEX_COLUMNS)).join(datasets, how="cross")

    hits = associations.group_by(["canonical_rsid", "eqtl_alt_allele", "dataset_id"]).agg(
        pl.col("fdr_bh").min().alias("min_fdr"),
        (pl.col("fdr_bh") < FDR_THRESHOLD).any().alias("any_fdr05"),
        pl.len().alias("n_hits"),
    )

    joined = grid.join(
        hits,
        left_on=["canonical_rsid", "alt_allele", "dataset_id"],
        right_on=["canonical_rsid", "eqtl_alt_allele", "dataset_id"],
        how="left",
    )

    unknown_status = joined.filter(
        ~pl.col("eqtl_query_status").is_in(
            [STATUS_TESTED, STATUS_ALLELE_ABSENT, STATUS_VARIANT_ABSENT]
        )
    )
    if unknown_status.height:
        _fail(
            "Unexpected eqtl_query_status values: "
            f"{unknown_status['eqtl_query_status'].unique().to_list()}"
        )

    cell = (
        pl.when(pl.col("eqtl_query_status") == STATUS_VARIANT_ABSENT)
        .then(pl.lit(CELL_VARIANT_ABSENT))
        .when(pl.col("eqtl_query_status") == STATUS_ALLELE_ABSENT)
        .then(pl.lit(CELL_ALLELE_ABSENT))
        .when(pl.col("eqtl_query_status") == STATUS_TESTED)
        .then(
            pl.when(pl.col("n_hits").is_null())
            .then(pl.lit(CELL_VARIANT_ABSENT))
            .when(pl.col("any_fdr05"))
            .then(pl.lit(CELL_ASSOC))
            .otherwise(pl.lit(CELL_TESTED_NS))
        )
        .otherwise(pl.lit("UNCLASSIFIED"))
    )
    classified = joined.with_columns(cell.alias("coverage_cell"))
    bad = classified.filter(pl.col("coverage_cell") == "UNCLASSIFIED")
    if bad.height:
        _fail(f"{bad.height} cells could not be classified.")

    null_cells = classified.filter(pl.col("coverage_cell").is_null())
    if null_cells.height:
        _fail(
            f"{null_cells.height} cells are null; absence must be an explicit "
            "category, never NA-meaning-no-effect."
        )
    return classified


def pivot_matrix(long_cells: pl.DataFrame, datasets: pl.DataFrame) -> pl.DataFrame:
    """Wide 73 × 14 matrix with index columns first."""
    wide = long_cells.pivot(
        on="column_label",
        index=list(INDEX_COLUMNS),
        values="coverage_cell",
        aggregate_function="first",
    )
    label_order = datasets["column_label"].to_list()
    wide = wide.select([*INDEX_COLUMNS, *label_order])
    wide = wide.sort(["gene", "canonical_rsid", "alt_allele"])
    return wide


def assert_matrix(wide: pl.DataFrame, associations: pl.DataFrame, annotation: pl.DataFrame) -> None:
    """Fail loudly if shape or FDR-significant row count is wrong."""
    dataset_cols = [c for c in wide.columns if c not in INDEX_COLUMNS]
    if wide.height != EXPECTED_ALLELE_ROWS:
        _fail(f"Coverage matrix has {wide.height} rows; expected {EXPECTED_ALLELE_ROWS}.")
    if len(dataset_cols) != EXPECTED_ASSOC_DATASETS:
        _fail(
            f"Coverage matrix has {len(dataset_cols)} dataset columns; "
            f"expected {EXPECTED_ASSOC_DATASETS}. A 13-column tissue matrix is wrong."
        )
    if any(wide[c].null_count() for c in dataset_cols):
        _fail("Coverage matrix contains null cells; absence must stay categorical.")

    sig_pairs = (
        associations.filter(pl.col("fdr_bh") < FDR_THRESHOLD)
        .select(["canonical_rsid", "eqtl_alt_allele"])
        .unique()
    )
    n_sig_variants = associations.filter(pl.col("fdr_bh") < FDR_THRESHOLD)[
        "canonical_rsid"
    ].n_unique()
    if n_sig_variants != EXPECTED_FDR_VARIANTS:
        _fail(
            f"Found {n_sig_variants} FDR<0.05 variants in associations; "
            f"expected {EXPECTED_FDR_VARIANTS}."
        )

    # Allele-keyed: a row "belongs" to the 16 FDR-significant variants when its
    # (rsid, alt) pair is the Catalogue allele that carries the FDR hit.
    # Sibling alts of those rsIDs are allele_not_represented and cannot have an
    # association cell. Matching on rsid alone would count 27 rows and mix
    # absence with association.
    belonging = annotation.join(
        sig_pairs,
        left_on=["canonical_rsid", "alt_allele"],
        right_on=["canonical_rsid", "eqtl_alt_allele"],
        how="inner",
    )
    n_belonging = belonging.height

    has_assoc = pl.any_horizontal([pl.col(c) == CELL_ASSOC for c in dataset_cols]).alias(
        "has_assoc"
    )
    n_rows_with_assoc = wide.select(has_assoc.alias("has_assoc")).filter(pl.col("has_assoc")).height

    if n_rows_with_assoc != n_belonging:
        _fail(
            "Rows with at least one 'association at FDR<0.05' cell "
            f"({n_rows_with_assoc}) != allele rows belonging to the "
            f"{n_sig_variants} FDR-significant variants ({n_belonging})."
        )


def build_provenance(
    *,
    associations_path: Path,
    annotation_path: Path,
    matrix_path: Path,
    figure_path: Path,
) -> str:
    """Format a provenance block (tool versions, checksums, timestamp)."""
    timestamp = datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")
    lines = [
        "provenance",
        f"timestamp_utc={timestamp}",
        f"python={sys.version.split()[0]}",
        f"polars={_package_version('polars')}",
        f"matplotlib={_package_version('matplotlib')}",
        f"numpy={_package_version('numpy')}",
        f"typer={_package_version('typer')}",
        f"source={associations_path.name} sha256={_sha256_hex(associations_path)}",
        f"source={annotation_path.name} sha256={_sha256_hex(annotation_path)}",
        f"output={matrix_path.name}",
        f"output={figure_path.name}",
        "fdr_method=Benjamini-Hochberg as stored in fdr_bh; threshold=0.05",
        "fdr_scope=within-query correction across the 1650 returned association rows; not genome-wide",
        "grouping=dataset_id (14 columns; cerebellum=QTD000161 and QTD000166)",
        "absence_encoding=categorical statements; never 0/NA/imputed",
    ]
    return "\n".join(lines)


def write_matrix_csv(wide: pl.DataFrame, path: Path, provenance: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    comment = "\n".join(f"# {line}" for line in provenance.splitlines()) + "\n"
    with path.open("w", encoding="utf-8", newline="") as handle:
        handle.write(comment)
        wide.write_csv(handle)


def render_heatmap(
    wide: pl.DataFrame,
    path: Path,
    provenance: str,
) -> None:
    """Draw a categorical heatmap; legend is discrete, not a continuous scale."""
    dataset_cols = [c for c in wide.columns if c not in INDEX_COLUMNS]
    code_of = {name: idx for idx, name in enumerate(CELL_ORDER)}
    values = wide.select(dataset_cols).to_numpy()
    codes = np.vectorize(code_of.__getitem__)(values).astype(int)

    cmap = ListedColormap([CELL_COLORS[name] for name in CELL_ORDER], name="eqtl_cov")
    cmap.set_extremes(under=CELL_COLORS[CELL_TESTED_NS], over=CELL_COLORS[CELL_VARIANT_ABSENT])
    bounds = np.arange(len(CELL_ORDER) + 1) - 0.5
    norm = BoundaryNorm(bounds, cmap.N)

    row_labels = [
        f"{gene}  {rsid}:{alt}"
        for gene, rsid, alt in wide.select(["gene", "canonical_rsid", "alt_allele"]).iter_rows()
    ]
    col_labels = [
        label.replace(" | ", "\n").replace("brain (", "").replace(")", "") for label in dataset_cols
    ]

    fig_h = max(16.0, 0.22 * len(row_labels) + 3.5)
    fig, ax = plt.subplots(figsize=(13.5, fig_h), layout="constrained")
    mesh = ax.imshow(
        codes,
        aspect="auto",
        interpolation="nearest",
        cmap=cmap,
        norm=norm,
    )
    # No colorbar: imshow+BoundaryNorm is discrete; the legend names the statements.
    mesh.set_clim(bounds[0], bounds[-1])

    ax.set_xticks(np.arange(len(col_labels)))
    ax.set_xticklabels(col_labels, fontsize=8, rotation=0, ha="center")
    ax.set_yticks(np.arange(len(row_labels)))
    ax.set_yticklabels(row_labels, fontsize=7, fontfamily="monospace")
    ax.set_xlim(-0.5, len(col_labels) - 0.5)
    ax.set_ylim(len(row_labels) - 0.5, -0.5)
    ax.set_xlabel("GTEx eQTL Catalogue dataset (grouped by dataset_id)")
    ax.set_ylabel("Allele row (gene, canonical_rsid, alt_allele)")
    ax.set_title(
        "eQTL Catalogue coverage — 73 alleles × 14 datasets\n"
        "Cerebellum is two studies (QTD000161 n=175, QTD000166 n=209). "
        "FDR<0.05 is a within-query BH correction across 1650 rows, not genome-wide."
    )
    ax.tick_params(length=0)
    ax.set_xticks(np.arange(-0.5, len(col_labels), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(row_labels), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.4)
    ax.tick_params(which="minor", bottom=False, left=False)

    legend_handles = [
        Patch(facecolor=CELL_COLORS[name], edgecolor="0.25", linewidth=0.6, label=name)
        for name in CELL_ORDER
    ]
    ax.legend(
        handles=legend_handles,
        title="Cell statement (categorical; not a numeric scale)",
        loc="upper left",
        bbox_to_anchor=(0.0, -0.04),
        ncol=2,
        frameon=True,
        fontsize=8,
        title_fontsize=8,
        borderaxespad=0.0,
    )
    fig.text(
        0.01,
        0.005,
        provenance.replace("\n", "  |  "),
        fontsize=6.5,
        family="monospace",
        va="bottom",
        wrap=True,
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def build_coverage(
    associations_path: Path,
    annotation_path: Path,
    matrix_path: Path,
    figure_path: Path,
) -> dict[str, Any]:
    """Load, classify, assert, and write matrix + figure."""
    associations = load_associations(associations_path)
    annotation = load_annotation(annotation_path)
    long_cells = classify_cells(annotation, associations)
    datasets = dataset_columns(associations)
    wide = pivot_matrix(long_cells, datasets)
    assert_matrix(wide, associations, annotation)

    provenance = build_provenance(
        associations_path=associations_path,
        annotation_path=annotation_path,
        matrix_path=matrix_path,
        figure_path=figure_path,
    )
    write_matrix_csv(wide, matrix_path, provenance)
    render_heatmap(wide, figure_path, provenance)
    return {
        "n_rows": wide.height,
        "n_dataset_columns": wide.width - len(INDEX_COLUMNS),
        "matrix_path": str(matrix_path),
        "figure_path": str(figure_path),
        "provenance": provenance,
    }


def main(
    associations: Path = typer.Option(
        DEFAULT_ASSOCIATIONS,
        "--associations",
        help="eQTL Catalogue association table (1650 rows, 14 dataset_id values).",
    ),
    annotation: Path = typer.Option(
        DEFAULT_ANNOTATION,
        "--annotation",
        help="Allele annotation layer (73 rows, eqtl_query_status).",
    ),
    matrix_out: Path = typer.Option(
        DEFAULT_MATRIX_OUT,
        "--matrix-out",
        help="Output coverage matrix CSV.",
    ),
    figure_out: Path = typer.Option(
        DEFAULT_FIGURE_OUT,
        "--figure-out",
        help="Output categorical heatmap PNG.",
    ),
) -> None:
    """Build the allele × dataset eQTL coverage matrix and heatmap."""
    result = build_coverage(associations, annotation, matrix_out, figure_out)
    typer.echo(
        f"Wrote {result['n_rows']} × {result['n_dataset_columns']} coverage matrix "
        f"to {result['matrix_path']}"
    )
    typer.echo(f"Wrote heatmap to {result['figure_path']}")
    typer.echo("")
    typer.echo(result["provenance"])


if __name__ == "__main__":
    typer.run(main)
