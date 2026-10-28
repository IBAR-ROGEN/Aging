#!/usr/bin/env python3
"""One evidence table per prioritized variant (GWAS × eQTL × panel frequencies).

Join keys are ``rsid`` + panel alt. Evidence labels are defined in
``docs/INTEGRATIVE_PIPELINE.md`` (not only here):

- ``GWAS support``: BH q < 0.05 in at least one longevity/lifespan trait
- ``eQTL support``: FDR < 0.05 in at least one eQTL Catalogue dataset for the
  panel alt
- ``population difference``: TSI vs NFE BH q < 0.05

No label implies causality.

Example:
    uv run python scripts/longevity/integrated_evidence.py
"""

from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import polars as pl  # noqa: E402
import typer  # noqa: E402

from rogen_aging.config import find_repo_root  # noqa: E402
from rogen_aging.ukb.gnomad import normalize_rsid  # noqa: E402

REPO_ROOT = find_repo_root()
DEFAULT_GWAS = REPO_ROOT / "analysis" / "longevity" / "longevity_gwas_lookup.csv"
DEFAULT_ASSOCIATIONS = REPO_ROOT / "results" / "eqtl_associations_full.csv"
DEFAULT_DATASETS = REPO_ROOT / "analysis" / "eqtl" / "eqtl_datasets_used.tsv"
DEFAULT_PANEL = REPO_ROOT / "analysis" / "panel" / "ro_longevity_panel.csv"
DEFAULT_POWER = REPO_ROOT / "analysis" / "panel" / "ro_longevity_power_grid.csv"
DEFAULT_OUT_DIR = REPO_ROOT / "analysis" / "longevity"

LABEL_GWAS = "GWAS support"
LABEL_EQTL = "eQTL support"
LABEL_POP = "population difference"
EVIDENCE_LABELS: tuple[str, ...] = (LABEL_GWAS, LABEL_EQTL, LABEL_POP)

POWER_RANK_N = 1000
POWER_RANK_OR = 1.5
POWER_RANK_ALPHA = "bonferroni_0.05_over_47"
POWER_HIGH = 0.8
POWER_MODERATE = 0.5
FDR_THRESHOLD = 0.05
RELEASE_CHECK_KIND = "summary"

DISCLAIMER = (
    "Integrated public-evidence labels for the prioritized longevity panel. "
    "Labels summarize GWAS, eQTL, and reference-population frequency tests. "
    "No label implies causality, Romanian enrichment, or a ROGEN association."
)


def _stop(message: str) -> None:
    typer.echo(message, err=True)
    raise typer.Exit(code=1)


def _repo_relative(path: Path) -> str:
    resolved = path.resolve()
    root = REPO_ROOT.resolve()
    try:
        return resolved.relative_to(root).as_posix()
    except ValueError:
        return str(path)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _require_file(path: Path, role: str) -> None:
    if not path.is_file():
        _stop(f"STOP: missing {role}: {path}")


def load_panel(path: Path) -> pl.DataFrame:
    """Panel frequencies and TSI vs NFE test (one row per variant)."""
    _require_file(path, "panel table")
    frame = pl.read_csv(path, schema_overrides={"rsid": pl.Utf8, "alt": pl.Utf8})
    required = {
        "rsid",
        "alt",
        "gene_symbol",
        "tsi_differs_fdr05",
        "q_tsi_vs_nfe",
        "testability_rank",
        "power_n1000_or1_5_alpha_bonferroni",
    }
    missing = sorted(required - set(frame.columns))
    if missing:
        _stop(f"STOP: panel missing columns: {missing}")
    return frame.with_columns(
        pl.col("rsid").map_elements(normalize_rsid, return_dtype=pl.Utf8),
        pl.col("alt").cast(pl.Utf8).str.to_uppercase().alias("panel_alt"),
    )


def load_gwas(path: Path) -> pl.DataFrame:
    """Variant × trait GWAS lookup (Prompt 4 output)."""
    _require_file(
        path,
        "GWAS lookup CSV (run scripts/longevity/longevity_gwas_lookup.py first)",
    )
    frame = pl.read_csv(
        path,
        schema_overrides={
            "rsid": pl.Utf8,
            "panel_alt": pl.Utf8,
            "trait": pl.Utf8,
            "direction": pl.Utf8,
            "status": pl.Utf8,
        },
    )
    required = {"rsid", "panel_alt", "trait", "direction", "q_bh", "status"}
    missing = sorted(required - set(frame.columns))
    if missing:
        _stop(f"STOP: GWAS lookup missing columns: {missing}")
    return frame.with_columns(
        pl.col("rsid").map_elements(normalize_rsid, return_dtype=pl.Utf8),
        pl.col("panel_alt").cast(pl.Utf8).str.to_uppercase(),
    )


def load_eqtl_associations(path: Path, datasets_path: Path) -> tuple[pl.DataFrame, list[str]]:
    """September eQTL Catalogue associations and the locked dataset_id list."""
    _require_file(path, "eQTL associations table")
    _require_file(datasets_path, "eQTL datasets_used TSV")
    datasets = (
        pl.read_csv(datasets_path, separator="\t", schema_overrides={"dataset_id": pl.Utf8})
        .get_column("dataset_id")
        .to_list()
    )
    if not datasets:
        _stop(f"STOP: {datasets_path} has no dataset_id rows.")
    frame = pl.read_csv(
        path,
        schema_overrides={
            "canonical_rsid": pl.Utf8,
            "eqtl_alt_allele": pl.Utf8,
            "dataset_id": pl.Utf8,
        },
    )
    required = {"canonical_rsid", "eqtl_alt_allele", "dataset_id", "beta", "fdr_bh"}
    missing = sorted(required - set(frame.columns))
    if missing:
        _stop(f"STOP: eQTL associations missing columns: {missing}")
    frame = frame.filter(pl.col("dataset_id").is_in(datasets)).with_columns(
        pl.col("canonical_rsid").map_elements(normalize_rsid, return_dtype=pl.Utf8),
        pl.col("eqtl_alt_allele").cast(pl.Utf8).str.to_uppercase(),
    )
    return frame, datasets


def load_power_grid(path: Path) -> pl.DataFrame:
    """Illustrative power grid; used for the rank-cell power tier."""
    _require_file(path, "power grid")
    frame = pl.read_csv(path, schema_overrides={"rsid": pl.Utf8})
    required = {"rsid", "n_cases_per_arm", "allelic_or", "alpha_name", "power"}
    missing = sorted(required - set(frame.columns))
    if missing:
        _stop(f"STOP: power grid missing columns: {missing}")
    return frame.with_columns(pl.col("rsid").map_elements(normalize_rsid, return_dtype=pl.Utf8))


def power_tier(power: float | None) -> str:
    """Tier at the illustrative rank cell (1000 cases, OR 1.5, Bonferroni)."""
    if power is None:
        return "unavailable"
    if power >= POWER_HIGH:
        return "high"
    if power >= POWER_MODERATE:
        return "moderate"
    return "low"


def summarize_gwas(gwas: pl.DataFrame) -> pl.DataFrame:
    """One row per rsid+panel_alt with per-trait direction/q and GWAS label."""
    traits = sorted(gwas.get_column("trait").unique().to_list())
    rows: list[dict[str, Any]] = []
    keys = gwas.select(["rsid", "panel_alt"]).unique()
    for key in keys.iter_rows(named=True):
        subset = gwas.filter(
            (pl.col("rsid") == key["rsid"]) & (pl.col("panel_alt") == key["panel_alt"])
        )
        out: dict[str, Any] = {
            "rsid": key["rsid"],
            "panel_alt": key["panel_alt"],
            "gwas_n_traits_q_lt_0_05": 0,
            "gwas_traits_q_lt_0_05": "",
        }
        n_sig = 0
        sig_traits: list[str] = []
        for trait in traits:
            trait_rows = subset.filter(pl.col("trait") == trait)
            direction = ""
            q_value: float | None = None
            if trait_rows.height:
                row = trait_rows.row(0, named=True)
                direction = str(row.get("direction") or "")
                raw_q = row.get("q_bh")
                if raw_q is not None and row.get("status") == "present":
                    q_value = float(raw_q)
                    if q_value < FDR_THRESHOLD:
                        n_sig += 1
                        sig_traits.append(trait)
            safe = trait.replace(" ", "_")
            out[f"gwas_direction_{safe}"] = direction
            out[f"gwas_q_{safe}"] = q_value
        out["gwas_n_traits_q_lt_0_05"] = n_sig
        out["gwas_traits_q_lt_0_05"] = ";".join(sig_traits)
        rows.append(out)
    if not rows:
        return pl.DataFrame(
            schema={
                "rsid": pl.Utf8,
                "panel_alt": pl.Utf8,
                "gwas_n_traits_q_lt_0_05": pl.Int64,
                "gwas_traits_q_lt_0_05": pl.Utf8,
            }
        )
    return pl.DataFrame(rows)


def summarize_eqtl(assoc: pl.DataFrame, panel: pl.DataFrame) -> pl.DataFrame:
    """Count significant datasets and alt-effect direction consistency per panel alt."""
    rows: list[dict[str, Any]] = []
    for panel_row in panel.iter_rows(named=True):
        rsid = panel_row["rsid"]
        alt = panel_row["panel_alt"]
        hits = assoc.filter(
            (pl.col("canonical_rsid") == rsid)
            & (pl.col("eqtl_alt_allele") == alt)
            & (pl.col("fdr_bh").is_not_null())
            & (pl.col("fdr_bh") < FDR_THRESHOLD)
        )
        # One sign per dataset (mean beta across eGenes in that dataset).
        if hits.height:
            by_dataset = hits.group_by("dataset_id").agg(pl.col("beta").mean().alias("beta_mean"))
            n_sig = by_dataset.height
            signs = [
                1 if float(beta) > 0 else (-1 if float(beta) < 0 else 0)
                for beta in by_dataset.get_column("beta_mean").to_list()
            ]
            nonzero = [sign for sign in signs if sign != 0]
            if not nonzero:
                direction = "none"
            elif len(set(nonzero)) == 1:
                direction = "consistent"
            else:
                direction = "mixed"
            datasets = ";".join(sorted(by_dataset.get_column("dataset_id").to_list()))
        else:
            n_sig = 0
            direction = "none"
            datasets = ""
        rows.append(
            {
                "rsid": rsid,
                "panel_alt": alt,
                "eqtl_n_datasets_sig": n_sig,
                "eqtl_alt_effect_direction": direction,
                "eqtl_datasets_sig": datasets,
            }
        )
    return pl.DataFrame(rows)


def summarize_power(grid: pl.DataFrame, panel: pl.DataFrame) -> pl.DataFrame:
    """Power tier from the illustrative rank cell; also keep rank-cell power."""
    rank_cell = grid.filter(
        (pl.col("n_cases_per_arm") == POWER_RANK_N)
        & (pl.col("allelic_or") == POWER_RANK_OR)
        & (pl.col("alpha_name") == POWER_RANK_ALPHA)
    ).select(["rsid", pl.col("power").alias("power_rank_cell")])
    joined = panel.select(["rsid", "panel_alt", "power_n1000_or1_5_alpha_bonferroni"]).join(
        rank_cell, on="rsid", how="left"
    )
    rows: list[dict[str, Any]] = []
    for row in joined.iter_rows(named=True):
        power = row.get("power_rank_cell")
        if power is None:
            power = row.get("power_n1000_or1_5_alpha_bonferroni")
        power_f = None if power is None else float(power)
        rows.append(
            {
                "rsid": row["rsid"],
                "panel_alt": row["panel_alt"],
                "power_rank_cell": power_f,
                "power_tier": power_tier(power_f),
            }
        )
    return pl.DataFrame(rows)


def assign_labels(frame: pl.DataFrame) -> pl.DataFrame:
    """Attach evidence labels and counts. Definitions live in the docs too."""
    labels_col: list[str] = []
    n_labels: list[int] = []
    has_gwas: list[bool] = []
    has_eqtl: list[bool] = []
    has_pop: list[bool] = []
    for row in frame.iter_rows(named=True):
        labels: list[str] = []
        gwas = int(row.get("gwas_n_traits_q_lt_0_05") or 0) > 0
        eqtl = int(row.get("eqtl_n_datasets_sig") or 0) > 0
        pop_flag = row.get("tsi_differs_fdr05")
        pop = bool(pop_flag) if pop_flag is not None else False
        if gwas:
            labels.append(LABEL_GWAS)
        if eqtl:
            labels.append(LABEL_EQTL)
        if pop:
            labels.append(LABEL_POP)
        labels_col.append("|".join(labels))
        n_labels.append(len(labels))
        has_gwas.append(gwas)
        has_eqtl.append(eqtl)
        has_pop.append(pop)
    return frame.with_columns(
        pl.Series("evidence_labels", labels_col),
        pl.Series("n_evidence_labels", n_labels),
        pl.Series("label_gwas_support", has_gwas),
        pl.Series("label_eqtl_support", has_eqtl),
        pl.Series("label_population_difference", has_pop),
    )


def build_table(
    *,
    panel: pl.DataFrame,
    gwas: pl.DataFrame,
    assoc: pl.DataFrame,
    grid: pl.DataFrame,
) -> pl.DataFrame:
    """Join all layers on rsid + panel alt."""
    gwas_sum = summarize_gwas(gwas)
    eqtl_sum = summarize_eqtl(assoc, panel)
    power_sum = summarize_power(grid, panel)
    base = panel.select(
        [
            "rsid",
            "panel_alt",
            pl.col("gene_symbol").alias("gene"),
            "tsi_differs_fdr05",
            "q_tsi_vs_nfe",
            "p_tsi_vs_nfe",
            "testability_rank",
            "af_alt_gnomad_nfe",
            "af_alt_1kg_tsi",
        ]
    )
    joined = (
        base.join(gwas_sum, on=["rsid", "panel_alt"], how="left")
        .join(eqtl_sum, on=["rsid", "panel_alt"], how="left")
        .join(power_sum, on=["rsid", "panel_alt"], how="left")
    )
    joined = joined.with_columns(
        pl.col("gwas_n_traits_q_lt_0_05").fill_null(0),
        pl.col("eqtl_n_datasets_sig").fill_null(0),
        pl.col("eqtl_alt_effect_direction").fill_null("none"),
        pl.col("gwas_traits_q_lt_0_05").fill_null(""),
        pl.col("eqtl_datasets_sig").fill_null(""),
    )
    labeled = assign_labels(joined)
    return labeled.sort(
        ["n_evidence_labels", "testability_rank", "rsid"],
        descending=[True, False, False],
    )


def render_dot_matrix(frame: pl.DataFrame, path: Path) -> None:
    """Variants × evidence columns, sorted by number of labels."""
    ordered = frame.sort(
        ["n_evidence_labels", "testability_rank", "rsid"],
        descending=[True, False, False],
    )
    labels = list(EVIDENCE_LABELS)
    n_var = ordered.height
    matrix = np.zeros((n_var, len(labels)), dtype=float)
    for row_index, row in enumerate(ordered.iter_rows(named=True)):
        matrix[row_index, 0] = 1.0 if row["label_gwas_support"] else 0.0
        matrix[row_index, 1] = 1.0 if row["label_eqtl_support"] else 0.0
        matrix[row_index, 2] = 1.0 if row["label_population_difference"] else 0.0
    y_labels = [f"{row['gene']} {row['rsid']}" for row in ordered.iter_rows(named=True)]
    fig_h = max(8.0, 0.28 * n_var + 2.0)
    fig, axis = plt.subplots(figsize=(6.5, fig_h), layout="constrained")
    axis.imshow(matrix, aspect="auto", cmap="Greys", vmin=0.0, vmax=1.0, interpolation="nearest")
    axis.set_xticks(np.arange(len(labels)))
    axis.set_xticklabels(labels, rotation=0, ha="center", fontsize=9)
    axis.set_yticks(np.arange(n_var))
    axis.set_yticklabels(y_labels, fontsize=7)
    axis.set_xlabel("Evidence label")
    axis.set_ylabel("Variant (sorted by number of labels)")
    axis.set_title("Integrated public evidence (dot = label present)")
    for row_index in range(n_var):
        for col_index in range(len(labels)):
            if matrix[row_index, col_index] > 0.5:
                axis.plot(col_index, row_index, "o", color="#1b4f72", markersize=5)
    axis.set_xticks(np.arange(-0.5, len(labels), 1), minor=True)
    axis.set_yticks(np.arange(-0.5, n_var, 1), minor=True)
    axis.grid(which="minor", color="0.85", linewidth=0.4)
    axis.tick_params(which="minor", bottom=False, left=False)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def release_check(path: Path, report: Path) -> dict[str, Any]:
    """Run rogen-release-check summary on the integrated table."""
    argv = [
        "uv",
        "run",
        "rogen-release-check",
        "--table",
        _repo_relative(path),
        "--table-kind",
        RELEASE_CHECK_KIND,
        "--report",
        _repo_relative(report),
    ]
    if shutil.which("uv") is None:
        return {
            "path": _repo_relative(path),
            "kind": RELEASE_CHECK_KIND,
            "report": _repo_relative(report),
            "argv": argv,
            "returncode": None,
            "status": "not_run",
            "stderr": "uv is not on PATH",
        }
    completed = subprocess.run(argv, cwd=REPO_ROOT, check=False, capture_output=True, text=True)
    return {
        "path": _repo_relative(path),
        "kind": RELEASE_CHECK_KIND,
        "report": _repo_relative(report),
        "argv": argv,
        "returncode": completed.returncode,
        "stdout": completed.stdout,
        "stderr": completed.stderr,
        "status": "passed" if completed.returncode == 0 else "failed",
    }


def main(
    gwas: Path = typer.Option(DEFAULT_GWAS, "--gwas"),
    associations: Path = typer.Option(DEFAULT_ASSOCIATIONS, "--associations"),
    datasets: Path = typer.Option(DEFAULT_DATASETS, "--datasets"),
    panel: Path = typer.Option(DEFAULT_PANEL, "--panel"),
    power_grid: Path = typer.Option(DEFAULT_POWER, "--power-grid"),
    out_dir: Path = typer.Option(DEFAULT_OUT_DIR, "--out-dir"),
) -> None:
    """Build the integrated evidence table, figure, and provenance."""
    started = datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")
    panel_frame = load_panel(panel)
    gwas_frame = load_gwas(gwas)
    assoc_frame, dataset_ids = load_eqtl_associations(associations, datasets)
    grid_frame = load_power_grid(power_grid)
    table = build_table(
        panel=panel_frame,
        gwas=gwas_frame,
        assoc=assoc_frame,
        grid=grid_frame,
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    table_path = out_dir / "integrated_evidence.csv"
    figure_path = out_dir / "integrated_evidence_dot_matrix.png"
    provenance_path = out_dir / "integrated_evidence_provenance.json"
    report_path = out_dir / "release_checks" / "integrated_evidence.json"
    table.write_csv(table_path)
    render_dot_matrix(table, figure_path)
    release = release_check(table_path, report_path)
    label_counts = {
        LABEL_GWAS: int(table.filter(pl.col("label_gwas_support")).height),
        LABEL_EQTL: int(table.filter(pl.col("label_eqtl_support")).height),
        LABEL_POP: int(table.filter(pl.col("label_population_difference")).height),
    }
    provenance = {
        "generated_at_utc": started,
        "disclaimer": DISCLAIMER,
        "join_keys": ["rsid", "panel_alt"],
        "evidence_label_definitions": {
            LABEL_GWAS: "BH q < 0.05 in at least one longevity/lifespan GWAS trait",
            LABEL_EQTL: (
                "FDR < 0.05 in at least one eQTL Catalogue dataset for the panel alt "
                f"(datasets from {_repo_relative(datasets)})"
            ),
            LABEL_POP: "TSI vs NFE BH q < 0.05 on the panel alt (public reference panels)",
        },
        "power_tier": {
            "cell": {
                "n_cases_per_arm": POWER_RANK_N,
                "allelic_or": POWER_RANK_OR,
                "alpha_name": POWER_RANK_ALPHA,
            },
            "high": f">= {POWER_HIGH}",
            "moderate": f">= {POWER_MODERATE} and < {POWER_HIGH}",
            "low": f"< {POWER_MODERATE}",
            "unavailable": "missing power",
        },
        "inputs": {
            "gwas_lookup": {
                "path": _repo_relative(gwas),
                "sha256": _sha256(gwas),
            },
            "eqtl_associations": {
                "path": _repo_relative(associations),
                "sha256": _sha256(associations),
            },
            "eqtl_datasets_used": {
                "path": _repo_relative(datasets),
                "dataset_ids": dataset_ids,
                "n_datasets": len(dataset_ids),
            },
            "panel": {
                "path": _repo_relative(panel),
                "sha256": _sha256(panel),
                "n_rows": panel_frame.height,
            },
            "power_grid": {
                "path": _repo_relative(power_grid),
                "sha256": _sha256(power_grid),
            },
        },
        "outputs": {
            "table": _repo_relative(table_path),
            "figure": _repo_relative(figure_path),
            "provenance": _repo_relative(provenance_path),
        },
        "n_variants": table.height,
        "label_counts": label_counts,
        "n_with_any_label": int(table.filter(pl.col("n_evidence_labels") > 0).height),
        "release_check": release,
    }
    provenance_path.write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    typer.echo(f"Wrote {table_path}")
    typer.echo(f"Wrote {figure_path}")
    typer.echo(f"Wrote {provenance_path}")
    typer.echo(
        "labels: "
        + ", ".join(f"{name}={count}" for name, count in label_counts.items())
        + f"; any={provenance['n_with_any_label']}/{table.height}"
    )
    if release["status"] != "passed":
        _stop("STOP: release check did not pass.")


if __name__ == "__main__":
    typer.run(main)
