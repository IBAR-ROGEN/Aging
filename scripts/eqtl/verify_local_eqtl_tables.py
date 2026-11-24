#!/usr/bin/env python3
"""Verify local eQTL Catalogue tables and re-run integrated evidence off-tree.

Inventories ``results/eqtl_associations_full.csv``,
``results/annotation_layer.csv``, and ``analysis/eqtl/eqtl_datasets_used.tsv``
against the 47-variant panel. Invokes
``scripts/longevity/integrated_evidence.py`` with ``--out-dir`` under a
temporary directory so committed ``analysis/longevity/`` artefacts are not
touched. Compares evidence-label counts to the ground-truth GWAS 1 /
eQTL 16 / population 5 and stops if they differ.

Does not query the eQTL Catalogue. Does not rewrite the local tables.

Example:
    uv run python scripts/eqtl/verify_local_eqtl_tables.py
"""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
import tempfile
from datetime import UTC, datetime
from importlib import metadata
from pathlib import Path
from typing import Any

import polars as pl
import typer

from rogen_aging.config import find_repo_root

REPO_ROOT = find_repo_root()

DEFAULT_PANEL = REPO_ROOT / "data" / "processed" / "variants_47_input_v2.csv"
DEFAULT_DATASETS = REPO_ROOT / "analysis" / "eqtl" / "eqtl_datasets_used.tsv"
DEFAULT_ASSOCIATIONS = REPO_ROOT / "results" / "eqtl_associations_full.csv"
DEFAULT_ANNOTATION = REPO_ROOT / "results" / "annotation_layer.csv"
DEFAULT_INTEGRATED_SCRIPT = REPO_ROOT / "scripts" / "longevity" / "integrated_evidence.py"
DEFAULT_PROVENANCE = REPO_ROOT / "analysis" / "eqtl" / "verify_local_eqtl_tables_provenance.json"
DEFAULT_DOC = REPO_ROOT / "docs" / "EQTL_TABLES_PROVENANCE.md"

EXPECTED_PANEL_N = 47
EXPECTED_DATASET_N = 14
EXPECTED_LABELS: dict[str, int] = {
    "GWAS support": 1,
    "eQTL support": 16,
    "population difference": 5,
}

RSID_COLUMN_CANDIDATES = ("canonical_rsid", "rsid")


def _stop(message: str) -> None:
    typer.echo(message, err=True)
    raise typer.Exit(code=1)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _mtime_utc(path: Path) -> str:
    return datetime.fromtimestamp(path.stat().st_mtime, tz=UTC).strftime("%Y-%m-%dT%H:%M:%SZ")


def _mtime_date_utc(path: Path) -> str:
    return datetime.fromtimestamp(path.stat().st_mtime, tz=UTC).strftime("%Y-%m-%d")


def _git_commit() -> str:
    completed = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=REPO_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    return completed.stdout.strip()


def _package_version(name: str) -> str:
    return metadata.version(name)


def _read_table(path: Path) -> pl.DataFrame:
    if path.suffix.lower() == ".tsv":
        return pl.read_csv(path, separator="\t", infer_schema_length=0)
    return pl.read_csv(path, infer_schema_length=0)


def _rsid_column(frame: pl.DataFrame) -> str | None:
    for name in RSID_COLUMN_CANDIDATES:
        if name in frame.columns:
            return name
    return None


def _load_panel_rsids(panel_path: Path) -> list[str]:
    assert panel_path.is_file(), f"MISSING INPUT: {panel_path}"
    frame = pl.read_csv(
        panel_path,
        schema_overrides={"rsid": pl.Utf8, "alt_status": pl.Utf8},
    )
    assert frame.height == EXPECTED_PANEL_N, (
        f"PANEL ROW COUNT DISAGREES: computed={frame.height} " f"ground_truth={EXPECTED_PANEL_N}"
    )
    assert "rsid" in frame.columns, f"panel missing rsid column: {panel_path}"
    if "alt_status" in frame.columns:
        n_resolved = frame.filter(pl.col("alt_status") == "resolved").height
        assert n_resolved == EXPECTED_PANEL_N, (
            f"PANEL RESOLVED COUNT DISAGREES: computed={n_resolved} "
            f"ground_truth={EXPECTED_PANEL_N}"
        )
    rsids = [str(value) for value in frame.get_column("rsid").to_list()]
    assert len(set(rsids)) == EXPECTED_PANEL_N, (
        f"PANEL UNIQUE RSID COUNT DISAGREES: computed={len(set(rsids))} "
        f"ground_truth={EXPECTED_PANEL_N}"
    )
    return rsids


def inventory_table(path: Path, panel_rsids: list[str]) -> dict[str, Any]:
    assert path.is_file(), f"MISSING INPUT: {path}"
    frame = _read_table(path)
    rsid_col = _rsid_column(frame)
    panel_set = set(panel_rsids)
    if rsid_col is None:
        n_rows_panel = None
        n_panel_rsids_present = None
        panel_note = "no rsid column; panel restriction not applicable"
    else:
        panel_rows = frame.filter(pl.col(rsid_col).is_in(list(panel_set)))
        n_rows_panel = int(panel_rows.height)
        n_panel_rsids_present = int(panel_rows.get_column(rsid_col).n_unique())
        panel_note = f"restricted on column {rsid_col}"
    return {
        "path": str(path.relative_to(REPO_ROOT)),
        "sha256": _sha256(path),
        "n_rows": int(frame.height),
        "columns": list(frame.columns),
        "mtime_utc": _mtime_utc(path),
        "mtime_date_utc": _mtime_date_utc(path),
        "rsid_column": rsid_col,
        "n_rows_panel_variants": n_rows_panel,
        "n_panel_rsids_present": n_panel_rsids_present,
        "panel_restriction_note": panel_note,
    }


def run_integrated_evidence(
    *,
    associations: Path,
    annotation: Path,
    datasets: Path,
    panel_v2: Path,
    integrated_script: Path,
    out_dir: Path,
) -> dict[str, Any]:
    assert integrated_script.is_file(), f"MISSING INPUT: {integrated_script}"
    assert associations.is_file(), f"MISSING INPUT: {associations}"
    assert annotation.is_file(), f"MISSING INPUT: {annotation}"
    assert datasets.is_file(), f"MISSING INPUT: {datasets}"
    committed_out = REPO_ROOT / "analysis" / "longevity"
    assert (
        out_dir.resolve() != committed_out.resolve()
    ), "REFUSING TO WRITE: out_dir must not be the committed analysis/longevity path"
    argv = [
        "uv",
        "run",
        "python",
        str(integrated_script),
        "--associations",
        str(associations),
        "--annotation",
        str(annotation),
        "--datasets",
        str(datasets),
        "--out-dir",
        str(out_dir),
    ]
    env = os.environ.copy()
    env.setdefault("MPLCONFIGDIR", str(out_dir / ".mplconfig"))
    completed = subprocess.run(
        argv,
        cwd=REPO_ROOT,
        check=False,
        capture_output=True,
        text=True,
        env=env,
    )
    assert completed.returncode == 0, (
        "integrated_evidence.py FAILED\n"
        f"argv={argv}\n"
        f"stdout:\n{completed.stdout}\n"
        f"stderr:\n{completed.stderr}"
    )
    table_path = out_dir / "integrated_evidence.csv"
    assert table_path.is_file(), f"MISSING OUTPUT: {table_path}"
    table = pl.read_csv(
        table_path,
        schema_overrides={
            "rsid": pl.Utf8,
            "label_gwas_support": pl.Boolean,
            "label_eqtl_support": pl.Utf8,
            "label_population_difference": pl.Boolean,
        },
    )
    n_variants = int(table.height)
    assert n_variants == EXPECTED_PANEL_N, (
        f"INTEGRATED VARIANT COUNT DISAGREES: computed={n_variants} "
        f"ground_truth={EXPECTED_PANEL_N}"
    )
    # label_eqtl_support is bool when present; may be Utf8 "true"/"false"/"not_available"
    if table.schema["label_eqtl_support"] == pl.Boolean:
        n_eqtl = int(table.filter(pl.col("label_eqtl_support") == True).height)  # noqa: E712
    else:
        n_eqtl = int(table.filter(pl.col("label_eqtl_support") == "true").height)
    observed_labels = {
        "GWAS support": int(table.filter(pl.col("label_gwas_support")).height),
        "eQTL support": n_eqtl,
        "population difference": int(table.filter(pl.col("label_population_difference")).height),
    }
    return {
        "argv": argv,
        "out_dir": str(out_dir),
        "table_path": str(table_path),
        "n_variants": n_variants,
        "observed_labels": observed_labels,
        "stdout": completed.stdout,
        "stderr": completed.stderr,
        "panel_v2_used_for_inventory": str(panel_v2.relative_to(REPO_ROOT)),
    }


def compare_labels(observed: dict[str, int], expected: dict[str, int]) -> None:
    typer.echo("Evidence labels (expected vs observed):")
    for name in expected:
        typer.echo(f"  {name}: expected={expected[name]} observed={observed[name]}")
    if observed != expected:
        _stop(
            "LABEL COUNT DISAGREES.\n"
            f"  expected={expected}\n"
            f"  observed={observed}\n"
            "Stopping without adjusting counts or inputs."
        )


def render_doc(payload: dict[str, Any]) -> str:
    tables = payload["tables"]
    datasets = tables["eqtl_datasets_used"]
    assoc = tables["eqtl_associations_full"]
    annotation = tables["annotation_layer"]
    panel = tables["variants_47_input_v2"]
    labels = payload["integrated_evidence_temp"]["observed_labels"]
    dataset_ids = payload["dataset_ids"]
    lines = [
        "# Local eQTL Catalogue tables — provenance",
        "",
        "Verification only. These tables were built outside this repo",
        "(see [INTEGRATIVE_PIPELINE.md](INTEGRATIVE_PIPELINE.md): there is no",
        "in-repo writer for `results/eqtl_associations_full.csv` /",
        "`results/annotation_layer.csv`). This document records what is known",
        "about the local copies and the result of running",
        "`scripts/eqtl/verify_local_eqtl_tables.py`.",
        "",
        "## What is known",
        "",
        "- Source: tables were built outside this repo (see",
        "  `docs/INTEGRATIVE_PIPELINE.md`).",
        "- File mtime on `results/eqtl_associations_full.csv` and",
        "  `results/annotation_layer.csv` is **2026-08-30** and is **not",
        "  authoritative** (filesystem timestamp only).",
        f"- Locked dataset list: **{len(dataset_ids)}** `dataset_id` values in",
        "  `analysis/eqtl/eqtl_datasets_used.tsv`.",
        "- Catalogue release version: **unknown, not recoverable from the repo**.",
        "",
        "## Dataset IDs",
        "",
    ]
    for dataset_id in dataset_ids:
        lines.append(f"- `{dataset_id}`")
    lines.extend(
        [
            "",
            "## Inventory (from verification)",
            "",
            "| Path | sha256 | n_rows | mtime (UTC) | rows on 47 panel rsIDs |"
            " panel rsIDs present |",
            "|------|--------|--------|-------------|------------------------|"
            "---------------------|",
        ]
    )
    for key in (
        "variants_47_input_v2",
        "eqtl_datasets_used",
        "eqtl_associations_full",
        "annotation_layer",
    ):
        row = tables[key]
        n_panel = row["n_rows_panel_variants"]
        n_rsids = row["n_panel_rsids_present"]
        lines.append(
            f"| `{row['path']}` | `{row['sha256']}` | {row['n_rows']} |"
            f" {row['mtime_utc']} |"
            f" {n_panel if n_panel is not None else 'n/a'} |"
            f" {n_rsids if n_rsids is not None else 'n/a'} |"
        )
    lines.extend(
        [
            "",
            "### Columns",
            "",
            f"- `{panel['path']}`: {', '.join(f'`{c}`' for c in panel['columns'])}",
            f"- `{datasets['path']}`: " + ", ".join(f"`{c}`" for c in datasets["columns"]),
            f"- `{assoc['path']}`: " + ", ".join(f"`{c}`" for c in assoc["columns"]),
            f"- `{annotation['path']}`: " + ", ".join(f"`{c}`" for c in annotation["columns"]),
            "",
            "## Integrated evidence check (temporary output)",
            "",
            "`scripts/longevity/integrated_evidence.py` was invoked with the local",
            "associations, annotation, and dataset list, writing under a temporary",
            "directory (committed `analysis/longevity/` outputs were not written).",
            "",
            f"- Variants in temp table: **{payload['integrated_evidence_temp']['n_variants']}**",
            f"- Labels: GWAS **{labels['GWAS support']}** /"
            f" eQTL **{labels['eQTL support']}** /"
            f" population **{labels['population difference']}**",
            f"- Ground truth compared: GWAS **{EXPECTED_LABELS['GWAS support']}** /"
            f" eQTL **{EXPECTED_LABELS['eQTL support']}** /"
            f" population **{EXPECTED_LABELS['population difference']}**",
            "",
            "## Verification provenance",
            "",
            f"- Git commit: `{payload['git_commit']}`",
            f"- Timestamp (UTC): `{payload['generated_at_utc']}`",
            f"- Python: `{payload['tool_versions']['python']}`",
            f"- Polars: `{payload['tool_versions']['polars']}`",
            f"- Typer: `{payload['tool_versions']['typer']}`",
            f"- Machine-readable sidecar: `{payload['provenance_path']}`",
            "",
        ]
    )
    return "\n".join(lines)


def main(
    panel: Path = typer.Option(DEFAULT_PANEL, "--panel"),
    datasets: Path = typer.Option(DEFAULT_DATASETS, "--datasets"),
    associations: Path = typer.Option(DEFAULT_ASSOCIATIONS, "--associations"),
    annotation: Path = typer.Option(DEFAULT_ANNOTATION, "--annotation"),
    integrated_script: Path = typer.Option(DEFAULT_INTEGRATED_SCRIPT, "--integrated-script"),
    provenance_out: Path = typer.Option(DEFAULT_PROVENANCE, "--provenance-out"),
    doc_out: Path = typer.Option(DEFAULT_DOC, "--doc-out"),
) -> None:
    """Inventory local eQTL tables and check integrated-evidence labels."""
    started = datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")
    panel_rsids = _load_panel_rsids(panel)

    panel_info = inventory_table(panel, panel_rsids)
    datasets_info = inventory_table(datasets, panel_rsids)
    assoc_info = inventory_table(associations, panel_rsids)
    annotation_info = inventory_table(annotation, panel_rsids)

    assert datasets_info["n_rows"] == EXPECTED_DATASET_N, (
        f"DATASET COUNT DISAGREES: computed={datasets_info['n_rows']} "
        f"ground_truth={EXPECTED_DATASET_N}"
    )
    dataset_frame = _read_table(datasets)
    assert "dataset_id" in dataset_frame.columns, "datasets file missing dataset_id"
    dataset_ids = [str(value) for value in dataset_frame.get_column("dataset_id").to_list()]
    assert len(dataset_ids) == EXPECTED_DATASET_N
    assert len(set(dataset_ids)) == EXPECTED_DATASET_N, f"DATASET IDS NOT UNIQUE: {dataset_ids}"

    for info in (panel_info, datasets_info, assoc_info, annotation_info):
        typer.echo(
            f"{info['path']}: n_rows={info['n_rows']} "
            f"n_rows_panel={info['n_rows_panel_variants']} "
            f"mtime={info['mtime_utc']} sha256={info['sha256']}"
        )
        typer.echo(f"  columns={info['columns']}")

    with tempfile.TemporaryDirectory(prefix="eqtl_verify_") as tmp:
        out_dir = Path(tmp) / "integrated"
        out_dir.mkdir(parents=True, exist_ok=True)
        integrated = run_integrated_evidence(
            associations=associations,
            annotation=annotation,
            datasets=datasets,
            panel_v2=panel,
            integrated_script=integrated_script,
            out_dir=out_dir,
        )
        compare_labels(integrated["observed_labels"], EXPECTED_LABELS)
        # Snapshot stdout before the temp dir is removed.
        integrated_record = {
            "argv": integrated["argv"],
            "n_variants": integrated["n_variants"],
            "observed_labels": integrated["observed_labels"],
            "expected_labels": dict(EXPECTED_LABELS),
            "labels_match": True,
            "wrote_committed_analysis_longevity": False,
            "panel_v2_used_for_inventory": integrated["panel_v2_used_for_inventory"],
            "stdout_tail": integrated["stdout"][-2000:],
        }

    payload: dict[str, Any] = {
        "generated_at_utc": started,
        "git_commit": _git_commit(),
        "tool_versions": {
            "python": sys.version.split()[0],
            "polars": _package_version("polars"),
            "typer": _package_version("typer"),
        },
        "ground_truth_labels": dict(EXPECTED_LABELS),
        "dataset_ids": dataset_ids,
        "catalogue_release_version": "unknown, not recoverable from the repo",
        "notes": [
            "Tables were built outside this repo (see docs/INTEGRATIVE_PIPELINE.md).",
            "File mtime 2026-08-30 on associations/annotation is not authoritative.",
            "Verification only; local tables were not rewritten.",
        ],
        "tables": {
            "variants_47_input_v2": panel_info,
            "eqtl_datasets_used": datasets_info,
            "eqtl_associations_full": assoc_info,
            "annotation_layer": annotation_info,
        },
        "integrated_evidence_temp": integrated_record,
        "provenance_path": str(provenance_out.relative_to(REPO_ROOT)),
        "doc_path": str(doc_out.relative_to(REPO_ROOT)),
    }

    provenance_out.parent.mkdir(parents=True, exist_ok=True)
    provenance_out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    doc_out.parent.mkdir(parents=True, exist_ok=True)
    doc_out.write_text(render_doc(payload), encoding="utf-8")
    typer.echo(f"Wrote {provenance_out}")
    typer.echo(f"Wrote {doc_out}")
    typer.echo(
        "OK: labels match ground truth "
        f"(GWAS {EXPECTED_LABELS['GWAS support']} / "
        f"eQTL {EXPECTED_LABELS['eQTL support']} / "
        f"population {EXPECTED_LABELS['population difference']})"
    )


if __name__ == "__main__":
    typer.run(main)
