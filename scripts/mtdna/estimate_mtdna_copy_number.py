#!/usr/bin/env python3
"""Estimate relative mtDNA copy number from BAM depth via ``samtools coverage``.

Formula (diploid nuclear genome):

    estimate = 2 * (mean_chrM_depth / mean_depth_chr1_to_chr22)

Autosomal mean depth is the length-weighted mean of ``meandepth`` over
``chr1``…``chr22`` from ``samtools coverage``. Contig names must be GRCh38-style
(``chrM``, not ``MT``).

Every output header states that blood-derived estimates are confounded by
platelets and cell composition and are relative, not absolute.

Example:

    uv run python scripts/mtdna/estimate_mtdna_copy_number.py
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
import sys
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import polars as pl
import typer

from rogen_aging.config import find_repo_root

REPO_ROOT = find_repo_root()
DEFAULT_SAMPLES = REPO_ROOT / "data" / "external" / "mtdna_validation_samples.csv"
DEFAULT_OUT_CSV = REPO_ROOT / "analysis" / "mtdna" / "mtdna_copy_number_estimates.csv"
DEFAULT_OUT_PROVENANCE = (
    REPO_ROOT / "analysis" / "mtdna" / "mtdna_copy_number_estimates_provenance.json"
)

WARNING = (
    "Blood-derived estimates are confounded by platelets and cell composition; "
    "values are relative, not absolute."
)

REQUIRED_SAMPLE_COLUMNS = ("sample_id", "bam_path")
AUTOSOMES = tuple(f"chr{i}" for i in range(1, 23))
MT_CONTIG = "chrM"
REQUIRED_CONTIGS = AUTOSOMES + (MT_CONTIG,)

app = typer.Typer(add_completion=False, no_args_is_help=False)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _git_commit() -> str:
    completed = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=REPO_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    return completed.stdout.strip()


def _candidate_bin_dirs() -> list[Path]:
    dirs: list[Path] = []
    which_root = shutil.which("samtools")
    if which_root:
        dirs.append(Path(which_root).resolve().parent)
    micro = Path.home() / "micromamba" / "envs" / "rogen-tools" / "bin"
    if micro.is_dir():
        dirs.append(micro)
    tools = REPO_ROOT / ".tools"
    if tools.is_dir():
        dirs.append(tools)
        if (tools / "bin").is_dir():
            dirs.append(tools / "bin")
    path_dirs = [Path(part) for part in os.environ.get("PATH", "").split(os.pathsep) if part]
    return dirs + path_dirs


def find_samtools() -> str:
    """Resolve ``samtools`` on PATH, micromamba ``rogen-tools``, or ``.tools``."""
    for directory in _candidate_bin_dirs():
        candidate = directory / "samtools"
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return str(candidate)
    found = shutil.which("samtools")
    if found:
        return found
    raise AssertionError(
        "MISSING TOOL: samtools. Install into micromamba env rogen-tools "
        "(bioconda samtools=1.24) or place a binary under .tools/."
    )


def samtools_version(samtools: str) -> str:
    """Return the first line of ``samtools --version``."""
    completed = subprocess.run(
        [samtools, "--version"],
        check=True,
        capture_output=True,
        text=True,
    )
    line = completed.stdout.strip().splitlines()[0]
    assert line, "samtools --version produced empty stdout"
    return line


def resolve_bam_path(raw: str, samples_csv: Path) -> Path:
    """Resolve a BAM path relative to the samples CSV directory, then repo root."""
    path = Path(raw)
    if path.is_file():
        return path.resolve()
    beside = (samples_csv.parent / path).resolve()
    if beside.is_file():
        return beside
    under_repo = (REPO_ROOT / path).resolve()
    if under_repo.is_file():
        return under_repo
    raise AssertionError(f"MISSING INPUT: BAM not found for path {raw!r}")


def load_samples(path: Path) -> pl.DataFrame:
    """Load the validation sample table; stop if the file or columns are missing."""
    assert path.is_file(), f"MISSING INPUT: {path}"
    frame = pl.read_csv(
        path,
        schema_overrides={"sample_id": pl.Utf8, "bam_path": pl.Utf8},
    )
    missing = [name for name in REQUIRED_SAMPLE_COLUMNS if name not in frame.columns]
    assert not missing, f"STOP: samples CSV missing columns: {missing}"
    assert frame.height > 0, f"STOP: samples CSV has no data rows: {path}"
    empty = frame.filter(
        (pl.col("sample_id").cast(pl.Utf8).str.strip_chars() == "")
        | (pl.col("bam_path").cast(pl.Utf8).str.strip_chars() == "")
    )
    assert empty.height == 0, "STOP: sample_id and bam_path must be non-empty"
    return frame


def run_samtools_coverage(samtools: str, bam: Path) -> pl.DataFrame:
    """Run ``samtools coverage`` and return the tabular contig summary."""
    assert bam.is_file(), f"MISSING INPUT: {bam}"
    index = Path(str(bam) + ".bai")
    csi = Path(str(bam) + ".csi")
    assert (
        index.is_file() or csi.is_file()
    ), f"MISSING INPUT: BAM index (.bai or .csi) required next to {bam}"
    completed = subprocess.run(
        [samtools, "coverage", str(bam)],
        check=True,
        capture_output=True,
        text=True,
    )
    lines = [line for line in completed.stdout.splitlines() if line.strip()]
    assert lines, f"STOP: samtools coverage returned no rows for {bam}"
    header = lines[0].lstrip("#").split("\t")
    required = {
        "rname",
        "startpos",
        "endpos",
        "numreads",
        "covbases",
        "coverage",
        "meandepth",
        "meanbaseq",
        "meanmapq",
    }
    assert (
        set(header) == required
    ), f"STOP: unexpected samtools coverage columns: {header} (expected {sorted(required)})"
    rows: list[dict[str, Any]] = []
    for line in lines[1:]:
        parts = line.split("\t")
        assert len(parts) == len(header), f"STOP: bad coverage row: {line!r}"
        row = dict(zip(header, parts, strict=True))
        rows.append(
            {
                "rname": row["rname"],
                "startpos": int(row["startpos"]),
                "endpos": int(row["endpos"]),
                "numreads": int(row["numreads"]),
                "covbases": int(row["covbases"]),
                "coverage": float(row["coverage"]),
                "meandepth": float(row["meandepth"]),
                "meanbaseq": float(row["meanbaseq"]),
                "meanmapq": float(row["meanmapq"]),
            }
        )
    return pl.DataFrame(rows)


def mean_autosomal_depth(coverage: pl.DataFrame) -> float:
    """Length-weighted mean ``meandepth`` over ``chr1``…``chr22``."""
    present = set(coverage["rname"].to_list())
    missing = [name for name in AUTOSOMES if name not in present]
    assert not missing, f"STOP: BAM coverage missing autosomal contigs: {missing}"
    auto = coverage.filter(pl.col("rname").is_in(list(AUTOSOMES))).with_columns(
        (pl.col("endpos") - pl.col("startpos") + 1).alias("length")
    )
    assert auto.height == 22, f"STOP: expected 22 autosomal rows, got {auto.height}"
    total_length = float(auto["length"].sum())
    assert total_length > 0, "STOP: autosomal contig lengths sum to 0"
    weighted = float((auto["meandepth"] * auto["length"]).sum())
    mean_depth = weighted / total_length
    assert mean_depth > 0, "STOP: mean autosomal depth is 0; cannot estimate mtDNA CN"
    return mean_depth


def mean_chrM_depth(coverage: pl.DataFrame) -> float:
    """Return ``meandepth`` for ``chrM``."""
    mt = coverage.filter(pl.col("rname") == MT_CONTIG)
    assert mt.height == 1, f"STOP: expected exactly one {MT_CONTIG} coverage row, got {mt.height}"
    depth = float(mt["meandepth"][0])
    assert depth > 0, f"STOP: {MT_CONTIG} meandepth is 0"
    return depth


def estimate_copy_number(mean_mt: float, mean_auto: float) -> float:
    """``2 * (mean chrM depth / mean depth over chr1-22)``."""
    assert mean_auto > 0, "STOP: autosomal mean depth must be positive"
    return 2.0 * (mean_mt / mean_auto)


def estimate_one_bam(samtools: str, bam: Path) -> dict[str, float]:
    """Run coverage and return depth + estimate fields for one BAM."""
    coverage = run_samtools_coverage(samtools, bam)
    present = set(coverage["rname"].to_list())
    missing = [name for name in REQUIRED_CONTIGS if name not in present]
    assert not missing, f"STOP: BAM missing required contigs {missing}: {bam}"
    mean_mt = mean_chrM_depth(coverage)
    mean_auto = mean_autosomal_depth(coverage)
    return {
        "mean_chrM_depth": mean_mt,
        "mean_depth_chr1_22": mean_auto,
        "mtdna_copy_number_estimate": estimate_copy_number(mean_mt, mean_auto),
    }


def write_estimates_csv(path: Path, frame: pl.DataFrame) -> None:
    """Write estimates CSV with the required warning as the first header line."""
    path.parent.mkdir(parents=True, exist_ok=True)
    body = frame.write_csv()
    path.write_text(f"# {WARNING}\n{body}", encoding="utf-8")


def write_provenance(path: Path, payload: dict[str, Any]) -> None:
    """Write provenance JSON next to the estimates output."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def run_estimates(
    samples_csv: Path,
    out_csv: Path,
    out_provenance: Path,
) -> pl.DataFrame:
    """Estimate mtDNA CN for every sample in ``samples_csv`` and write outputs."""
    samtools = find_samtools()
    version = samtools_version(samtools)
    samples = load_samples(samples_csv)
    rows: list[dict[str, Any]] = []
    per_sample: list[dict[str, Any]] = []
    for record in samples.iter_rows(named=True):
        sample_id = str(record["sample_id"]).strip()
        bam = resolve_bam_path(str(record["bam_path"]).strip(), samples_csv)
        metrics = estimate_one_bam(samtools, bam)
        row = {
            "sample_id": sample_id,
            "bam_path": str(bam),
            **metrics,
        }
        rows.append(row)
        per_sample.append(
            {
                "sample_id": sample_id,
                "bam_path": str(bam),
                "bam_sha256": _sha256(bam),
                **metrics,
            }
        )
    result = pl.DataFrame(rows)
    write_estimates_csv(out_csv, result)
    provenance = {
        "warning": WARNING,
        "generated_at_utc": datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "git_commit": _git_commit(),
        "formula": "2 * (mean_chrM_depth / length_weighted_mean_meandepth_chr1_to_chr22)",
        "tool": {
            "name": "samtools",
            "path": samtools,
            "version": version,
            "command": "samtools coverage <bam>",
        },
        "python": sys.version.split()[0],
        "inputs": {
            "samples_csv": {
                "path": str(samples_csv),
                "sha256": _sha256(samples_csv),
                "n_rows": samples.height,
            }
        },
        "outputs": {
            "estimates_csv": {
                "path": str(out_csv),
                "sha256": _sha256(out_csv),
            }
        },
        "samples": per_sample,
    }
    write_provenance(out_provenance, provenance)
    return result


@app.command()
def main(
    samples_csv: Path = typer.Option(
        DEFAULT_SAMPLES,
        "--samples-csv",
        help="CSV with sample_id,bam_path columns.",
    ),
    out_csv: Path = typer.Option(
        DEFAULT_OUT_CSV,
        "--out-csv",
        help="Estimates CSV (warning line prefixed).",
    ),
    out_provenance: Path = typer.Option(
        DEFAULT_OUT_PROVENANCE,
        "--out-provenance",
        help="Provenance JSON written next to the estimates.",
    ),
) -> None:
    """Estimate relative mtDNA copy number for the validation sample table."""
    frame = run_estimates(samples_csv, out_csv, out_provenance)
    typer.echo(f"# {WARNING}")
    typer.echo(frame.write_csv())
    typer.echo(f"wrote {out_csv}", err=True)
    typer.echo(f"wrote {out_provenance}", err=True)


if __name__ == "__main__":
    app()
