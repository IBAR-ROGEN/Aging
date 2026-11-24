#!/usr/bin/env python3
"""Extract chrM reads, call with Mutect2 mitochondria mode, HaploGrep 3.

Public GRCh38 alignments only. No Romanian samples. No population-heteroplasmy
or aging interpretation.

Example:
    uv run python scripts/mtdna/extract_and_call_mtdna.py \\
        --samples data/external/mtdna_validation_samples.csv \\
        --reference /path/to/GRCh38.fa
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import re
import shutil
import subprocess
import urllib.request
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, TextIO

import polars as pl
import typer
from omegaconf import OmegaConf

from rogen_aging.config import find_repo_root

REPO_ROOT = find_repo_root()
DEFAULT_SAMPLES = REPO_ROOT / "data" / "external" / "mtdna_validation_samples.csv"
DEFAULT_CONFIG = REPO_ROOT / "config" / "mtdna_calling.yaml"
DEFAULT_OUT_DIR = REPO_ROOT / "analysis" / "mtdna"
REQUIRED_TOOLS = ("samtools", "gatk", "haplogrep3")
SAMPLE_COLUMNS = (
    "sample_id",
    "path_or_url",
    "published_haplogroup",
    "source_citation",
)
NUMT_ABSENT_NOTE = "NUMT filtering not applied"
BAM_SQ_RE = re.compile(r"^@SQ\tSN:([^\t]+)")


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


def _git_commit() -> str:
    proc = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode == 0, f"git rev-parse failed: {proc.stderr.strip()}"
    return proc.stdout.strip()


def require_tools_via_which() -> dict[str, str]:
    """Resolve samtools, gatk, haplogrep3 with ``which``. Stop if any is missing."""
    found: dict[str, str] = {}
    missing: list[str] = []
    which_bin = shutil.which("which")
    assert which_bin is not None, "STOP: shell builtin/helper `which` is not available."
    for name in REQUIRED_TOOLS:
        proc = subprocess.run(
            [which_bin, name],
            capture_output=True,
            text=True,
            check=False,
        )
        path = proc.stdout.strip()
        if proc.returncode != 0 or not path:
            missing.append(name)
        else:
            found[name] = path
    if missing:
        _stop(f"STOP: missing required tools (which): {', '.join(missing)}")
    return found


def tool_version(exe: str, name: str) -> str:
    """Best-effort version string; assert the tool runs."""
    candidates = (
        [exe, "--version"],
        [exe, "version"],
        [exe, "-version"],
    )
    for args in candidates:
        proc = subprocess.run(args, capture_output=True, text=True, check=False)
        text = (proc.stdout or proc.stderr or "").strip()
        if proc.returncode == 0 and text:
            return text.splitlines()[0]
    assert False, f"STOP: could not obtain version from {name} ({exe})"


def load_config(path: Path) -> dict[str, Any]:
    assert path.is_file(), f"STOP: missing config: {path}"
    payload = OmegaConf.to_container(OmegaConf.load(path), resolve=True)
    assert isinstance(payload, dict), f"STOP: config is not a mapping: {path}"
    thresholds = payload.get("thresholds")
    assert isinstance(thresholds, dict), "STOP: config.thresholds missing"
    assert (
        "min_depth" in thresholds and "min_vaf" in thresholds
    ), "STOP: config.thresholds must define min_depth and min_vaf"
    min_depth = int(thresholds["min_depth"])
    min_vaf = float(thresholds["min_vaf"])
    assert min_depth > 0, f"STOP: min_depth must be > 0, got {min_depth}"
    assert 0.0 < min_vaf <= 1.0, f"STOP: min_vaf must be in (0, 1], got {min_vaf}"
    contig_name = str(payload.get("contig_name") or "").strip()
    assert contig_name, "STOP: config.contig_name is empty"
    haplogrep = payload.get("haplogrep")
    assert isinstance(haplogrep, dict) and haplogrep.get(
        "tree"
    ), "STOP: config.haplogrep.tree is required"
    numt_raw = payload.get("numt_mask")
    numt_path = Path(str(numt_raw)) if numt_raw else None
    if numt_path is not None and not numt_path.is_absolute():
        numt_path = REPO_ROOT / numt_path
    return {
        "min_depth": min_depth,
        "min_vaf": min_vaf,
        "contig_name": contig_name,
        "haplogrep_tree": str(haplogrep["tree"]),
        "numt_mask": numt_path,
    }


def load_samples(path: Path) -> pl.DataFrame:
    assert path.is_file(), f"STOP: missing samples CSV: {path}"
    frame = pl.read_csv(path, infer_schema_length=0)
    missing = [col for col in SAMPLE_COLUMNS if col not in frame.columns]
    assert not missing, f"STOP: samples CSV missing columns: {missing}"
    assert frame.height > 0, f"STOP: samples CSV has no rows: {path}"
    for col in SAMPLE_COLUMNS:
        empty = frame.filter(pl.col(col).cast(pl.Utf8).str.strip_chars() == "")
        assert empty.height == 0, f"STOP: empty values in column {col}"
    ids = frame["sample_id"].to_list()
    assert len(ids) == len(set(ids)), "STOP: duplicate sample_id in samples CSV"
    return frame.select(list(SAMPLE_COLUMNS))


def major_clade(haplogroup: str) -> str:
    """Leading alphabetic prefix before the first digit (H1a → H, HV1 → HV)."""
    text = haplogroup.strip()
    assert text, "STOP: empty haplogroup"
    idx = 0
    while idx < len(text) and text[idx].isalpha():
        idx += 1
    assert idx > 0, f"STOP: cannot parse major clade from {haplogroup!r}"
    return text[:idx]


def reference_has_contig(reference: Path, contig: str) -> bool:
    fai = Path(str(reference) + ".fai")
    assert fai.is_file(), f"STOP: reference index missing: {fai}. Run: samtools faidx {reference}"
    with fai.open(encoding="utf-8") as handle:
        for line in handle:
            name = line.split("\t", 1)[0]
            if name == contig:
                return True
    return False


def bam_contig_names(samtools: str, bam: Path) -> list[str]:
    proc = subprocess.run(
        [samtools, "view", "-H", str(bam)],
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode == 0, f"STOP: samtools view -H failed for {bam}: {proc.stderr.strip()}"
    names: list[str] = []
    for line in proc.stdout.splitlines():
        match = BAM_SQ_RE.match(line)
        if match:
            names.append(match.group(1))
    return names


def resolve_alignment(
    path_or_url: str,
    sample_id: str,
    cache_dir: Path,
) -> Path:
    """Return a local BAM/CRAM path. Download HTTP(S) URLs into cache_dir."""
    text = path_or_url.strip()
    if text.startswith("http://") or text.startswith("https://"):
        cache_dir.mkdir(parents=True, exist_ok=True)
        suffix = ".cram" if ".cram" in Path(text).name.lower() else ".bam"
        dest = cache_dir / f"{sample_id}.source{suffix}"
        if not dest.is_file():
            typer.echo(f"download {text} -> {dest}")
            urllib.request.urlretrieve(text, dest)
        assert dest.is_file() and dest.stat().st_size > 0, f"STOP: download failed or empty: {text}"
        return dest
    local = Path(text)
    if not local.is_absolute():
        local = REPO_ROOT / local
    assert local.is_file(), f"STOP: alignment not found for {sample_id}: {local}"
    return local


def extract_chrm(
    samtools: str,
    alignment: Path,
    contig: str,
    out_bam: Path,
    reference: Path,
) -> tuple[int, float]:
    """Write contig-only BAM; return (read_count, mean_depth)."""
    out_bam.parent.mkdir(parents=True, exist_ok=True)
    cmd = [samtools, "view", "-b", "-o", str(out_bam)]
    if alignment.suffix.lower() == ".cram":
        cmd.extend(["-T", str(reference)])
    cmd.extend([str(alignment), contig])
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    assert proc.returncode == 0, f"STOP: samtools view extract failed: {proc.stderr.strip()}"
    idx = subprocess.run(
        [samtools, "index", str(out_bam)],
        capture_output=True,
        text=True,
        check=False,
    )
    assert idx.returncode == 0, f"STOP: samtools index failed: {idx.stderr.strip()}"

    count_proc = subprocess.run(
        [samtools, "view", "-c", str(out_bam)],
        capture_output=True,
        text=True,
        check=False,
    )
    assert count_proc.returncode == 0, f"STOP: samtools view -c failed: {count_proc.stderr.strip()}"
    read_count = int(count_proc.stdout.strip())
    if read_count == 0:
        return 0, 0.0

    depth_proc = subprocess.run(
        [samtools, "depth", "-a", "-r", contig, str(out_bam)],
        capture_output=True,
        text=True,
        check=False,
    )
    assert depth_proc.returncode == 0, f"STOP: samtools depth failed: {depth_proc.stderr.strip()}"
    depths = [int(line.split("\t")[2]) for line in depth_proc.stdout.splitlines() if line]
    assert depths, f"STOP: samtools depth returned no rows for {contig}"
    mean_depth = float(sum(depths)) / float(len(depths))
    return read_count, mean_depth


def run_mutect2(
    gatk: str,
    reference: Path,
    bam: Path,
    contig: str,
    out_vcf: Path,
) -> None:
    out_vcf.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        gatk,
        "Mutect2",
        "-R",
        str(reference),
        "-I",
        str(bam),
        "-L",
        contig,
        "--mitochondria-mode",
        "-O",
        str(out_vcf),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    assert proc.returncode == 0, f"STOP: gatk Mutect2 failed:\n{proc.stderr[-4000:]}"
    assert out_vcf.is_file(), f"STOP: Mutect2 did not write {out_vcf}"


def load_numt_intervals(path: Path | None, contig: str) -> list[tuple[int, int]] | None:
    """Return 1-based inclusive intervals on ``contig``, or None if mask absent."""
    if path is None or not path.is_file():
        return None
    intervals: list[tuple[int, int]] = []
    with path.open(encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#") or line.startswith("track"):
                continue
            parts = line.rstrip("\n").split("\t")
            assert len(parts) >= 3, f"STOP: bad BED line {line_no} in {path}"
            if parts[0] != contig:
                continue
            start0 = int(parts[1])
            end = int(parts[2])
            assert end > start0, f"STOP: empty BED interval at line {line_no}"
            intervals.append((start0 + 1, end))
    return intervals


def overlaps_numt(pos: int, intervals: list[tuple[int, int]]) -> bool:
    for start, end in intervals:
        if start <= pos <= end:
            return True
    return False


def _open_vcf_text(path: Path) -> TextIO:
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")  # type: ignore[return-value]
    return path.open("r", encoding="utf-8")


def parse_vcf_to_rows(
    vcf_path: Path,
    min_depth: int,
    min_vaf: float,
    numt_intervals: list[tuple[int, int]] | None,
) -> list[dict[str, Any]]:
    """Parse Mutect2 VCF into depth/VAF rows (one row per ALT allele)."""
    rows: list[dict[str, Any]] = []
    with _open_vcf_text(vcf_path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            assert len(parts) >= 10, f"STOP: short VCF line in {vcf_path}"
            chrom = parts[0]
            pos_s = parts[1]
            ref = parts[3]
            alt_field = parts[4]
            fmt = parts[8]
            sample = parts[9]
            pos = int(pos_s)
            alts = alt_field.split(",")
            fmt_keys = fmt.split(":")
            fmt_vals = sample.split(":")
            assert len(fmt_keys) == len(
                fmt_vals
            ), f"STOP: FORMAT/sample length mismatch at {chrom}:{pos}"
            fields = dict(zip(fmt_keys, fmt_vals, strict=True))
            depth = None
            if "DP" in fields and fields["DP"] not in {".", ""}:
                depth = int(fields["DP"])
            ad_values: list[int] | None = None
            if "AD" in fields and fields["AD"] not in {".", ""}:
                ad_values = [int(x) for x in fields["AD"].split(",")]
            af_values: list[float] | None = None
            if "AF" in fields and fields["AF"] not in {".", ""}:
                af_values = [float(x) for x in fields["AF"].split(",")]

            for alt_i, alt in enumerate(alts):
                if ad_values is not None:
                    assert (
                        len(ad_values) >= alt_i + 2
                    ), f"STOP: AD length too short at {chrom}:{pos}"
                    alt_ad = ad_values[alt_i + 1]
                    ref_ad = ad_values[0]
                    if depth is None:
                        depth = sum(ad_values)
                    vaf = float(alt_ad) / float(depth) if depth else 0.0
                elif af_values is not None:
                    assert len(af_values) > alt_i, f"STOP: AF length too short at {chrom}:{pos}"
                    vaf = float(af_values[alt_i])
                    alt_ad = None
                    ref_ad = None
                else:
                    assert False, f"STOP: VCF site {chrom}:{pos} has neither AD nor AF"
                assert depth is not None
                passes = depth >= min_depth and vaf >= min_vaf
                if numt_intervals is None:
                    numt_flag = ""
                else:
                    numt_flag = "true" if overlaps_numt(pos, numt_intervals) else "false"
                rows.append(
                    {
                        "chrom": chrom,
                        "pos": pos,
                        "ref": ref,
                        "alt": alt,
                        "depth": depth,
                        "ref_ad": ref_ad if ref_ad is not None else "",
                        "alt_ad": alt_ad if alt_ad is not None else "",
                        "vaf": round(vaf, 6),
                        "passes_thresholds": passes,
                        "numt_overlap": numt_flag,
                    }
                )
    return rows


def write_variant_tsv(
    path: Path,
    rows: list[dict[str, Any]],
    numt_note: str,
    provenance_line: str,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "chrom",
        "pos",
        "ref",
        "alt",
        "depth",
        "ref_ad",
        "alt_ad",
        "vaf",
        "passes_thresholds",
        "numt_overlap",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        handle.write(f"# {numt_note}\n")
        handle.write(f"# {provenance_line}\n")
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def annotate_vcf_header(vcf_path: Path, numt_note: str, provenance_line: str) -> None:
    """Rewrite VCF header to record NUMT status and provenance."""
    assert not str(vcf_path).endswith(
        ".gz"
    ), f"STOP: annotate_vcf_header expects uncompressed VCF, got {vcf_path}"
    lines = vcf_path.read_text(encoding="utf-8").splitlines(keepends=True)
    insert_at = 0
    for i, line in enumerate(lines):
        if line.startswith("##"):
            insert_at = i + 1
        if line.startswith("#CHROM"):
            break
    meta = [
        f'##NUMT_FILTERING="{numt_note}"\n',
        f'##rogen_mtdna_provenance="{provenance_line}"\n',
    ]
    new_lines = lines[:insert_at] + meta + lines[insert_at:]
    tmp = vcf_path.with_name(vcf_path.name + ".tmp")
    tmp.write_text("".join(new_lines), encoding="utf-8")
    tmp.replace(vcf_path)


def write_haplogrep_vcf(
    dest_vcf: Path,
    passing_rows: list[dict[str, Any]],
    contig: str,
    numt_note: str,
    provenance_line: str,
) -> None:
    """Write threshold-passing alleles as 1/1 for HaploGrep (already filtered)."""
    dest_vcf.parent.mkdir(parents=True, exist_ok=True)
    with dest_vcf.open("w", encoding="utf-8") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write(f'##NUMT_FILTERING="{numt_note}"\n')
        handle.write(f'##rogen_mtdna_provenance="{provenance_line}"\n')
        handle.write(f"##contig=<ID={contig},length=16569>\n")
        handle.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
        for row in sorted(passing_rows, key=lambda r: (r["chrom"], int(r["pos"]), r["alt"])):
            handle.write(
                f"{row['chrom']}\t{row['pos']}\t.\t{row['ref']}\t{row['alt']}\t"
                f".\tPASS\t.\tGT\t1/1\n"
            )


def run_haplogrep3(
    haplogrep3: str,
    tree: str,
    in_vcf: Path,
    out_path: Path,
    het_level: float,
) -> str:
    """Classify haplogroup; return the top haplogroup string."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        haplogrep3,
        "classify",
        "--tree",
        tree,
        "--in",
        str(in_vcf),
        "--out",
        str(out_path),
        f"--hetLevel={het_level}",
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    assert proc.returncode == 0, f"STOP: haplogrep3 classify failed:\n{proc.stderr[-4000:]}"
    assert out_path.is_file(), f"STOP: haplogrep3 did not write {out_path}"
    with out_path.open(encoding="utf-8") as handle:
        lines = [ln for ln in handle.read().splitlines() if ln and not ln.startswith("#")]
    assert lines, f"STOP: haplogrep3 output empty: {out_path}"
    header = lines[0].split("\t")
    assert len(lines) >= 2, f"STOP: haplogrep3 output has no data rows: {out_path}"
    row = lines[1].split("\t")
    assert len(row) >= len(header), f"STOP: haplogrep3 row shorter than header in {out_path}"
    record = dict(zip(header, row[: len(header)], strict=True))
    for key in ("Haplogroup", "haplogroup", "Hg", "hg"):
        if key in record and record[key].strip():
            return record[key].strip()
    assert len(row) >= 2, f"STOP: cannot parse haplogroup from {out_path}"
    return row[1].strip()


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def main(
    samples: Path = typer.Option(
        DEFAULT_SAMPLES,
        "--samples",
        help="CSV: sample_id, path_or_url, published_haplogroup, source_citation",
    ),
    reference: Path = typer.Option(
        ...,
        "--reference",
        help="GRCh38 FASTA used for the alignments (indexed .fai required).",
    ),
    config: Path = typer.Option(
        DEFAULT_CONFIG,
        "--config",
        help="Thresholds and optional NUMT mask path.",
    ),
    out_dir: Path = typer.Option(
        DEFAULT_OUT_DIR,
        "--out-dir",
        help="Output directory (per-sample subfolders + concordance).",
    ),
) -> None:
    """Extract chrM, Mutect2 --mitochondria-mode, HaploGrep 3, major-clade n/N."""
    tools = require_tools_via_which()
    versions = {name: tool_version(tools[name], name) for name in REQUIRED_TOOLS}
    cfg = load_config(config)
    sample_table = load_samples(samples)
    assert reference.is_file(), f"STOP: missing reference FASTA: {reference}"
    contig = str(cfg["contig_name"])
    assert reference_has_contig(reference, contig), (
        f"STOP: reference does not contain contig {contig!r} "
        f"(reference contig name must match alignments)"
    )

    numt_mask_path = cfg["numt_mask"]
    assert numt_mask_path is None or isinstance(numt_mask_path, Path)
    numt_intervals = load_numt_intervals(numt_mask_path, contig)
    if numt_intervals is None:
        numt_note = NUMT_ABSENT_NOTE
        numt_status = "absent"
    else:
        assert numt_mask_path is not None
        numt_note = (
            f"NUMT mask applied from {_repo_relative(numt_mask_path)}; "
            "overlapping variants flagged in numt_overlap"
        )
        numt_status = "applied"

    git_commit = _git_commit()
    timestamp = datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")
    samples_sha = _sha256(samples)
    config_sha = _sha256(config)
    provenance_line = (
        f"git={git_commit}; timestamp_utc={timestamp}; "
        f"samples_sha256={samples_sha}; config_sha256={config_sha}; "
        f"tools={{samtools:{versions['samtools']}; "
        f"gatk:{versions['gatk']}; haplogrep3:{versions['haplogrep3']}}}"
    )

    out_dir.mkdir(parents=True, exist_ok=True)
    cache_dir = out_dir / "_downloads"
    summary_rows: list[dict[str, Any]] = []

    for row in sample_table.iter_rows(named=True):
        sample_id = str(row["sample_id"]).strip()
        published = str(row["published_haplogroup"]).strip()
        citation = str(row["source_citation"]).strip()
        sample_dir = out_dir / sample_id
        sample_dir.mkdir(parents=True, exist_ok=True)
        typer.echo(f"== {sample_id} ==")

        alignment = resolve_alignment(str(row["path_or_url"]), sample_id, cache_dir)
        header_contigs = bam_contig_names(tools["samtools"], alignment)
        assert contig in header_contigs, (
            f"STOP: alignment for {sample_id} has no contig {contig!r}. "
            f"Header contigs include: {header_contigs[:20]}"
        )

        chr_bam = sample_dir / f"{sample_id}.chrM.bam"
        read_count, mean_depth = extract_chrm(
            tools["samtools"],
            alignment,
            contig,
            chr_bam,
            reference,
        )
        metrics: dict[str, Any] = {
            "sample_id": sample_id,
            "alignment": str(alignment),
            "contig": contig,
            "read_count": read_count,
            "mean_depth": mean_depth,
            "published_haplogroup": published,
            "published_major_clade": major_clade(published),
            "source_citation": citation,
            "numt_filtering": numt_note,
            "status": "",
        }

        if read_count == 0:
            metrics["status"] = "skipped_zero_chrm_reads"
            write_json(
                sample_dir / "metrics.json",
                {**metrics, "provenance": provenance_line},
            )
            summary_rows.append(
                {
                    "sample_id": sample_id,
                    "status": metrics["status"],
                    "read_count": 0,
                    "mean_depth": 0.0,
                    "published_haplogroup": published,
                    "published_major_clade": major_clade(published),
                    "called_haplogroup": "",
                    "called_major_clade": "",
                    "major_clade_concordant": "",
                    "source_citation": citation,
                }
            )
            typer.echo(f"SKIP {sample_id}: 0 chrM reads")
            continue

        raw_vcf = sample_dir / f"{sample_id}.mutect2.vcf"
        run_mutect2(
            tools["gatk"],
            reference,
            chr_bam,
            contig,
            raw_vcf,
        )
        annotate_vcf_header(raw_vcf, numt_note, provenance_line)

        variant_rows = parse_vcf_to_rows(
            raw_vcf,
            int(cfg["min_depth"]),
            float(cfg["min_vaf"]),
            numt_intervals,
        )
        tsv_path = sample_dir / f"{sample_id}.variants.tsv"
        write_variant_tsv(tsv_path, variant_rows, numt_note, provenance_line)

        passing_rows = [r for r in variant_rows if r["passes_thresholds"]]
        assert passing_rows, (
            f"STOP: no variants passed thresholds for {sample_id} "
            f"(min_depth={cfg['min_depth']}, min_vaf={cfg['min_vaf']})"
        )
        haplo_vcf = sample_dir / f"{sample_id}.haplogrep_input.vcf"
        write_haplogrep_vcf(
            haplo_vcf,
            passing_rows,
            contig,
            numt_note,
            provenance_line,
        )

        haplo_raw = sample_dir / f"{sample_id}.haplogrep.raw.tsv"
        called = run_haplogrep3(
            tools["haplogrep3"],
            str(cfg["haplogrep_tree"]),
            haplo_vcf,
            haplo_raw,
            float(cfg["min_vaf"]),
        )
        haplo_out = sample_dir / f"{sample_id}.haplogrep.tsv"
        haplo_out.write_text(
            f"# {numt_note}\n# {provenance_line}\n{haplo_raw.read_text(encoding='utf-8')}",
            encoding="utf-8",
        )

        called_major = major_clade(called)
        published_major = major_clade(published)
        concordant = called_major == published_major
        metrics.update(
            {
                "status": "ok",
                "called_haplogroup": called,
                "called_major_clade": called_major,
                "major_clade_concordant": concordant,
                "n_variant_rows": len(variant_rows),
                "n_passing_thresholds": len(passing_rows),
                "vcf": _repo_relative(raw_vcf),
                "variants_tsv": _repo_relative(tsv_path),
            }
        )
        write_json(
            sample_dir / "metrics.json",
            {**metrics, "provenance": provenance_line},
        )
        summary_rows.append(
            {
                "sample_id": sample_id,
                "status": "ok",
                "read_count": read_count,
                "mean_depth": round(mean_depth, 4),
                "published_haplogroup": published,
                "published_major_clade": published_major,
                "called_haplogroup": called,
                "called_major_clade": called_major,
                "major_clade_concordant": concordant,
                "source_citation": citation,
            }
        )

    scored = [r for r in summary_rows if r["status"] == "ok"]
    n_concordant = sum(1 for r in scored if r["major_clade_concordant"] is True)
    n_scored = len(scored)
    concordance_str = f"{n_concordant}/{n_scored}"

    summary_path = out_dir / "haplogroup_concordance.tsv"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        handle.write(f"# {numt_note}\n")
        handle.write(f"# {provenance_line}\n")
        handle.write(f"# major_clade_concordance={concordance_str}\n")
        fieldnames = [
            "sample_id",
            "status",
            "read_count",
            "mean_depth",
            "published_haplogroup",
            "published_major_clade",
            "called_haplogroup",
            "called_major_clade",
            "major_clade_concordant",
            "source_citation",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for summary_row in summary_rows:
            writer.writerow(summary_row)

    provenance = {
        "generated_at_utc": timestamp,
        "git_commit": git_commit,
        "tool_versions": versions,
        "tool_paths": tools,
        "inputs": {
            "samples": {
                "path": _repo_relative(samples),
                "sha256": samples_sha,
                "n_rows": sample_table.height,
            },
            "config": {
                "path": _repo_relative(config),
                "sha256": config_sha,
                "min_depth": cfg["min_depth"],
                "min_vaf": cfg["min_vaf"],
                "contig_name": contig,
                "haplogrep_tree": cfg["haplogrep_tree"],
            },
            "reference": str(reference),
            "numt_mask": {
                "path": (_repo_relative(numt_mask_path) if numt_mask_path is not None else None),
                "status": numt_status,
                "note": numt_note,
            },
        },
        "outputs": {
            "out_dir": _repo_relative(out_dir),
            "concordance_tsv": _repo_relative(summary_path),
        },
        "major_clade_concordance": {
            "n_concordant": n_concordant,
            "n_scored": n_scored,
            "n_over_N": concordance_str,
            "definition": (
                "alphabetic prefix before first digit "
                "(H1a→H, HV1→HV); skipped_zero_chrm_reads excluded from N"
            ),
        },
        "prohibited": [
            "Romanian sequencing data",
            "population heteroplasmy claims",
            "aging interpretation",
        ],
    }
    provenance_path = out_dir / "run_provenance.json"
    write_json(provenance_path, provenance)

    typer.echo(f"major_clade_concordance={concordance_str}")
    typer.echo(f"concordance_tsv={_repo_relative(summary_path)}")
    typer.echo(f"provenance={_repo_relative(provenance_path)}")
    typer.echo(f"numt={numt_note}")


if __name__ == "__main__":
    typer.run(main)
