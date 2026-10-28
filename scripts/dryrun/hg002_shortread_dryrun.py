#!/usr/bin/env python3
"""HG002 short-read dry run: download → intake → align → call → GIAB benchmark.

HG002 = NA24385 = GM24385 (public GIAB sample). Uses only public URLs.
bcftools is the baseline caller; numbers are reported as produced (not tuned).

Example:
    uv run python scripts/dryrun/hg002_shortread_dryrun.py check-urls
    uv run python scripts/dryrun/hg002_shortread_dryrun.py fetch
"""

from __future__ import annotations

import gzip
import hashlib
import importlib.util
import json
import os
import re
import shutil
import subprocess
import sys
import time
import urllib.error
import urllib.request
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any
from urllib.parse import urljoin, urlparse

import typer

from rogen_aging.config import find_repo_root
from rogen_aging.intake.errors import IntakeError
from rogen_aging.intake.key import KEY_ENV_VAR, make_key
from rogen_aging.intake.run import run_intake

REPO_ROOT = find_repo_root()

# ---------------------------------------------------------------------------
# Constants (task-specified)
# ---------------------------------------------------------------------------

REGION = "chr20:10000000-20000000"
DOWNSAMPLE = "0.1"  # samtools view -s seed.fraction; seed 0 → about 30x
THREADS = 8
ALIGNER = "minimap2"  # or "bwa-mem2"

REF_URL = (
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/"
    "GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/"
    "GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
)
GIAB_DIR = (
    "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/"
    "AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/"
)
GIAB_VCF_NAME = "HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
GIAB_TBI_NAME = "HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi"
GIAB_BED_NAME = "HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"

WORKDIR = REPO_ROOT / "data" / "dryrun" / "hg002_shortread"
OUTDIR = REPO_ROOT / "analysis" / "dryrun"
CUSTODIAN_DIR = Path.home() / ".rogen" / "custodian" / "hg002_shortread"
KEY_FILE = Path.home() / ".rogen" / "keys" / "hg002_shortread.key"

SAMPLE_ID = "HG002"
SAMPLE_ALIASES = ("NA24385", "GM24385")
DELIVERY_R1 = "delivery_R1.fastq.gz"
DELIVERY_R2 = "delivery_R2.fastq.gz"

_CLOCK_SCRIPT = REPO_ROOT / "scripts" / "clock" / "ont_clock_site_coverage.py"


def _load_illumina_bam_url() -> str:
    """Import ILLUMINA_CHR20_BAM_URL from the ONT clock script (do not copy)."""
    name = "ont_clock_site_coverage_for_dryrun"
    spec = importlib.util.spec_from_file_location(name, _CLOCK_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load {_CLOCK_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    url = getattr(module, "ILLUMINA_CHR20_BAM_URL", None)
    if not isinstance(url, str) or not url.strip():
        raise RuntimeError("ILLUMINA_CHR20_BAM_URL missing from ont_clock_site_coverage.py")
    return url.strip()


ILLUMINA_BAM_URL = _load_illumina_bam_url()

app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    help="HG002 short-read dry run against GIAB v4.2.1 (chr20 window).",
)


# ---------------------------------------------------------------------------
# Paths under WORKDIR / OUTDIR
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Paths:
    """Resolved artifact paths for one dry-run workspace."""

    work: Path
    out: Path
    ref_gz: Path
    ref_fa: Path
    chr20_fa: Path
    chr20_fai: Path
    mm2_index: Path
    bwa_index_prefix: Path
    bam_slice: Path
    delivery_dir: Path
    delivery_r1: Path
    delivery_r2: Path
    manifest: Path
    intake_dir: Path
    giab_vcf: Path
    giab_tbi: Path
    giab_bed: Path
    eval_bed: Path
    aligned_bam: Path
    flagstat: Path
    depth_txt: Path
    metrics_align: Path
    calls_vcf: Path
    sdf_dir: Path
    vcfeval_snv: Path
    vcfeval_indel: Path
    timings: Path
    state: Path
    metrics_json: Path
    report_md: Path


def paths_for(work: Path = WORKDIR, out: Path = OUTDIR) -> Paths:
    """Build the path map for ``work`` / ``out``."""
    delivery = work / "delivery"
    return Paths(
        work=work,
        out=out,
        ref_gz=work / "GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz",
        ref_fa=work / "GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
        chr20_fa=work / "chr20.fa",
        chr20_fai=work / "chr20.fa.fai",
        mm2_index=work / "chr20.fa.mmi",
        bwa_index_prefix=work / "chr20.fa",
        bam_slice=work / "region.downsampled.bam",
        delivery_dir=delivery,
        delivery_r1=delivery / DELIVERY_R1,
        delivery_r2=delivery / DELIVERY_R2,
        manifest=work / "delivery_manifest.tsv",
        intake_dir=work / "intake_out",
        giab_vcf=work / GIAB_VCF_NAME,
        giab_tbi=work / GIAB_TBI_NAME,
        giab_bed=work / GIAB_BED_NAME,
        eval_bed=work / "eval_region.bed",
        aligned_bam=work / "aligned.markdup.bam",
        flagstat=work / "aligned.flagstat.txt",
        depth_txt=work / "aligned.depth.txt",
        metrics_align=work / "align_metrics.json",
        calls_vcf=work / "calls.norm.vcf.gz",
        sdf_dir=work / "chr20.sdf",
        vcfeval_snv=work / "vcfeval_snv",
        vcfeval_indel=work / "vcfeval_indel",
        timings=work / "timings.json",
        state=work / "state.json",
        metrics_json=out / "hg002_shortread_metrics.json",
        report_md=REPO_ROOT / "docs" / "HG002_DRYRUN.md",
    )


# ---------------------------------------------------------------------------
# Region / BED helpers (unit-tested)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class GenomicRegion:
    """1-based inclusive samtools-style region."""

    chrom: str
    start: int | None
    end: int | None

    def as_samtools(self) -> str:
        if self.start is None or self.end is None:
            return self.chrom
        return f"{self.chrom}:{self.start}-{self.end}"

    def bed_start0(self) -> int:
        if self.start is None:
            return 0
        return self.start - 1

    def bed_end(self) -> int | None:
        return self.end


def parse_region(region: str) -> GenomicRegion:
    """Parse ``chr20`` or ``chr20:10000000-20000000``."""
    text = region.strip()
    if not text:
        raise ValueError("region is empty")
    if ":" not in text:
        return GenomicRegion(chrom=text, start=None, end=None)
    chrom, span = text.split(":", 1)
    if "-" not in span:
        raise ValueError(f"region span must be start-end, got {region!r}")
    start_s, end_s = span.split("-", 1)
    start = int(start_s.replace(",", ""))
    end = int(end_s.replace(",", ""))
    if start < 1 or end < start:
        raise ValueError(f"invalid region bounds: {region!r}")
    return GenomicRegion(chrom=chrom, start=start, end=end)


def intersect_bed_interval(
    chrom: str,
    bed_start0: int,
    bed_end: int,
    region: GenomicRegion,
) -> tuple[str, int, int] | None:
    """Intersect one BED interval with ``region``. Returns None if empty."""
    if chrom != region.chrom:
        return None
    r0 = region.bed_start0()
    r1 = region.bed_end()
    start = max(bed_start0, r0)
    end = bed_end if r1 is None else min(bed_end, r1)
    if end <= start:
        return None
    return chrom, start, end


def intersect_bed_with_region(bed_path: Path, region: str, out_path: Path) -> Path:
    """Write BED intervals intersected with ``region`` (0-based half-open)."""
    parsed = parse_region(region)
    lines: list[str] = []
    with bed_path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            chrom, start_s, end_s = parts[0], parts[1], parts[2]
            hit = intersect_bed_interval(chrom, int(start_s), int(end_s), parsed)
            if hit is None:
                continue
            out_chrom, out_start, out_end = hit
            rest = parts[3:]
            lines.append("\t".join([out_chrom, str(out_start), str(out_end), *rest]))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(lines) + ("\n" if lines else ""), encoding="utf-8")
    return out_path


# ---------------------------------------------------------------------------
# vcfeval summary parsing (unit-tested)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class VcfevalSummary:
    """One row from ``rtg vcfeval`` ``summary.txt``."""

    true_pos_baseline: int
    true_pos_call: int
    false_pos: int
    false_neg: int
    precision: float
    sensitivity: float
    f_measure: float

    @property
    def recall(self) -> float:
        return self.sensitivity

    @property
    def f1(self) -> float:
        return self.f_measure


def parse_vcfeval_summary(text: str) -> VcfevalSummary:
    """Parse ``summary.txt`` from ``rtg vcfeval``. Prefer the ``None`` threshold row."""
    rows: list[list[str]] = []
    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        rows.append(re.split(r"\s+", line))
    if not rows:
        raise ValueError("vcfeval summary.txt has no data rows")
    chosen = rows[0]
    for row in rows:
        if row and row[0].lower() == "none":
            chosen = row
            break
    else:
        chosen = rows[-1]
    if len(chosen) < 8:
        raise ValueError(f"unexpected vcfeval summary row: {chosen!r}")
    return VcfevalSummary(
        true_pos_baseline=int(float(chosen[1])),
        true_pos_call=int(float(chosen[2])),
        false_pos=int(float(chosen[3])),
        false_neg=int(float(chosen[4])),
        precision=float(chosen[5]),
        sensitivity=float(chosen[6]),
        f_measure=float(chosen[7]),
    )


# ---------------------------------------------------------------------------
# Command builders (unit-tested)
# ---------------------------------------------------------------------------


def downsample_view_command(
    bam_url: str,
    region: str,
    downsample: str,
    threads: int,
) -> list[str]:
    """``samtools view`` that streams ``region`` and downsamples."""
    return [
        "samtools",
        "view",
        "-b",
        "-@",
        str(threads),
        "-s",
        downsample,
        bam_url,
        region,
    ]


def collate_fastq_commands(
    bam_path: Path,
    r1: Path,
    r2: Path,
    threads: int,
) -> tuple[list[str], list[str]]:
    """Return (collate argv, fastq argv) for paired FASTQ extraction."""
    collate = [
        "samtools",
        "collate",
        "-@",
        str(threads),
        "-O",
        "-u",
        str(bam_path),
        str(bam_path.parent / "collate_tmp"),
    ]
    fastq = [
        "samtools",
        "fastq",
        "-@",
        str(threads),
        "-1",
        str(r1),
        "-2",
        str(r2),
        "-0",
        "/dev/null",
        "-s",
        "/dev/null",
        "-n",
        "-",
    ]
    return collate, fastq


def minimap2_align_command(
    ref_fa: Path,
    index: Path,
    r1: Path,
    r2: Path,
    pseudonym: str,
    threads: int,
) -> list[str]:
    """minimap2 short-read align argv with an ``@RG`` using ``pseudonym``."""
    rg = f"@RG\\tID:{pseudonym}\\tSM:{pseudonym}\\tLB:{pseudonym}\\tPL:ILLUMINA"
    target = str(index) if index.is_file() else str(ref_fa)
    return [
        "minimap2",
        "-ax",
        "sr",
        "-t",
        str(threads),
        "-R",
        rg,
        target,
        str(r1),
        str(r2),
    ]


def bwa_mem2_align_command(
    ref_fa: Path,
    r1: Path,
    r2: Path,
    pseudonym: str,
    threads: int,
) -> list[str]:
    """bwa-mem2 mem argv with an ``@RG`` using ``pseudonym``."""
    rg = f"@RG\\tID:{pseudonym}\\tSM:{pseudonym}\\tLB:{pseudonym}\\tPL:ILLUMINA"
    return [
        "bwa-mem2",
        "mem",
        "-t",
        str(threads),
        "-R",
        rg,
        str(ref_fa),
        str(r1),
        str(r2),
    ]


def markdup_pipeline_commands(threads: int) -> list[list[str]]:
    """fixmate → sort → markdup stages reading SAM/BAM from stdin where applicable."""
    return [
        ["samtools", "collate", "-@", str(threads), "-O", "-u", "-", "markdup_collate"],
        ["samtools", "fixmate", "-@", str(threads), "-m", "-u", "-", "-"],
        ["samtools", "sort", "-@", str(threads), "-u", "-"],
        ["samtools", "markdup", "-@", str(threads), "-", "ALIGNED_OUT"],
    ]


def call_pipeline_commands(ref_fa: Path, region: str, threads: int) -> list[list[str]]:
    """bcftools mpileup | call | norm command stages."""
    return [
        [
            "bcftools",
            "mpileup",
            "-f",
            str(ref_fa),
            "-r",
            region,
            "-a",
            "AD,DP",
            "-Ou",
            "--threads",
            str(threads),
            "BAM",
        ],
        ["bcftools", "call", "-mv", "-Ou", "--threads", str(threads)],
        ["bcftools", "norm", "-f", str(ref_fa), "-m", "-any", "-Oz", "-o", "OUT"],
    ]


def vcfeval_command(
    sdf: Path,
    baseline_vcf: Path,
    calls_vcf: Path,
    bed: Path,
    output_dir: Path,
) -> list[str]:
    """``rtg vcfeval`` argv restricted to ``bed``."""
    return [
        "rtg",
        "vcfeval",
        "-t",
        str(sdf),
        "-b",
        str(baseline_vcf),
        "-c",
        str(calls_vcf),
        "-e",
        str(bed),
        "-o",
        str(output_dir),
        "--vcf-score-field",
        "QUAL",
    ]


# ---------------------------------------------------------------------------
# Tool discovery / process helpers
# ---------------------------------------------------------------------------


def _candidate_bin_dirs() -> list[Path]:
    dirs: list[Path] = []
    which_root = shutil.which("samtools")
    if which_root:
        dirs.append(Path(which_root).resolve().parent)
    micro = Path.home() / "micromamba" / "envs" / "rogen-tools" / "bin"
    if micro.is_dir():
        dirs.append(micro)
    tools = REPO_ROOT / ".tools"
    if (tools / "rtg").is_file() or (tools / "bin").is_dir():
        dirs.append(tools)
        if (tools / "bin").is_dir():
            dirs.append(tools / "bin")
    path_dirs = [Path(part) for part in os.environ.get("PATH", "").split(os.pathsep) if part]
    return dirs + path_dirs


def find_tool(name: str) -> str:
    """Resolve ``name`` on PATH, micromamba ``rogen-tools``, or repo ``.tools``."""
    if name == "rtg":
        wrapped = REPO_ROOT / ".tools" / "rtg"
        if wrapped.is_file() and os.access(wrapped, os.X_OK):
            return str(wrapped)
        nested = REPO_ROOT / ".tools" / "rtg-tools" / "rtg"
        if nested.is_file() and os.access(nested, os.X_OK):
            return str(nested)
    for directory in _candidate_bin_dirs():
        candidate = directory / name
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return str(candidate)
    found = shutil.which(name)
    if found:
        return found
    raise FileNotFoundError(
        f"{name} not found. Install into micromamba env rogen-tools "
        f"(bioconda) or place a binary under .tools/."
    )


def tool_version(name: str) -> str:
    """Best-effort version string for metrics."""
    exe = find_tool(name)
    for args in ([exe, "--version"], [exe, "version"], [exe, "-V"]):
        try:
            proc = subprocess.run(
                args,
                capture_output=True,
                text=True,
                check=False,
                timeout=60,
            )
        except (OSError, subprocess.TimeoutExpired):
            continue
        text = (proc.stdout or proc.stderr or "").strip()
        if text:
            return text.splitlines()[0][:200]
    return exe


def run_cmd(
    argv: list[str],
    *,
    cwd: Path | None = None,
    env: dict[str, str] | None = None,
    check: bool = True,
) -> subprocess.CompletedProcess[str]:
    """Run a command with the resolved first-token tool path when needed."""
    resolved = list(argv)
    if resolved:
        try:
            if "/" not in resolved[0] and not Path(resolved[0]).is_file():
                resolved[0] = find_tool(resolved[0])
        except FileNotFoundError:
            pass
    typer.echo("$ " + " ".join(resolved))
    merged = os.environ.copy()
    # Prefer micromamba + .tools on PATH for pipelines that spawn helpers.
    extras = [
        str(Path.home() / "micromamba" / "envs" / "rogen-tools" / "bin"),
        str(REPO_ROOT / ".tools"),
        str(REPO_ROOT / ".tools" / "bin"),
        str(REPO_ROOT / ".tools" / "jdk-17" / "Contents" / "Home" / "bin"),
    ]
    merged["PATH"] = os.pathsep.join([*extras, merged.get("PATH", "")])
    java_home = REPO_ROOT / ".tools" / "jdk-17" / "Contents" / "Home"
    if java_home.is_dir():
        merged.setdefault("JAVA_HOME", str(java_home))
    if env:
        merged.update(env)
    proc = subprocess.run(
        resolved,
        cwd=str(cwd) if cwd else None,
        env=merged,
        text=True,
        capture_output=True,
        check=False,
    )
    if proc.stdout:
        typer.echo(proc.stdout.rstrip())
    if proc.stderr:
        typer.echo(proc.stderr.rstrip(), err=True)
    if check and proc.returncode != 0:
        raise RuntimeError(
            f"command failed ({proc.returncode}): {' '.join(resolved)}\n{proc.stderr}"
        )
    return proc


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download(url: str, dest: Path) -> None:
    """Download ``url`` to ``dest`` if missing (curl preferred for TLS on macOS)."""
    if dest.is_file() and dest.stat().st_size > 0:
        typer.echo(f"skip download (exists): {dest}")
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    partial = dest.with_suffix(dest.suffix + ".partial")
    typer.echo(f"download {url} -> {dest}")
    curl = shutil.which("curl")
    if curl:
        proc = subprocess.run(
            [curl, "-fL", "--retry", "3", "--connect-timeout", "60", "-o", str(partial), url],
            capture_output=True,
            text=True,
            check=False,
        )
        if proc.returncode != 0:
            partial.unlink(missing_ok=True)
            raise RuntimeError(f"curl download failed for {url}: {proc.stderr.strip()}")
        partial.replace(dest)
        return
    with urllib.request.urlopen(url, timeout=120) as response, partial.open("wb") as handle:
        shutil.copyfileobj(response, handle)
    partial.replace(dest)


def head_ok(url: str) -> tuple[bool, int | None, str]:
    """HEAD ``url``. Returns (ok, status, detail). Prefer curl for system TLS trust."""
    curl = shutil.which("curl")
    if curl:
        proc = subprocess.run(
            [
                curl,
                "-sI",
                "-L",
                "--connect-timeout",
                "60",
                "-o",
                "/dev/null",
                "-w",
                "%{http_code}",
                url,
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        code_s = (proc.stdout or "").strip()
        try:
            code = int(code_s)
        except ValueError:
            return False, None, f"curl HEAD failed: {proc.stderr.strip() or code_s}"
        if 200 <= code < 400:
            return True, code, "ok"
        # Some mirrors reject HEAD; try a 1-byte ranged GET.
        ranged = subprocess.run(
            [
                curl,
                "-sI",
                "-L",
                "--connect-timeout",
                "60",
                "-H",
                "Range: bytes=0-0",
                "-o",
                "/dev/null",
                "-w",
                "%{http_code}",
                url,
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        try:
            ranged_code = int((ranged.stdout or "").strip())
        except ValueError:
            return False, code, f"HEAD {code}; Range probe failed"
        if 200 <= ranged_code < 400:
            return True, ranged_code, "ok-via-range"
        return False, code, f"HTTP {code}"

    request = urllib.request.Request(url, method="HEAD")
    try:
        with urllib.request.urlopen(request, timeout=60) as response:
            return True, getattr(response, "status", 200), "ok"
    except urllib.error.HTTPError as exc:
        if exc.code in {403, 405, 501}:
            get_req = urllib.request.Request(
                url,
                headers={"Range": "bytes=0-0"},
                method="GET",
            )
            try:
                with urllib.request.urlopen(get_req, timeout=60) as response:
                    return True, getattr(response, "status", 200), "ok-via-range"
            except Exception as retry_exc:  # noqa: BLE001 — report and stop upstream
                return False, None, f"HEAD {exc.code}; Range GET failed: {retry_exc}"
        return False, exc.code, str(exc)
    except Exception as exc:  # noqa: BLE001
        return False, None, str(exc)


def parent_listing_url(url: str) -> str:
    parsed = urlparse(url)
    path = parsed.path
    if path.endswith("/"):
        return url
    parent = path.rsplit("/", 1)[0] + "/"
    return urljoin(url, parent)


def load_json(path: Path) -> dict[str, Any]:
    if not path.is_file():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def save_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def record_timing(paths: Paths, step: str, seconds: float) -> None:
    data = load_json(paths.timings)
    data[step] = round(seconds, 3)
    save_json(paths.timings, data)


# ---------------------------------------------------------------------------
# Subcommands
# ---------------------------------------------------------------------------


def all_urls() -> list[tuple[str, str]]:
    """Named URLs checked by ``check-urls``."""
    return [
        ("ILLUMINA_BAM", ILLUMINA_BAM_URL),
        ("ILLUMINA_BAI", ILLUMINA_BAM_URL + ".bai"),
        ("REF", REF_URL),
        ("GIAB_VCF", urljoin(GIAB_DIR, GIAB_VCF_NAME)),
        ("GIAB_TBI", urljoin(GIAB_DIR, GIAB_TBI_NAME)),
        ("GIAB_BED", urljoin(GIAB_DIR, GIAB_BED_NAME)),
    ]


@app.command("check-urls")
def check_urls_cmd() -> None:
    """HEAD every URL. On failure, list the parent folder and stop."""
    failures: list[tuple[str, str, str]] = []
    for name, url in all_urls():
        ok, status, detail = head_ok(url)
        status_s = "?" if status is None else str(status)
        typer.echo(f"{name}: {status_s} {url} ({detail})")
        if not ok:
            failures.append((name, url, detail))
    if not failures:
        typer.echo("all URLs reachable")
        return
    for name, url, detail in failures:
        parent = parent_listing_url(url)
        typer.echo(f"FAILED {name}: {url}", err=True)
        typer.echo(f"  detail: {detail}", err=True)
        typer.echo(f"  parent folder: {parent}", err=True)
        curl = shutil.which("curl")
        try:
            if curl:
                proc = subprocess.run(
                    [curl, "-fsSL", "--connect-timeout", "60", parent],
                    capture_output=True,
                    text=True,
                    check=False,
                )
                body = proc.stdout if proc.returncode == 0 else ""
                if proc.returncode != 0:
                    raise RuntimeError(proc.stderr.strip() or f"curl exit {proc.returncode}")
            else:
                with urllib.request.urlopen(parent, timeout=60) as response:
                    body = response.read(8000).decode("utf-8", errors="replace")
            typer.echo("  parent listing (truncated):", err=True)
            typer.echo(body[:2000], err=True)
        except Exception as exc:  # noqa: BLE001
            typer.echo(f"  could not list parent: {exc}", err=True)
    raise typer.Exit(code=1)


@app.command("fetch")
def fetch_cmd(
    region: str = typer.Option(REGION, "--region"),
    downsample: str = typer.Option(DOWNSAMPLE, "--downsample"),
    threads: int = typer.Option(THREADS, "--threads"),
) -> None:
    """Download reference, GIAB truth, and region FASTQs (skips existing outputs)."""
    paths = paths_for()
    paths.work.mkdir(parents=True, exist_ok=True)
    paths.delivery_dir.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()

    find_tool("samtools")
    find_tool("minimap2" if ALIGNER == "minimap2" else "bwa-mem2")

    download(REF_URL, paths.ref_gz)
    if not paths.ref_fa.is_file():
        typer.echo("decompress reference")
        with gzip.open(paths.ref_gz, "rb") as src, paths.ref_fa.open("wb") as dst:
            shutil.copyfileobj(src, dst)

    ref_fai = Path(str(paths.ref_fa) + ".fai")
    if not ref_fai.is_file():
        run_cmd(["samtools", "faidx", str(paths.ref_fa)])

    if not paths.chr20_fa.is_file():
        proc = subprocess.run(
            [find_tool("samtools"), "faidx", str(paths.ref_fa), "chr20"],
            capture_output=True,
            check=True,
        )
        paths.chr20_fa.write_bytes(proc.stdout)
    if not paths.chr20_fai.is_file():
        run_cmd(["samtools", "faidx", str(paths.chr20_fa)])

    if ALIGNER == "minimap2" and not paths.mm2_index.is_file():
        run_cmd(["minimap2", "-d", str(paths.mm2_index), str(paths.chr20_fa)])
    if ALIGNER == "bwa-mem2" and not Path(str(paths.bwa_index_prefix) + ".bwt.2bit.64").is_file():
        run_cmd(["bwa-mem2", "index", str(paths.chr20_fa)])

    download(urljoin(GIAB_DIR, GIAB_VCF_NAME), paths.giab_vcf)
    download(urljoin(GIAB_DIR, GIAB_TBI_NAME), paths.giab_tbi)
    download(urljoin(GIAB_DIR, GIAB_BED_NAME), paths.giab_bed)
    if not paths.giab_tbi.is_file():
        run_cmd(["bcftools", "index", "-t", "-f", str(paths.giab_vcf)])

    if not (paths.delivery_r1.is_file() and paths.delivery_r2.is_file()):
        if not paths.bam_slice.is_file():
            view = downsample_view_command(ILLUMINA_BAM_URL, region, downsample, threads)
            view[0] = find_tool("samtools")
            typer.echo("$ " + " ".join(view) + f" > {paths.bam_slice}")
            with paths.bam_slice.open("wb") as handle:
                proc = subprocess.run(
                    view,
                    stdout=handle,
                    stderr=subprocess.PIPE,
                    check=False,
                    env={
                        **os.environ,
                        "PATH": str(Path.home() / "micromamba" / "envs" / "rogen-tools" / "bin")
                        + os.pathsep
                        + os.environ.get("PATH", ""),
                    },
                )
            if proc.returncode != 0:
                raise RuntimeError(
                    f"samtools view failed: {proc.stderr.decode('utf-8', errors='replace')}"
                )
            run_cmd(["samtools", "index", str(paths.bam_slice)])
        collate, fastq = collate_fastq_commands(
            paths.bam_slice, paths.delivery_r1, paths.delivery_r2, threads
        )
        collate[0] = find_tool("samtools")
        fastq[0] = find_tool("samtools")
        typer.echo("$ " + " ".join(collate) + " | " + " ".join(fastq))
        env = {
            **os.environ,
            "PATH": str(Path.home() / "micromamba" / "envs" / "rogen-tools" / "bin")
            + os.pathsep
            + os.environ.get("PATH", ""),
        }
        p1 = subprocess.Popen(collate, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=env)
        assert p1.stdout is not None
        p2 = subprocess.run(fastq, stdin=p1.stdout, capture_output=True, env=env, check=False)
        p1.stdout.close()
        c1 = p1.wait()
        if c1 != 0 or p2.returncode != 0:
            err = (p1.stderr.read().decode() if p1.stderr else "") + (
                p2.stderr.decode() if p2.stderr else ""
            )
            raise RuntimeError(f"collate|fastq failed: {err}")

    record_timing(paths, "fetch", time.perf_counter() - started)
    typer.echo(f"fetch complete under {paths.work}")


@app.command("intake")
def intake_cmd() -> None:
    """Write the delivery manifest and run ``rogen-intake`` (pseudonymize)."""
    paths = paths_for()
    started = time.perf_counter()
    if not paths.delivery_r1.is_file() or not paths.delivery_r2.is_file():
        raise typer.Exit("FASTQs missing; run fetch first")

    if paths.intake_dir.is_dir() and any(paths.intake_dir.iterdir()):
        typer.echo(f"skip intake (exists): {paths.intake_dir}")
        record_timing(paths, "intake", 0.0)
        return

    digest_r1 = sha256_file(paths.delivery_r1)
    digest_r2 = sha256_file(paths.delivery_r2)
    aliases = ",".join(SAMPLE_ALIASES)
    lines = [
        "original_sample_id\tfile_name\tfile_type\tchecksum_type\tchecksum\taliases",
        f"{SAMPLE_ID}\t{DELIVERY_R1}\tfastq\tsha256\t{digest_r1}\t{aliases}",
        f"{SAMPLE_ID}\t{DELIVERY_R2}\tfastq\tsha256\t{digest_r2}\t{aliases}",
        "",
    ]
    paths.manifest.write_text("\n".join(lines), encoding="utf-8")

    KEY_FILE.parent.mkdir(parents=True, exist_ok=True)
    if not KEY_FILE.is_file():
        make_key(KEY_FILE)
    os.environ[KEY_ENV_VAR] = str(KEY_FILE)
    CUSTODIAN_DIR.mkdir(parents=True, exist_ok=True)

    try:
        run_intake(paths.manifest, paths.delivery_dir, paths.intake_dir, CUSTODIAN_DIR)
    except IntakeError as exc:
        typer.echo(f"INTAKE LEAK/ERROR: {exc}", err=True)
        typer.echo(
            "Leak check failed. The identifier was not stripped silently. "
            "Inspect the message above for which surface leaked.",
            err=True,
        )
        raise typer.Exit(code=1) from exc

    # Persist pseudonym for later steps (never write original ids into OUTDIR metrics alone).
    linkage = CUSTODIAN_DIR / "linkage.tsv"
    save_json(
        paths.state,
        {
            "intake_dir": str(paths.intake_dir),
            "linkage": str(linkage),
            "manifest": str(paths.manifest),
        },
    )
    record_timing(paths, "intake", time.perf_counter() - started)
    typer.echo(f"intake complete: {paths.intake_dir}")


def _pseudonym_and_fastqs(paths: Paths) -> tuple[str, Path, Path]:
    fastqs = sorted(paths.intake_dir.glob("*.fastq.gz")) + sorted(paths.intake_dir.glob("*.fq.gz"))
    if len(fastqs) < 2:
        # uncompressed
        fastqs = sorted(paths.intake_dir.glob("*.fastq")) + sorted(paths.intake_dir.glob("*.fq"))
    if len(fastqs) < 2:
        raise FileNotFoundError(f"need two FASTQs under {paths.intake_dir}")
    r1 = next((p for p in fastqs if re.search(r"[._]R?1\.", p.name, re.I)), fastqs[0])
    r2 = next((p for p in fastqs if p != r1 and re.search(r"[._]R?2\.", p.name, re.I)), None)
    if r2 is None:
        r2 = next(p for p in fastqs if p != r1)
    # Pseudonym is the filename stem without pair marker / extensions.
    stem = r1.name
    for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if stem.lower().endswith(suffix):
            stem = stem[: -len(suffix)]
            break
    stem = re.sub(r"[._]R?[12]$", "", stem, flags=re.I)
    return stem, r1, r2


@app.command("align")
def align_cmd(
    region: str = typer.Option(REGION, "--region"),
    threads: int = typer.Option(THREADS, "--threads"),
    aligner: str = typer.Option(ALIGNER, "--aligner"),
) -> None:
    """Align pseudonymized FASTQs; sort, markdup, index; write depth metrics."""
    paths = paths_for()
    started = time.perf_counter()
    if paths.aligned_bam.is_file() and paths.metrics_align.is_file():
        typer.echo(f"skip align (exists): {paths.aligned_bam}")
        return
    if not paths.intake_dir.is_dir():
        raise typer.Exit("intake_out missing; run intake first")

    find_tool("samtools")
    find_tool(aligner if aligner != "bwa-mem2" else "bwa-mem2")
    pseudonym, r1, r2 = _pseudonym_and_fastqs(paths)

    if aligner == "bwa-mem2":
        align = bwa_mem2_align_command(paths.chr20_fa, r1, r2, pseudonym, threads)
    else:
        align = minimap2_align_command(paths.chr20_fa, paths.mm2_index, r1, r2, pseudonym, threads)
    align[0] = find_tool(align[0])
    # minimap2/bwa expect literal backslash-t in -R, not real tab characters.

    env = {
        **os.environ,
        "PATH": str(Path.home() / "micromamba" / "envs" / "rogen-tools" / "bin")
        + os.pathsep
        + os.environ.get("PATH", ""),
    }
    tmp_bam = paths.work / "aligned.tmpdir"
    # minimap2 → collate → fixmate → sort → markdup
    cmd_align = align
    cmd_collate = [
        find_tool("samtools"),
        "collate",
        "-@",
        str(threads),
        "-O",
        "-u",
        "-",
        str(paths.work / "aln_collate"),
    ]
    cmd_fixmate = [find_tool("samtools"), "fixmate", "-@", str(threads), "-m", "-u", "-", "-"]
    cmd_sort = [find_tool("samtools"), "sort", "-@", str(threads), "-u", "-"]
    cmd_markdup = [
        find_tool("samtools"),
        "markdup",
        "-@",
        str(threads),
        "-",
        str(paths.aligned_bam),
    ]

    typer.echo(
        "$ "
        + " ".join(cmd_align)
        + " | "
        + " ".join(cmd_collate)
        + " | "
        + " ".join(cmd_fixmate)
        + " | "
        + " ".join(cmd_sort)
        + " | "
        + " ".join(cmd_markdup)
    )
    p_align = subprocess.Popen(cmd_align, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=env)
    assert p_align.stdout is not None
    p_collate = subprocess.Popen(
        cmd_collate, stdin=p_align.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=env
    )
    p_align.stdout.close()
    assert p_collate.stdout is not None
    p_fix = subprocess.Popen(
        cmd_fixmate,
        stdin=p_collate.stdout,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        env=env,
    )
    p_collate.stdout.close()
    assert p_fix.stdout is not None
    p_sort = subprocess.Popen(
        cmd_sort, stdin=p_fix.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=env
    )
    p_fix.stdout.close()
    assert p_sort.stdout is not None
    p_dup = subprocess.run(
        cmd_markdup, stdin=p_sort.stdout, capture_output=True, env=env, check=False
    )
    p_sort.stdout.close()
    codes = [p_align.wait(), p_collate.wait(), p_fix.wait(), p_sort.wait(), p_dup.returncode]
    if any(code != 0 for code in codes):
        err = b"".join(
            filter(
                None,
                [
                    p_align.stderr.read() if p_align.stderr else b"",
                    p_collate.stderr.read() if p_collate.stderr else b"",
                    p_fix.stderr.read() if p_fix.stderr else b"",
                    p_sort.stderr.read() if p_sort.stderr else b"",
                    p_dup.stderr or b"",
                ],
            )
        ).decode("utf-8", errors="replace")
        raise RuntimeError(f"align pipeline failed codes={codes}: {err}")

    run_cmd(["samtools", "index", str(paths.aligned_bam)])
    flag = run_cmd(["samtools", "flagstat", str(paths.aligned_bam)])
    paths.flagstat.write_text(flag.stdout, encoding="utf-8")

    # Mean depth over REGION
    depth = run_cmd(
        [
            "samtools",
            "depth",
            "-a",
            "-r",
            region,
            str(paths.aligned_bam),
        ]
    )
    paths.depth_txt.write_text(depth.stdout, encoding="utf-8")
    depths = [int(line.split("\t")[2]) for line in depth.stdout.splitlines() if line.strip()]
    mean_depth = (sum(depths) / len(depths)) if depths else 0.0

    mapped_pct, dup_pct = _parse_flagstat(flag.stdout)
    metrics = {
        "pseudonym": pseudonym,
        "region": region,
        "aligner": aligner,
        "mean_depth": mean_depth,
        "percent_mapped": mapped_pct,
        "percent_duplicates": dup_pct,
        "bam": str(paths.aligned_bam),
    }
    save_json(paths.metrics_align, metrics)
    if tmp_bam.exists():
        shutil.rmtree(tmp_bam, ignore_errors=True)
    record_timing(paths, "align", time.perf_counter() - started)
    typer.echo(json.dumps(metrics, indent=2))


def _parse_flagstat(text: str) -> tuple[float | None, float | None]:
    """Return (% mapped, % duplicates) from samtools flagstat text."""
    mapped = None
    dups = None
    total = None
    for line in text.splitlines():
        if "in total" in line:
            total = float(line.split()[0])
        if "duplicates" in line and "mate mapped" not in line:
            dups_n = float(line.split()[0])
            if total and total > 0:
                dups = 100.0 * dups_n / total
        if "mapped (" in line and "primary" not in line and "mate mapped" not in line:
            match = re.search(r"mapped \(([0-9.]+)%", line)
            if match:
                mapped = float(match.group(1))
    return mapped, dups


@app.command("call")
def call_cmd(
    region: str = typer.Option(REGION, "--region"),
    threads: int = typer.Option(THREADS, "--threads"),
) -> None:
    """bcftools mpileup | call | norm on REGION (untuned baseline)."""
    paths = paths_for()
    started = time.perf_counter()
    if paths.calls_vcf.is_file():
        typer.echo(f"skip call (exists): {paths.calls_vcf}")
        return
    if not paths.aligned_bam.is_file():
        raise typer.Exit("aligned BAM missing; run align first")
    find_tool("bcftools")
    stages = call_pipeline_commands(paths.chr20_fa, region, threads)
    mpileup, call, norm = stages
    mpileup[0] = find_tool("bcftools")
    call[0] = find_tool("bcftools")
    norm[0] = find_tool("bcftools")
    mpileup[-1] = str(paths.aligned_bam)
    # norm writes to stdout compressed when -Oz without -o? use -o
    if "-o" in norm:
        norm[norm.index("-o") + 1] = str(paths.calls_vcf)
    else:
        norm.extend(["-o", str(paths.calls_vcf)])

    env = {
        **os.environ,
        "PATH": str(Path.home() / "micromamba" / "envs" / "rogen-tools" / "bin")
        + os.pathsep
        + os.environ.get("PATH", ""),
    }
    typer.echo("$ " + " | ".join(" ".join(s) for s in (mpileup, call, norm)))
    p1 = subprocess.Popen(mpileup, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=env)
    assert p1.stdout is not None
    p2 = subprocess.Popen(
        call, stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=env
    )
    p1.stdout.close()
    assert p2.stdout is not None
    p3 = subprocess.run(norm, stdin=p2.stdout, capture_output=True, env=env, check=False)
    p2.stdout.close()
    codes = [p1.wait(), p2.wait(), p3.returncode]
    if any(c != 0 for c in codes):
        err = b"".join(
            filter(
                None,
                [
                    p1.stderr.read() if p1.stderr else b"",
                    p2.stderr.read() if p2.stderr else b"",
                    p3.stderr or b"",
                ],
            )
        ).decode("utf-8", errors="replace")
        raise RuntimeError(f"call pipeline failed codes={codes}: {err}")
    run_cmd(["bcftools", "index", "-t", "-f", str(paths.calls_vcf)])
    record_timing(paths, "call", time.perf_counter() - started)
    typer.echo(f"wrote {paths.calls_vcf}")


def _filter_vcf_by_type(src: Path, dest: Path, variant_type: str) -> None:
    """Write SNVs or indels from ``src`` into ``dest`` using bcftools."""
    if dest.is_file():
        return
    run_cmd(
        [
            "bcftools",
            "view",
            "-v",
            variant_type,
            "-Oz",
            "-o",
            str(dest),
            str(src),
        ]
    )
    run_cmd(["bcftools", "index", "-t", "-f", str(dest)])


@app.command("benchmark")
def benchmark_cmd(
    region: str = typer.Option(REGION, "--region"),
) -> None:
    """rtg vcfeval vs GIAB v4.2.1 for SNVs and indels on BED ∩ REGION."""
    paths = paths_for()
    started = time.perf_counter()
    find_tool("rtg")
    find_tool("bcftools")
    if not paths.calls_vcf.is_file():
        raise typer.Exit("calls missing; run call first")

    intersect_bed_with_region(paths.giab_bed, region, paths.eval_bed)
    if not paths.sdf_dir.is_dir():
        run_cmd(["rtg", "format", "-o", str(paths.sdf_dir), str(paths.chr20_fa)])

    # Restrict truth + calls to REGION for speed, then split by type.
    truth_region = paths.work / "truth.region.vcf.gz"
    calls_region = paths.work / "calls.region.vcf.gz"
    if not truth_region.is_file():
        run_cmd(
            [
                "bcftools",
                "view",
                "-r",
                region,
                "-Oz",
                "-o",
                str(truth_region),
                str(paths.giab_vcf),
            ]
        )
        run_cmd(["bcftools", "index", "-t", "-f", str(truth_region)])
    if not calls_region.is_file():
        run_cmd(
            [
                "bcftools",
                "view",
                "-r",
                region,
                "-Oz",
                "-o",
                str(calls_region),
                str(paths.calls_vcf),
            ]
        )
        run_cmd(["bcftools", "index", "-t", "-f", str(calls_region)])

    def _summary_dict(parsed: VcfevalSummary) -> dict[str, Any]:
        return {
            "true_pos_baseline": parsed.true_pos_baseline,
            "true_pos_call": parsed.true_pos_call,
            "false_pos": parsed.false_pos,
            "false_neg": parsed.false_neg,
            "precision": parsed.precision,
            "recall": parsed.recall,
            "f1": parsed.f1,
        }

    summaries: dict[str, Any] = {}
    for label, vtype in (("snv", "snps"), ("indel", "indels")):
        out_dir = paths.vcfeval_snv if label == "snv" else paths.vcfeval_indel
        summary_path = out_dir / "summary.txt"
        if summary_path.is_file():
            summaries[label] = _summary_dict(
                parse_vcfeval_summary(summary_path.read_text(encoding="utf-8"))
            )
            continue
        truth_f = paths.work / f"truth.{label}.vcf.gz"
        calls_f = paths.work / f"calls.{label}.vcf.gz"
        _filter_vcf_by_type(truth_region, truth_f, vtype)
        _filter_vcf_by_type(calls_region, calls_f, vtype)
        if out_dir.exists():
            shutil.rmtree(out_dir)
        cmd = vcfeval_command(paths.sdf_dir, truth_f, calls_f, paths.eval_bed, out_dir)
        run_cmd(cmd)
        summaries[label] = _summary_dict(
            parse_vcfeval_summary(summary_path.read_text(encoding="utf-8"))
        )

    save_json(paths.work / "benchmark_summaries.json", summaries)
    record_timing(paths, "benchmark", time.perf_counter() - started)
    typer.echo(json.dumps(summaries, indent=2))


@app.command("report")
def report_cmd(
    region: str = typer.Option(REGION, "--region"),
    downsample: str = typer.Option(DOWNSAMPLE, "--downsample"),
) -> None:
    """Write metrics JSON and ``docs/HG002_DRYRUN.md``."""
    paths = paths_for()
    paths.out.mkdir(parents=True, exist_ok=True)
    align_metrics = load_json(paths.metrics_align)
    summaries = load_json(paths.work / "benchmark_summaries.json")
    timings = load_json(paths.timings)

    def _summ(label: str) -> dict[str, Any]:
        block = summaries.get(label) or {}
        return {
            "TP": block.get("true_pos_baseline"),
            "TP_call": block.get("true_pos_call"),
            "FP": block.get("false_pos"),
            "FN": block.get("false_neg"),
            "precision": block.get("precision"),
            "recall": block.get("recall"),
            "f1": block.get("f1"),
        }

    versions: dict[str, str] = {}
    for name in ("samtools", "bcftools", "minimap2", "rtg"):
        try:
            versions[name] = tool_version(name)
        except FileNotFoundError:
            versions[name] = "not found"

    metrics = {
        "sample": SAMPLE_ID,
        "aliases": list(SAMPLE_ALIASES),
        "region": region,
        "downsample_spec": downsample,
        "downsample_seed": 0 if downsample.startswith(".") or downsample.startswith("0.") else None,
        "aligner": ALIGNER,
        "mean_depth": align_metrics.get("mean_depth"),
        "percent_mapped": align_metrics.get("percent_mapped"),
        "percent_duplicates": align_metrics.get("percent_duplicates"),
        "snv": _summ("snv"),
        "indel": _summ("indel"),
        "tool_versions": versions,
        "wall_time_seconds": timings,
        "generated_at_utc": datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "notes": (
            "bcftools mpileup/call/norm is an untuned baseline caller. "
            "F1 is reported as produced; not optimized."
        ),
    }
    # Fix seed parsing: DOWNSAMPLE "0.1" → seed 0
    whole, _, _frac = downsample.partition(".")
    metrics["downsample_seed"] = int(whole) if whole else 0

    save_json(paths.metrics_json, metrics)

    snv = metrics["snv"]
    indel = metrics["indel"]
    md = f"""# HG002 short-read dry run

Public GIAB sample **HG002** (= NA24385 = GM24385). Pipeline:
download → `rogen-intake` → align → bcftools call → `rtg vcfeval` vs NISTv4.2.1.

bcftools is a baseline caller. Numbers below are **not** tuned for F1.

## Window

| Item | Value |
|------|-------|
| Region | `{region}` |
| Downsample | samtools `-s {downsample}` (seed {metrics["downsample_seed"]}, ~30x target) |
| Aligner | `{ALIGNER}` |
| Mean depth | {metrics.get("mean_depth")} |
| % mapped | {metrics.get("percent_mapped")} |
| % duplicates | {metrics.get("percent_duplicates")} |

## Benchmark (GIAB v4.2.1, BED ∩ region)

### SNVs

| TP | FP | FN | Precision | Recall | F1 |
|----|----|----|-----------|--------|----|
| {snv.get("TP")} | {snv.get("FP")} | {snv.get("FN")} | {snv.get("precision")} | {snv.get("recall")} | {snv.get("f1")} |

### Indels

| TP | FP | FN | Precision | Recall | F1 |
|----|----|----|-----------|--------|----|
| {indel.get("TP")} | {indel.get("FP")} | {indel.get("FN")} | {indel.get("precision")} | {indel.get("recall")} | {indel.get("f1")} |

## Tool versions

```text
{json.dumps(versions, indent=2)}
```

## Wall time (seconds)

```text
{json.dumps(timings, indent=2)}
```

## Reproduce

```bash
uv run python scripts/dryrun/hg002_shortread_dryrun.py check-urls
uv run python scripts/dryrun/hg002_shortread_dryrun.py fetch
uv run python scripts/dryrun/hg002_shortread_dryrun.py intake
uv run python scripts/dryrun/hg002_shortread_dryrun.py align
uv run python scripts/dryrun/hg002_shortread_dryrun.py call
uv run python scripts/dryrun/hg002_shortread_dryrun.py benchmark
uv run python scripts/dryrun/hg002_shortread_dryrun.py report
```

Metrics JSON: [`analysis/dryrun/hg002_shortread_metrics.json`](../analysis/dryrun/hg002_shortread_metrics.json).
"""
    paths.report_md.write_text(md, encoding="utf-8")
    typer.echo(f"Wrote {paths.metrics_json}")
    typer.echo(f"Wrote {paths.report_md}")


if __name__ == "__main__":
    app()
