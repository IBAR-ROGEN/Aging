#!/usr/bin/env python3
"""Look up prioritized longevity variants in large lifespan/longevity GWAS.

Replication check only. Most candidates are not expected to replicate.

Example:
    uv run python scripts/longevity/longevity_gwas_lookup.py
"""

from __future__ import annotations

import hashlib
import json
import math
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
from omegaconf import OmegaConf  # noqa: E402

from rogen_aging.config import find_repo_root  # noqa: E402
from rogen_aging.longevity.allele_harmonize import (  # noqa: E402
    harmonize_to_panel_alt,
    maf_from_af,
)
from rogen_aging.ukb.gnomad import normalize_rsid  # noqa: E402

REPO_ROOT = find_repo_root()
DEFAULT_VARIANTS = REPO_ROOT / "data" / "processed" / "variants_47_input_v2.csv"
DEFAULT_ALLELE_EVIDENCE = REPO_ROOT / "data" / "processed" / "allele_evidence_47.csv"
DEFAULT_ALLELE_AUDIT = REPO_ROOT / "analysis" / "panel" / "allele_audit_47.csv"
DEFAULT_SOURCES = REPO_ROOT / "config" / "gwas_sources.yaml"
DEFAULT_OUT_DIR = REPO_ROOT / "analysis" / "longevity"
RELEASE_CHECK_KIND = "summary"


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


def load_sources(path: Path) -> list[dict[str, Any]]:
    """Load GWAS source blocks from YAML."""
    payload = OmegaConf.to_container(OmegaConf.load(path), resolve=True)
    if not isinstance(payload, dict):
        _stop(f"STOP: {path} did not load as a mapping.")
    sources = payload.get("sources")
    if not isinstance(sources, list) or not sources:
        _stop(f"STOP: {path} has no sources list.")
    return [dict(item) for item in sources if isinstance(item, dict)]


def load_panel(path: Path) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Return (resolved rows, not-tested rows)."""
    frame = pl.read_csv(path, schema_overrides={"rsid": pl.Utf8, "alt_status": pl.Utf8})
    required = {"chrom", "pos", "ref", "alt", "rsid", "gene_symbol", "alt_status"}
    missing = sorted(required - set(frame.columns))
    if missing:
        _stop(f"STOP: panel is missing columns: {missing}")
    frame = frame.with_columns(pl.col("rsid").map_elements(normalize_rsid, return_dtype=pl.Utf8))
    resolved = frame.filter(pl.col("alt_status") == "resolved")
    held = frame.filter(pl.col("alt_status") != "resolved")
    if resolved.height == 0:
        _stop("STOP: no alt_status==resolved rows in the panel.")
    return resolved, held


def load_longevitymap_alleles(path: Path) -> dict[str, str]:
    """rsID → LongevityMap allele from allele_evidence, if present."""
    if not path.is_file():
        return {}
    frame = pl.read_csv(
        path, schema_overrides={"rsid": pl.Utf8, "longevitymap_alt_allele": pl.Utf8}
    )
    out: dict[str, str] = {}
    for row in frame.iter_rows(named=True):
        rsid = normalize_rsid(str(row.get("rsid") or ""))
        allele = str(row.get("longevitymap_alt_allele") or "").strip().upper()
        if rsid and allele:
            out[rsid] = allele
    return out


def load_panel_maf(path: Path) -> dict[str, float]:
    """rsID → gnomAD NFE AF of the audit panel alt (for palindrome MAF)."""
    if not path.is_file():
        return {}
    frame = pl.read_csv(path, schema_overrides={"rsid": pl.Utf8})
    out: dict[str, float] = {}
    for row in frame.iter_rows(named=True):
        rsid = normalize_rsid(str(row.get("rsid") or ""))
        raw = row.get("gnomad_nfe_af")
        if not rsid or raw is None:
            continue
        try:
            af = float(raw)
        except (TypeError, ValueError):
            continue
        maf = maf_from_af(af)
        if maf is not None:
            out[rsid] = maf
    return out


def _as_float(value: Any) -> float | None:
    if value is None:
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(number):
        return None
    return number


def scan_source_hits(
    source: dict[str, Any],
    *,
    rsids: set[str],
    chrom_pos_keys: set[str],
) -> pl.DataFrame:
    """Lazily scan one GWAS file and keep panel-matching rows only."""
    raw_path = Path(str(source["path"]))
    path = raw_path if raw_path.is_absolute() else REPO_ROOT / raw_path
    if not path.is_file():
        _stop(
            f"STOP: missing GWAS file for {source.get('id')}: {path}. "
            "Download into data/external/gwas/ (see docs/LONGEVITY_GWAS_LOOKUP.md)."
        )
    separator = str(source.get("separator") or "\t")
    lazy = pl.scan_csv(
        path,
        separator=separator,
        infer_schema_length=10_000,
        ignore_errors=False,
        truncate_ragged_lines=False,
    )
    rsid_col = source.get("rsid")
    match = str(source.get("match") or ("rsid" if rsid_col else "chrom_pos"))
    if match == "rsid":
        if not rsid_col:
            _stop(f"STOP: source {source.get('id')} match=rsid but rsid column is null.")
        col = str(rsid_col)
        hits = (
            lazy.with_columns(pl.col(col).cast(pl.Utf8).alias("_rsid_raw"))
            .filter(pl.col("_rsid_raw").is_in(sorted(rsids)))
            .collect()
        )
    elif match == "chrom_pos":
        chrom_col = str(source["chrom"])
        pos_col = str(source["pos"])
        hits = (
            lazy.with_columns(
                (
                    pl.col(chrom_col).cast(pl.Utf8).str.replace(r"(?i)^chr", "")
                    + pl.lit(":")
                    + pl.col(pos_col).cast(pl.Utf8)
                ).alias("_chrom_pos")
            )
            .filter(pl.col("_chrom_pos").is_in(sorted(chrom_pos_keys)))
            .collect()
        )
    else:
        _stop(f"STOP: unknown match mode {match!r} for {source.get('id')}")
    return hits


def grch37_keys_from_timmers(
    timmers_path: Path,
    rsids: set[str],
) -> dict[str, str]:
    """Map panel rsID → GRCh37 chrom:pos using the Timmers table when present."""
    if not timmers_path.is_file():
        return {}
    hits = (
        pl.scan_csv(timmers_path, separator="\t", infer_schema_length=5_000)
        .filter(pl.col("rsid").cast(pl.Utf8).is_in(sorted(rsids)))
        .select(["rsid", "chr", "pos"])
        .collect()
    )
    out: dict[str, str] = {}
    for row in hits.iter_rows(named=True):
        rsid = normalize_rsid(str(row["rsid"]))
        chrom = str(row["chr"]).removeprefix("chr").removeprefix("CHR")
        out[rsid] = f"{chrom}:{int(row['pos'])}"
    return out


def effect_direction(beta: float | None) -> str:
    if beta is None or not math.isfinite(beta) or beta == 0.0:
        return "zero_or_missing"
    return "longevity_increasing" if beta > 0.0 else "longevity_decreasing"


def longevitymap_concordance(
    *,
    panel_ref: str,
    panel_alt: str,
    longevitymap_allele: str | None,
    beta_panel_alt: float | None,
) -> str:
    """Compare LongevityMap allele with GWAS sign on the panel alt.

    Does not infer the LongevityMap allele from frequency.
    """
    if not longevitymap_allele:
        return "needs_curation"
    if beta_panel_alt is None or not math.isfinite(beta_panel_alt) or beta_panel_alt == 0.0:
        return "needs_curation"
    lm = longevitymap_allele.strip().upper()
    ref = panel_ref.strip().upper()
    alt = panel_alt.strip().upper()
    gwas_favors_alt = beta_panel_alt > 0.0
    if lm == alt:
        return "concordant" if gwas_favors_alt else "discordant"
    if lm == ref:
        return "concordant" if not gwas_favors_alt else "discordant"
    return "needs_curation"


def benjamini_hochberg(p_values: list[float | None]) -> list[float | None]:
    """Return BH q-values aligned with ``p_values`` (None stays None)."""
    indexed = [
        (index, float(value))
        for index, value in enumerate(p_values)
        if value is not None and math.isfinite(float(value))
    ]
    n = len(indexed)
    q_values: list[float | None] = [None] * len(p_values)
    if n == 0:
        return q_values
    indexed.sort(key=lambda item: item[1])
    min_q = 1.0
    ordered_q = [0.0] * n
    for i in range(n - 1, -1, -1):
        rank = i + 1
        min_q = min(min_q, indexed[i][1] * n / rank)
        ordered_q[i] = min(1.0, min_q)
    for (index, _p), q_value in zip(indexed, ordered_q, strict=True):
        q_values[index] = q_value
    return q_values


def build_rows(
    *,
    resolved: pl.DataFrame,
    held: pl.DataFrame,
    sources: list[dict[str, Any]],
    longevitymap: dict[str, str],
    panel_maf: dict[str, float],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Run all lookups and return result rows plus source digests."""
    rsids = {normalize_rsid(str(r)) for r in resolved["rsid"].to_list()}
    timmers = next((s for s in sources if s.get("id") == "timmers2019_parental_lifespan"), None)
    grch37_by_rsid: dict[str, str] = {}
    if timmers is not None:
        timmers_raw = Path(str(timmers["path"]))
        timmers_path = timmers_raw if timmers_raw.is_absolute() else REPO_ROOT / timmers_raw
        grch37_by_rsid = grch37_keys_from_timmers(timmers_path, rsids)

    source_meta: dict[str, Any] = {}
    rows: list[dict[str, Any]] = []

    for source in sources:
        source_id = str(source["id"])
        trait = str(source["trait"])
        raw_path = Path(str(source["path"]))
        path = raw_path if raw_path.is_absolute() else REPO_ROOT / raw_path
        chrom_pos_keys = (
            set(grch37_by_rsid.values()) if source.get("match") == "chrom_pos" else set()
        )
        if source.get("match") == "chrom_pos" and not chrom_pos_keys:
            _stop(
                f"STOP: source {source_id} matches on GRCh37 chrom:pos but no "
                "coordinates are available. Place the Timmers file first so "
                "rsID→GRCh37 keys can be built, or add coordinates another way."
            )
        hits = scan_source_hits(source, rsids=rsids, chrom_pos_keys=chrom_pos_keys)
        source_meta[source_id] = {
            "path": _repo_relative(path),
            "sha256": _sha256(path),
            "citation": source.get("citation"),
            "effect_units": source.get("effect_units"),
            "genome_build": source.get("genome_build"),
            "n_rows_kept": hits.height,
            "readme": source.get("readme"),
            "notes": source.get("notes"),
        }

        by_rsid: dict[str, dict[str, Any]] = {}
        rsid_col = source.get("rsid")
        if rsid_col:
            for row in hits.iter_rows(named=True):
                rsid = normalize_rsid(str(row[str(rsid_col)]))
                by_rsid[rsid] = row
        else:
            # Invert GRCh37 key → rsid for chrom_pos match.
            key_to_rsid = {key: rsid for rsid, key in grch37_by_rsid.items()}
            for row in hits.iter_rows(named=True):
                chrom = str(row[str(source["chrom"])]).removeprefix("chr").removeprefix("CHR")
                key = f"{chrom}:{int(row[str(source['pos'])])}"
                rsid = key_to_rsid.get(key)
                if rsid:
                    by_rsid[rsid] = row

        for panel_row in resolved.iter_rows(named=True):
            rsid = normalize_rsid(str(panel_row["rsid"]))
            hit = by_rsid.get(rsid)
            base = {
                "rsid": rsid,
                "gene_symbol": panel_row["gene_symbol"],
                "chrom_grch38": panel_row["chrom"],
                "pos_grch38": panel_row["pos"],
                "panel_ref": panel_row["ref"],
                "panel_alt": panel_row["alt"],
                "alt_status": panel_row["alt_status"],
                "trait": trait,
                "source_id": source_id,
                "effect_units": source.get("effect_units"),
                "longevitymap_allele": longevitymap.get(rsid) or "",
            }
            if hit is None:
                rows.append(
                    {
                        **base,
                        "status": "absent",
                        "harmonize_status": "",
                        "beta_panel_alt": None,
                        "se": None,
                        "p": None,
                        "eaf_panel_alt": None,
                        "direction": "",
                        "longevitymap_direction": (
                            "needs_curation" if not longevitymap.get(rsid) else ""
                        ),
                        "q_bh": None,
                        "bonferroni_sig": False,
                        "nominal_p_lt_0_05": False,
                    }
                )
                continue

            ea = str(hit[str(source["effect_allele"])])
            oa = str(hit[str(source["other_allele"])])
            beta_raw = _as_float(hit[str(source["beta"])])
            se = _as_float(hit[str(source["se"])])
            p_value = _as_float(hit[str(source["p"])])
            eaf_col = source.get("effect_allele_frequency")
            eaf = _as_float(hit[str(eaf_col)]) if eaf_col else None
            if beta_raw is None:
                rows.append(
                    {
                        **base,
                        "status": "absent",
                        "harmonize_status": "missing_beta",
                        "beta_panel_alt": None,
                        "se": se,
                        "p": p_value,
                        "eaf_panel_alt": eaf,
                        "direction": "",
                        "longevitymap_direction": "needs_curation",
                        "q_bh": None,
                        "bonferroni_sig": False,
                        "nominal_p_lt_0_05": False,
                    }
                )
                continue

            harm = harmonize_to_panel_alt(
                panel_ref=str(panel_row["ref"]),
                panel_alt=str(panel_row["alt"]),
                effect_allele=ea,
                other_allele=oa,
                beta=beta_raw,
                eaf=eaf,
                panel_maf=panel_maf.get(rsid),
            )
            if harm.status == "ambiguous_palindromic":
                status = "ambiguous"
            elif harm.status == "allele_mismatch":
                status = "absent"
            else:
                status = "present"
            lm_dir = longevitymap_concordance(
                panel_ref=str(panel_row["ref"]),
                panel_alt=str(panel_row["alt"]),
                longevitymap_allele=longevitymap.get(rsid),
                beta_panel_alt=harm.beta_panel_alt,
            )
            rows.append(
                {
                    **base,
                    "status": status,
                    "harmonize_status": harm.status,
                    "beta_panel_alt": harm.beta_panel_alt,
                    "se": se,
                    "p": p_value,
                    "eaf_panel_alt": harm.eaf_panel_alt,
                    "direction": effect_direction(harm.beta_panel_alt),
                    "longevitymap_direction": lm_dir,
                    "q_bh": None,
                    "bonferroni_sig": False,
                    "nominal_p_lt_0_05": bool(p_value is not None and p_value < 0.05),
                }
            )

    for panel_row in held.iter_rows(named=True):
        rsid = normalize_rsid(str(panel_row["rsid"]))
        for source in sources:
            rows.append(
                {
                    "rsid": rsid,
                    "gene_symbol": panel_row["gene_symbol"],
                    "chrom_grch38": panel_row["chrom"],
                    "pos_grch38": panel_row["pos"],
                    "panel_ref": panel_row["ref"],
                    "panel_alt": panel_row["alt"],
                    "alt_status": panel_row["alt_status"],
                    "trait": source["trait"],
                    "source_id": source["id"],
                    "effect_units": source.get("effect_units"),
                    "longevitymap_allele": longevitymap.get(rsid) or "",
                    "status": "not_tested",
                    "harmonize_status": "",
                    "beta_panel_alt": None,
                    "se": None,
                    "p": None,
                    "eaf_panel_alt": None,
                    "direction": "",
                    "longevitymap_direction": "needs_curation",
                    "q_bh": None,
                    "bonferroni_sig": False,
                    "nominal_p_lt_0_05": False,
                }
            )
    return rows, source_meta


def apply_multiple_testing(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """BH and Bonferroni over present+ambiguous tests with a p-value."""
    test_indices = [
        index
        for index, row in enumerate(rows)
        if row["status"] == "present" and row["p"] is not None
    ]
    n_tests = len(test_indices)
    alpha_bonf = 0.05 / n_tests if n_tests else None
    p_list: list[float | None] = [None] * len(rows)
    for index in test_indices:
        p_list[index] = float(rows[index]["p"])
    q_list = benjamini_hochberg(p_list)
    n_nominal = 0
    n_bh = 0
    n_bonf = 0
    for index in test_indices:
        rows[index]["q_bh"] = q_list[index]
        p_value = float(rows[index]["p"])
        rows[index]["nominal_p_lt_0_05"] = p_value < 0.05
        if p_value < 0.05:
            n_nominal += 1
        if q_list[index] is not None and float(q_list[index]) < 0.05:
            n_bh += 1
        if alpha_bonf is not None and p_value < alpha_bonf:
            rows[index]["bonferroni_sig"] = True
            n_bonf += 1
    expected_null = 0.05 * n_tests
    return {
        "n_tests": n_tests,
        "alpha_bonferroni": alpha_bonf,
        "n_p_lt_0_05": n_nominal,
        "n_expected_p_lt_0_05_under_null": expected_null,
        "n_q_bh_lt_0_05": n_bh,
        "n_bonferroni_sig": n_bonf,
    }


def render_forest(rows: list[dict[str, Any]], path: Path) -> dict[str, Any]:
    """Forest plot of q<0.05 present hits, else top 10 by p."""
    present = [
        row
        for row in rows
        if row["status"] == "present" and row["beta_panel_alt"] is not None and row["p"] is not None
    ]
    sig = [row for row in present if row["q_bh"] is not None and float(row["q_bh"]) < 0.05]
    selected = sig if sig else sorted(present, key=lambda row: float(row["p"]))[:10]
    selected = list(reversed(selected))
    path.parent.mkdir(parents=True, exist_ok=True)
    if not selected:
        fig, axis = plt.subplots(figsize=(7, 3))
        axis.text(0.5, 0.5, "No present GWAS hits to plot", ha="center", va="center")
        axis.set_axis_off()
        fig.savefig(path, dpi=200, bbox_inches="tight")
        plt.close(fig)
        return {"n_plotted": 0, "selection": "empty"}

    labels = [f"{row['gene_symbol']} {row['rsid']} ({row['trait']})" for row in selected]
    betas = np.array([float(row["beta_panel_alt"]) for row in selected], dtype=float)
    ses = np.array(
        [float(row["se"]) if row["se"] is not None else float("nan") for row in selected],
        dtype=float,
    )
    lo = betas - 1.96 * ses
    hi = betas + 1.96 * ses
    y = np.arange(len(selected))
    fig_h = max(3.5, 0.45 * len(selected) + 1.5)
    fig, axis = plt.subplots(figsize=(8.5, fig_h))
    axis.errorbar(betas, y, xerr=[betas - lo, hi - betas], fmt="o", color="#1b4f72", capsize=3)
    axis.axvline(0.0, color="0.5", linewidth=1.0)
    axis.set_yticks(y)
    axis.set_yticklabels(labels, fontsize=8)
    axis.set_xlabel("Effect on panel alt (source units; see provenance)")
    title = (
        "Variants with BH q < 0.05" if sig else "Top 10 present associations by p (no BH q < 0.05)"
    )
    axis.set_title(title)
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return {"n_plotted": len(selected), "selection": "q_lt_0_05" if sig else "top10_by_p"}


def release_check(path: Path, report: Path) -> dict[str, Any]:
    """Run rogen-release-check summary on the result table."""
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
    status = "passed" if completed.returncode == 0 else "failed"
    return {
        "path": _repo_relative(path),
        "kind": RELEASE_CHECK_KIND,
        "report": _repo_relative(report),
        "argv": argv,
        "returncode": completed.returncode,
        "stdout": completed.stdout,
        "stderr": completed.stderr,
        "status": status,
    }


def main(
    variants: Path = typer.Option(DEFAULT_VARIANTS, "--variants"),
    allele_evidence: Path = typer.Option(DEFAULT_ALLELE_EVIDENCE, "--allele-evidence"),
    allele_audit: Path = typer.Option(DEFAULT_ALLELE_AUDIT, "--allele-audit"),
    sources_yaml: Path = typer.Option(DEFAULT_SOURCES, "--sources"),
    out_dir: Path = typer.Option(DEFAULT_OUT_DIR, "--out-dir"),
) -> None:
    """Look up resolved panel variants in configured longevity GWAS files."""
    started = datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")
    sources = load_sources(sources_yaml)
    resolved, held = load_panel(variants)
    longevitymap = load_longevitymap_alleles(allele_evidence)
    panel_maf = load_panel_maf(allele_audit)
    rows, source_meta = build_rows(
        resolved=resolved,
        held=held,
        sources=sources,
        longevitymap=longevitymap,
        panel_maf=panel_maf,
    )
    mt = apply_multiple_testing(rows)
    out_dir.mkdir(parents=True, exist_ok=True)
    table_path = out_dir / "longevity_gwas_lookup.csv"
    figure_path = out_dir / "longevity_gwas_lookup_forest.png"
    provenance_path = out_dir / "longevity_gwas_lookup_provenance.json"
    report_path = out_dir / "release_checks" / "longevity_gwas_lookup.json"

    frame = pl.DataFrame(rows)
    frame.write_csv(table_path)
    forest_meta = render_forest(rows, figure_path)
    release = release_check(table_path, report_path)
    provenance = {
        "generated_at_utc": started,
        "panel": {
            "path": _repo_relative(variants),
            "n_resolved": resolved.height,
            "n_not_tested": held.height,
            "not_tested_rsids": [normalize_rsid(str(r)) for r in held["rsid"].to_list()],
        },
        "allele_evidence": _repo_relative(allele_evidence),
        "allele_audit": _repo_relative(allele_audit),
        "sources": source_meta,
        "multiple_testing": mt,
        "forest": forest_meta,
        "outputs": {
            "table": _repo_relative(table_path),
            "figure": _repo_relative(figure_path),
            "provenance": _repo_relative(provenance_path),
        },
        "release_check": release,
        "disclaimer": (
            "Candidate-gene replication lookup in public lifespan/longevity GWAS. "
            "Most variants are not expected to replicate; absence of association is "
            "a valid result. Effects are oriented to the panel alt allele."
        ),
    }
    provenance_path.write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    typer.echo(f"Wrote {table_path}")
    typer.echo(f"Wrote {figure_path}")
    typer.echo(f"Wrote {provenance_path}")
    typer.echo(
        f"tests={mt['n_tests']} p<0.05={mt['n_p_lt_0_05']} "
        f"expected_under_null={mt['n_expected_p_lt_0_05_under_null']:.2f} "
        f"q<0.05={mt['n_q_bh_lt_0_05']}"
    )
    if release["status"] != "passed":
        _stop("STOP: release check did not pass.")


if __name__ == "__main__":
    typer.run(main)
