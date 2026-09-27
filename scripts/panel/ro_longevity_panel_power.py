#!/usr/bin/env python3
"""Rank prioritized variants by illustrative detectability in public European frequencies.

This is a detectability calculation. It does not estimate a longevity effect in
any Romanian or ROGEN sample. Control allele frequency is the gnomAD v4
non-Finnish European minor-allele frequency. 1000 Genomes phase 3 EUR
subpopulation frequencies are reported separately and are never averaged into
a Romanian frequency.

Example:
    uv run python scripts/panel/ro_longevity_panel_power.py
"""

from __future__ import annotations

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
from scipy.stats import norm  # noqa: E402

from rogen_aging.config import find_repo_root, get_config  # noqa: E402
from rogen_aging.ensembl.cache import cache_key_for, open_cache  # noqa: E402
from rogen_aging.ensembl.client import (  # noqa: E402
    ENSEMBL_ASSEMBLY,
    ENSEMBL_RELEASE,
    EnsemblApiError,
    EnsemblClient,
)
from rogen_aging.ukb.gnomad import (  # noqa: E402
    GNOMAD_API_URL,
    GNOMAD_DATASET,
    GNOMAD_POPULATION,
    LookupPlan,
    load_cache,
    normalize_chromosome_label,
    normalize_rsid,
    to_gnomad_variant_id,
)

REPO_ROOT = find_repo_root()
EXPECTED_VARIANT_ROWS = 47
EUR_SUBPOPS: tuple[str, ...] = ("TSI", "IBS", "CEU", "GBR", "FIN")
N_CASES_PER_ARM: tuple[int, ...] = (250, 500, 1000, 2500)
ALLELIC_ODDS_RATIOS: tuple[float, ...] = (1.1, 1.2, 1.5, 2.0)
ALPHA_BONFERRONI = 0.05 / EXPECTED_VARIANT_ROWS
ALPHA_GENOMEWIDE = 5e-8
RANK_N_CASES = 1000
RANK_ODDS_RATIO = 1.5
TSI_NFE_MAF_GAP = 0.05
RARE_MAF = 0.01

DISCLAIMER = (
    "Cohort sizes and odds ratios are illustrative. ROGEN sample size "
    "and design are not fixed. Frequencies are from public reference populations, "
    "not from Romanian data."
)

DEFAULT_VARIANTS = REPO_ROOT / "data" / "processed" / "variants_47_input.csv"
DEFAULT_GENE_LIST = REPO_ROOT / "manuscript" / "tables" / "41_gene_candidate_list.csv"
DEFAULT_GNOMAD_CACHE = REPO_ROOT / "data" / "geo" / "gnomad_r4_nfe_cache.json"
DEFAULT_ENSEMBL_CACHE = REPO_ROOT / "data" / "geo" / "ensembl_r116_variation_cache"
DEFAULT_OUT_DIR = REPO_ROOT / "analysis" / "panel"

VARIANT_COLUMNS = ("chrom", "pos", "ref", "alt", "rsid", "gene_symbol")
RELEASE_CHECK_KIND = "summary"


def _stop(message: str) -> None:
    typer.echo(message, err=True)
    raise typer.Exit(code=1)


def maf_from_af(allele_frequency: float | None) -> float | None:
    """Return the minor-allele frequency of a diallelic frequency, or None."""
    if allele_frequency is None:
        return None
    if not math.isfinite(allele_frequency) or allele_frequency < 0.0 or allele_frequency > 1.0:
        return None
    return float(min(allele_frequency, 1.0 - allele_frequency))


def case_allele_frequency(control_maf: float, odds_ratio: float) -> float:
    """Allele frequency in cases under an allelic odds ratio.

    Args:
        control_maf: Control minor-allele frequency, in (0, 1).
        odds_ratio: Allelic odds ratio for the minor allele.

    Returns:
        Case minor-allele frequency.
    """
    return (odds_ratio * control_maf) / (1.0 - control_maf + odds_ratio * control_maf)


def allelic_power(
    control_maf: float | None,
    odds_ratio: float,
    n_per_arm: int,
    alpha: float,
) -> float | None:
    """Two-sided allelic z-test power, 1:1 case-control, normal approximation.

    Allele counts are ``2 * n_per_arm`` in each arm. The critical value uses
    the null standard error of the pooled allele-frequency difference. The
    rejection probability uses the alternative standard error. Power is
    undefined when the control minor-allele frequency is missing, 0, or 1.
    """
    if control_maf is None or not math.isfinite(control_maf):
        return None
    if control_maf <= 0.0 or control_maf >= 1.0:
        return None
    if n_per_arm < 1 or odds_ratio <= 0.0 or not (0.0 < alpha < 1.0):
        return None
    p0 = control_maf
    p1 = case_allele_frequency(p0, odds_ratio)
    n0 = 2 * n_per_arm
    n1 = 2 * n_per_arm
    pooled = (n1 * p1 + n0 * p0) / (n1 + n0)
    se_null = math.sqrt(pooled * (1.0 - pooled) * (1.0 / n1 + 1.0 / n0))
    se_alt = math.sqrt(p1 * (1.0 - p1) / n1 + p0 * (1.0 - p0) / n0)
    if se_null == 0.0 or se_alt == 0.0:
        return None
    z_crit = float(norm.ppf(1.0 - alpha / 2.0))
    delta = p1 - p0
    upper = z_crit * se_null / se_alt - delta / se_alt
    lower = -z_crit * se_null / se_alt - delta / se_alt
    power = float(norm.sf(upper) + norm.cdf(lower))
    if power < 0.0:
        return 0.0
    if power > 1.0:
        return 1.0
    return power


def _population_code(label: str) -> str | None:
    """Return a phase 3 EUR subpopulation code, or None for other panels."""
    upper = label.strip().upper()
    if "1000GENOMES" not in upper or "PHASE_3" not in upper:
        return None
    code = upper.rsplit(":", 1)[-1]
    if code in EUR_SUBPOPS:
        return code
    return None


def alt_allele_frequency(
    populations: list[dict[str, Any]] | None,
    subpopulation: str,
    ref: str,
    alt: str,
) -> tuple[float | None, str | None]:
    """Frequency of ``alt`` in one 1000 Genomes phase 3 subpopulation.

    If Ensembl lists only the reference allele, the alternate frequency is
    the complement. A missing population, or alleles that match neither ref
    nor alt, stays empty. Nothing is imputed from another population.
    """
    if not populations:
        return (
            None,
            f"no population frequencies in the Ensembl variation record for {subpopulation}",
        )
    ref_u = ref.strip().upper()
    alt_u = alt.strip().upper()
    listed: list[tuple[str, float]] = []
    for row in populations:
        if not isinstance(row, dict):
            continue
        if _population_code(str(row.get("population", ""))) != subpopulation:
            continue
        frequency = row.get("frequency")
        allele = row.get("allele")
        if frequency is None or allele is None:
            continue
        listed.append((str(allele).strip().upper(), float(frequency)))
    if not listed:
        return None, f"1000 Genomes phase 3 {subpopulation} frequency not in the Ensembl record"
    for allele, frequency in listed:
        if allele == alt_u:
            return frequency, None
    for allele, frequency in listed:
        if allele == ref_u:
            return 1.0 - frequency, None
    return None, (
        f"1000 Genomes phase 3 {subpopulation} alleles do not match ref/alt; frequency left empty"
    )


def load_variants(path: Path) -> pl.DataFrame:
    """Load the prioritized variant table and require exactly 47 data rows."""
    if not path.is_file():
        _stop(
            f"STOP: prioritized variant file not found: {path}. "
            "The 47-variant list was not rebuilt from other files."
        )
    frame = pl.read_csv(path, schema_overrides={"chrom": pl.Utf8, "rsid": pl.Utf8})
    missing = [column for column in VARIANT_COLUMNS if column not in frame.columns]
    if missing:
        _stop(f"STOP: {path} is missing columns {missing}. Found {frame.columns}.")
    if frame.height != EXPECTED_VARIANT_ROWS:
        _stop(
            f"STOP: {path} has {frame.height} data rows; expected {EXPECTED_VARIANT_ROWS}. "
            "The variant list was not rebuilt."
        )
    return frame.select(list(VARIANT_COLUMNS)).with_columns(
        pl.col("pos").cast(pl.Int64),
        pl.col("ref").cast(pl.Utf8).str.strip_chars(),
        pl.col("alt").cast(pl.Utf8).str.strip_chars(),
        pl.col("rsid").cast(pl.Utf8).str.strip_chars(),
        pl.col("gene_symbol").cast(pl.Utf8).str.strip_chars(),
    )


def load_longevitymap_classes(path: Path) -> pl.DataFrame:
    """Load gene-level LongevityMap classes from the manuscript candidate list.

    The file is an annotation table. This function does not query LongevityMap.
    """
    if not path.is_file():
        _stop(f"STOP: LongevityMap annotation table not found: {path}")
    frame = pl.read_csv(path)
    required = {"Gene_Symbol", "Longevity_Class"}
    missing = required - set(frame.columns)
    if missing:
        _stop(f"STOP: {path} is missing columns {sorted(missing)}.")
    classes = frame.select(
        pl.col("Gene_Symbol").cast(pl.Utf8).str.strip_chars().alias("gene_symbol"),
        pl.col("Longevity_Class").cast(pl.Utf8).str.strip_chars().alias("longevitymap_class"),
    )
    if classes["gene_symbol"].n_unique() != classes.height:
        _stop(f"STOP: {path} has duplicate Gene_Symbol values; join would duplicate variants.")
    return classes


def _gnomad_plans(variants: pl.DataFrame) -> list[LookupPlan]:
    plans: list[LookupPlan] = []
    for row in variants.iter_rows(named=True):
        rsid = normalize_rsid(str(row["rsid"]))
        chrom = normalize_chromosome_label(str(row["chrom"]))
        plans.append(
            LookupPlan(
                rsid=rsid,
                af_1kg=None,
                variant_id=to_gnomad_variant_id(
                    chrom, int(row["pos"]), str(row["ref"]), str(row["alt"])
                ),
                chromosome=chrom,
                position=int(row["pos"]),
                ref=str(row["ref"]),
                alt=str(row["alt"]),
            )
        )
    return plans


def fetch_nfe_frequencies(
    variants: pl.DataFrame,
    cache_path: Path,
) -> tuple[dict[str, dict[str, Any]], dict[str, int]]:
    """Read the gnomAD NFE cache. Do not query the API and do not rewrite the file."""
    plans = _gnomad_plans(variants)
    existing = load_cache(cache_path) if cache_path.is_file() else {}
    hits = sum(1 for plan in plans if plan.rsid in existing)
    return existing, {
        "cache_hits": hits,
        "new_queries": 0,
        "absent_from_cache": len(plans) - hits,
    }


def release_number_from_info(payload: dict[str, Any]) -> int | None:
    """Return the single release reported by ``GET /info/data``, or None."""
    release = payload.get("release")
    if isinstance(release, int):
        return release
    if isinstance(release, str) and release.isdigit():
        return int(release)
    releases = payload.get("releases")
    if isinstance(releases, list) and len(releases) == 1:
        item = releases[0]
        if isinstance(item, int):
            return item
        if isinstance(item, str) and item.isdigit():
            return int(item)
    return None


def ensembl_client_from_config() -> EnsemblClient:
    """Build the variation client from ``apis.ensembl_rest_dated`` and the UKB timeout."""
    cfg = get_config()
    host = cfg.apis.ensembl_rest_dated
    if host is None or not str(host).strip():
        _stop("STOP: config apis.ensembl_rest_dated is empty.")
    return EnsemblClient(
        base_url=str(host).strip(),
        min_interval_sec=float(cfg.ukb.ensembl.min_interval_sec),
        timeout_sec=float(cfg.ukb.ensembl.timeout_sec),
        max_retries=int(cfg.ukb.ensembl.max_retries),
    )


def assert_ensembl_release_116(client: EnsemblClient) -> tuple[int, str]:
    """Stop unless ``GET /info/data`` reports release 116 on this host."""
    host = client.base_url.rstrip("/")
    try:
        info = client.info_data()
    except EnsemblApiError as exc:
        _stop(f"STOP: Ensembl /info/data failed on {host}: {exc}")
    release = release_number_from_info(info)
    if release is None:
        _stop(f"STOP: Ensembl /info/data on {host} did not report a single release: {info}")
    if release != ENSEMBL_RELEASE:
        _stop(f"STOP: Ensembl host {host} reports release {release}, expected {ENSEMBL_RELEASE}.")
    return release, host


def fetch_ensembl_variations(
    rsids: list[str],
    cache_path: Path,
) -> tuple[dict[str, Any | None], dict[str, int], int, str]:
    """Fetch Ensembl variation records. One failed rsID does not stop the run."""
    client = ensembl_client_from_config()
    release, host = assert_ensembl_release_116(client)
    cache = open_cache(cache_path, backend="json")
    payloads: dict[str, Any | None] = {}
    hits = 0
    queries = 0
    failures = 0
    try:
        for rsid in rsids:
            key = cache_key_for("variation", rsid, params={"pops": 1})
            cached = cache.get(key)
            if cached is not None:
                hits += 1
                if isinstance(cached, dict) and cached.get("__ensembl_missing__") is True:
                    payloads[rsid] = None
                else:
                    payloads[rsid] = cached
                continue
            queries += 1
            try:
                payload = client.get_variation(rsid, pops=True)
            except EnsemblApiError as exc:
                failures += 1
                payloads[rsid] = {
                    "__ensembl_request_failed__": True,
                    "error": str(exc),
                }
                continue
            if payload is None:
                cache.set(key, {"__ensembl_missing__": True, "id": rsid})
                payloads[rsid] = None
            else:
                cache.set(key, payload)
                payloads[rsid] = payload
    finally:
        cache.close()
        client.close()
    counts = {
        "cache_hits": hits,
        "new_queries": queries,
        "request_failures": failures,
    }
    return payloads, counts, release, host


def _join_notes(notes: list[str]) -> str:
    unique: list[str] = []
    for note in notes:
        if note and note not in unique:
            unique.append(note)
    return "; ".join(unique)


def build_panel_rows(
    variants: pl.DataFrame,
    classes: pl.DataFrame,
    gnomad_cache: dict[str, dict[str, Any]],
    ensembl_payloads: dict[str, Any | None],
) -> list[dict[str, Any]]:
    """One row per prioritized variant. Missing frequencies stay empty."""
    class_by_gene = {
        str(row["gene_symbol"]).upper(): row["longevitymap_class"]
        for row in classes.iter_rows(named=True)
    }
    rows: list[dict[str, Any]] = []
    for variant in variants.iter_rows(named=True):
        rsid = normalize_rsid(str(variant["rsid"]))
        ref = str(variant["ref"])
        alt = str(variant["alt"])
        notes: list[str] = []
        cached = gnomad_cache.get(rsid)
        af_nfe: float | None = None
        if cached is None:
            notes.append(
                "gnomAD v4 NFE cache has no entry; frequency left empty and gnomAD was not re-queried"
            )
        else:
            raw_af = cached.get("af_gnomad_nfe")
            if raw_af is None:
                notes.append(
                    "variant not found in gnomAD v4 NFE or NFE allele count was empty; "
                    "frequency left empty"
                )
            else:
                af_nfe = float(raw_af)
        maf_nfe = maf_from_af(af_nfe)
        if af_nfe is not None and maf_nfe is None:
            notes.append("gnomAD NFE allele frequency was outside [0, 1]; MAF left empty")

        payload = ensembl_payloads.get(rsid)
        populations: list[dict[str, Any]] | None = None
        if isinstance(payload, dict) and payload.get("__ensembl_request_failed__") is True:
            notes.append(
                "Ensembl variation request failed after retries; "
                "1000 Genomes TSI/IBS/CEU/GBR/FIN frequencies left empty"
            )
            payload = None
        record_missing = payload is None
        if record_missing:
            if not any("1000 Genomes TSI/IBS/CEU/GBR/FIN" in note for note in notes):
                notes.append(
                    "variant not found in Ensembl variation; 1000 Genomes frequencies left empty"
                )
        elif isinstance(payload, dict):
            raw_pops = payload.get("populations")
            if isinstance(raw_pops, list):
                populations = [item for item in raw_pops if isinstance(item, dict)]
            else:
                notes.append("Ensembl variation record has no populations array")
        else:
            notes.append("Ensembl variation payload was not an object; frequencies left empty")

        af_alt: dict[str, float | None] = {}
        maf_sub: dict[str, float | None] = {}
        for code in EUR_SUBPOPS:
            frequency, note = alt_allele_frequency(populations, code, ref, alt)
            if note and not record_missing:
                notes.append(note)
            af_alt[code] = frequency
            maf_sub[code] = maf_from_af(frequency)

        present = [value for value in maf_sub.values() if value is not None]
        maf_min = min(present) if present else None
        maf_max = max(present) if present else None
        maf_range = (
            (maf_max - maf_min)
            if maf_min is not None and maf_max is not None and len(present) >= 2
            else None
        )
        if len(present) < 2:
            notes.append(
                "EUR subpopulation MAF range left empty because fewer than two subpopulations had a frequency"
            )

        tsi_gap: float | None = None
        tsi_flag: bool | None = None
        if maf_sub["TSI"] is None or maf_nfe is None:
            notes.append(
                "TSI versus NFE MAF comparison left empty because one frequency is missing"
            )
        else:
            tsi_gap = maf_sub["TSI"] - maf_nfe
            tsi_flag = abs(tsi_gap) > TSI_NFE_MAF_GAP

        gene = str(variant["gene_symbol"])
        longevity_class = class_by_gene.get(gene.upper())
        if longevity_class is None:
            notes.append(
                "gene_symbol is absent from manuscript/tables/41_gene_candidate_list.csv; "
                "LongevityMap class left empty"
            )

        rare_flag: bool | None = None if maf_nfe is None else maf_nfe < RARE_MAF
        reference_power = allelic_power(maf_nfe, RANK_ODDS_RATIO, RANK_N_CASES, ALPHA_BONFERRONI)
        if maf_nfe is None:
            notes.append(
                "illustrative power left empty because gnomAD v4 NFE MAF is missing; "
                "control frequency was not imputed"
            )
        elif maf_nfe == 0.0 or maf_nfe == 1.0:
            notes.append("illustrative power left empty because NFE MAF is 0 or 1")

        rows.append(
            {
                "chrom": str(variant["chrom"]),
                "pos": int(variant["pos"]),
                "ref": ref,
                "alt": alt,
                "rsid": rsid,
                "gene_symbol": gene,
                "longevitymap_class": longevity_class,
                "af_alt_gnomad_nfe": af_nfe,
                "maf_gnomad_nfe": maf_nfe,
                "nfe_maf_below_0_01": rare_flag,
                "af_alt_1kg_tsi": af_alt["TSI"],
                "maf_1kg_tsi": maf_sub["TSI"],
                "af_alt_1kg_ibs": af_alt["IBS"],
                "maf_1kg_ibs": maf_sub["IBS"],
                "af_alt_1kg_ceu": af_alt["CEU"],
                "maf_1kg_ceu": maf_sub["CEU"],
                "af_alt_1kg_gbr": af_alt["GBR"],
                "maf_1kg_gbr": maf_sub["GBR"],
                "af_alt_1kg_fin": af_alt["FIN"],
                "maf_1kg_fin": maf_sub["FIN"],
                "eur_subpop_maf_min": maf_min,
                "eur_subpop_maf_max": maf_max,
                "eur_subpop_maf_range": maf_range,
                "tsi_minus_nfe_maf": tsi_gap,
                "tsi_nfe_maf_abs_diff_gt_0_05": tsi_flag,
                "power_n1000_or1_5_alpha_bonferroni": reference_power,
                "notes": _join_notes(notes),
            }
        )
    return rows


def assign_testability_rank(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Rank by illustrative power at 1000 cases/arm, OR 1.5, alpha 0.05/47.

    Higher power ranks first. Missing power ranks last. Ties break on higher
    NFE MAF, then rsID. This is a detectability order, not a biological one.
    """

    def sort_key(row: dict[str, Any]) -> tuple[int, float, float, str]:
        power = row["power_n1000_or1_5_alpha_bonferroni"]
        maf = row["maf_gnomad_nfe"]
        missing_power = power is None
        return (
            1 if missing_power else 0,
            0.0 if power is None else -float(power),
            0.0 if maf is None else -float(maf),
            str(row["rsid"]),
        )

    ordered = sorted(rows, key=sort_key)
    for index, row in enumerate(ordered, start=1):
        row["testability_rank"] = index
    return ordered


def build_power_grid(rows: list[dict[str, Any]]) -> pl.DataFrame:
    """Long grid of illustrative power. Control frequency is NFE MAF."""
    records: list[dict[str, Any]] = []
    alphas = (
        ("bonferroni_0.05_over_47", ALPHA_BONFERRONI),
        ("genomewide_5e-8", ALPHA_GENOMEWIDE),
    )
    for row in rows:
        for n_cases in N_CASES_PER_ARM:
            for odds_ratio in ALLELIC_ODDS_RATIOS:
                for alpha_name, alpha in alphas:
                    records.append(
                        {
                            "rsid": row["rsid"],
                            "gene_symbol": row["gene_symbol"],
                            "testability_rank": row["testability_rank"],
                            "n_cases_per_arm": n_cases,
                            "allelic_or": odds_ratio,
                            "alpha_name": alpha_name,
                            "alpha": alpha,
                            "maf_control_nfe": row["maf_gnomad_nfe"],
                            "power": allelic_power(
                                row["maf_gnomad_nfe"],
                                odds_ratio,
                                n_cases,
                                alpha,
                            ),
                        }
                    )
    return pl.DataFrame(records)


def panel_frame(rows: list[dict[str, Any]]) -> pl.DataFrame:
    """Column order for the variant-level detectability table."""
    columns = [
        "testability_rank",
        "chrom",
        "pos",
        "ref",
        "alt",
        "rsid",
        "gene_symbol",
        "longevitymap_class",
        "af_alt_gnomad_nfe",
        "maf_gnomad_nfe",
        "nfe_maf_below_0_01",
        "af_alt_1kg_tsi",
        "maf_1kg_tsi",
        "af_alt_1kg_ibs",
        "maf_1kg_ibs",
        "af_alt_1kg_ceu",
        "maf_1kg_ceu",
        "af_alt_1kg_gbr",
        "maf_1kg_gbr",
        "af_alt_1kg_fin",
        "maf_1kg_fin",
        "eur_subpop_maf_min",
        "eur_subpop_maf_max",
        "eur_subpop_maf_range",
        "tsi_minus_nfe_maf",
        "tsi_nfe_maf_abs_diff_gt_0_05",
        "power_n1000_or1_5_alpha_bonferroni",
        "notes",
    ]
    return pl.DataFrame(rows).select(columns).sort("testability_rank")


def render_power_heatmap(grid: pl.DataFrame, panel: pl.DataFrame, path: Path) -> None:
    """Variants by cohort size, one panel per allelic OR, Bonferroni alpha only."""
    bonferroni = grid.filter(pl.col("alpha_name") == "bonferroni_0.05_over_47")
    labels = [
        f"{row['gene_symbol']} {row['rsid']}"
        for row in panel.sort("testability_rank").iter_rows(named=True)
    ]
    rsids = panel.sort("testability_rank")["rsid"].to_list()
    n_variants = len(rsids)
    fig_h = max(10.0, 0.28 * n_variants + 2.2)
    fig, axes = plt.subplots(
        1,
        len(ALLELIC_ODDS_RATIOS),
        figsize=(2.4 * len(ALLELIC_ODDS_RATIOS) + 3.2, fig_h),
        sharey=True,
        layout="constrained",
    )
    image = None
    for axis, odds_ratio in zip(axes, ALLELIC_ODDS_RATIOS, strict=True):
        slice_frame = bonferroni.filter(pl.col("allelic_or") == odds_ratio)
        matrix = np.full((n_variants, len(N_CASES_PER_ARM)), np.nan)
        lookup = {
            (row["rsid"], int(row["n_cases_per_arm"])): row["power"]
            for row in slice_frame.iter_rows(named=True)
        }
        for row_index, rsid in enumerate(rsids):
            for col_index, n_cases in enumerate(N_CASES_PER_ARM):
                value = lookup.get((rsid, n_cases))
                if value is not None:
                    matrix[row_index, col_index] = float(value)
        image = axis.imshow(
            matrix,
            aspect="auto",
            interpolation="nearest",
            cmap="viridis",
            vmin=0.0,
            vmax=1.0,
        )
        axis.set_title(f"allelic OR {odds_ratio:g}")
        axis.set_xticks(np.arange(len(N_CASES_PER_ARM)))
        axis.set_xticklabels([str(n) for n in N_CASES_PER_ARM], rotation=45, ha="right")
        axis.set_xlabel("Cases per arm")
    axes[0].set_yticks(np.arange(n_variants))
    axes[0].set_yticklabels(labels, fontsize=7)
    axes[0].set_ylabel("Variant (detectability rank, top = higher illustrative power)")
    if image is not None:
        colorbar = fig.colorbar(image, ax=axes, fraction=0.03, pad=0.02)
        colorbar.set_label("Illustrative power")
    fig.suptitle(
        "Illustrative two-sided allelic-test power at alpha = 0.05/47\n"
        "Control frequency: gnomAD v4 non-Finnish European MAF",
        fontsize=11,
    )
    fig.text(0.01, 0.005, DISCLAIMER, fontsize=8, ha="left", va="bottom", wrap=True)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def release_check_argv(table: Path, report: Path) -> list[str]:
    """Argv for the installed ``rogen-release-check`` summary check."""
    return [
        "uv",
        "run",
        "rogen-release-check",
        "--table",
        str(table),
        "--table-kind",
        RELEASE_CHECK_KIND,
        "--report",
        str(report),
    ]


def run_release_check(path: Path, report: Path) -> dict[str, Any]:
    """Run ``rogen-release-check --table-kind summary`` and record the outcome."""
    argv = release_check_argv(path, report)
    if shutil.which("uv") is None:
        return {
            "path": str(path),
            "kind": RELEASE_CHECK_KIND,
            "report": str(report),
            "argv": argv,
            "returncode": None,
            "stdout": "",
            "stderr": "uv is not on PATH",
            "status": "not_run",
        }
    completed = subprocess.run(
        argv,
        cwd=REPO_ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    status = "passed" if completed.returncode == 0 else "failed"
    if (
        completed.returncode == 127
        or "No such file" in completed.stderr
        or "Failed to spawn" in completed.stderr
    ):
        status = "command_not_found"
    return {
        "path": str(path),
        "kind": RELEASE_CHECK_KIND,
        "report": str(report),
        "argv": argv,
        "returncode": completed.returncode,
        "stdout": completed.stdout,
        "stderr": completed.stderr,
        "status": status,
    }


def write_provenance(
    path: Path,
    *,
    started_at: str,
    gnomad_counts: dict[str, int],
    ensembl_counts: dict[str, int],
    ensembl_release: int,
    ensembl_host: str,
    gnomad_cache: Path,
    ensembl_cache: Path,
    variants_path: Path,
    gene_list_path: Path,
    release_checks: list[dict[str, Any]],
) -> None:
    """Record data sources, cache use, and the release-check result."""
    payload = {
        "generated_at_utc": started_at,
        "disclaimer": DISCLAIMER,
        "input_variants": {"path": str(variants_path), "n_rows": EXPECTED_VARIANT_ROWS},
        "longevitymap_annotation": {
            "database": "LongevityMap",
            "table": str(gene_list_path),
            "join": "gene_symbol to Gene_Symbol",
            "queried_live": False,
            "role": "gene-level class annotation only",
        },
        "gnomad": {
            "version": "gnomAD v4",
            "dataset": GNOMAD_DATASET,
            "population": GNOMAD_POPULATION,
            "api": GNOMAD_API_URL,
            "cache_path": str(gnomad_cache),
            "cache_hits": gnomad_counts["cache_hits"],
            "new_queries": gnomad_counts["new_queries"],
            "absent_from_cache": gnomad_counts["absent_from_cache"],
            "requery": False,
            "query_date_utc": None,
        },
        "ensembl": {
            "release": ensembl_release,
            "assembly": ENSEMBL_ASSEMBLY,
            "host": ensembl_host,
            "info_data": "GET /info/data",
            "endpoint": "GET /variation/human/{id}?pops=1",
            "populations": [f"1000GENOMES:phase_3:{code}" for code in EUR_SUBPOPS],
            "cache_path": str(ensembl_cache),
            "cache_hits": ensembl_counts["cache_hits"],
            "new_queries": ensembl_counts["new_queries"],
            "request_failures": ensembl_counts["request_failures"],
            "query_date_utc": started_at if ensembl_counts["new_queries"] else None,
        },
        "frequency_policy": (
            "gnomAD v4 NFE and each 1000 Genomes phase 3 EUR subpopulation are stored "
            "separately. They are not averaged, and none is labeled as a Romanian frequency."
        ),
        "power": {
            "design": "1:1 case-control, two-sided allelic z-test, normal approximation",
            "control_frequency": "gnomAD v4 NFE minor-allele frequency",
            "n_cases_per_arm": list(N_CASES_PER_ARM),
            "allelic_odds_ratios": list(ALLELIC_ODDS_RATIOS),
            "alpha_bonferroni_0.05_over_47": ALPHA_BONFERRONI,
            "alpha_genomewide": ALPHA_GENOMEWIDE,
            "rank_cell": {
                "n_cases_per_arm": RANK_N_CASES,
                "allelic_or": RANK_ODDS_RATIO,
                "alpha": ALPHA_BONFERRONI,
            },
        },
        "release_check": release_checks,
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def main(
    variants: Path = typer.Option(
        DEFAULT_VARIANTS,
        "--variants",
        help="47 prioritized GRCh38 variants (chrom, pos, ref, alt, rsid, gene_symbol).",
    ),
    gene_list: Path = typer.Option(
        DEFAULT_GENE_LIST,
        "--gene-list",
        help="Manuscript gene table with LongevityMap Longevity_Class. Annotation only.",
    ),
    gnomad_cache: Path = typer.Option(
        DEFAULT_GNOMAD_CACHE,
        "--gnomad-cache",
        help="gnomAD v4 NFE JSON cache. Read only; the API is not queried.",
    ),
    ensembl_cache: Path = typer.Option(
        DEFAULT_ENSEMBL_CACHE,
        "--ensembl-cache",
        help="JSON directory cache for Ensembl variation responses.",
    ),
    out_dir: Path = typer.Option(
        DEFAULT_OUT_DIR,
        "--out-dir",
        help="Directory for the panel table, power grid, heatmap, and provenance JSON.",
    ),
) -> None:
    """Build the detectability table, power grid, heatmap, and provenance record."""
    started_at = datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")
    variant_frame = load_variants(variants)
    classes = load_longevitymap_classes(gene_list)
    gnomad_cache_map, gnomad_counts = fetch_nfe_frequencies(variant_frame, gnomad_cache)
    rsids = [normalize_rsid(str(rsid)) for rsid in variant_frame["rsid"].to_list()]
    ensembl_payloads, ensembl_counts, ensembl_release, ensembl_host = fetch_ensembl_variations(
        rsids, ensembl_cache
    )
    rows = assign_testability_rank(
        build_panel_rows(variant_frame, classes, gnomad_cache_map, ensembl_payloads)
    )
    panel = panel_frame(rows)
    grid = build_power_grid(rows)
    out_dir.mkdir(parents=True, exist_ok=True)
    panel_path = out_dir / "ro_longevity_panel.csv"
    grid_path = out_dir / "ro_longevity_power_grid.csv"
    figure_path = out_dir / "ro_longevity_power_heatmap.png"
    provenance_path = out_dir / "ro_longevity_panel_provenance.json"
    panel.write_csv(panel_path)
    grid.write_csv(grid_path)
    render_power_heatmap(grid, panel, figure_path)
    report_dir = out_dir / "release_checks"
    release_checks = [
        run_release_check(panel_path, report_dir / "ro_longevity_panel.json"),
        run_release_check(grid_path, report_dir / "ro_longevity_power_grid.json"),
    ]
    write_provenance(
        provenance_path,
        started_at=started_at,
        gnomad_counts=gnomad_counts,
        ensembl_counts=ensembl_counts,
        ensembl_release=ensembl_release,
        ensembl_host=ensembl_host,
        gnomad_cache=gnomad_cache,
        ensembl_cache=ensembl_cache,
        variants_path=variants,
        gene_list_path=gene_list,
        release_checks=release_checks,
    )
    typer.echo(f"Wrote {panel_path}")
    typer.echo(f"Wrote {grid_path}")
    typer.echo(f"Wrote {figure_path}")
    typer.echo(f"Wrote {provenance_path}")
    for check in release_checks:
        typer.echo(
            f"rogen-release-check kind={RELEASE_CHECK_KIND} {check['path']}: {check['status']} "
            f"returncode={check['returncode']}"
        )
    if any(check["status"] != "passed" for check in release_checks):
        _stop(
            "STOP: rogen-release-check did not pass for every output table. "
            "The result is recorded in the provenance JSON."
        )


if __name__ == "__main__":
    typer.run(main)
