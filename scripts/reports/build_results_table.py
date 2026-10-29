#!/usr/bin/env python3
"""Build the Results markdown table for the October technical report.

Each row is ``(label, file, locator)``. Locators:

- JSON: dotted path, e.g. ``snv.f1`` or ``sha256.reference``.
  List elements: ``levels[level=hannum_nonzero].pearson_r``.
- CSV: ``filter_col=value,...:column``, e.g.
  ``clock=horvath2013,cohort=GSE87571:mae``.
- TSV: same as CSV (detected by suffix).
- Special: ``row_count``, ``count_where:col=value``, ``sum:column``.

Exits non-zero if any file, key, or filtered row is missing.
"""

from __future__ import annotations

import csv
import json
import re
from pathlib import Path
from typing import Any

import typer

from rogen_aging.config import find_repo_root

REPO_ROOT = find_repo_root()

# (label, relative path, locator)
SPECS: list[tuple[str, str, str]] = [
    ("HG002 region", "analysis/dryrun/hg002_shortread_metrics.json", "region"),
    ("Mean depth", "analysis/dryrun/hg002_shortread_metrics.json", "mean_depth"),
    ("Percent mapped", "analysis/dryrun/hg002_shortread_metrics.json", "percent_mapped"),
    (
        "Percent duplicates",
        "analysis/dryrun/hg002_shortread_metrics.json",
        "percent_duplicates",
    ),
    ("SNV TP", "analysis/dryrun/hg002_shortread_metrics.json", "snv.TP"),
    ("SNV FP", "analysis/dryrun/hg002_shortread_metrics.json", "snv.FP"),
    ("SNV FN", "analysis/dryrun/hg002_shortread_metrics.json", "snv.FN"),
    ("SNV precision", "analysis/dryrun/hg002_shortread_metrics.json", "snv.precision"),
    ("SNV recall", "analysis/dryrun/hg002_shortread_metrics.json", "snv.recall"),
    ("SNV F1", "analysis/dryrun/hg002_shortread_metrics.json", "snv.f1"),
    ("Indel TP", "analysis/dryrun/hg002_shortread_metrics.json", "indel.TP"),
    ("Indel FP", "analysis/dryrun/hg002_shortread_metrics.json", "indel.FP"),
    ("Indel FN", "analysis/dryrun/hg002_shortread_metrics.json", "indel.FN"),
    ("Indel precision", "analysis/dryrun/hg002_shortread_metrics.json", "indel.precision"),
    ("Indel recall", "analysis/dryrun/hg002_shortread_metrics.json", "indel.recall"),
    ("Indel F1", "analysis/dryrun/hg002_shortread_metrics.json", "indel.f1"),
    (
        "Short-read metrics time",
        "analysis/dryrun/hg002_shortread_metrics.json",
        "generated_at_utc",
    ),
    (
        "Short-read metrics git commit",
        "analysis/dryrun/hg002_shortread_metrics.json",
        "git_commit",
    ),
    ("Benchmark VCF", "analysis/dryrun/hg002_shortread_metrics.json", "benchmark_vcf"),
    ("Benchmark BED", "analysis/dryrun/hg002_shortread_metrics.json", "benchmark_bed"),
    (
        "Reference sha256",
        "analysis/dryrun/hg002_shortread_metrics.json",
        "sha256.reference",
    ),
    (
        "Benchmark VCF sha256",
        "analysis/dryrun/hg002_shortread_metrics.json",
        "sha256.benchmark_vcf",
    ),
    (
        "Benchmark BED sha256",
        "analysis/dryrun/hg002_shortread_metrics.json",
        "sha256.benchmark_bed",
    ),
    (
        "Wall time fetch (s)",
        "analysis/dryrun/hg002_shortread_metrics.json",
        "wall_time_seconds.fetch",
    ),
    (
        "Wall time intake (s)",
        "analysis/dryrun/hg002_shortread_metrics.json",
        "wall_time_seconds.intake",
    ),
    (
        "Wall time align (s)",
        "analysis/dryrun/hg002_shortread_metrics.json",
        "wall_time_seconds.align",
    ),
    (
        "Wall time call (s)",
        "analysis/dryrun/hg002_shortread_metrics.json",
        "wall_time_seconds.call",
    ),
    (
        "Wall time benchmark (s)",
        "analysis/dryrun/hg002_shortread_metrics.json",
        "wall_time_seconds.benchmark",
    ),
    (
        "ONT all_cpgs n",
        "analysis/dryrun/hg002_ont_methylation_metrics.json",
        "levels[level=all_cpgs].n_cpgs",
    ),
    (
        "ONT all_cpgs depth≥10",
        "analysis/dryrun/hg002_ont_methylation_metrics.json",
        "levels[level=all_cpgs].n_ont_depth_ge_min",
    ),
    (
        "ONT all_cpgs Pearson r",
        "analysis/dryrun/hg002_ont_methylation_metrics.json",
        "levels[level=all_cpgs].pearson_r",
    ),
    (
        "ONT all_cpgs mean ONT−bisulfite",
        "analysis/dryrun/hg002_ont_methylation_metrics.json",
        "levels[level=all_cpgs].mean_ont_minus_bisulfite",
    ),
    (
        "ONT hm450 n",
        "analysis/dryrun/hg002_ont_methylation_metrics.json",
        "levels[level=hm450].n_cpgs",
    ),
    (
        "ONT hm450 depth≥10",
        "analysis/dryrun/hg002_ont_methylation_metrics.json",
        "levels[level=hm450].n_ont_depth_ge_min",
    ),
    (
        "ONT hm450 Pearson r",
        "analysis/dryrun/hg002_ont_methylation_metrics.json",
        "levels[level=hm450].pearson_r",
    ),
    (
        "ONT hm450 mean diff",
        "analysis/dryrun/hg002_ont_methylation_metrics.json",
        "levels[level=hm450].mean_ont_minus_bisulfite",
    ),
    (
        "ONT hannum_nonzero n",
        "analysis/dryrun/hg002_ont_methylation_metrics.json",
        "levels[level=hannum_nonzero].n_cpgs",
    ),
    (
        "ONT hannum_nonzero depth≥10",
        "analysis/dryrun/hg002_ont_methylation_metrics.json",
        "levels[level=hannum_nonzero].n_ont_depth_ge_min",
    ),
    (
        "ONT hannum_nonzero Pearson r",
        "analysis/dryrun/hg002_ont_methylation_metrics.json",
        "levels[level=hannum_nonzero].pearson_r",
    ),
    (
        "ONT hannum_nonzero mean diff",
        "analysis/dryrun/hg002_ont_methylation_metrics.json",
        "levels[level=hannum_nonzero].mean_ont_minus_bisulfite",
    ),
    ("modkit", "analysis/dryrun/hg002_ont_methylation_metrics.json", "modkit"),
    (
        "Panel variants",
        "analysis/panel/ro_longevity_panel_provenance.json",
        "input_variants.n_rows",
    ),
    (
        "Panel provenance time",
        "analysis/panel/ro_longevity_panel_provenance.json",
        "generated_at_utc",
    ),
    (
        "Ensembl release",
        "analysis/panel/ro_longevity_panel_provenance.json",
        "ensembl.release",
    ),
    (
        "Ensembl host",
        "analysis/panel/ro_longevity_panel_provenance.json",
        "ensembl.host",
    ),
    (
        "Ensembl cache hits",
        "analysis/panel/ro_longevity_panel_provenance.json",
        "ensembl.cache_hits",
    ),
    (
        "Ensembl new queries",
        "analysis/panel/ro_longevity_panel_provenance.json",
        "ensembl.new_queries",
    ),
    (
        "Ensembl request failures",
        "analysis/panel/ro_longevity_panel_provenance.json",
        "ensembl.request_failures",
    ),
    (
        "gnomAD cache hits",
        "analysis/panel/ro_longevity_panel_provenance.json",
        "gnomad.cache_hits",
    ),
    (
        "gnomAD new queries",
        "analysis/panel/ro_longevity_panel_provenance.json",
        "gnomad.new_queries",
    ),
    (
        "TSI vs NFE n tested",
        "analysis/panel/ro_longevity_panel_provenance.json",
        "tsi_vs_nfe.n_tests",
    ),
    (
        "TSI vs NFE n_chi_square",
        "analysis/panel/ro_longevity_panel_provenance.json",
        "tsi_vs_nfe.n_chi_square",
    ),
    (
        "TSI vs NFE n_fisher",
        "analysis/panel/ro_longevity_panel_provenance.json",
        "tsi_vs_nfe.n_fisher",
    ),
    (
        "TSI vs NFE q_below_fdr05 count",
        "analysis/panel/ro_longevity_panel_provenance.json",
        "len:tsi_vs_nfe.q_below_fdr05",
    ),
    (
        "Panel release check status (first)",
        "analysis/panel/ro_longevity_panel_provenance.json",
        "release_check[0].status",
    ),
    (
        "Sanity flag_count",
        "analysis/panel/ro_longevity_panel_provenance.json",
        "sanity_check.flag_count",
    ),
    (
        "GWAS n_resolved",
        "analysis/longevity/longevity_gwas_lookup_provenance.json",
        "panel.n_resolved",
    ),
    (
        "GWAS n_not_tested",
        "analysis/longevity/longevity_gwas_lookup_provenance.json",
        "panel.n_not_tested",
    ),
    (
        "GWAS n_tests",
        "analysis/longevity/longevity_gwas_lookup_provenance.json",
        "multiple_testing.n_tests",
    ),
    (
        "GWAS n_p_lt_0_05",
        "analysis/longevity/longevity_gwas_lookup_provenance.json",
        "multiple_testing.n_p_lt_0_05",
    ),
    (
        "GWAS n_expected_p_lt_0_05_under_null",
        "analysis/longevity/longevity_gwas_lookup_provenance.json",
        "multiple_testing.n_expected_p_lt_0_05_under_null",
    ),
    (
        "GWAS n_q_bh_lt_0_05",
        "analysis/longevity/longevity_gwas_lookup_provenance.json",
        "multiple_testing.n_q_bh_lt_0_05",
    ),
    (
        "GWAS n_bonferroni_sig",
        "analysis/longevity/longevity_gwas_lookup_provenance.json",
        "multiple_testing.n_bonferroni_sig",
    ),
    (
        "GWAS lookup time",
        "analysis/longevity/longevity_gwas_lookup_provenance.json",
        "generated_at_utc",
    ),
    (
        "GWAS release check",
        "analysis/longevity/longevity_gwas_lookup_provenance.json",
        "release_check.status",
    ),
    (
        "Integrated n_variants",
        "analysis/longevity/integrated_evidence_provenance.json",
        "n_variants",
    ),
    (
        "Integrated GWAS support count",
        "analysis/longevity/integrated_evidence_provenance.json",
        "label_counts.GWAS support",
    ),
    (
        "Integrated eQTL support count",
        "analysis/longevity/integrated_evidence_provenance.json",
        "label_counts.eQTL support",
    ),
    (
        "Integrated population difference count",
        "analysis/longevity/integrated_evidence_provenance.json",
        "label_counts.population difference",
    ),
    (
        "Integrated eQTL input status",
        "analysis/longevity/integrated_evidence_provenance.json",
        "inputs.eqtl_input_status",
    ),
    (
        "eQTL datasets locked",
        "analysis/longevity/integrated_evidence_provenance.json",
        "inputs.eqtl_datasets_used.n_datasets",
    ),
    (
        "Integrated eqtl_n_ambiguous (sum)",
        "analysis/longevity/integrated_evidence.csv",
        "sum:eqtl_n_ambiguous",
    ),
    (
        "Integrated eqtl_n_flipped (sum)",
        "analysis/longevity/integrated_evidence.csv",
        "sum:eqtl_n_flipped",
    ),
    (
        "Integrated release check",
        "analysis/longevity/integrated_evidence_provenance.json",
        "release_check.status",
    ),
    (
        "Integrated time",
        "analysis/longevity/integrated_evidence_provenance.json",
        "generated_at_utc",
    ),
    (
        "Clocks benchmark generated_at_utc",
        "analysis/clocks_benchmark/provenance.json",
        "generated_at_utc",
    ),
    (
        "Clocks benchmark finished_at_utc",
        "analysis/clocks_benchmark/provenance.json",
        "finished_at_utc",
    ),
    (
        "pyaging version",
        "analysis/clocks_benchmark/provenance.json",
        "pyaging_version",
    ),
    (
        "CpG columns used",
        "analysis/clocks_benchmark/provenance.json",
        "cpg_restriction.n_probes",
    ),
    (
        "GSE87571 horvath2013 MAE",
        "analysis/clocks_benchmark/metrics_by_clock_cohort.csv",
        "clock=horvath2013,cohort=GSE87571:mae",
    ),
    (
        "GSE87571 horvath2013 r",
        "analysis/clocks_benchmark/metrics_by_clock_cohort.csv",
        "clock=horvath2013,cohort=GSE87571:pearson_r",
    ),
    (
        "GSE87571 horvath2013 mean signed error",
        "analysis/clocks_benchmark/metrics_by_clock_cohort.csv",
        "clock=horvath2013,cohort=GSE87571:mean_signed_error",
    ),
    (
        "GSE87571 hannum MAE",
        "analysis/clocks_benchmark/metrics_by_clock_cohort.csv",
        "clock=hannum,cohort=GSE87571:mae",
    ),
    (
        "GSE87571 hannum r",
        "analysis/clocks_benchmark/metrics_by_clock_cohort.csv",
        "clock=hannum,cohort=GSE87571:pearson_r",
    ),
    (
        "GSE87571 hannum mean signed error",
        "analysis/clocks_benchmark/metrics_by_clock_cohort.csv",
        "clock=hannum,cohort=GSE87571:mean_signed_error",
    ),
    (
        "GSE87571 dnamphenoage MAE",
        "analysis/clocks_benchmark/metrics_by_clock_cohort.csv",
        "clock=dnamphenoage,cohort=GSE87571:mae",
    ),
    (
        "GSE87571 dnamphenoage r",
        "analysis/clocks_benchmark/metrics_by_clock_cohort.csv",
        "clock=dnamphenoage,cohort=GSE87571:pearson_r",
    ),
    (
        "GSE87571 dnamphenoage mean signed error",
        "analysis/clocks_benchmark/metrics_by_clock_cohort.csv",
        "clock=dnamphenoage,cohort=GSE87571:mean_signed_error",
    ),
    (
        "GSE87571 dunedinpace r",
        "analysis/clocks_benchmark/metrics_by_clock_cohort.csv",
        "clock=dunedinpace,cohort=GSE87571:pearson_r",
    ),
    (
        "GSE87571 grimage2 MAE",
        "analysis/clocks_benchmark/metrics_by_clock_cohort.csv",
        "clock=grimage2,cohort=GSE87571:mae",
    ),
    (
        "GSE87571 grimage2 r",
        "analysis/clocks_benchmark/metrics_by_clock_cohort.csv",
        "clock=grimage2,cohort=GSE87571:pearson_r",
    ),
    (
        "GSE87571 grimage2 mean signed error",
        "analysis/clocks_benchmark/metrics_by_clock_cohort.csv",
        "clock=grimage2,cohort=GSE87571:mean_signed_error",
    ),
    (
        "GSE87571 rogen_hannum_retrain MAE",
        "analysis/clocks_benchmark/metrics_by_clock_cohort.csv",
        "clock=rogen_hannum_retrain,cohort=GSE87571:mae",
    ),
    (
        "GSE87571 rogen_hannum_retrain r",
        "analysis/clocks_benchmark/metrics_by_clock_cohort.csv",
        "clock=rogen_hannum_retrain,cohort=GSE87571:pearson_r",
    ),
    (
        "GSE87571 rogen_hannum_retrain mean signed error",
        "analysis/clocks_benchmark/metrics_by_clock_cohort.csv",
        "clock=rogen_hannum_retrain,cohort=GSE87571:mean_signed_error",
    ),
    (
        "GSE40279 held-out rogen_hannum_retrain MAE",
        "analysis/clocks_benchmark/metrics_by_clock_cohort.csv",
        "clock=rogen_hannum_retrain,cohort=GSE40279_heldout:mae",
    ),
    (
        "GSE40279 held-out rogen_hannum_retrain r",
        "analysis/clocks_benchmark/metrics_by_clock_cohort.csv",
        "clock=rogen_hannum_retrain,cohort=GSE40279_heldout:pearson_r",
    ),
    (
        "GSE40279 held-out rogen_hannum_retrain mean signed error",
        "analysis/clocks_benchmark/metrics_by_clock_cohort.csv",
        "clock=rogen_hannum_retrain,cohort=GSE40279_heldout:mean_signed_error",
    ),
    (
        "GSE87571 horvath2013 sex coef",
        "analysis/clocks_benchmark/sex_acceleration_models.csv",
        "clock=horvath2013,cohort=GSE87571:sex_coef_female_minus_male",
    ),
    (
        "GSE87571 horvath2013 sex Bonf p",
        "analysis/clocks_benchmark/sex_acceleration_models.csv",
        "clock=horvath2013,cohort=GSE87571:sex_p_bonferroni",
    ),
    (
        "GSE87571 hannum sex coef",
        "analysis/clocks_benchmark/sex_acceleration_models.csv",
        "clock=hannum,cohort=GSE87571:sex_coef_female_minus_male",
    ),
    (
        "GSE87571 hannum sex Bonf p",
        "analysis/clocks_benchmark/sex_acceleration_models.csv",
        "clock=hannum,cohort=GSE87571:sex_p_bonferroni",
    ),
    (
        "GSE87571 dnamphenoage sex coef",
        "analysis/clocks_benchmark/sex_acceleration_models.csv",
        "clock=dnamphenoage,cohort=GSE87571:sex_coef_female_minus_male",
    ),
    (
        "GSE87571 dnamphenoage sex Bonf p",
        "analysis/clocks_benchmark/sex_acceleration_models.csv",
        "clock=dnamphenoage,cohort=GSE87571:sex_p_bonferroni",
    ),
    (
        "GSE87571 grimage2 sex coef",
        "analysis/clocks_benchmark/sex_acceleration_models.csv",
        "clock=grimage2,cohort=GSE87571:sex_coef_female_minus_male",
    ),
    (
        "GSE87571 grimage2 sex Bonf p",
        "analysis/clocks_benchmark/sex_acceleration_models.csv",
        "clock=grimage2,cohort=GSE87571:sex_p_bonferroni",
    ),
    (
        "GSE87571 rogen_hannum_retrain sex coef",
        "analysis/clocks_benchmark/sex_acceleration_models.csv",
        "clock=rogen_hannum_retrain,cohort=GSE87571:sex_coef_female_minus_male",
    ),
    (
        "GSE87571 rogen_hannum_retrain sex Bonf p",
        "analysis/clocks_benchmark/sex_acceleration_models.csv",
        "clock=rogen_hannum_retrain,cohort=GSE87571:sex_p_bonferroni",
    ),
    (
        "v2 alt_status=resolved count",
        "data/processed/variants_47_input_v2.csv",
        "count_where:alt_status=resolved",
    ),
    (
        "Allele audit keep count",
        "analysis/panel/tsi_nfe_allele_audit.tsv",
        "count_where:action=keep",
    ),
    (
        "Sanity TSV row count",
        "analysis/panel/ro_longevity_panel_sanity.tsv",
        "row_count",
    ),
    (
        "Sanity complement_skipped count",
        "analysis/panel/ro_longevity_panel_sanity.tsv",
        "count_where:complement_skipped=true",
    ),
]

LIST_FILTER_RE = re.compile(r"^([^.\[]+)\[([^=\]]+)=([^\]]+)\](?:\.(.+))?$")
LIST_INDEX_RE = re.compile(r"^([^.\[]+)\[(\d+)\](?:\.(.+))?$")


def _format_value(value: Any) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        return repr(value) if value != int(value) else str(int(value))
    if value is None:
        return "null"
    return str(value)


def _walk_json(node: Any, parts: list[str]) -> Any:
    cur = node
    i = 0
    while i < len(parts):
        part = parts[i]
        if part == "#len":
            if not isinstance(cur, list):
                raise KeyError("#len requires a list")
            return len(cur)
        filt = LIST_FILTER_RE.match(".".join(parts[i:]))
        if filt and filt.group(0).startswith(part + "["):
            # handled below via full-locator match; keep simple walk for dotted keys
            pass
        idx = LIST_INDEX_RE.match(part) if False else None
        _ = idx
        if "[" in part and part.endswith("]"):
            raise KeyError(f"use full-locator list syntax for {part}")
        if isinstance(cur, dict):
            if part not in cur:
                raise KeyError(part)
            cur = cur[part]
        else:
            raise KeyError(part)
        i += 1
    return cur


def _resolve_json(data: Any, locator: str) -> Any:
    if locator.startswith("len:"):
        node = _resolve_json(data, locator.removeprefix("len:"))
        if not isinstance(node, list):
            raise KeyError("len: requires a list")
        return len(node)
    match = LIST_FILTER_RE.match(locator)
    if match:
        list_key, filter_key, filter_val, rest = match.groups()
        items = data[list_key]
        if not isinstance(items, list):
            raise KeyError(f"{list_key} is not a list")
        selected = [item for item in items if str(item.get(filter_key)) == filter_val]
        if not selected:
            raise KeyError(f"no list element with {filter_key}={filter_val}")
        if len(selected) > 1:
            raise KeyError(f"multiple list elements with {filter_key}={filter_val}")
        node: Any = selected[0]
        if rest:
            return _resolve_json(node, rest) if "[" in rest else _walk_json(node, rest.split("."))
        return node
    idx_match = LIST_INDEX_RE.match(locator)
    if idx_match:
        list_key, index_s, rest = idx_match.groups()
        items = data[list_key]
        if not isinstance(items, list):
            raise KeyError(f"{list_key} is not a list")
        index = int(index_s)
        if index >= len(items):
            raise KeyError(f"{list_key}[{index}] out of range")
        node = items[index]
        if rest:
            return _resolve_json(node, rest) if "[" in rest else _walk_json(node, rest.split("."))
        return node
    return _walk_json(data, locator.split("."))


def _read_table(path: Path) -> list[dict[str, str]]:
    sep = "\t" if path.suffix.lower() == ".tsv" else ","
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter=sep))


def _resolve_table(rows: list[dict[str, str]], locator: str) -> Any:
    if locator == "row_count":
        return len(rows)
    if locator.startswith("sum:"):
        col = locator.removeprefix("sum:")
        total = 0.0
        for row in rows:
            raw = row.get(col, "")
            if raw == "" or raw is None:
                continue
            total += float(raw)
        return total
    if locator.startswith("count_where:"):
        clause = locator.removeprefix("count_where:")
        key, _, want = clause.partition("=")
        return sum(1 for row in rows if str(row.get(key, "")).lower() == want.lower())
    if ":" not in locator:
        raise KeyError(f"bad CSV locator: {locator}")
    filter_part, _, column = locator.rpartition(":")
    filters: dict[str, str] = {}
    for piece in filter_part.split(","):
        k, _, v = piece.partition("=")
        filters[k] = v
    matched = [row for row in rows if all(str(row.get(k, "")) == v for k, v in filters.items())]
    if not matched:
        raise KeyError(f"no CSV row matching {filters}")
    if len(matched) > 1:
        raise KeyError(f"multiple CSV rows matching {filters}")
    if column not in matched[0]:
        raise KeyError(column)
    return matched[0][column]


def resolve_one(rel_path: str, locator: str) -> Any:
    """Resolve one locator against a repo-relative file."""
    path = REPO_ROOT / rel_path
    if not path.is_file():
        raise FileNotFoundError(rel_path)
    if path.suffix.lower() == ".json":
        return _resolve_json(json.loads(path.read_text(encoding="utf-8")), locator)
    if path.suffix.lower() in {".csv", ".tsv"}:
        return _resolve_table(_read_table(path), locator)
    raise ValueError(f"unsupported file type: {rel_path}")


def collect() -> tuple[list[tuple[str, str, str, str]], list[str]]:
    """Return resolved rows and error messages."""
    rows: list[tuple[str, str, str, str]] = []
    errors: list[str] = []
    for label, rel, locator in SPECS:
        try:
            value = resolve_one(rel, locator)
            rows.append((label, _format_value(value), rel, locator))
        except (FileNotFoundError, KeyError, ValueError, TypeError, json.JSONDecodeError) as exc:
            errors.append(f"{label} | {rel} | {locator} -> {exc}")
    return rows, errors


# Report display overrides (raw value still validated from the file/locator).
DISPLAY_OVERRIDES: dict[str, str] = {
    "ONT hannum_nonzero Pearson r": "n = 4, too few for a correlation",
}


def render_markdown(rows: list[tuple[str, str, str, str]]) -> str:
    """Markdown Results table with Source = file + locator."""
    lines = [
        "| Item | Value | Source |",
        "|------|-------|--------|",
    ]
    for label, value, rel, locator in rows:
        shown = DISPLAY_OVERRIDES.get(label, value)
        source = f"`{rel}` `{locator}`"
        lines.append(f"| {label} | {shown} | {source} |")
    return "\n".join(lines) + "\n"


def main(
    out: Path | None = typer.Option(
        None,
        "--out",
        help="Write the markdown table here (default: stdout).",
    ),
    list_values: bool = typer.Option(
        False,
        "--list-values",
        help="Print label, value, file, locator (one per line) for spot-checks.",
    ),
) -> None:
    """Build the Results table; exit 1 if any locator fails."""
    rows, errors = collect()
    if errors:
        for err in errors:
            typer.echo(err, err=True)
        raise typer.Exit(code=1)
    if list_values:
        for label, value, rel, locator in rows:
            typer.echo(f"{value}\t{rel}\t{locator}\t{label}")
        return
    text = render_markdown(rows)
    if out is None:
        typer.echo(text, nl=False)
    else:
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(text, encoding="utf-8")
        typer.echo(f"Wrote {out}")


if __name__ == "__main__":
    typer.run(main)
