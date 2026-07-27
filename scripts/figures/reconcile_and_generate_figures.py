#!/usr/bin/env python3
"""Audit legacy gene nomenclature and generate manuscript AF / network figures.

Loads Supplementary Table 3 and manuscript text, reconciles historical CETP and
HLA locus labels to HGNC symbols, cross-checks the 41-gene candidate list, then
writes:

* ``outputs/nomenclature_audit.log``
* ``outputs/figures/Figure_AF_Scatter.{pdf,png}``
* ``outputs/figures/Figure_41_Gene_Network.{pdf,png}``

Example:
    uv run python reconcile_and_generate_figures.py
    uv run python reconcile_and_generate_figures.py --af-csv analysis/rogen_vs_gnomad_af.csv
"""

from __future__ import annotations

import logging
import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Annotated

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import typer
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

REPO_ROOT = Path(__file__).resolve().parent.parent.parent

DEFAULT_SUPP_TABLE = REPO_ROOT / "manuscript" / "tables" / "Supplementary_Table_3.xlsx"
DEFAULT_CANDIDATE_LIST = REPO_ROOT / "manuscript" / "tables" / "41_gene_candidate_list.csv"
DEFAULT_TEXT_DIR = REPO_ROOT / "manuscript" / "text"
DEFAULT_AF_CSV = REPO_ROOT / "analysis" / "rogen_vs_gnomad_af.csv"
DEFAULT_NETWORK_CSV = REPO_ROOT / "data" / "network" / "41_gene_interactions.csv"
DEFAULT_NODES_CSV = REPO_ROOT / "data" / "network" / "41_gene_nodes.csv"
DEFAULT_AUDIT_LOG = REPO_ROOT / "outputs" / "nomenclature_audit.log"
DEFAULT_FIG_DIR = REPO_ROOT / "outputs" / "figures"

DELTA_AF_THRESHOLD = 0.05
FIGURE_DPI = 300
STRING_SCORE_MIN = 0.15

# Historical CETP / HLA (and related HSP) locus identifiers → HGNC symbols.
LEGACY_TO_HGNC: dict[str, str] = {
    "cetp": "CETP",
    "cholesteryl ester transfer protein": "CETP",
    "cholesteryl ester transfer protein (cetp)": "CETP",
    "cetp gene": "CETP",
    "cetp locus": "CETP",
    "dqb1": "HLA-DQB1",
    "dqb103": "HLA-DQB1",
    "dqb105": "HLA-DQB1",
    "dqb1*03": "HLA-DQB1",
    "dqb1*05": "HLA-DQB1",
    "hla-dqb1*03": "HLA-DQB1",
    "hla-dqb1*05": "HLA-DQB1",
    "hla-dq": "HLA-DQB1",
    "hla class ii dq beta": "HLA-DQB1",
    "hla class ii dq beta chain": "HLA-DQB1",
    "hsp70-1": "HSPA1A",
    "hsp70-1a": "HSPA1A",
    "hsp70-1b": "HSPA1B",
    "hsp70-2": "HSPA1B",
    "hsp70-hom": "HSPA1L",
}

# Phrase / token patterns scanned in manuscript text (longest first).
LEGACY_TEXT_PATTERNS: tuple[tuple[str, str], ...] = tuple(
    sorted(
        (
            (r"cholesteryl ester transfer protein\s*\(CETP\)", "CETP"),
            (r"cholesteryl ester transfer protein", "CETP"),
            (r"HLA class II DQ beta chain", "HLA-DQB1"),
            (r"HLA class II DQ beta", "HLA-DQB1"),
            (r"(?<!HLA-)(?<![\w-])DQB1(?![\w*])", "HLA-DQB1"),
            (r"\bDQB103\b", "HLA-DQB1"),
            (r"\bDQB105\b", "HLA-DQB1"),
            (r"DQB1\*03", "HLA-DQB1"),
            (r"DQB1\*05", "HLA-DQB1"),
            (r"HLA-DQB1\*03", "HLA-DQB1"),
            (r"HLA-DQB1\*05", "HLA-DQB1"),
            (r"\bHLA-DQ\b", "HLA-DQB1"),
            (r"\bCETP\b", "CETP"),
            (r"\bHSP70-1A\b", "HSPA1A"),
            (r"\bHSP70-1B\b", "HSPA1B"),
            (r"\bHSP70-Hom\b", "HSPA1L"),
            (r"\bHSP70-1\b", "HSPA1A"),
            (r"\bHSP70-2\b", "HSPA1B"),
        ),
        key=lambda item: -len(item[0]),
    )
)

LONGEVITY_COLORS: dict[str, str] = {
    "Pro-Longevity": "#2A6F97",
    "Anti-Longevity": "#C45C3E",
    "Context-Dependent": "#6B7280",
}

CLUSTER_LAYOUT_ORDER: tuple[str, ...] = (
    "Protein Folding/HSP",
    "Lipid Metabolism",
    "Neuroinflammation",
    "Mitochondrial Function",
    "Neuronal Signaling",
    "Transport/Metabolism",
    "Cell Adhesion/Structure",
    "Growth Signaling",
    "Other/Context",
)

app = typer.Typer(
    add_completion=False,
    help="Reconcile legacy gene nomenclature and generate manuscript figures.",
)


@dataclass
class AuditFinding:
    """Single audit finding recorded during nomenclature / concordance checks.

    Attributes:
        severity: Finding level such as ``OK``, ``WARNING``, ``DISCREPANCY``,
            or ``ERROR``.
        category: Short machine-readable category (e.g. ``gene_symbols``).
        detail: Human-readable description of the finding.
    """

    severity: str
    category: str
    detail: str


@dataclass
class AuditReport:
    """Mutable container for reconciliations and concordance findings.

    Attributes:
        findings: Ordered list of ``AuditFinding`` objects.
        reconciliations: Human-readable legacy → HGNC mapping notes.
    """

    findings: list[AuditFinding] = field(default_factory=list)
    reconciliations: list[str] = field(default_factory=list)

    def add(self, severity: str, category: str, detail: str) -> None:
        """Append a structured finding.

        Args:
            severity: Finding level (``OK``, ``WARNING``, ``DISCREPANCY``, …).
            category: Short category key for grouping in the audit log.
            detail: Free-text explanation.
        """
        self.findings.append(AuditFinding(severity=severity, category=category, detail=detail))

    def add_reconciliation(self, message: str) -> None:
        """Record a legacy → HGNC reconciliation line for the audit log.

        Args:
            message: Already-formatted reconciliation string.
        """
        self.reconciliations.append(message)


def normalize_symbol(value: object) -> str:
    """Normalize a gene-like token for lookups.

    Args:
        value: Raw cell / token value (may be ``NaN`` or non-string).

    Returns:
        Stripped string, or ``""`` when empty / missing.
    """
    text = str(value).strip()
    if not text or text.lower() in {"nan", "none"}:
        return ""
    return text


def map_legacy_symbol(raw: object) -> tuple[str, bool]:
    """Map a gene or locus identifier to an HGNC symbol when possible.

    Args:
        raw: Raw identifier from a table or text token.

    Returns:
        Tuple of ``(hgnc_symbol, was_legacy_mapped)``. ``was_legacy_mapped`` is
        ``True`` only when ``LEGACY_TO_HGNC`` changed the token.
    """
    text = normalize_symbol(raw)
    if not text:
        return "", False
    key = text.lower()
    if key in LEGACY_TO_HGNC:
        hgnc = LEGACY_TO_HGNC[key]
        return hgnc, hgnc != text
    # Preserve caller casing for symbols already treated as HGNC.
    return text, False


def load_supp_table(path: Path) -> pd.DataFrame:
    """Load Supplementary Table 3 and normalize column names.

    Args:
        path: Excel workbook with the 41-gene candidate summary.

    Returns:
        DataFrame with at least ``Gene_Symbol``, ``Variant_Count``, and
        ``Functional_Cluster``. Adds ``Gene_Symbol_Raw`` and
        ``Mapped_From_Legacy`` after HGNC reconciliation.

    Raises:
        FileNotFoundError: If ``path`` does not exist.
        ValueError: If required columns are missing after renaming.
    """
    if not path.is_file():
        raise FileNotFoundError(f"Supplementary Table 3 not found: {path}")

    df = pd.read_excel(path)
    rename = {
        "Gene Symbol": "Gene_Symbol",
        "gene_symbol": "Gene_Symbol",
        "Gene": "Gene_Symbol",
        "Variant Count": "Variant_Count",
        "variant_count": "Variant_Count",
        "n_variants": "Variant_Count",
        "Functional Cluster": "Functional_Cluster",
        "functional_cluster": "Functional_Cluster",
        "Pathway": "Functional_Cluster",
        "Longevity Class": "Longevity_Class",
        "longevity_class": "Longevity_Class",
        "Direction": "Longevity_Class",
    }
    df = df.rename(columns={c: rename.get(c, c) for c in df.columns})
    required = {"Gene_Symbol", "Variant_Count", "Functional_Cluster"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"Supplementary Table 3 missing columns {sorted(missing)}; "
            f"found {list(df.columns)}"
        )

    out = df.copy()
    mapped: list[str] = []
    legacy_flags: list[bool] = []
    for value in out["Gene_Symbol"]:
        symbol, was_legacy = map_legacy_symbol(value)
        mapped.append(symbol)
        legacy_flags.append(was_legacy)
    out["Gene_Symbol_Raw"] = out["Gene_Symbol"].map(normalize_symbol)
    out["Gene_Symbol"] = mapped
    out["Mapped_From_Legacy"] = legacy_flags
    out["Variant_Count"] = pd.to_numeric(out["Variant_Count"], errors="coerce")
    out["Functional_Cluster"] = out["Functional_Cluster"].astype(str).str.strip()
    if "Longevity_Class" in out.columns:
        out["Longevity_Class"] = out["Longevity_Class"].astype(str).str.strip()
    return out


def load_candidate_list(path: Path) -> pd.DataFrame:
    """Load the 41-gene candidate list used for concordance checks.

    Args:
        path: CSV with gene symbols, variant counts, and functional clusters.

    Returns:
        Normalized DataFrame aligned to Supplementary Table 3 column names.

    Raises:
        FileNotFoundError: If ``path`` does not exist.
        ValueError: If required columns are missing after renaming.
    """
    if not path.is_file():
        raise FileNotFoundError(f"Candidate list not found: {path}")

    df = pd.read_csv(path)
    rename = {
        "Gene Symbol": "Gene_Symbol",
        "gene_symbol": "Gene_Symbol",
        "Gene": "Gene_Symbol",
        "Variant Count": "Variant_Count",
        "variant_count": "Variant_Count",
        "Functional Cluster": "Functional_Cluster",
        "functional_cluster": "Functional_Cluster",
        "Longevity Class": "Longevity_Class",
        "longevity_class": "Longevity_Class",
    }
    df = df.rename(columns={c: rename.get(c, c) for c in df.columns})
    required = {"Gene_Symbol", "Variant_Count", "Functional_Cluster"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"Candidate list missing columns {sorted(missing)}; found {list(df.columns)}"
        )

    out = df.copy()
    out["Gene_Symbol"] = [map_legacy_symbol(v)[0] for v in out["Gene_Symbol"]]
    out["Variant_Count"] = pd.to_numeric(out["Variant_Count"], errors="coerce")
    out["Functional_Cluster"] = out["Functional_Cluster"].astype(str).str.strip()
    if "Longevity_Class" in out.columns:
        out["Longevity_Class"] = out["Longevity_Class"].astype(str).str.strip()
    return out


def scan_manuscript_text(text_dir: Path, report: AuditReport) -> dict[str, str]:
    """Find legacy locus identifiers in manuscript text and map them to HGNC.

    Patterns are applied longest-first. Overlapping matches are skipped so that
    e.g. ``HSP70-1A`` is not also counted as ``HSP70-1``. Whitespace is
    collapsed so phrase patterns survive line wraps.

    Args:
        text_dir: Directory containing ``.md`` / ``.txt`` manuscript drafts.
        report: Audit report that receives reconciliation notes and errors.

    Returns:
        Mapping of ``\"{path}:{matched_text}\"`` → HGNC symbol. Empty when the
        directory is missing or contains no text files.
    """
    if not text_dir.is_dir():
        report.add("ERROR", "text_dir", f"Manuscript text directory missing: {text_dir}")
        return {}

    replacements: dict[str, str] = {}
    files = sorted(
        p for p in text_dir.rglob("*") if p.is_file() and p.suffix.lower() in {".md", ".txt"}
    )
    if not files:
        report.add("WARNING", "text_dir", f"No .md/.txt files under {text_dir}")
        return {}

    for path in files:
        content = path.read_text(encoding="utf-8")
        # Collapse whitespace so phrase patterns match across line wraps.
        flat = re.sub(r"\s+", " ", content)
        claimed: list[tuple[int, int]] = []
        for pattern, hgnc in LEGACY_TEXT_PATTERNS:
            for match in re.finditer(pattern, flat, flags=re.IGNORECASE):
                raw = match.group(0).strip()
                display = re.sub(r"\s+", " ", raw)
                if display in {hgnc, "CETP", "HLA-DQB1"}:
                    continue
                start, end = match.span()
                # Skip spans already claimed by a longer / earlier pattern.
                if any(start < c_end and end > c_start for c_start, c_end in claimed):
                    continue
                claimed.append((start, end))
                key = f"{path}:{display}"
                if key not in replacements:
                    replacements[key] = hgnc
                    report.add_reconciliation(
                        f"{path.relative_to(REPO_ROOT)}: '{display}' → {hgnc}"
                    )
    return replacements


def crosscheck_candidate_vs_supp(
    candidates: pd.DataFrame,
    supp: pd.DataFrame,
    report: AuditReport,
) -> None:
    """Require concordance for symbols, variant counts, and clusters.

    Compares the candidate list against Supplementary Table 3 gene-by-gene.
    Set differences and per-gene mismatches are recorded as ``DISCREPANCY``
    findings. Perfect agreement adds an ``OK`` concordance finding.

    Args:
        candidates: Normalized 41-gene candidate list.
        supp: Normalized Supplementary Table 3.
        report: Audit report mutated in place.
    """
    cand = candidates.drop_duplicates(subset=["Gene_Symbol"], keep="first").set_index(
        "Gene_Symbol", drop=False
    )
    table = supp.drop_duplicates(subset=["Gene_Symbol"], keep="first").set_index(
        "Gene_Symbol", drop=False
    )

    cand_genes = set(cand.index.astype(str))
    supp_genes = set(table.index.astype(str))

    only_cand = sorted(cand_genes - supp_genes)
    only_supp = sorted(supp_genes - cand_genes)
    if only_cand:
        report.add(
            "DISCREPANCY",
            "gene_symbols",
            f"In candidate list only ({len(only_cand)}): {', '.join(only_cand)}",
        )
    if only_supp:
        report.add(
            "DISCREPANCY",
            "gene_symbols",
            f"In Supplementary Table 3 only ({len(only_supp)}): {', '.join(only_supp)}",
        )

    shared = sorted(cand_genes & supp_genes)
    for gene in shared:
        c_count = cand.at[gene, "Variant_Count"]
        s_count = table.at[gene, "Variant_Count"]
        if pd.isna(c_count) or pd.isna(s_count) or int(c_count) != int(s_count):
            report.add(
                "DISCREPANCY",
                "variant_counts",
                f"{gene}: candidate={c_count} vs SuppTable3={s_count}",
            )

        c_cluster = str(cand.at[gene, "Functional_Cluster"]).strip()
        s_cluster = str(table.at[gene, "Functional_Cluster"]).strip()
        if c_cluster != s_cluster:
            report.add(
                "DISCREPANCY",
                "functional_clusters",
                f"{gene}: candidate='{c_cluster}' vs SuppTable3='{s_cluster}'",
            )

        if "Longevity_Class" in cand.columns and "Longevity_Class" in table.columns:
            c_class = str(cand.at[gene, "Longevity_Class"]).strip()
            s_class = str(table.at[gene, "Longevity_Class"]).strip()
            if c_class != s_class:
                report.add(
                    "DISCREPANCY",
                    "longevity_class",
                    f"{gene}: candidate='{c_class}' vs SuppTable3='{s_class}'",
                )

    n_legacy = int(supp["Mapped_From_Legacy"].sum()) if "Mapped_From_Legacy" in supp.columns else 0
    if n_legacy:
        legacy_rows = supp.loc[supp["Mapped_From_Legacy"], ["Gene_Symbol_Raw", "Gene_Symbol"]]
        for raw, hgnc in legacy_rows.itertuples(index=False):
            report.add_reconciliation(f"Supplementary Table 3 gene column: '{raw}' → {hgnc}")

    if not only_cand and not only_supp:
        n_shared = len(shared)
        n_disc = sum(1 for f in report.findings if f.severity == "DISCREPANCY")
        if n_disc == 0:
            report.add(
                "OK",
                "concordance",
                (
                    f"100% concordance for {n_shared} genes "
                    "(symbols, variant counts, functional clusters)."
                ),
            )


def write_audit_log(path: Path, report: AuditReport) -> None:
    """Write a human-readable nomenclature audit log.

    Args:
        path: Destination log path (parent directories are created).
        report: Populated audit report.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    lines: list[str] = [
        "ROGEN Aging — gene nomenclature audit",
        f"Repository root: {REPO_ROOT}",
        "",
        "== Legacy → HGNC reconciliations ==",
    ]
    if report.reconciliations:
        lines.extend(f"  - {msg}" for msg in report.reconciliations)
    else:
        lines.append("  (none detected)")

    lines.extend(["", "== Findings =="])
    if not report.findings:
        lines.append("  (none)")
    else:
        for finding in report.findings:
            lines.append(f"  [{finding.severity}] {finding.category}: {finding.detail}")

    n_disc = sum(1 for f in report.findings if f.severity == "DISCREPANCY")
    lines.extend(
        [
            "",
            "== Summary ==",
            f"  reconciliations: {len(report.reconciliations)}",
            f"  discrepancies:   {n_disc}",
            f"  other findings:  {len(report.findings) - n_disc}",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def load_af_table(path: Path) -> pd.DataFrame:
    """Load ROGEN vs gnomAD v4 NFE allele frequencies.

    Accepts either manuscript column names (``ROGEN_AF``, ``gnomAD_AF``) or
    legacy pipeline aliases (``AF_1kg``, ``AF_gnomad_nfe``). ``AF_1kg`` is
    renamed to ``ROGEN_AF`` for plotting compatibility; the original source is
    recorded in ``af_source`` (``\"AF_1kg\"`` or ``\"ROGEN\"``).

    Args:
        path: Allele-frequency comparison CSV.

    Returns:
        DataFrame with ``rsID``, ``ROGEN_AF``, ``gnomAD_AF``, ``delta_af``,
        ``abs_delta_af``, ``af_source``, and boolean ``outlier``
        (|ΔAF| ≥ threshold).

    Raises:
        FileNotFoundError: If ``path`` does not exist.
        ValueError: If required columns are missing after renaming.
    """
    if not path.is_file():
        raise FileNotFoundError(f"AF comparison CSV not found: {path}")

    df = pd.read_csv(path)
    had_rogen = any(c in df.columns for c in ("ROGEN_AF", "AF_rogen", "rogen_af"))
    had_1kg = "AF_1kg" in df.columns
    rename = {
        "AF_gnomad_nfe": "gnomAD_AF",
        "gnomad_af": "gnomAD_AF",
        "AF_rogen": "ROGEN_AF",
        "rogen_af": "ROGEN_AF",
        "AF_1kg": "ROGEN_AF",
        "rsid": "rsID",
        "SNP_rsID": "rsID",
        "gene": "Gene",
    }
    df = df.rename(columns={c: rename.get(c, c) for c in df.columns})
    required = {"rsID", "ROGEN_AF", "gnomAD_AF"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"AF CSV missing columns {sorted(missing)}; found {list(df.columns)}")

    out = df.copy()
    out["rsID"] = out["rsID"].astype(str).str.strip()
    out["ROGEN_AF"] = pd.to_numeric(out["ROGEN_AF"], errors="coerce")
    out["gnomAD_AF"] = pd.to_numeric(out["gnomAD_AF"], errors="coerce")
    out["delta_af"] = out["gnomAD_AF"] - out["ROGEN_AF"]
    out["abs_delta_af"] = out["delta_af"].abs()
    out["outlier"] = out["abs_delta_af"] >= DELTA_AF_THRESHOLD
    if "af_source" not in out.columns:
        if had_rogen:
            out["af_source"] = "ROGEN"
        elif had_1kg:
            out["af_source"] = "AF_1kg"
        else:
            out["af_source"] = "ROGEN"
    if "Gene" in out.columns:
        out["Gene"] = out["Gene"].map(lambda v: map_legacy_symbol(v)[0] if pd.notna(v) else "")
    return out


def _configure_matplotlib() -> None:
    """Apply publication defaults and quiet noisy font subsetting logs."""
    logging.getLogger("fontTools").setLevel(logging.WARNING)
    logging.getLogger("matplotlib.font_manager").setLevel(logging.WARNING)
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica Neue", "Helvetica", "DejaVu Sans"],
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.linewidth": 0.8,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )


def _save_figure(fig: plt.Figure, stem: Path) -> tuple[Path, Path]:
    """Save a figure as both vector PDF and raster PNG.

    Args:
        fig: Matplotlib figure to write.
        stem: Output path without or with any suffix; ``.pdf`` / ``.png`` are
            applied explicitly.

    Returns:
        Tuple of ``(pdf_path, png_path)``.
    """
    stem.parent.mkdir(parents=True, exist_ok=True)
    pdf_path = stem.with_suffix(".pdf")
    png_path = stem.with_suffix(".png")
    for path, fmt in ((pdf_path, "pdf"), (png_path, "png")):
        fig.savefig(
            path,
            format=fmt,
            dpi=FIGURE_DPI if fmt == "png" else None,
            bbox_inches="tight",
            facecolor="white",
            edgecolor="none",
        )
    plt.close(fig)
    return pdf_path, png_path


def plot_af_scatter(af: pd.DataFrame, output_stem: Path) -> tuple[Path, Path]:
    """Draw a log-scale scatter of gnomAD_AF versus ROGEN_AF.

    Concordant points (|ΔAF| < threshold) are grey; outliers (|ΔAF| ≥ threshold)
    are highlighted and labelled by ``rsID``. Axes use a small epsilon floor so
    AF = 0 does not break the log scale (display-only; values are not rewritten
    on disk).

    Args:
        af: Output of :func:`load_af_table`.
        output_stem: Basename for ``Figure_AF_Scatter`` PDF/PNG.

    Returns:
        Paths to the written PDF and PNG.

    Raises:
        ValueError: If no paired AF values remain after dropping nulls.
    """
    _configure_matplotlib()
    work = af.dropna(subset=["ROGEN_AF", "gnomAD_AF"]).copy()
    if work.empty:
        raise ValueError("No paired ROGEN_AF / gnomAD_AF values to plot")

    # Avoid log(0): clamp tiny AFs for display only.
    eps = 1e-4
    x = work["ROGEN_AF"].clip(lower=eps)
    y = work["gnomAD_AF"].clip(lower=eps)
    outliers = work["outlier"].fillna(False).astype(bool)

    fig, ax = plt.subplots(figsize=(6.2, 6.0), dpi=FIGURE_DPI)
    ax.scatter(
        x[~outliers],
        y[~outliers],
        s=28,
        c="#9CA3AF",
        alpha=0.85,
        edgecolors="none",
        label=f"|ΔAF| < {DELTA_AF_THRESHOLD:.2f}",
        zorder=2,
    )
    ax.scatter(
        x[outliers],
        y[outliers],
        s=48,
        c="#C45C3E",
        alpha=0.95,
        edgecolors="#7F1D1D",
        linewidths=0.6,
        label=f"|ΔAF| ≥ {DELTA_AF_THRESHOLD:.2f}",
        zorder=3,
    )

    lim_lo = min(float(x.min()), float(y.min())) * 0.8
    lim_hi = max(float(x.max()), float(y.max())) * 1.2
    lim_lo = max(lim_lo, eps)
    ax.plot([lim_lo, lim_hi], [lim_lo, lim_hi], ls="--", lw=1.0, color="#4B5563", zorder=1)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lim_lo, lim_hi)
    ax.set_ylim(lim_lo, lim_hi)
    ax.set_xlabel("ROGEN AF (log scale)")
    ax.set_ylabel("gnomAD v4 NFE AF (log scale)")
    ax.set_title("Allele-frequency concordance: ROGEN vs gnomAD v4 NFE")
    ax.grid(True, which="both", ls=":", lw=0.5, color="#E5E7EB")
    ax.legend(frameon=False, loc="upper left")

    label_df = work.loc[outliers].sort_values("abs_delta_af", ascending=False)
    texts: list[plt.Text] = []
    for row in label_df.itertuples(index=False):
        label = str(row.rsID)
        text = ax.annotate(
            label,
            (max(float(row.ROGEN_AF), eps), max(float(row.gnomAD_AF), eps)),
            fontsize=8,
            color="#7F1D1D",
            xytext=(5, 5),
            textcoords="offset points",
        )
        texts.append(text)

    try:
        from adjustText import adjust_text

        if texts:
            adjust_text(texts, ax=ax, arrowprops={"arrowstyle": "-", "color": "#9CA3AF", "lw": 0.5})
    except ImportError:
        pass

    n_out = int(outliers.sum())
    ax.text(
        0.98,
        0.02,
        f"n = {len(work)} SNPs; outliers = {n_out}",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=8,
        color="#374151",
    )

    fig.tight_layout()
    return _save_figure(fig, output_stem)


def load_network_tables(
    interactions_path: Path,
    nodes_path: Path | None,
    gene_meta: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load STRING edges and node metadata for the 41-gene network.

    Args:
        interactions_path: CSV of pairwise interactions with a confidence score.
        nodes_path: Optional node metadata CSV; when missing, ``gene_meta``
            (typically Supplementary Table 3) is used.
        gene_meta: Fallback node table with gene / cluster / longevity columns.

    Returns:
        Tuple of ``(edges, nodes)`` DataFrames with reconciled HGNC symbols.

    Raises:
        FileNotFoundError: If the interactions file is missing.
        ValueError: If required edge or node columns are absent.
    """
    if not interactions_path.is_file():
        raise FileNotFoundError(f"Network interactions CSV not found: {interactions_path}")

    edges = pd.read_csv(interactions_path)
    rename = {
        "Gene_A": "gene_a",
        "Gene_B": "gene_b",
        "protein1": "gene_a",
        "protein2": "gene_b",
        "score": "string_score",
        "combined_score": "string_score",
    }
    edges = edges.rename(columns={c: rename.get(c, c) for c in edges.columns})
    required = {"gene_a", "gene_b", "string_score"}
    missing = required - set(edges.columns)
    if missing:
        raise ValueError(
            f"Interactions CSV missing columns {sorted(missing)}; found {list(edges.columns)}"
        )

    edges = edges.copy()
    edges["gene_a"] = [map_legacy_symbol(v)[0] for v in edges["gene_a"]]
    edges["gene_b"] = [map_legacy_symbol(v)[0] for v in edges["gene_b"]]
    edges["string_score"] = pd.to_numeric(edges["string_score"], errors="coerce")
    edges = edges.dropna(subset=["gene_a", "gene_b", "string_score"])
    edges = edges.loc[edges["string_score"] >= STRING_SCORE_MIN]

    if nodes_path is not None and nodes_path.is_file():
        nodes = pd.read_csv(nodes_path)
        nodes = nodes.rename(
            columns={
                "Gene_Symbol": "gene",
                "Gene": "gene",
                "Functional_Cluster": "functional_cluster",
                "Longevity_Class": "longevity_class",
                "Variant_Count": "variant_count",
            }
        )
    else:
        nodes = gene_meta.rename(
            columns={
                "Gene_Symbol": "gene",
                "Functional_Cluster": "functional_cluster",
                "Longevity_Class": "longevity_class",
                "Variant_Count": "variant_count",
            }
        ).copy()

    if "gene" not in nodes.columns:
        raise ValueError(f"Node table missing 'gene' column; found {list(nodes.columns)}")

    nodes = nodes.copy()
    nodes["gene"] = [map_legacy_symbol(v)[0] for v in nodes["gene"]]
    if "functional_cluster" not in nodes.columns:
        nodes["functional_cluster"] = "Other/Context"
    if "longevity_class" not in nodes.columns:
        nodes["longevity_class"] = "Context-Dependent"
    return edges, nodes


def _cluster_positions(
    genes_by_cluster: Mapping[str, Sequence[str]],
) -> dict[str, np.ndarray]:
    """Place each functional cluster on a circle with genes fanned on an arc.

    Virtual ``__hub__{cluster}`` coordinates are included for pathway label
    placement inside the ring.

    Args:
        genes_by_cluster: Mapping of pathway name → ordered gene symbols.

    Returns:
        Mapping of node id → ``(x, y)`` position arrays.
    """
    clusters = [c for c in CLUSTER_LAYOUT_ORDER if c in genes_by_cluster]
    clusters.extend(sorted(set(genes_by_cluster) - set(clusters)))
    n = max(len(clusters), 1)
    pos: dict[str, np.ndarray] = {}
    hub_radius = 2.15
    gene_radius = 3.55

    for i, cluster in enumerate(clusters):
        # Equally spaced cluster angles, starting at the top of the circle.
        theta_c = np.pi / 2.0 - (2.0 * np.pi * i) / n
        genes = list(genes_by_cluster[cluster])
        k = len(genes)
        if k == 1:
            offsets = [0.0]
        else:
            spread = min(0.95, 0.18 * k)
            offsets = list(np.linspace(-spread / 2.0, spread / 2.0, k))
        for off, gene in zip(offsets, genes, strict=True):
            theta = theta_c + off
            # Slight radial stagger so labels collide less within a cluster.
            r = gene_radius + (0.18 if abs(off) > 1e-9 else 0.0)
            pos[gene] = np.array([r * np.cos(theta), r * np.sin(theta)], dtype=float)
        pos[f"__hub__{cluster}"] = np.array(
            [hub_radius * np.cos(theta_c), hub_radius * np.sin(theta_c)],
            dtype=float,
        )
    return pos


def plot_gene_network(
    edges: pd.DataFrame,
    nodes: pd.DataFrame,
    output_stem: Path,
) -> tuple[Path, Path]:
    """Render the 41-gene STRING network with pathway grouping.

    Nodes are coloured by Pro-/Anti-/Context-Dependent longevity class, sized
    by degree centrality, and laid out by functional cluster. Edge width and
    alpha scale with STRING confidence (scores > 1 are treated as 0–1000).

    Args:
        edges: Interaction table from :func:`load_network_tables`.
        nodes: Node metadata from :func:`load_network_tables`.
        output_stem: Basename for ``Figure_41_Gene_Network`` PDF/PNG.

    Returns:
        Paths to the written PDF and PNG.
    """
    _configure_matplotlib()

    graph = nx.Graph()
    for row in nodes.itertuples(index=False):
        graph.add_node(
            str(row.gene),
            functional_cluster=str(row.functional_cluster),
            longevity_class=str(getattr(row, "longevity_class", "Context-Dependent")),
        )

    for row in edges.itertuples(index=False):
        a, b = str(row.gene_a), str(row.gene_b)
        if a == b:
            continue
        if a not in graph or b not in graph:
            continue
        score = float(row.string_score)
        # STRING scores may be 0–1 or 0–1000.
        if score > 1.0:
            score = score / 1000.0
        if graph.has_edge(a, b):
            graph[a][b]["weight"] = max(graph[a][b]["weight"], score)
        else:
            graph.add_edge(a, b, weight=score)

    genes_by_cluster: dict[str, list[str]] = {}
    for gene, data in graph.nodes(data=True):
        cluster = str(data.get("functional_cluster", "Other/Context"))
        genes_by_cluster.setdefault(cluster, []).append(gene)
    for cluster in genes_by_cluster:
        genes_by_cluster[cluster].sort()

    raw_pos = _cluster_positions(genes_by_cluster)
    hub_pos = {k: v for k, v in raw_pos.items() if k.startswith("__hub__")}
    pos = {k: v for k, v in raw_pos.items() if not k.startswith("__hub__")}

    centrality = nx.degree_centrality(graph)
    node_list = list(graph.nodes())
    sizes = [320.0 + 2600.0 * centrality.get(n, 0.0) for n in node_list]
    colors = [
        LONGEVITY_COLORS.get(
            str(graph.nodes[n].get("longevity_class", "Context-Dependent")),
            LONGEVITY_COLORS["Context-Dependent"],
        )
        for n in node_list
    ]

    edge_weights = [graph[u][v].get("weight", 0.2) for u, v in graph.edges()]
    edge_widths = [0.35 + 2.8 * w for w in edge_weights]
    edge_alphas = [0.20 + 0.55 * w for w in edge_weights]

    fig, ax = plt.subplots(figsize=(11.0, 10.0), dpi=FIGURE_DPI)
    ax.set_aspect("equal")
    ax.axis("off")

    for (u, v), width, alpha in zip(graph.edges(), edge_widths, edge_alphas, strict=True):
        nx.draw_networkx_edges(
            graph,
            pos,
            edgelist=[(u, v)],
            width=width,
            alpha=alpha,
            edge_color="#94A3B8",
            ax=ax,
        )

    nx.draw_networkx_nodes(
        graph,
        pos,
        nodelist=node_list,
        node_size=sizes,
        node_color=colors,
        linewidths=0.8,
        edgecolors="#1F2937",
        alpha=0.95,
        ax=ax,
    )
    nx.draw_networkx_labels(
        graph,
        pos,
        labels={n: n for n in node_list},
        font_size=6.8,
        font_color="#111827",
        ax=ax,
    )

    for cluster in genes_by_cluster:
        hub_key = f"__hub__{cluster}"
        if hub_key not in hub_pos:
            continue
        cx, cy = hub_pos[hub_key]
        ax.text(
            cx,
            cy,
            cluster.replace("/", "/" + "\n"),
            ha="center",
            va="center",
            fontsize=7.5,
            fontweight="600",
            color="#334155",
            alpha=0.95,
            bbox={
                "boxstyle": "round,pad=0.25",
                "facecolor": "white",
                "edgecolor": "#CBD5E1",
                "linewidth": 0.6,
                "alpha": 0.92,
            },
        )

    longevity_handles = [
        Patch(facecolor=color, edgecolor="#1F2937", label=label)
        for label, color in LONGEVITY_COLORS.items()
    ]
    size_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor="#6B7280",
            markeredgecolor="#1F2937",
            markersize=6,
            label="Low degree centrality",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor="#6B7280",
            markeredgecolor="#1F2937",
            markersize=14,
            label="High degree centrality",
        ),
    ]
    edge_handle = Line2D([0], [0], color="#94A3B8", lw=2.5, label="STRING confidence")
    ax.legend(
        handles=longevity_handles + size_handles + [edge_handle],
        loc="lower left",
        frameon=False,
        fontsize=8,
        title="Node / edge encoding",
        title_fontsize=8,
    )
    ax.set_title(
        "41-gene longevity network (STRING interactions; nodes by longevity class)",
        fontsize=11,
        fontweight="600",
        pad=12,
    )

    # Keep a small margin around the circular layout.
    ax.set_xlim(-5.1, 5.1)
    ax.set_ylim(-5.1, 5.1)

    fig.tight_layout()
    return _save_figure(fig, output_stem)


def run_pipeline(
    supp_table: Path,
    candidate_list: Path,
    text_dir: Path,
    af_csv: Path,
    network_csv: Path,
    nodes_csv: Path,
    audit_log: Path,
    fig_dir: Path,
) -> None:
    """Execute nomenclature audit and figure generation.

    Args:
        supp_table: Supplementary Table 3 Excel path.
        candidate_list: 41-gene candidate list CSV.
        text_dir: Manuscript text directory for legacy-name scanning.
        af_csv: ROGEN vs gnomAD AF comparison CSV.
        network_csv: STRING interaction edges CSV.
        nodes_csv: Optional node metadata CSV.
        audit_log: Destination for ``nomenclature_audit.log``.
        fig_dir: Directory for PDF/PNG figure exports.

    Raises:
        typer.Exit: Exit code ``1`` when concordance discrepancies remain.
    """
    report = AuditReport()
    typer.echo("1/3  Auditing nomenclature and concordance…")

    supp = load_supp_table(supp_table)
    candidates = load_candidate_list(candidate_list)
    scan_manuscript_text(text_dir, report)
    crosscheck_candidate_vs_supp(candidates, supp, report)
    write_audit_log(audit_log, report)
    typer.echo(f"     Audit log → {audit_log}")

    typer.echo("2/3  Allele-frequency scatter (ROGEN vs gnomAD v4 NFE)…")
    af = load_af_table(af_csv)
    af_pdf, af_png = plot_af_scatter(af, fig_dir / "Figure_AF_Scatter")
    typer.echo(f"     {af_pdf}")
    typer.echo(f"     {af_png}")
    n_out = int(af["outlier"].fillna(False).sum())
    typer.echo(f"     Outliers (|ΔAF| ≥ {DELTA_AF_THRESHOLD}): {n_out}")

    typer.echo("3/3  41-gene functional network…")
    edges, nodes = load_network_tables(network_csv, nodes_csv, supp)
    net_pdf, net_png = plot_gene_network(edges, nodes, fig_dir / "Figure_41_Gene_Network")
    typer.echo(f"     {net_pdf}")
    typer.echo(f"     {net_png}")

    n_disc = sum(1 for f in report.findings if f.severity == "DISCREPANCY")
    if n_disc:
        typer.echo(f"Done with {n_disc} discrepancy(ies) — see audit log.")
        raise typer.Exit(code=1)
    typer.echo("Done — full concordance; figures written.")


@app.command()
def main(
    supp_table: Annotated[
        Path,
        typer.Option(help="Supplementary Table 3 Excel path"),
    ] = DEFAULT_SUPP_TABLE,
    candidate_list: Annotated[
        Path,
        typer.Option(help="41-gene candidate list CSV"),
    ] = DEFAULT_CANDIDATE_LIST,
    text_dir: Annotated[
        Path,
        typer.Option(help="Manuscript text directory"),
    ] = DEFAULT_TEXT_DIR,
    af_csv: Annotated[
        Path,
        typer.Option(help="ROGEN_AF vs gnomAD_AF comparison CSV"),
    ] = DEFAULT_AF_CSV,
    network_csv: Annotated[
        Path,
        typer.Option(help="STRING interaction edges for the 41 genes"),
    ] = DEFAULT_NETWORK_CSV,
    nodes_csv: Annotated[
        Path,
        typer.Option(help="Optional node metadata CSV"),
    ] = DEFAULT_NODES_CSV,
    audit_log: Annotated[
        Path,
        typer.Option(help="Nomenclature audit log path"),
    ] = DEFAULT_AUDIT_LOG,
    fig_dir: Annotated[
        Path,
        typer.Option(help="Output directory for PDF/PNG figures"),
    ] = DEFAULT_FIG_DIR,
) -> None:
    """Reconcile legacy gene names and build AF + network publication figures.

    Args:
        supp_table: Supplementary Table 3 Excel path.
        candidate_list: 41-gene candidate list CSV.
        text_dir: Manuscript text directory.
        af_csv: ROGEN vs gnomAD AF comparison CSV.
        network_csv: STRING interaction edges CSV.
        nodes_csv: Optional node metadata CSV.
        audit_log: Nomenclature audit log path.
        fig_dir: Output directory for PDF/PNG figures.
    """
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    run_pipeline(
        supp_table=supp_table,
        candidate_list=candidate_list,
        text_dir=text_dir,
        af_csv=af_csv,
        network_csv=network_csv,
        nodes_csv=nodes_csv,
        audit_log=audit_log,
        fig_dir=fig_dir,
    )


if __name__ == "__main__":
    app()
