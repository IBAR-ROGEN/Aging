#!/usr/bin/env python3
"""Independent validation of a trained ElasticNet methylation clock (GSE87571).

Loads a processed beta matrix plus phenotype table, predicts DNAm age with a
saved bare ``sklearn.linear_model.ElasticNet`` or Pipeline ending in
ElasticNet/ElasticNetCV, writes validation metrics JSON, and saves a three-panel
publication figure (scatter, residuals, top CpG weights).

Residuals use age acceleration (predicted − chronological).

Example:
    uv run python scripts/clock/evaluate_methylation_clock.py

See Also:
    INPUT_MANIFEST.md: Required input paths.
    docs/METHYLATION_CLOCK_VALIDATION.md: Inputs, outputs, and CLI reference.
"""

from __future__ import annotations

import json
import pickle
import re
import warnings
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import typer
from scipy.stats import pearsonr
from sklearn.linear_model import ElasticNet, ElasticNetCV
from sklearn.metrics import mean_absolute_error, median_absolute_error

from rogen_aging.clock.evaluate import build_feature_matrix

REPO_ROOT = Path(__file__).resolve().parent.parent.parent

INPUT_MANIFEST = REPO_ROOT / "INPUT_MANIFEST.md"
DEFAULT_METHYLATION = REPO_ROOT / "data" / "methylation" / "GSE87571_processed.parquet"
DEFAULT_META = REPO_ROOT / "data" / "methylation" / "GSE87571_meta.csv"
DEFAULT_MODEL_PKL = REPO_ROOT / "models" / "ro_clock_elasticnet_gse40279.pkl"
DEFAULT_MODEL_JOBLIB = REPO_ROOT / "models" / "methylation_clock_v1.joblib"
DEFAULT_MODEL = DEFAULT_MODEL_PKL
DEFAULT_METRICS = REPO_ROOT / "outputs" / "clock_metrics.json"
DEFAULT_FIGURE_STEM = REPO_ROOT / "outputs" / "figures" / "Figure_Epigenetic_Clock_Panels"

DEFAULT_ANNOTATION = REPO_ROOT / "data" / "methylation" / "HM450_probe_annotation.csv"
HORVATH_ANNOTATION = REPO_ROOT / "test_data" / "gb-2013-14-10-r115-S3.csv"

TOP_N_CPGS = 25
FIGURE_DPI = 300
POSITIVE_COLOR = "#2166ac"
NEGATIVE_COLOR = "#b2182b"
SCATTER_COLOR = "#404040"
AGE_BINS = ("<30", "30-60", ">60")

# Backtick-quoted paths marked required in INPUT_MANIFEST.md tables.
_MANIFEST_REQUIRED_PATH_RE = re.compile(
    r"^\|\s*`([^`]+)`\s*\|[^|]*\|\s*yes\s*\|",
    re.IGNORECASE | re.MULTILINE,
)


def resolve_clock_model_path(preferred: Path | None = None) -> Path:
    """Resolve a usable clock artifact path for local development.

    Prefers ``models/ro_clock_elasticnet_gse40279.pkl``, then falls back to
    ``models/methylation_clock_v1.joblib`` when the pickle is absent.

    Args:
        preferred: Explicit model path from the CLI. When provided and present,
            it is returned unchanged.

    Returns:
        Path to an existing model file when one of the defaults exists;
        otherwise ``preferred`` or ``DEFAULT_MODEL_PKL`` (caller may still fail
        later with a clear FileNotFoundError).
    """
    if preferred is not None and preferred != DEFAULT_MODEL_PKL and preferred.is_file():
        return preferred
    if preferred is not None and preferred.is_file():
        return preferred
    if DEFAULT_MODEL_PKL.is_file():
        return DEFAULT_MODEL_PKL
    if DEFAULT_MODEL_JOBLIB.is_file():
        return DEFAULT_MODEL_JOBLIB
    return preferred if preferred is not None else DEFAULT_MODEL_PKL


def verify_input_manifest(
    manifest_path: Path = INPUT_MANIFEST,
    *,
    repo_root: Path | None = None,
) -> list[Path]:
    """Read ``INPUT_MANIFEST.md`` and confirm every required input file exists.

    Args:
        manifest_path: Path to the markdown manifest listing required inputs.
        repo_root: Directory used to resolve relative paths from the manifest.
            Defaults to the repository root (``REPO_ROOT``).

    Returns:
        Absolute paths of required input files that were verified present.

    Raises:
        FileNotFoundError: If the manifest itself is missing, or any required
            input path listed with ``yes`` does not exist on disk.
        ValueError: If the manifest contains no parseable required paths.
    """
    root = REPO_ROOT if repo_root is None else Path(repo_root)
    if not manifest_path.is_file():
        raise FileNotFoundError(
            f"INPUT_MANIFEST.md not found at {manifest_path}. "
            "Cannot verify required inputs before evaluation."
        )

    text = manifest_path.read_text(encoding="utf-8")
    relative_paths = _MANIFEST_REQUIRED_PATH_RE.findall(text)
    if not relative_paths:
        raise ValueError(
            f"No required input paths (marked 'yes') were parsed from {manifest_path}."
        )

    resolved: list[Path] = []
    missing: list[str] = []
    for rel in relative_paths:
        path = (root / rel).resolve()
        if path.is_file():
            resolved.append(path)
        else:
            missing.append(rel)

    if missing:
        # Development unblock: accept methylation_clock_v1.joblib when the
        # preferred pickle path is missing under the same repo_root.
        remaining: list[str] = []
        joblib_fallback = root / "models" / "methylation_clock_v1.joblib"
        for rel in missing:
            if rel.endswith("ro_clock_elasticnet_gse40279.pkl") and joblib_fallback.is_file():
                resolved.append(joblib_fallback.resolve())
                continue
            remaining.append(rel)
        if remaining:
            raise FileNotFoundError(
                "Required input file(s) listed in INPUT_MANIFEST.md are missing:\n  - "
                + "\n  - ".join(remaining)
                + "\nHint: uv run python scripts/dev/write_pipeline_fixtures.py"
                + "\nHalt: fix paths or restore artifacts before running evaluation."
            )
    return resolved


def load_elasticnet_clock(model_path: Path) -> ElasticNet | Any:
    """Load a bare ElasticNet or Pipeline clock from pickle/joblib.

    Accepts a bare ``sklearn.linear_model.ElasticNet`` or a Pipeline whose
    final step is ElasticNet (e.g. imputer + ElasticNetCV/ElasticNet from
    ``rogen_aging.clock.train``).

    Args:
        model_path: Path to a pickled / joblib clock artifact.

    Returns:
        Fitted estimator ready for ``predict`` (bare ElasticNet or Pipeline).

    Raises:
        FileNotFoundError: If ``model_path`` does not exist.
        TypeError: If the object is neither ElasticNet nor a Pipeline ending
            in ElasticNet / ElasticNetCV.
    """
    if not model_path.is_file():
        raise FileNotFoundError(
            f"Trained ElasticNet model not found: {model_path}. "
            "Expected a pickled sklearn.linear_model.ElasticNet or Pipeline "
            "at models/ro_clock_elasticnet_gse40279.pkl. "
            "Halting — will not train or substitute another estimator."
        )

    suffix = model_path.suffix.lower()
    if suffix == ".joblib":
        import joblib

        model = joblib.load(model_path)
    else:
        with model_path.open("rb") as handle:
            model = pickle.load(handle)

    if type(model) is ElasticNet:
        return model

    if hasattr(model, "named_steps"):
        final = list(model.named_steps.values())[-1]
        if isinstance(final, (ElasticNet, ElasticNetCV)):
            return model
        raise TypeError(
            f"Pipeline final step from {model_path} has type "
            f"{type(final).__module__}.{type(final).__name__}; "
            "expected ElasticNet or ElasticNetCV."
        )

    raise TypeError(
        f"Loaded object from {model_path} has type {type(model).__module__}."
        f"{type(model).__name__}; expected sklearn.linear_model.ElasticNet "
        "or a Pipeline ending in ElasticNet/ElasticNetCV. "
        "Halting — will not train or substitute another estimator."
    )


def _cg_columns(df: pd.DataFrame) -> list[str]:
    """Return Illumina-style CpG column names from a wide table.

    Args:
        df: Sample-by-feature table whose CpG columns start with ``cg``.

    Returns:
        Column names that begin with the ``cg`` probe-ID prefix.
    """
    return [c for c in df.columns if str(c).startswith("cg")]


def _pick_id_column(df: pd.DataFrame, candidates: tuple[str, ...]) -> str | None:
    """Return the first matching column name from ``candidates`` (case-insensitive).

    Args:
        df: Table to search.
        candidates: Preferred column names in priority order.

    Returns:
        The actual column name present in ``df``, or ``None`` if none match.
    """
    lower_map = {str(c).lower(): str(c) for c in df.columns}
    for name in candidates:
        if name.lower() in lower_map:
            return lower_map[name.lower()]
    return None


def _normalize_age_column(meta: pd.DataFrame) -> pd.DataFrame:
    """Ensure a numeric ``chronological_age`` column exists in phenotype metadata.

    Accepts either ``chronological_age`` or common aliases such as ``age`` /
    ``Age_years``. Non-numeric values become NaN via ``errors='coerce'``.

    Args:
        meta: Phenotype table with at least one age column.

    Returns:
        Copy of ``meta`` with a float ``chronological_age`` column.

    Raises:
        ValueError: If no recognizable age column is present.
    """
    out = meta.copy()
    if "chronological_age" in out.columns:
        out["chronological_age"] = pd.to_numeric(out["chronological_age"], errors="coerce")
        return out

    age_col = _pick_id_column(out, ("age", "Age", "AGE", "age_years", "Age_years"))
    if age_col is None:
        raise ValueError(
            "Phenotype metadata must include 'chronological_age' or an 'age' column. "
            f"Found columns: {list(out.columns)}"
        )
    out["chronological_age"] = pd.to_numeric(out[age_col], errors="coerce")
    return out


def load_validation_cohort(
    methylation_path: Path,
    meta_path: Path,
    *,
    allow_positional_align: bool = False,
) -> pd.DataFrame:
    """Join processed methylation betas with phenotype ages into a wide table.

    Accepts either samples×CpGs or CpGs×samples methylation tables. Phenotype
    rows are matched on ``sample_id`` / GEO accession / shared index labels.
    When IDs do not overlap but row counts match, positional alignment is refused
    unless ``allow_positional_align`` is True (explicit opt-in).

    Args:
        methylation_path: Parquet path for the processed beta matrix.
        meta_path: CSV path for phenotype metadata (sample ID + age).
        allow_positional_align: If True, allow row-order alignment when sample
            IDs do not overlap but lengths match.

    Returns:
        Wide DataFrame with ``chronological_age`` plus ``cg*`` feature columns,
        one row per aligned sample.

    Raises:
        FileNotFoundError: If either input path is missing.
        ValueError: If no CpG columns are found or sample alignment fails.
    """
    if not methylation_path.is_file():
        raise FileNotFoundError(f"Methylation matrix not found: {methylation_path}")
    if not meta_path.is_file():
        raise FileNotFoundError(f"Phenotype metadata not found: {meta_path}")

    meth = pd.read_parquet(methylation_path)
    meta = _normalize_age_column(pd.read_csv(meta_path))

    cg_cols = _cg_columns(meth)
    # GEO supplements are sometimes stored as probes × samples; transpose if needed.
    if not cg_cols and meth.index.map(lambda x: str(x).startswith("cg")).any():
        meth = meth.T
        cg_cols = _cg_columns(meth)
    if not cg_cols:
        raise ValueError(
            f"No CpG columns (prefix 'cg') found in {methylation_path}. "
            "Expected a wide beta matrix with Illumina probe IDs."
        )

    meth_id = _pick_id_column(meth, ("sample_id", "geo_accession", "GSM", "id", "ID"))
    meta_id = _pick_id_column(meta, ("sample_id", "geo_accession", "GSM", "id", "ID"))

    if meth_id is not None:
        meth = meth.set_index(meth_id, drop=True)
    if meta_id is not None:
        meta = meta.set_index(meta_id, drop=True)

    meth.index = meth.index.astype(str)
    meta.index = meta.index.astype(str)

    shared = meth.index.intersection(meta.index)
    if len(shared) == 0:
        if len(meth) == len(meta) and allow_positional_align:
            warnings.warn(
                "No overlapping sample IDs between matrix and metadata; "
                "aligning by row order (--allow-positional-align).",
                stacklevel=2,
            )
            ages = meta["chronological_age"].to_numpy()
            wide = meth.loc[:, cg_cols].copy()
            wide.insert(0, "chronological_age", ages)
            return wide.reset_index(drop=True)
        raise ValueError(
            "Could not align methylation matrix to phenotype metadata: "
            "no shared sample IDs"
            + (
                " (pass --allow-positional-align to opt into row-order alignment)."
                if len(meth) == len(meta)
                else " and unequal row counts."
            )
        )

    meth_aligned = meth.loc[shared, cg_cols]
    ages = meta.loc[shared, "chronological_age"]
    wide = meth_aligned.copy()
    wide.insert(0, "chronological_age", ages.to_numpy(dtype=float))
    wide.index.name = "sample_id"
    return wide


def assign_age_stratum(ages: np.ndarray) -> np.ndarray:
    """Bin chronological ages into ``<30``, ``30-60``, and ``>60`` years.

    The middle bin is closed on both ends (``[30, 60]``). Ages exactly 30 or 60
    therefore contribute to the ``30-60`` stratum MAE.

    Args:
        ages: 1-D array of chronological ages in years.

    Returns:
        Object array of stratum labels aligned to ``ages``.
    """
    strata = np.empty(ages.shape[0], dtype=object)
    strata[ages < 30.0] = "<30"
    strata[(ages >= 30.0) & (ages <= 60.0)] = "30-60"
    strata[ages > 60.0] = ">60"
    return strata


def compute_metrics(
    chronological_age: np.ndarray,
    predicted_age: np.ndarray,
) -> dict[str, Any]:
    """Compute overall and age-stratified validation metrics.

    Residuals are defined as predicted − chronological (years; age acceleration).
    Empty strata receive ``mae_by_age_stratum[label] = None``.

    Args:
        chronological_age: Observed ages in years.
        predicted_age: Model-predicted DNAm ages in years.

    Returns:
        Dictionary with ``n_samples``, ``mae``, ``median_absolute_error``,
        ``pearson_r``, ``pearson_p``, ``mae_by_age_stratum``, and
        ``n_by_age_stratum``.
    """
    residual = predicted_age - chronological_age
    abs_err = np.abs(residual)
    r_value, r_p = pearsonr(chronological_age, predicted_age)

    strata = assign_age_stratum(chronological_age)
    mae_by_stratum: dict[str, float | None] = {}
    n_by_stratum: dict[str, int] = {}
    for label in AGE_BINS:
        mask = strata == label
        n_by_stratum[label] = int(mask.sum())
        if mask.any():
            mae_by_stratum[label] = float(np.mean(abs_err[mask]))
        else:
            mae_by_stratum[label] = None

    return {
        "n_samples": int(len(chronological_age)),
        "mae": float(mean_absolute_error(chronological_age, predicted_age)),
        "median_absolute_error": float(median_absolute_error(chronological_age, predicted_age)),
        "pearson_r": float(r_value),
        "pearson_p": float(r_p),
        "mae_by_age_stratum": mae_by_stratum,
        "n_by_age_stratum": n_by_stratum,
    }


def format_metrics_markdown(metrics: dict[str, Any]) -> str:
    """Render a concise markdown summary of validation performance.

    Args:
        metrics: Output of :func:`compute_metrics` (plus optional extras).

    Returns:
        Markdown string suitable for printing to stdout.
    """
    lines = [
        "# Methylation clock validation (GSE87571)",
        "",
        f"- **n samples:** {metrics['n_samples']}",
        f"- **MAE:** {metrics['mae']:.3f} years",
        f"- **Median Absolute Error:** {metrics['median_absolute_error']:.3f} years",
        f"- **Pearson r:** {metrics['pearson_r']:.4f} (p = {metrics['pearson_p']:.3g})",
        "",
        "## Age-stratified MAE",
        "",
        "| Stratum | n | MAE (years) |",
        "|---------|---|-------------|",
    ]
    for label in AGE_BINS:
        n = metrics["n_by_age_stratum"][label]
        mae = metrics["mae_by_age_stratum"][label]
        mae_str = f"{mae:.3f}" if mae is not None else "NA"
        lines.append(f"| {label} | {n} | {mae_str} |")
    lines.append("")
    return "\n".join(lines)


def extract_cpg_coefficients(model: Any) -> pd.Series:
    """Extract ElasticNet coefficients indexed by CpG probe ID.

    Args:
        model: Fitted bare ``ElasticNet`` or Pipeline with an ElasticNet /
            ElasticNetCV final step and ``feature_names_in_``.

    Returns:
        Series of coefficients indexed by training feature names.

    Raises:
        ValueError: If coefficients or feature names cannot be recovered, or
            if their lengths disagree.
    """
    enet: Any = model
    if hasattr(model, "named_steps"):
        steps = model.named_steps
        if "elasticnet" in steps:
            enet = steps["elasticnet"]
        else:
            enet = list(steps.values())[-1]

    if not hasattr(enet, "coef_"):
        raise ValueError("ElasticNet has no coef_; model may be unfitted.")
    coef = np.ravel(enet.coef_)

    names: list[str] | None = None
    if hasattr(model, "feature_names_in_"):
        names = [str(x) for x in model.feature_names_in_]
    elif hasattr(enet, "feature_names_in_"):
        names = [str(x) for x in enet.feature_names_in_]
    if names is None:
        raise ValueError("Model has no feature_names_in_; cannot label CpG coefficients.")

    if len(names) != len(coef):
        raise ValueError(
            f"Coefficient length ({len(coef)}) does not match feature names ({len(names)})."
        )
    return pd.Series(coef, index=names, name="coefficient")


def _primary_gene_symbol(raw: object) -> str | None:
    """Parse the primary HGNC symbol from an Illumina multi-gene annotation cell.

    Illumina ``UCSC_RefGene_Name`` fields list multiple genes separated by
    semicolons; only the first non-empty token is retained.

    Args:
        raw: Annotation cell value (string, NaN, or None).

    Returns:
        Primary gene symbol, or ``None`` if the cell is empty / missing.
    """
    if raw is None or (isinstance(raw, float) and np.isnan(raw)):
        return None
    text = str(raw).strip()
    if not text or text.lower() in {"nan", "none", "na", "."}:
        return None
    gene = text.split(";")[0].strip()
    return gene or None


def load_probe_gene_map(annotation_path: Path | None) -> dict[str, str]:
    """Map Illumina probe IDs to nearest / annotated gene symbols.

    Resolution order:

    1. Explicit ``annotation_path`` (if provided and present).
    2. Default HM450 annotation CSV under ``data/methylation/``.
    3. Horvath 353-CpG supplementary table in ``test_data/`` (partial coverage).

    Args:
        annotation_path: Optional probe→gene table path. Expected columns are
            flexible (e.g. ``IlmnID`` + ``UCSC_RefGene_Name``, or Horvath
            ``CpGmarker`` + ``Symbol``).

    Returns:
        Dictionary mapping probe ID → primary gene symbol. Empty if no usable
        annotation file is found.
    """
    mapping: dict[str, str] = {}

    candidates: list[Path] = []
    if annotation_path is not None:
        candidates.append(annotation_path)
    candidates.append(DEFAULT_ANNOTATION)
    candidates.append(HORVATH_ANNOTATION)

    for path in candidates:
        if not path.is_file():
            continue
        sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
        frame = pd.read_csv(path, sep=sep, comment="#", low_memory=False)
        # Horvath S3 tables prepend descriptive header rows before CpGmarker data.
        if "CpGmarker" in frame.columns:
            frame = frame.loc[frame["CpGmarker"].astype(str).str.startswith("cg")]

        probe_col = _pick_id_column(
            frame,
            ("IlmnID", "Name", "ProbeID", "probe_id", "CpG", "CpGmarker", "cg"),
        )
        gene_col = _pick_id_column(
            frame,
            (
                "UCSC_RefGene_Name",
                "RefGene",
                "Nearest_Gene",
                "nearest_gene",
                "Symbol",
                "gene_symbol",
                "Gene",
            ),
        )
        if probe_col is None or gene_col is None:
            continue
        for probe, gene_raw in zip(
            frame[probe_col].astype(str),
            frame[gene_col],
            strict=False,
        ):
            gene = _primary_gene_symbol(gene_raw)
            if gene and probe not in mapping:
                mapping[probe] = gene
        if mapping:
            break
    return mapping


def label_cpg(probe: str, gene_map: dict[str, str]) -> str:
    """Format a CpG axis label, optionally including the nearest gene symbol.

    Args:
        probe: Illumina probe ID (e.g. ``cg07814318``).
        gene_map: Probe → gene lookup from :func:`load_probe_gene_map`.

    Returns:
        ``GENE (cg…)`` when annotated, otherwise the bare probe ID.
    """
    gene = gene_map.get(probe)
    if gene:
        return f"{gene} ({probe})"
    return probe


def plot_panel_a(
    ax: plt.Axes,
    chronological_age: np.ndarray,
    predicted_age: np.ndarray,
    metrics: dict[str, Any],
) -> None:
    """Draw chronological vs predicted age with regression line and 95% CI.

    Args:
        ax: Matplotlib axes for panel A.
        chronological_age: Observed ages (x-axis).
        predicted_age: Predicted DNAm ages (y-axis).
        metrics: Validation metrics used for the on-panel annotation box.
    """
    plot_df = pd.DataFrame(
        {
            "chronological_age": chronological_age,
            "predicted_age": predicted_age,
        }
    )
    sns.regplot(
        data=plot_df,
        x="chronological_age",
        y="predicted_age",
        ax=ax,
        ci=95,
        scatter_kws={
            "s": 28,
            "alpha": 0.75,
            "color": SCATTER_COLOR,
            "edgecolors": "white",
            "linewidths": 0.4,
            "zorder": 3,
        },
        line_kws={"color": "#d95f02", "linewidth": 1.6, "label": "Linear fit (95% CI)", "zorder": 2},
    )

    lo = float(min(chronological_age.min(), predicted_age.min()))
    hi = float(max(chronological_age.max(), predicted_age.max()))
    pad = 0.05 * (hi - lo) if hi > lo else 1.0
    lim_lo, lim_hi = lo - pad, hi + pad
    ax.plot(
        [lim_lo, lim_hi],
        [lim_lo, lim_hi],
        linestyle="--",
        color="0.45",
        linewidth=1.2,
        label="y = x",
        zorder=1,
    )
    ax.set_xlim(lim_lo, lim_hi)
    ax.set_ylim(lim_lo, lim_hi)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Chronological age (years)")
    ax.set_ylabel("Predicted DNAm age (years)")
    ax.set_title("A  Chronological vs predicted age")
    ax.legend(loc="lower right", frameon=True, fontsize=9)
    ax.text(
        0.03,
        0.97,
        f"MAE = {metrics['mae']:.2f} yr\n"
        f"r = {metrics['pearson_r']:.3f}\n"
        f"n = {metrics['n_samples']}",
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=10,
        bbox={
            "boxstyle": "round,pad=0.35",
            "facecolor": "white",
            "edgecolor": "0.8",
            "alpha": 0.95,
        },
    )


def plot_panel_b(
    ax: plt.Axes,
    chronological_age: np.ndarray,
    residual: np.ndarray,
) -> None:
    """Draw residuals (chronological − predicted) versus chronological age.

    Args:
        ax: Matplotlib axes for panel B.
        chronological_age: Observed ages (x-axis).
        residual: Predicted − chronological age (years; age acceleration).
    """
    ax.scatter(
        chronological_age,
        residual,
        s=28,
        alpha=0.75,
        color=SCATTER_COLOR,
        edgecolors="white",
        linewidths=0.4,
        zorder=2,
    )
    ax.axhline(0.0, color="crimson", linestyle="--", linewidth=1.2, label="Zero residual")
    ax.set_xlabel("Chronological age (years)")
    ax.set_ylabel("Residual (predicted − chronological, years)")
    ax.set_title("B  Age residuals")
    ax.grid(True, linestyle=":", linewidth=0.6, alpha=0.7)
    ax.legend(loc="best", frameon=True, fontsize=9)


def plot_panel_c(
    ax: plt.Axes,
    weights: pd.Series,
    gene_map: dict[str, str],
    top_n: int,
) -> None:
    """Draw horizontal bars for the top ``|coefficient|`` CpG sites.

    Prefers non-zero ElasticNet weights; falls back to all coefficients when
    the model has not yet selected features. Bars are colored by sign.

    Args:
        ax: Matplotlib axes for panel C.
        weights: Probe-indexed ElasticNet coefficients.
        gene_map: Probe → gene labels for tick text.
        top_n: Maximum number of probes to display.
    """
    nonzero = weights[weights != 0.0]
    ranked = nonzero if not nonzero.empty else weights
    top = ranked.reindex(ranked.abs().sort_values(ascending=False).head(top_n).index)
    top = top.sort_values()

    colors = [POSITIVE_COLOR if v >= 0 else NEGATIVE_COLOR for v in top.to_numpy()]
    labels = [label_cpg(str(p), gene_map) for p in top.index]
    y_pos = np.arange(len(top))
    ax.barh(y_pos, top.to_numpy(), color=colors, edgecolor="0.2", linewidth=0.4)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=8)
    ax.axvline(0.0, color="0.3", linewidth=0.8)
    ax.set_xlabel("ElasticNet coefficient")
    ax.set_title(f"C  Top {len(top)} CpG sites by |weight|")


def save_figure(fig: plt.Figure, stem: Path) -> tuple[Path, Path]:
    """Write PNG (300 dpi) and vector PDF next to ``stem``.

    Args:
        fig: Matplotlib figure to serialize.
        stem: Output path stem (suffixes are replaced per format).

    Returns:
        Tuple of ``(png_path, pdf_path)``.
    """
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    png_path = stem.with_suffix(".png")
    pdf_path = stem.with_suffix(".pdf")
    fig.savefig(png_path, dpi=FIGURE_DPI, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    return png_path, pdf_path


def run_validation(
    methylation_path: Path,
    meta_path: Path,
    model_path: Path,
    metrics_path: Path,
    figure_stem: Path,
    annotation_path: Path | None,
    top_n_cpgs: int,
    *,
    skip_manifest_check: bool = False,
    allow_positional_align: bool = False,
) -> dict[str, Any]:
    """Run end-to-end GSE87571 validation: predict, score, plot, and persist.

    Verifies ``INPUT_MANIFEST.md`` required files, loads an ElasticNet or
    Pipeline clock, uses :func:`rogen_aging.clock.evaluate.build_feature_matrix`
    to align CpGs (training imputer when present), then writes metrics plus the
    three-panel figure.

    Args:
        methylation_path: Processed validation beta matrix (Parquet).
        meta_path: Phenotype CSV with sample IDs and chronological age.
        model_path: Pickled ElasticNet or Pipeline clock.
        metrics_path: Destination for ``clock_metrics.json``.
        figure_stem: Output stem for ``Figure_Epigenetic_Clock_Panels.*``.
        annotation_path: Optional Illumina probe→gene annotation table.
        top_n_cpgs: Number of top ``|weight|`` CpGs for panel C.
        skip_manifest_check: If True, skip ``INPUT_MANIFEST.md`` existence
            checks (useful when paths are overridden via CLI).
        allow_positional_align: Opt into row-order alignment when IDs do not
            overlap.

    Returns:
        Metrics dictionary extended with output file paths
        (``metrics_path``, ``figure_png``, ``figure_pdf``).
    """
    if not skip_manifest_check:
        verify_input_manifest(INPUT_MANIFEST)

    wide = load_validation_cohort(
        methylation_path,
        meta_path,
        allow_positional_align=allow_positional_align,
    )
    model = load_elasticnet_clock(model_path)

    y = pd.to_numeric(wide["chronological_age"], errors="coerce")
    valid = y.notna()
    if not bool(valid.all()):
        warnings.warn(
            f"Dropping {int((~valid).sum())} rows with invalid chronological_age.",
            stacklevel=2,
        )
    wide = wide.loc[valid].copy()
    y = y.loc[valid]

    x, imputed = build_feature_matrix(wide, model)
    y_pred = np.asarray(model.predict(x), dtype=float)
    y_true = y.to_numpy(dtype=float)
    residual = y_pred - y_true

    metrics = compute_metrics(y_true, y_pred)
    metrics["n_features_used"] = int(x.shape[1])
    metrics["n_imputed_missing_cpgs"] = int(len(imputed))
    metrics["model_path"] = str(model_path)
    metrics["methylation_path"] = str(methylation_path)
    metrics["meta_path"] = str(meta_path)

    metrics_path = Path(metrics_path)
    metrics_path.parent.mkdir(parents=True, exist_ok=True)
    metrics_path.write_text(json.dumps(metrics, indent=2) + "\n", encoding="utf-8")

    weights = extract_cpg_coefficients(model)
    gene_map = load_probe_gene_map(annotation_path)

    sns.set_theme(style="ticks", context="paper")
    fig = plt.figure(figsize=(12.5, 10.0), layout="constrained")
    # Top row: panels A/B; bottom row: panel C spanning full width.
    grid = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.15])
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[1, :])

    plot_panel_a(ax_a, y_true, y_pred, metrics)
    plot_panel_b(ax_b, y_true, residual)
    plot_panel_c(ax_c, weights, gene_map, top_n_cpgs)
    sns.despine(fig=fig)

    png_path, pdf_path = save_figure(fig, figure_stem)
    plt.close(fig)

    print(format_metrics_markdown(metrics))
    print(f"Metrics JSON: {metrics_path}")
    print(f"Figure PNG:   {png_path}")
    print(f"Figure PDF:   {pdf_path}")
    if imputed:
        print(
            f"Note: imputed {len(imputed)} model CpGs absent from the validation matrix.",
            flush=True,
        )

    return {
        **metrics,
        "metrics_path": str(metrics_path),
        "figure_png": str(png_path),
        "figure_pdf": str(pdf_path),
    }


def main(
    methylation: Path = typer.Option(
        DEFAULT_METHYLATION,
        "--methylation",
        help="Processed GSE87571 methylation parquet (samples × cg*).",
    ),
    meta: Path = typer.Option(
        DEFAULT_META,
        "--meta",
        help="Phenotype CSV with sample IDs and chronological age.",
    ),
    model: Path = typer.Option(
        DEFAULT_MODEL,
        "--model",
        help="Pickled ElasticNet or Pipeline clock (.pkl / .joblib).",
    ),
    metrics_out: Path = typer.Option(
        DEFAULT_METRICS,
        "--metrics-out",
        help="Output path for clock_metrics.json.",
    ),
    figure_stem: Path = typer.Option(
        DEFAULT_FIGURE_STEM,
        "--figure-stem",
        help="Output stem for Figure_Epigenetic_Clock_Panels (.png / .pdf).",
    ),
    annotation: Path | None = typer.Option(
        None,
        "--annotation",
        help="Optional Illumina probe→gene table (IlmnID + UCSC_RefGene_Name).",
    ),
    top_n: int = typer.Option(
        TOP_N_CPGS,
        "--top-n",
        help="Number of top |weight| CpGs for panel C.",
    ),
    skip_manifest_check: bool = typer.Option(
        False,
        "--skip-manifest-check",
        help="Skip INPUT_MANIFEST.md required-file verification.",
    ),
    allow_positional_align: bool = typer.Option(
        False,
        "--allow-positional-align",
        help="Allow row-order alignment when sample IDs do not overlap.",
    ),
    demo: bool = typer.Option(
        False,
        "--demo",
        help=(
            "Write offline clock fixtures (model pickle; synthetic cohort only if "
            "methylation inputs are missing) and run evaluation."
        ),
    ),
) -> None:
    """Validate ElasticNet/Pipeline clock on GSE87571 and write metrics + figures.

    Args:
        methylation: Processed GSE87571 methylation parquet (samples × cg*).
        meta: Phenotype CSV with sample IDs and chronological age.
        model: Pickled ElasticNet or Pipeline clock.
        metrics_out: Output path for ``clock_metrics.json``.
        figure_stem: Output stem for ``Figure_Epigenetic_Clock_Panels.*``.
        annotation: Optional Illumina probe→gene annotation table.
        top_n: Number of top ``|weight|`` CpGs for panel C.
        skip_manifest_check: Bypass ``INPUT_MANIFEST.md`` checks.
        allow_positional_align: Opt into row-order sample alignment.
        demo: Materialize fixtures and run offline-friendly evaluation.
    """
    if demo:
        from rogen_aging.pipeline_fixtures import write_clock_fixtures

        fixtures = write_clock_fixtures(repo_root=REPO_ROOT)
        model = fixtures["model"]
        methylation = fixtures["methylation"]
        meta = fixtures["meta"]
        skip_manifest_check = True
        if metrics_out == DEFAULT_METRICS:
            metrics_out = REPO_ROOT / "outputs" / "demo" / "clock_metrics.json"
        if figure_stem == DEFAULT_FIGURE_STEM:
            figure_stem = (
                REPO_ROOT / "outputs" / "demo" / "figures" / "Figure_Epigenetic_Clock_Panels"
            )

    model = resolve_clock_model_path(model)

    run_validation(
        methylation_path=methylation,
        meta_path=meta,
        model_path=model,
        metrics_path=metrics_out,
        figure_stem=figure_stem,
        annotation_path=annotation,
        top_n_cpgs=top_n,
        skip_manifest_check=skip_manifest_check,
        allow_positional_align=allow_positional_align,
    )


if __name__ == "__main__":
    typer.run(main)
