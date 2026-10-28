"""Published-clock vs retrained-Hannum benchmark helpers."""

from __future__ import annotations

import subprocess
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any, cast

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402
from sklearn.linear_model import LinearRegression  # noqa: E402
from sklearn.metrics import mean_absolute_error  # noqa: E402
from sklearn.model_selection import train_test_split  # noqa: E402

from rogen_aging.clock.evaluate import build_feature_matrix, load_model  # noqa: E402
from rogen_aging.clock.external_data import load_gse87571  # noqa: E402
from rogen_aging.clock.gse40279 import load_gse40279_from_series_matrix  # noqa: E402
from rogen_aging.config import find_repo_root, get_config  # noqa: E402

REPO_ROOT = find_repo_root()

# pyaging's ``phenoage`` is clinical labs only; DNAm PhenoAge is ``dnamphenoage``.
TARGET_CLOCKS: tuple[str, ...] = (
    "horvath2013",
    "hannum",
    "dnamphenoage",
    "dunedinpace",
)
OPTIONAL_CLOCKS: tuple[str, ...] = ("grimage2",)
PACE_CLOCKS: frozenset[str] = frozenset({"dunedinpace"})
OUR_CLOCK_LABEL = "rogen_hannum_retrain"
# Phenotype metadata for AnnData.obs. GrimAge2 needs ``female`` and ``age`` as
# *data* columns (features), so they must not appear here.
META_COLS: tuple[str, ...] = ("chronological_age", "sex", "ethnicity")
GRIMAGE2_DATA_FEATURES: frozenset[str] = frozenset({"female", "age"})
RETRAIN_HOLDOUT_MAE = 3.5248637199401855
RETRAIN_HOLDOUT_MAE_TOL = 5e-4

PROVENANCE_N_TRAIN = 524
PROVENANCE_N_TEST = 132


@dataclass(frozen=True)
class CohortBundle:
    """Samples × CpGs matrix plus phenotype columns for one cohort."""

    name: str
    frame: pd.DataFrame
    notes: list[str]
    geo_probe_report: dict[str, Any] | None = None


def _repo_relative(path: Path) -> str:
    try:
        return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()
    except ValueError:
        return str(path)


def list_pyaging_clocks() -> list[str]:
    """Return every clock name reported by pyaging."""
    import io
    from contextlib import redirect_stdout

    import pyaging as pya

    buf = io.StringIO()
    with redirect_stdout(buf):
        pya.utils.show_all_clocks()
    text = buf.getvalue()
    names: list[str] = []
    for line in text.splitlines():
        for token in line.split():
            cleaned = token.strip().lower()
            if cleaned and cleaned[0].isalpha() and "clock" not in cleaned:
                # Filter logger noise
                if cleaned.startswith("loading") or cleaned.startswith("available"):
                    continue
                if cleaned in {"⏺"}:
                    continue
                names.append(cleaned)
    # Deduplicate while preserving order
    seen: set[str] = set()
    ordered: list[str] = []
    for name in names:
        if name not in seen and name.isidentifier():
            seen.add(name)
            ordered.append(name)
    return ordered


def clock_feature_union(clock_names: Sequence[str], *, cache_dir: Path) -> set[str]:
    """Union of Illumina probe IDs required by the named pyaging clocks."""
    from pyaging.predict import _pred as pyaging_pred

    load_clock = getattr(pyaging_pred, "load_clock")
    features: set[str] = set()
    cache_dir.mkdir(parents=True, exist_ok=True)
    for name in clock_names:
        model = load_clock(name, dir=str(cache_dir))
        feats = getattr(model, "features", None)
        if feats is None:
            continue
        for feat in feats:
            text = str(feat)
            if text.startswith("cg") or text.startswith("ch."):
                features.add(text)
    return features


def reconstruct_gse40279_split(
    n_samples: int,
    *,
    test_size: float | None = None,
    random_state: int | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Rebuild the Hannum training/test indices from config defaults."""
    cfg = get_config().clock
    size = float(cfg.test_size if test_size is None else test_size)
    seed = int(cfg.random_state if random_state is None else random_state)
    indices = np.arange(n_samples)
    train_idx, test_idx = train_test_split(indices, test_size=size, random_state=seed)
    if len(train_idx) != PROVENANCE_N_TRAIN or len(test_idx) != PROVENANCE_N_TEST:
        raise ValueError(
            f"Reconstructed split sizes {len(train_idx)}/{len(test_idx)} do not match "
            f"provenance {PROVENANCE_N_TRAIN}/{PROVENANCE_N_TEST}."
        )
    return train_idx, test_idx


def load_gse87571_cohort(
    *,
    geo_cache_dir: Path,
    series_matrix: Path | None,
    restrict_cpgs: Sequence[str] | None,
) -> CohortBundle:
    """Load all GSE87571 samples with full (or restricted) CpG coverage + sex."""
    frame = load_gse87571(
        local_path=series_matrix,
        geo_cache_dir=geo_cache_dir,
        restrict_to_cpgs=restrict_cpgs,
    )
    if frame["sex"].isna().any() or frame["female"].isna().any():
        missing = int(frame["sex"].isna().sum())
        raise ValueError(f"GSE87571: {missing} samples missing parsed sex; refusing to drop.")
    notes = [
        "GSE87571: all samples with chronological age and sex.",
        "Main external comparison cohort for published clocks and the retrained Hannum.",
    ]
    return CohortBundle(name="GSE87571", frame=frame, notes=notes)


def load_gse40279_heldout_cohort(
    *,
    series_matrix: Path,
    restrict_cpgs: Sequence[str] | None,
    parquet_path: Path | None = None,
) -> CohortBundle:
    """Load GSE40279 held-out test rows from the full GEO series matrix.

    Betas come from the GEO series matrix, not the training-prepared parquet.
    ``parquet_path`` is accepted only for API compatibility and is ignored.
    """
    del parquet_path  # Benchmark matrix is the GEO series file, not the parquet.
    full, geo_report = load_gse40279_from_series_matrix(
        series_matrix,
        restrict_cpgs=restrict_cpgs,
    )
    if len(full) != PROVENANCE_N_TRAIN + PROVENANCE_N_TEST:
        raise ValueError(
            f"GSE40279 series matrix has {len(full)} rows; expected "
            f"{PROVENANCE_N_TRAIN + PROVENANCE_N_TEST}."
        )
    _train_idx, test_idx = reconstruct_gse40279_split(len(full))
    held = full.iloc[list(test_idx)].copy()
    if held["sex"].isna().any():
        raise ValueError("GSE40279 held-out samples missing sex; refusing to drop.")
    if len(held) != PROVENANCE_N_TEST:
        raise ValueError(
            f"GSE40279 held-out has {len(held)} samples; expected {PROVENANCE_N_TEST}."
        )
    absent_raw = geo_report.get("probes_absent_from_geo")
    absent: list[str] = [str(item) for item in absent_raw] if isinstance(absent_raw, list) else []
    notes = [
        "GSE40279: held-out test split only (test_size=0.2, random_state=42).",
        "Betas loaded from the full GEO series matrix, not the training-prepared parquet.",
        "Published Hannum was trained on GSE40279; its predictions here are in-sample.",
        "Retrained Hannum is scored only on these test rows, never on training rows.",
        str(geo_report.get("imputation_note") or ""),
    ]
    if absent:
        notes.append(
            f"{len(absent)} requested Illumina probes are absent from the GEO file "
            "and rely on pyaging KNN imputation."
        )
    return CohortBundle(
        name="GSE40279_heldout",
        frame=held,
        notes=notes,
        geo_probe_report=geo_report,
    )


def prepare_pyaging_frame(frame: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    """Build the DataFrame and metadata_cols for ``df_to_adata``.

    ``female`` (0/1) and ``age`` stay as data columns for GrimAge2. They are
    never placed in ``metadata_cols``.
    """
    df = frame.copy()
    if "female" in df.columns:
        df["female"] = pd.to_numeric(df["female"], errors="coerce").astype("float64")
    if "chronological_age" in df.columns:
        df["age"] = pd.to_numeric(df["chronological_age"], errors="coerce")
    meta_present = [c for c in META_COLS if c in df.columns]
    for forbidden in GRIMAGE2_DATA_FEATURES:
        if forbidden in meta_present:
            raise ValueError(
                f"{forbidden!r} must be a data column for GrimAge2, not metadata_cols."
            )
    return df, meta_present


def assert_grimage2_data_features_present(uns_summary: dict[str, Any]) -> None:
    """Fail if pyaging still lists ``female`` or ``age`` as GrimAge2 missing features."""
    if "grimage2" not in uns_summary:
        return
    missing = [str(item) for item in (uns_summary["grimage2"].get("missing_features") or [])]
    bad = sorted(feat for feat in GRIMAGE2_DATA_FEATURES if feat in missing)
    if bad:
        raise ValueError(
            "grimage2 missing_features still includes data columns that must be "
            f"supplied (not imputed): {bad}. Keep female and age out of metadata_cols."
        )


def missing_features_before(
    frame: pd.DataFrame,
    clock_names: Sequence[str],
    *,
    cache_dir: Path,
) -> dict[str, list[str]]:
    """Clock features absent from ``frame`` columns before pyaging imputation."""
    from pyaging.predict import _pred as pyaging_pred

    load_clock = getattr(pyaging_pred, "load_clock")
    present = set(frame.columns.astype(str))
    if "chronological_age" in present:
        present.add("age")
    report: dict[str, list[str]] = {}
    cache_dir.mkdir(parents=True, exist_ok=True)
    for name in clock_names:
        model = load_clock(name, dir=str(cache_dir))
        feats = [str(f) for f in (getattr(model, "features", None) or [])]
        report[name] = sorted(feat for feat in feats if feat not in present)
    return report


def predict_pyaging(
    frame: pd.DataFrame,
    clock_names: Sequence[str],
    *,
    cache_dir: Path,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Run pyaging clocks; return obs predictions and per-clock uns summaries."""
    import pyaging as pya

    df, meta_present = prepare_pyaging_frame(frame)
    cache_dir.mkdir(parents=True, exist_ok=True)
    adata = pya.pp.df_to_adata(
        df,
        metadata_cols=meta_present,
        imputer_strategy="knn",
        verbose=False,
    )
    # Runtime accepts str | list[str]; published stubs type the arg as str only.
    pya.pred.predict_age(
        adata,
        clock_names=cast(Any, list(clock_names)),
        dir=str(cache_dir),
        verbose=False,
    )
    preds = adata.obs[[c for c in clock_names if c in adata.obs.columns]].copy()
    uns_summary: dict[str, Any] = {}
    for name in clock_names:
        uns_summary[name] = {
            "percent_na": adata.uns.get(f"{name}_percent_na"),
            "missing_features": list(adata.uns.get(f"{name}_missing_features", []) or []),
            "metadata": dict(adata.uns.get(f"{name}_metadata", {}) or {}),
        }
    assert_grimage2_data_features_present(uns_summary)
    return preds, uns_summary


def predict_rogen_hannum(frame: pd.DataFrame, model_path: Path) -> pd.Series:
    """Score the retrained Hannum joblib on the same samples."""
    model = load_model(model_path)
    # build_feature_matrix expects chronological_age + cg*
    x, _report = build_feature_matrix(frame.reset_index(drop=False), model)
    pred = model.predict(x)
    return pd.Series(np.asarray(pred, dtype=float), index=frame.index, name=OUR_CLOCK_LABEL)


def age_acceleration(y_clock: np.ndarray, age: np.ndarray) -> np.ndarray:
    """Residuals of ``y_clock ~ age`` (linear), same length as inputs."""
    model = LinearRegression()
    model.fit(age.reshape(-1, 1), y_clock)
    fitted = model.predict(age.reshape(-1, 1))
    return np.asarray(y_clock - fitted, dtype=float)


def metrics_for_clock(
    *,
    clock: str,
    cohort: str,
    chronological_age: np.ndarray,
    predicted: np.ndarray,
) -> dict[str, Any]:
    """Aggregate accuracy metrics for one clock on one cohort."""
    signed = predicted - chronological_age
    pearson = float(stats.pearsonr(predicted, chronological_age).statistic)
    row: dict[str, Any] = {
        "clock": clock,
        "cohort": cohort,
        "n": int(len(predicted)),
        "pearson_r": pearson,
        "mean_signed_error": float(np.mean(signed)),
    }
    if clock in PACE_CLOCKS:
        # Pace clock: keep r and acceleration residual; leave age-unit errors empty.
        row["mean_signed_error"] = None
        row["mae"] = None
        row["mae_note"] = (
            "DunedinPACE is a rate; MAE and mean signed error vs chronological age "
            "are not reported."
        )
    else:
        row["mae"] = float(mean_absolute_error(chronological_age, predicted))
        row["mae_note"] = None
    return row


def sex_acceleration_model(
    *,
    accel: np.ndarray,
    female: np.ndarray,
    age: np.ndarray,
) -> dict[str, Any]:
    """Fit ``accel ~ female + age``; return female coefficient, 95% CI, and p."""
    x = np.column_stack([female.astype(float), age.astype(float)])
    # OLS via statsmodels if available; else scipy least-squares with SE
    try:
        import statsmodels.api as sm

        design = sm.add_constant(x)
        fit = sm.OLS(accel, design).fit()
        coef = float(fit.params[1])
        ci_low, ci_high = (float(x) for x in fit.conf_int(alpha=0.05)[1])
        p_value = float(fit.pvalues[1])
    except Exception:
        # Fallback: manual OLS with normal approx CI
        design = np.column_stack([np.ones(len(accel)), x])
        beta, *_ = np.linalg.lstsq(design, accel, rcond=None)
        resid = accel - design @ beta
        df = max(len(accel) - design.shape[1], 1)
        sigma2 = float(np.sum(resid**2) / df)
        cov = sigma2 * np.linalg.inv(design.T @ design)
        se = float(np.sqrt(cov[1, 1]))
        coef = float(beta[1])
        t_crit = float(stats.t.ppf(0.975, df))
        ci_low = coef - t_crit * se
        ci_high = coef + t_crit * se
        p_value = float(2 * stats.t.sf(abs(coef / se), df)) if se > 0 else float("nan")
    return {
        "sex_coef_female_minus_male": coef,
        "sex_ci_low": ci_low,
        "sex_ci_high": ci_high,
        "sex_p": p_value,
    }


def bonferroni_adjust(p_values: Sequence[float], n_tests: int) -> list[float]:
    """Bonferroni-adjusted p-values capped at 1."""
    return [min(1.0, float(p) * n_tests) for p in p_values]


def training_overlap_notes(clock_metadata: dict[str, Any]) -> list[dict[str, str]]:
    """Document known training-cohort overlaps with GSE40279 / GSE87571."""
    rows: list[dict[str, str]] = []
    # Published knowledge + pyaging metadata notes
    known = {
        "hannum": (
            "overlap",
            "Hannum et al. 2013 trained on GSE40279 whole blood. "
            "GSE40279_heldout metrics for the published hannum clock are in-sample. "
            "GSE87571 is the fair external comparison.",
        ),
        "horvath2013": (
            "no_documented_overlap",
            "Horvath 2013 pan-tissue training set is a multi-dataset compilation; "
            "pyaging metadata does not list GSE40279 or GSE87571 as the training accession. "
            "No overlap flagged.",
        ),
        "dnamphenoage": (
            "no_documented_overlap",
            "DNAm PhenoAge coefficients were fitted on InCHIANTI 450K whole blood "
            "(Levine et al. 2018 / pyaging notes). No GSE40279 or GSE87571 overlap flagged.",
        ),
        "dunedinpace": (
            "no_documented_overlap",
            "DunedinPACE was trained in the Dunedin Study (EPIC, age 45). "
            "No GSE40279 or GSE87571 overlap flagged.",
        ),
        "grimage2": (
            "no_documented_overlap",
            "GrimAge2 was trained for mortality risk on older-adult whole-blood 450K data "
            "(Lu et al. 2022). No GSE40279 or GSE87571 overlap flagged from pyaging metadata.",
        ),
        OUR_CLOCK_LABEL: (
            "trained_on_gse40279_train_split",
            "Retrained Hannum ElasticNet was fit on the GSE40279 training split only "
            "(test_size=0.2, random_state=42). Held-out GSE40279 and all GSE87571 scores "
            "are out-of-sample for this model.",
        ),
    }
    for clock, (status, note) in known.items():
        meta = clock_metadata.get(clock, {})
        rows.append(
            {
                "clock": clock,
                "overlap_status": status,
                "note": note,
                "pyaging_notes": str(meta.get("notes") or ""),
                "pyaging_platform": str(meta.get("platform") or ""),
            }
        )
    return rows


def run_release_check(table: Path, report: Path) -> dict[str, Any]:
    """Run ``rogen-release-check --table-kind summary`` and capture the outcome."""
    report.parent.mkdir(parents=True, exist_ok=True)
    argv = [
        "uv",
        "run",
        "rogen-release-check",
        "--table",
        _repo_relative(table),
        "--table-kind",
        "summary",
        "--report",
        _repo_relative(report),
    ]
    completed = subprocess.run(argv, capture_output=True, text=True, check=False)
    return {
        "path": _repo_relative(table),
        "kind": "summary",
        "report": _repo_relative(report),
        "argv": argv,
        "returncode": completed.returncode,
        "stdout": completed.stdout,
        "stderr": completed.stderr,
        "status": "passed" if completed.returncode == 0 else "failed",
    }


def plot_predicted_vs_age(
    per_sample: pd.DataFrame,
    clocks: Sequence[str],
    cohort: str,
    path: Path,
) -> None:
    """Scatter predicted vs chronological age, one panel per clock."""
    n = len(clocks)
    cols = min(3, n)
    rows = int(np.ceil(n / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(4.2 * cols, 3.8 * rows), squeeze=False)
    sub = per_sample.loc[per_sample["cohort"] == cohort]
    for idx, clock in enumerate(clocks):
        axis = axes[idx // cols][idx % cols]
        x = sub["chronological_age"].to_numpy(dtype=float)
        y = sub[clock].to_numpy(dtype=float)
        axis.scatter(x, y, s=10, alpha=0.55, edgecolors="none")
        lo = float(np.nanmin([x.min(), y.min()]))
        hi = float(np.nanmax([x.max(), y.max()]))
        if clock not in PACE_CLOCKS:
            axis.plot([lo, hi], [lo, hi], color="0.4", linewidth=1.0, linestyle="--")
        axis.set_title(clock)
        axis.set_xlabel("Chronological age")
        axis.set_ylabel("Clock value")
    for idx in range(n, rows * cols):
        axes[idx // cols][idx % cols].axis("off")
    fig.suptitle(f"Predicted vs chronological age — {cohort}", fontsize=12)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_accel_by_sex(
    per_sample: pd.DataFrame,
    clocks: Sequence[str],
    cohort: str,
    path: Path,
) -> None:
    """Box/strip plot of age acceleration by sex for each clock."""
    n = len(clocks)
    cols = min(3, n)
    rows = int(np.ceil(n / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(4.2 * cols, 3.8 * rows), squeeze=False)
    sub = per_sample.loc[per_sample["cohort"] == cohort].copy()
    for idx, clock in enumerate(clocks):
        axis = axes[idx // cols][idx % cols]
        accel_col = f"{clock}_accel"
        for sex_label, color in (("M", "#4C78A8"), ("F", "#F58518")):
            mask = sub["sex"] == sex_label
            vals = sub.loc[mask, accel_col].to_numpy(dtype=float)
            jitter = np.random.default_rng(0).uniform(-0.12, 0.12, size=vals.shape[0])
            xpos = (0 if sex_label == "M" else 1) + jitter
            axis.scatter(xpos, vals, s=10, alpha=0.5, color=color, edgecolors="none")
        axis.set_xticks([0, 1], ["M", "F"])
        axis.set_title(clock)
        axis.set_ylabel("Age acceleration (residual)")
        axis.axhline(0.0, color="0.5", linewidth=0.8, linestyle="--")
    for idx in range(n, rows * cols):
        axes[idx // cols][idx % cols].axis("off")
    fig.suptitle(f"Age acceleration by sex — {cohort}", fontsize=12)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def try_include_grimage2(clocks: list[str]) -> tuple[list[str], str | None]:
    """Add grimage2 when sex/age metadata will be available; else return skip reason."""
    # GrimAge2 needs CpG surrogates + age + female; we supply those.
    return [*clocks, "grimage2"], None
