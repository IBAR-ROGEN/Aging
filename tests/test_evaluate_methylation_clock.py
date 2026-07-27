"""Unit tests for ``evaluate_methylation_clock`` (GSE87571 ElasticNet validation)."""

from __future__ import annotations

import importlib.util
import json
import pickle
import sys
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
import pandas as pd
import pytest
from sklearn.linear_model import ElasticNet, ElasticNetCV
from sklearn.pipeline import Pipeline

_REPO_ROOT = Path(__file__).resolve().parent.parent
_SCRIPT = _REPO_ROOT / "scripts" / "clock" / "evaluate_methylation_clock.py"


def _load_script() -> ModuleType:
    """Import the repo-root CLI module (not on the default ``pythonpath``)."""
    name = "evaluate_methylation_clock"
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, _SCRIPT)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load {_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


emc = _load_script()


CPG_A = "cg00000001"
CPG_B = "cg00000002"
CPG_C = "cg00000003"


def _fit_elasticnet(rng: np.random.Generator | None = None) -> ElasticNet:
    """Fit a tiny bare ElasticNet with ``feature_names_in_`` set."""
    rng = rng or np.random.default_rng(0)
    x = pd.DataFrame(
        {
            CPG_A: rng.uniform(0.1, 0.9, size=40),
            CPG_B: rng.uniform(0.1, 0.9, size=40),
            CPG_C: rng.uniform(0.1, 0.9, size=40),
        }
    )
    y = 25.0 + 30.0 * x[CPG_A] + 10.0 * x[CPG_B] + rng.normal(0.0, 0.5, size=40)
    model = ElasticNet(alpha=0.01, l1_ratio=0.5, max_iter=10_000, random_state=0)
    model.fit(x, y)
    return model


def _write_pickle(path: Path, obj: Any) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as handle:
        pickle.dump(obj, handle)
    return path


def _write_cohort(
    tmp_path: Path,
    *,
    n: int = 12,
    include_missing_cpg: bool = False,
    age_col: str = "chronological_age",
    transpose_meth: bool = False,
    sample_ids: list[str] | None = None,
) -> tuple[Path, Path, ElasticNet]:
    """Write parquet + meta CSV fixtures and a matching ElasticNet pickle."""
    rng = np.random.default_rng(7)
    model = _fit_elasticnet(rng)
    ids = sample_ids or [f"GSM{i:04d}" for i in range(n)]
    ages = np.concatenate(
        [
            rng.uniform(18.0, 29.0, size=max(1, n // 3)),
            rng.uniform(30.0, 60.0, size=max(1, n // 3)),
            rng.uniform(61.0, 85.0, size=n - 2 * max(1, n // 3)),
        ]
    )[:n]

    cols = [CPG_A, CPG_B] if include_missing_cpg else [CPG_A, CPG_B, CPG_C]
    meth = pd.DataFrame(
        {c: rng.uniform(0.05, 0.95, size=n) for c in cols},
        index=pd.Index(ids, name="sample_id"),
    )
    meth = meth.reset_index()
    if transpose_meth:
        # probes × samples layout (no sample_id column; probes as index).
        probe_frame = meth.set_index("sample_id").T
        probe_frame.index.name = "probe_id"
        meth_path = tmp_path / "meth.parquet"
        probe_frame.to_parquet(meth_path)
    else:
        meth_path = tmp_path / "meth.parquet"
        meth.to_parquet(meth_path, index=False)

    meta = pd.DataFrame({"sample_id": ids, age_col: ages})
    meta_path = tmp_path / "meta.csv"
    meta.to_csv(meta_path, index=False)

    _write_pickle(tmp_path / "ro_clock.pkl", model)
    return meth_path, meta_path, model


# ---------------------------------------------------------------------------
# Manifest verification
# ---------------------------------------------------------------------------


def test_verify_input_manifest_ok(tmp_path: Path) -> None:
    required = ["a/one.parquet", "b/two.csv", "c/model.pkl"]
    for rel in required:
        path = tmp_path / rel
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("x", encoding="utf-8")

    manifest = tmp_path / "INPUT_MANIFEST.md"
    lines = [
        "| Path | Role | Required |",
        "|------|------|----------|",
        "| `a/one.parquet` | meth | yes |",
        "| `b/two.csv` | meta | yes |",
        "| `c/model.pkl` | model | yes |",
        "| `d/optional.csv` | annot | no |",
        "",
    ]
    manifest.write_text("\n".join(lines), encoding="utf-8")

    resolved = emc.verify_input_manifest(manifest, repo_root=tmp_path)
    assert len(resolved) == 3
    assert all(p.is_file() for p in resolved)


def test_verify_input_manifest_missing_file(tmp_path: Path) -> None:
    (tmp_path / "a").mkdir()
    (tmp_path / "a" / "one.parquet").write_text("x", encoding="utf-8")
    manifest = tmp_path / "INPUT_MANIFEST.md"
    manifest.write_text(
        "| Path | Role | Required |\n"
        "|------|------|----------|\n"
        "| `a/one.parquet` | meth | yes |\n"
        "| `missing/model.pkl` | model | yes |\n",
        encoding="utf-8",
    )
    with pytest.raises(FileNotFoundError, match="missing/model.pkl"):
        emc.verify_input_manifest(manifest, repo_root=tmp_path)


def test_verify_input_manifest_absent_manifest(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError, match="INPUT_MANIFEST.md not found"):
        emc.verify_input_manifest(tmp_path / "nope.md", repo_root=tmp_path)


def test_verify_input_manifest_no_required_rows(tmp_path: Path) -> None:
    manifest = tmp_path / "INPUT_MANIFEST.md"
    manifest.write_text("# empty\n\n| Path | Role | Required |\n|------|------|----------|\n", encoding="utf-8")
    with pytest.raises(ValueError, match="No required input paths"):
        emc.verify_input_manifest(manifest, repo_root=tmp_path)


# ---------------------------------------------------------------------------
# Model loading
# ---------------------------------------------------------------------------


def test_load_elasticnet_clock_ok(tmp_path: Path) -> None:
    model = _fit_elasticnet()
    path = _write_pickle(tmp_path / "clock.pkl", model)
    loaded = emc.load_elasticnet_clock(path)
    assert type(loaded) is ElasticNet
    assert hasattr(loaded, "coef_")
    np.testing.assert_allclose(loaded.coef_, model.coef_)


def test_load_elasticnet_clock_missing(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError, match="Trained ElasticNet model not found"):
        emc.load_elasticnet_clock(tmp_path / "absent.pkl")


def test_load_elasticnet_clock_accepts_pipeline(tmp_path: Path) -> None:
    from sklearn.impute import SimpleImputer

    pipe = Pipeline(
        [
            ("imputer", SimpleImputer(strategy="mean")),
            ("elasticnet", _fit_elasticnet()),
        ]
    )
    # Fit imputer on dummy data so statistics_ exist.
    rng = np.random.default_rng(0)
    x = pd.DataFrame(
        {
            CPG_A: rng.uniform(0.1, 0.9, size=20),
            CPG_B: rng.uniform(0.1, 0.9, size=20),
            CPG_C: rng.uniform(0.1, 0.9, size=20),
        }
    )
    y = 30.0 + 10.0 * x[CPG_A]
    pipe.fit(x, y)
    path = _write_pickle(tmp_path / "pipe.pkl", pipe)
    loaded = emc.load_elasticnet_clock(path)
    assert hasattr(loaded, "named_steps")
    assert "elasticnet" in loaded.named_steps


def test_load_elasticnet_clock_rejects_bare_elasticnetcv(tmp_path: Path) -> None:
    cv = ElasticNetCV(l1_ratio=[0.5], alphas=[0.1], cv=2)
    path = _write_pickle(tmp_path / "cv.pkl", cv)
    with pytest.raises(TypeError, match="ElasticNet"):
        emc.load_elasticnet_clock(path)


# ---------------------------------------------------------------------------
# Age strata and metrics
# ---------------------------------------------------------------------------


def test_assign_age_stratum_boundaries() -> None:
    ages = np.array([29.9, 30.0, 45.0, 60.0, 60.1])
    labels = emc.assign_age_stratum(ages)
    assert list(labels) == ["<30", "30-60", "30-60", "30-60", ">60"]


def test_compute_metrics_known_values() -> None:
    chrono = np.array([20.0, 40.0, 70.0])
    pred = np.array([22.0, 38.0, 65.0])
    metrics = emc.compute_metrics(chrono, pred)
    assert metrics["n_samples"] == 3
    assert metrics["mae"] == pytest.approx(3.0)
    assert metrics["median_absolute_error"] == pytest.approx(2.0)
    assert metrics["n_by_age_stratum"]["<30"] == 1
    assert metrics["n_by_age_stratum"]["30-60"] == 1
    assert metrics["n_by_age_stratum"][">60"] == 1
    assert metrics["mae_by_age_stratum"]["<30"] == pytest.approx(2.0)
    assert metrics["mae_by_age_stratum"]["30-60"] == pytest.approx(2.0)
    assert metrics["mae_by_age_stratum"][">60"] == pytest.approx(5.0)
    assert -1.0 <= metrics["pearson_r"] <= 1.0


def test_compute_metrics_empty_stratum() -> None:
    chrono = np.array([35.0, 40.0, 50.0])
    pred = np.array([34.0, 41.0, 49.0])
    metrics = emc.compute_metrics(chrono, pred)
    assert metrics["mae_by_age_stratum"]["<30"] is None
    assert metrics["mae_by_age_stratum"][">60"] is None
    assert metrics["n_by_age_stratum"]["<30"] == 0
    assert metrics["n_by_age_stratum"][">60"] == 0


def test_format_metrics_markdown_contains_strata() -> None:
    metrics = emc.compute_metrics(np.array([25.0, 50.0]), np.array([26.0, 48.0]))
    text = emc.format_metrics_markdown(metrics)
    assert "MAE:" in text
    assert "| <30 |" in text
    assert "| 30-60 |" in text
    assert "| >60 |" in text
    assert "NA" in text  # empty >60 stratum


# ---------------------------------------------------------------------------
# Coefficients and gene labels
# ---------------------------------------------------------------------------


def test_extract_cpg_coefficients_ok() -> None:
    model = _fit_elasticnet()
    weights = emc.extract_cpg_coefficients(model)
    assert list(weights.index) == [CPG_A, CPG_B, CPG_C]
    assert weights.name == "coefficient"


def test_extract_cpg_coefficients_unfitted() -> None:
    with pytest.raises(ValueError, match="no coef_"):
        emc.extract_cpg_coefficients(ElasticNet())


def test_extract_cpg_coefficients_length_mismatch() -> None:
    model = _fit_elasticnet()
    model.coef_ = np.array([1.0, 2.0])  # type: ignore[assignment]
    with pytest.raises(ValueError, match="does not match feature names"):
        emc.extract_cpg_coefficients(model)


@pytest.mark.parametrize(
    ("raw", "expected"),
    [
        ("TP53;MDM2", "TP53"),
        ("  BRCA1 ", "BRCA1"),
        (None, None),
        (float("nan"), None),
        ("nan", None),
        (".", None),
        ("", None),
    ],
)
def test_primary_gene_symbol(raw: object, expected: str | None) -> None:
    assert emc._primary_gene_symbol(raw) == expected


def test_label_cpg() -> None:
    assert emc.label_cpg(CPG_A, {CPG_A: "FOXO3"}) == f"FOXO3 ({CPG_A})"
    assert emc.label_cpg(CPG_A, {}) == CPG_A


def test_load_probe_gene_map_from_csv(tmp_path: Path) -> None:
    annot = tmp_path / "annot.csv"
    annot.write_text(
        "IlmnID,UCSC_RefGene_Name\n"
        f"{CPG_A},GENEA;GENEB\n"
        f"{CPG_B},\n",
        encoding="utf-8",
    )
    mapping = emc.load_probe_gene_map(annot)
    assert mapping[CPG_A] == "GENEA"
    assert CPG_B not in mapping


# ---------------------------------------------------------------------------
# Cohort loading
# ---------------------------------------------------------------------------


def test_load_validation_cohort_by_sample_id(tmp_path: Path) -> None:
    meth_path, meta_path, _model = _write_cohort(tmp_path, n=8)
    wide = emc.load_validation_cohort(meth_path, meta_path)
    assert "chronological_age" in wide.columns
    assert set(_cg_cols := emc._cg_columns(wide)) == {CPG_A, CPG_B, CPG_C}
    assert len(wide) == 8
    assert wide["chronological_age"].notna().all()
    del _cg_cols


def test_load_validation_cohort_age_alias(tmp_path: Path) -> None:
    meth_path, meta_path, _model = _write_cohort(tmp_path, n=6, age_col="Age_years")
    wide = emc.load_validation_cohort(meth_path, meta_path)
    assert "chronological_age" in wide.columns
    assert len(wide) == 6


def test_load_validation_cohort_transpose(tmp_path: Path) -> None:
    meth_path, meta_path, _model = _write_cohort(tmp_path, n=5, transpose_meth=True)
    wide = emc.load_validation_cohort(meth_path, meta_path)
    assert set(emc._cg_columns(wide)) == {CPG_A, CPG_B, CPG_C}
    assert len(wide) == 5


def test_load_validation_cohort_positional_fallback(
    tmp_path: Path,
) -> None:
    """Positional align requires explicit allow_positional_align=True."""
    rng = np.random.default_rng(1)
    n = 4
    # Methylation carries GSM IDs; metadata has only ages (RangeIndex → "0"…"3").
    meth = pd.DataFrame(
        {
            "sample_id": [f"GSM{i}" for i in range(n)],
            CPG_A: rng.random(n),
            CPG_B: rng.random(n),
            CPG_C: rng.random(n),
        }
    )
    meth_path = tmp_path / "meth.parquet"
    meth.to_parquet(meth_path, index=False)
    meta = pd.DataFrame({"chronological_age": [20.0, 40.0, 55.0, 70.0]})
    meta_path = tmp_path / "meta.csv"
    meta.to_csv(meta_path, index=False)

    with pytest.raises(ValueError, match="allow-positional-align"):
        emc.load_validation_cohort(meth_path, meta_path)

    with pytest.warns(UserWarning, match="aligning by row order"):
        wide = emc.load_validation_cohort(
            meth_path, meta_path, allow_positional_align=True
        )
    assert len(wide) == n
    np.testing.assert_allclose(wide["chronological_age"].to_numpy(), [20.0, 40.0, 55.0, 70.0])


def test_load_validation_cohort_no_cpg_raises(tmp_path: Path) -> None:
    meth = pd.DataFrame({"sample_id": ["s1"], "not_a_probe": [0.5]})
    meth_path = tmp_path / "meth.parquet"
    meth.to_parquet(meth_path, index=False)
    meta = pd.DataFrame({"sample_id": ["s1"], "chronological_age": [40.0]})
    meta_path = tmp_path / "meta.csv"
    meta.to_csv(meta_path, index=False)
    with pytest.raises(ValueError, match="No CpG columns"):
        emc.load_validation_cohort(meth_path, meta_path)


def test_load_validation_cohort_missing_paths(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError, match="Methylation matrix not found"):
        emc.load_validation_cohort(tmp_path / "no.parquet", tmp_path / "no.csv")


# ---------------------------------------------------------------------------
# End-to-end validation
# ---------------------------------------------------------------------------


def test_run_validation_writes_outputs(tmp_path: Path) -> None:
    meth_path, meta_path, model = _write_cohort(tmp_path, n=12, include_missing_cpg=True)
    model_path = tmp_path / "ro_clock.pkl"
    metrics_path = tmp_path / "out" / "clock_metrics.json"
    figure_stem = tmp_path / "out" / "Figure_Epigenetic_Clock_Panels"
    annot = tmp_path / "annot.csv"
    annot.write_text(
        f"IlmnID,UCSC_RefGene_Name\n{CPG_A},TESTGENE\n",
        encoding="utf-8",
    )

    result = emc.run_validation(
        methylation_path=meth_path,
        meta_path=meta_path,
        model_path=model_path,
        metrics_path=metrics_path,
        figure_stem=figure_stem,
        annotation_path=annot,
        top_n_cpgs=2,
        skip_manifest_check=True,
    )

    assert metrics_path.is_file()
    assert Path(result["figure_png"]).is_file()
    assert Path(result["figure_pdf"]).is_file()
    on_disk = json.loads(metrics_path.read_text(encoding="utf-8"))
    assert on_disk["n_samples"] == 12
    assert on_disk["n_features_used"] == 3
    assert on_disk["n_imputed_missing_cpgs"] == 1  # CPG_C absent from matrix
    assert "mae" in on_disk
    assert "pearson_r" in on_disk
    assert set(on_disk["mae_by_age_stratum"]) == {"<30", "30-60", ">60"}
    # Model was fitted; coefficients recoverable.
    weights = emc.extract_cpg_coefficients(model)
    assert len(weights) == 3


def test_run_validation_respects_manifest_check(tmp_path: Path) -> None:
    meth_path, meta_path, _model = _write_cohort(tmp_path, n=6)
    model_path = tmp_path / "ro_clock.pkl"
    # Force the default manifest check against the real INPUT_MANIFEST (likely missing
    # production data in CI). With skip_manifest_check=False this should raise when
    # required production artifacts are absent — skip only if they happen to exist.
    if not emc.DEFAULT_MODEL.is_file():
        with pytest.raises(FileNotFoundError):
            emc.run_validation(
                methylation_path=meth_path,
                meta_path=meta_path,
                model_path=model_path,
                metrics_path=tmp_path / "m.json",
                figure_stem=tmp_path / "fig",
                annotation_path=None,
                top_n_cpgs=3,
                skip_manifest_check=False,
            )


def test_save_figure(tmp_path: Path) -> None:
    import matplotlib.pyplot as plt

    fig, _ax = plt.subplots()
    png, pdf = emc.save_figure(fig, tmp_path / "panel" / "stem")
    plt.close(fig)
    assert png.is_file()
    assert pdf.is_file()
    assert png.suffix == ".png"
    assert pdf.suffix == ".pdf"
