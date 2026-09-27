"""GSE40279 preparation, train-only probe filters, and the fixture guard."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl
import pytest
from sklearn.impute import SimpleImputer
from sklearn.linear_model import ElasticNet
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline

from rogen_aging.clock.evaluate import evaluate_clock
from rogen_aging.clock.gse40279 import Gse40279SampleCountError, prepare_gse40279
from rogen_aging.clock.hannum import select_training_probes


def _series_matrix(path: Path, n_samples: int) -> None:
    gsms = [f"GSM{i}" for i in range(n_samples)]
    ages = [str(30 + i) for i in range(n_samples)]
    lines = [
        "!Sample_geo_accession\t" + "\t".join(f'"{gsm}"' for gsm in gsms),
        "!Sample_characteristics_ch1\t" + "\t".join(f'"age: {age}"' for age in ages),
        "!series_matrix_table_begin",
        "ID_REF\t" + "\t".join(gsms),
        "cg00000001\t" + "\t".join("0.2" for _ in gsms),
        "cg00000002\t" + "\t".join("0.4" for _ in gsms),
        "rs0000001\t" + "\t".join("0.9" for _ in gsms),
        "!series_matrix_table_end",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def test_prepare_writes_declared_sample_count(tmp_path: Path) -> None:
    matrix = tmp_path / "tiny_series_matrix.txt"
    _series_matrix(matrix, n_samples=4)
    output = tmp_path / "GSE40279_processed.parquet"
    counts = prepare_gse40279(
        output,
        geo_cache_dir=tmp_path,
        series_matrix_path=matrix,
        expected_n_samples=4,
    )
    assert counts == {"n_samples": 4, "n_probes": 2}
    frame = pd.read_parquet(output)
    assert list(frame.columns[:2]) == ["sample_id", "chronological_age"]
    assert set(frame.columns[2:]) == {"cg00000001", "cg00000002"}
    assert frame["chronological_age"].dtype == np.float32
    assert frame["cg00000001"].dtype == np.float32
    assert len(frame) == 4


def test_prepare_stops_when_sample_count_differs(tmp_path: Path) -> None:
    matrix = tmp_path / "tiny_series_matrix.txt"
    _series_matrix(matrix, n_samples=4)
    output = tmp_path / "GSE40279_processed.parquet"
    with pytest.raises(Gse40279SampleCountError, match="difference -652"):
        prepare_gse40279(
            output,
            geo_cache_dir=tmp_path,
            series_matrix_path=matrix,
            expected_n_samples=656,
        )
    assert not output.exists()


def test_age_correlation_filter_uses_training_rows_only() -> None:
    n = 20
    indices = np.arange(n)
    train_idx, test_idx = train_test_split(indices, test_size=0.2, random_state=42)
    age = np.linspace(30.0, 70.0, n)
    probe_keep = age + 0.01
    probe_leak = np.zeros(n)
    probe_leak[test_idx] = age[test_idx]
    probe_missing = probe_keep.copy()
    probe_missing[train_idx[:6]] = np.nan
    frame = pd.DataFrame(
        {
            "chronological_age": age,
            "cg00000001": probe_keep,
            "cg00000002": probe_leak,
            "cg00000003": probe_missing,
        }
    )
    report = select_training_probes(
        frame,
        {"cg00000001", "cg00000002", "cg00000003"},
        test_size=0.2,
        random_state=42,
    )
    assert report.columns == ["cg00000001"]
    assert report.n_probes_overlap_gse87571 == 3
    assert report.n_train_samples == len(train_idx)
    assert report.n_test_samples == len(test_idx)


def test_evaluate_refuses_cg_test_fixture_without_demo(tmp_path: Path) -> None:
    rng = np.random.default_rng(0)
    x = pd.DataFrame(
        {
            "cg_test_0000": rng.uniform(0.1, 0.9, size=12),
            "cg_test_0001": rng.uniform(0.1, 0.9, size=12),
        }
    )
    y = 40.0 + x["cg_test_0000"]
    pipe = Pipeline(
        [
            ("imputer", SimpleImputer()),
            ("elasticnet", ElasticNet(alpha=0.01, max_iter=5000, random_state=0)),
        ]
    )
    pipe.fit(x, y)
    model_path = tmp_path / "fixture.joblib"
    import joblib

    joblib.dump(pipe, model_path)
    table = x.copy()
    table.insert(0, "chronological_age", y.to_numpy())
    data_path = tmp_path / "betas.parquet"
    table.to_parquet(data_path, index=False)
    with pytest.raises(ValueError, match="fixture"):
        evaluate_clock(model_path, data_path, tmp_path / "out")
    result = evaluate_clock(model_path, data_path, tmp_path / "demo", demo=True)
    assert Path(result["metrics_path"]).is_file()
    assert (tmp_path / "demo" / "Fig_Clock_Predicted_vs_Chronological.png").is_file()


def test_evaluate_saves_predictions_and_signed_error(tmp_path: Path) -> None:
    rng = np.random.default_rng(1)
    ages = np.linspace(20.0, 80.0, 12)
    x = pd.DataFrame({"cg00000001": rng.uniform(0.1, 0.9, size=12)})
    y = ages
    pipe = Pipeline(
        [
            ("imputer", SimpleImputer()),
            ("elasticnet", ElasticNet(alpha=0.01, max_iter=5000, random_state=0)),
        ]
    )
    pipe.fit(x, y)
    model_path = tmp_path / "clock.joblib"
    import joblib

    joblib.dump(pipe, model_path)
    table = x.copy()
    table.insert(0, "sample_id", [f"S{index:02d}" for index in range(len(table))])
    table.insert(1, "chronological_age", ages)
    data_path = tmp_path / "betas.parquet"
    table.to_parquet(data_path, index=False)
    predictions_path = tmp_path / "predictions.csv"
    result = evaluate_clock(
        model_path,
        data_path,
        tmp_path / "out",
        predictions_path=predictions_path,
    )
    saved = pl.read_csv(predictions_path)
    assert saved.columns == ["sample_id", "chronological_age", "predicted_age"]
    assert saved.height == 12
    signed = saved["predicted_age"] - saved["chronological_age"]
    assert result["mean_signed_error"] == pytest.approx(float(signed.mean()))
    metrics = json.loads(Path(result["metrics_path"]).read_text(encoding="utf-8"))
    assert "mean_signed_error_by_decade" in metrics
