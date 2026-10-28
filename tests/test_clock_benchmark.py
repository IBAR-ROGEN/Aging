"""Unit tests for published-clock benchmark helpers."""

from __future__ import annotations

from typing import Any
from unittest.mock import MagicMock

import pandas as pd
import pytest

from rogen_aging.clock.benchmark import (
    GRIMAGE2_DATA_FEATURES,
    META_COLS,
    assert_grimage2_data_features_present,
    metrics_for_clock,
    prepare_pyaging_frame,
)


def test_prepare_pyaging_frame_keeps_female_and_age_as_data_columns() -> None:
    frame = pd.DataFrame(
        {
            "chronological_age": [40.0, 50.0],
            "sex": ["F", "M"],
            "female": [1, 0],
            "ethnicity": ["A", "B"],
            "cg00000001": [0.1, 0.2],
        },
        index=["s1", "s2"],
    )
    prepared, meta_cols = prepare_pyaging_frame(frame)
    assert "female" in prepared.columns
    assert "age" in prepared.columns
    assert list(prepared["female"]) == [1.0, 0.0]
    assert list(prepared["age"]) == [40.0, 50.0]
    assert "female" not in meta_cols
    assert "age" not in meta_cols
    for col in META_COLS:
        if col in frame.columns:
            assert col in meta_cols
    assert GRIMAGE2_DATA_FEATURES.isdisjoint(meta_cols)


def test_assert_grimage2_rejects_female_or_age_in_missing_features() -> None:
    with pytest.raises(ValueError, match="female"):
        assert_grimage2_data_features_present(
            {"grimage2": {"missing_features": ["female", "cg00000001"]}}
        )
    with pytest.raises(ValueError, match="age"):
        assert_grimage2_data_features_present({"grimage2": {"missing_features": ["age"]}})
    # CpG-only missing list is allowed.
    assert_grimage2_data_features_present(
        {"grimage2": {"missing_features": ["cg00000001", "cg00000002"]}}
    )


def test_predict_pyaging_asserts_grimage2_data_features(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Any,
) -> None:
    import sys
    from pathlib import Path

    import rogen_aging.clock.benchmark as bench

    frame = pd.DataFrame(
        {
            "chronological_age": [40.0, 50.0],
            "sex": ["F", "M"],
            "female": [1, 0],
            "cg00000001": [0.1, 0.2],
        },
        index=["s1", "s2"],
    )

    captured: dict[str, Any] = {}

    class _FakeAdata:
        def __init__(self) -> None:
            self.obs = pd.DataFrame({"grimage2": [45.0, 55.0]}, index=frame.index)
            self.uns = {
                "grimage2_percent_na": 0.0,
                "grimage2_missing_features": ["female", "age"],
                "grimage2_metadata": {"clock_name": "grimage2"},
            }

    def fake_df_to_adata(df: pd.DataFrame, metadata_cols: list[str], **_kwargs: Any) -> _FakeAdata:
        captured["metadata_cols"] = list(metadata_cols)
        captured["columns"] = list(df.columns)
        return _FakeAdata()

    fake_pya = MagicMock()
    fake_pya.pp.df_to_adata = fake_df_to_adata
    fake_pya.pred.predict_age = MagicMock()
    monkeypatch.setitem(sys.modules, "pyaging", fake_pya)

    with pytest.raises(ValueError, match="female"):
        bench.predict_pyaging(frame, ["grimage2"], cache_dir=Path(tmp_path))

    assert "female" in captured["columns"]
    assert "age" in captured["columns"]
    assert "female" not in captured["metadata_cols"]
    assert "age" not in captured["metadata_cols"]


def test_dunedinpace_metrics_leave_mae_and_signed_error_empty() -> None:
    row = metrics_for_clock(
        clock="dunedinpace",
        cohort="GSE87571",
        chronological_age=__import__("numpy").asarray([40.0, 50.0]),
        predicted=__import__("numpy").asarray([1.0, 1.1]),
    )
    assert row["mae"] is None
    assert row["mean_signed_error"] is None
    assert row["pearson_r"] is not None
    assert "rate" in (row["mae_note"] or "")
