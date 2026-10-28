"""Rebuild predictions from the Hannum coefficient CSV and match the joblib."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from sklearn.impute import SimpleImputer
from sklearn.linear_model import ElasticNet
from sklearn.pipeline import Pipeline

from rogen_aging.clock.hannum import (
    INTERCEPT_PROBE_ID,
    load_probe_chromosomes,
    predict_from_coefficient_table,
    write_reproducible_coefficient_table,
)

REPO_ROOT = Path(__file__).resolve().parent.parent
MAIN_MODEL = REPO_ROOT / "models" / "gse40279_hannum450k_elasticnet.joblib"
MAIN_COEFF = REPO_ROOT / "models" / "gse40279_hannum450k_elasticnet_coefficients.csv"


def _tiny_pipeline(tmp_path: Path) -> tuple[Path, pd.DataFrame]:
    """Fit a two-probe mean-imputer + ElasticNet and return model path + matrix."""
    rng = np.random.default_rng(0)
    n = 40
    age = np.linspace(20.0, 80.0, n)
    cg_a = age / 100.0 + rng.normal(0.0, 0.01, size=n)
    cg_b = 1.0 - age / 100.0 + rng.normal(0.0, 0.01, size=n)
    cg_a[0] = np.nan
    frame = pd.DataFrame({"cg00000001": cg_a, "cg00000002": cg_b})
    pipe = Pipeline(
        [
            ("imputer", SimpleImputer(strategy="mean")),
            ("elasticnet", ElasticNet(alpha=0.01, l1_ratio=0.5, max_iter=5000, random_state=0)),
        ]
    )
    pipe.fit(frame, age)
    model_path = tmp_path / "tiny.joblib"
    import joblib

    joblib.dump(pipe, model_path)
    return model_path, frame


def test_coefficient_csv_rebuilds_tiny_pipeline(tmp_path: Path) -> None:
    model_path, frame = _tiny_pipeline(tmp_path)
    csv_path = tmp_path / "coefficients.csv"
    nonzero = write_reproducible_coefficient_table(model_path, csv_path)
    assert nonzero >= 1
    table = pd.read_csv(csv_path)
    assert INTERCEPT_PROBE_ID in set(table["probe_id"])
    assert {"probe_id", "coefficient", "imputer_mean"} <= set(table.columns)

    import joblib

    model = joblib.load(model_path)
    expected = np.asarray(model.predict(frame), dtype=float)
    rebuilt = predict_from_coefficient_table(frame, csv_path)
    np.testing.assert_allclose(rebuilt, expected, rtol=0.0, atol=1e-6)


def test_load_probe_chromosomes_maps_probe_to_chr(tmp_path: Path) -> None:
    manifest = tmp_path / "zhou.tsv"
    manifest.write_text(
        "CpG_chrm\tprobeID\nchrX\tcg00000001\nchr1\tcg00000002\n",
        encoding="utf-8",
    )
    mapping = load_probe_chromosomes(manifest)
    assert mapping == {"cg00000001": "chrX", "cg00000002": "chr1"}


@pytest.mark.skipif(
    not MAIN_MODEL.is_file() or not MAIN_COEFF.is_file(), reason="main model CSV absent"
)
def test_main_model_coefficient_csv_matches_joblib_on_fixture_rows(tmp_path: Path) -> None:
    """Use a few synthetic rows on the recorded feature names (imputed means fill gaps)."""
    import joblib

    model = joblib.load(MAIN_MODEL)
    names = [str(n) for n in model.feature_names_in_]
    rng = np.random.default_rng(1)
    values = rng.random((5, len(names)))
    values[0, :10] = np.nan
    frame = pd.DataFrame(values, columns=names)
    expected = np.asarray(model.predict(frame), dtype=float)
    rebuilt = predict_from_coefficient_table(frame, MAIN_COEFF)
    np.testing.assert_allclose(rebuilt, expected, rtol=0.0, atol=1e-6)
