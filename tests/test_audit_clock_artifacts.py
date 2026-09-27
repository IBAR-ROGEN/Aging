"""Provenance rules for the clock-artifact audit."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path
from types import ModuleType

_REPO_ROOT = Path(__file__).resolve().parent.parent
_SCRIPT = _REPO_ROOT / "scripts" / "audit_clock_artifacts.py"


def _load_script() -> ModuleType:
    name = "audit_clock_artifacts"
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, _SCRIPT)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load {_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


audit = _load_script()

FIXTURE = {
    "path": "models/ro_clock_elasticnet_gse40279.pkl",
    "n_features": 12,
    "feature_name_pattern": r"^cg_test_\d+$ (all 12 names)",
    "feature_flag": "FIXTURE",
    "n_nonzero_coefficients": 3,
    "intercept": 1.0,
    "pipeline_steps": "(bare estimator, no Pipeline)",
    "sklearn_version_recorded": "1.6.0",
    "mtime": "2026-01-01T00:00:00",
    "feature_names_head": ["cg_test_0000"],
    "loader": "pickle",
}
OTHER_FIXTURE = {
    **FIXTURE,
    "path": "models/methylation_clock_v1.joblib",
    "loader": "joblib",
}


def _hit(path: str, value: str, snippet: str) -> dict[str, str]:
    return {
        "file": path,
        "line": "1",
        "value": value,
        "label": "mae",
        "source": "working-tree file",
        "gse40279_context": "yes",
        "snippet": snippet,
    }


def test_working_tree_scan_skips_its_own_report(tmp_path: Path) -> None:
    report = tmp_path / "CLOCK_ARTIFACT_AUDIT.md"
    report.write_text('GSE40279 "test_mae": 1.23\n', encoding="utf-8")
    notes = tmp_path / "notes.md"
    notes.write_text("GSE40279 MAE ~ 4.56 years\n", encoding="utf-8")
    nested = tmp_path / "tests" / "notes.md"
    nested.parent.mkdir()
    nested.write_text("GSE40279 MAE ~ 8.00 years\n", encoding="utf-8")
    hits = audit.find_mae_in_working_tree(tmp_path)
    files = {hit["file"] for hit in hits}
    assert "CLOCK_ARTIFACT_AUDIT.md" not in files
    assert "notes.md" in files
    assert "tests/notes.md" not in files


def test_feature_count_alone_is_not_traceable(tmp_path: Path) -> None:
    payload = {"mae": 1.5, "n_features_used": 12}
    (tmp_path / "clock_metrics.json").write_text(json.dumps(payload), encoding="utf-8")
    verdict, reason = audit.verdict_for_mae(
        _hit("clock_metrics.json", "1.5", '"mae": 1.5'),
        [FIXTURE, OTHER_FIXTURE],
        root=tmp_path,
    )
    assert verdict == "UNTRACEABLE"
    assert "feature count" in reason


def test_model_path_and_feature_count_trace_one_fixture(tmp_path: Path) -> None:
    payload = {
        "mae": 3.5,
        "n_features_used": 12,
        "model_path": "/tmp/models/ro_clock_elasticnet_gse40279.pkl",
    }
    (tmp_path / "clock_metrics.json").write_text(json.dumps(payload), encoding="utf-8")
    verdict, reason = audit.verdict_for_mae(
        _hit("clock_metrics.json", "3.5", '"mae": 3.5'),
        [FIXTURE, OTHER_FIXTURE],
        root=tmp_path,
    )
    assert verdict == "TRACEABLE"
    assert FIXTURE["path"] in reason
    assert OTHER_FIXTURE["path"] not in reason
    assert "FIXTURE" in reason


def test_selected_cg_test_names_trace_matching_fixtures(tmp_path: Path) -> None:
    payload = {
        "test_mae": 3.3,
        "n_cpgs_features": 12,
        "selected_cpgs": ["cg_test_0000", "cg_test_0001"],
    }
    (tmp_path / "gse40279_train_metrics.json").write_text(json.dumps(payload), encoding="utf-8")
    verdict, reason = audit.verdict_for_mae(
        _hit("gse40279_train_metrics.json", "3.3", '"test_mae": 3.3'),
        [FIXTURE, OTHER_FIXTURE],
        root=tmp_path,
    )
    assert verdict == "TRACEABLE"
    assert FIXTURE["path"] in reason
    assert OTHER_FIXTURE["path"] in reason


def test_simulated_target_stays_untraceable_when_a_real_clock_exists() -> None:
    real = {
        **FIXTURE,
        "path": "models/real_clock.pkl",
        "feature_flag": "CANDIDATE REAL",
        "feature_name_pattern": r"^cg\d{8}$ (all 12 names)",
        "n_features": 12,
    }
    verdict, reason = audit.verdict_for_mae(
        _hit("scripts/figures/generate_clock_validation.py", "2.1", "MAE ~ 2.1 years"),
        [real],
    )
    assert verdict == "UNTRACEABLE"
    assert "simulated" in reason


def test_illumina_names_trace_the_real_clock(tmp_path: Path) -> None:
    real = {
        **FIXTURE,
        "path": "models/real_clock.pkl",
        "feature_flag": "CANDIDATE REAL",
        "n_features": 2,
    }
    payload = {
        "mae": 4.0,
        "n_features_used": 2,
        "selected_cpgs": ["cg00000001", "cg00000002"],
    }
    (tmp_path / "clock_metrics.json").write_text(json.dumps(payload), encoding="utf-8")
    verdict, reason = audit.verdict_for_mae(
        _hit("clock_metrics.json", "4.0", '"mae": 4.0'),
        [FIXTURE, real],
        root=tmp_path,
    )
    assert verdict == "TRACEABLE"
    assert reason == "models/real_clock.pkl"


def test_report_follows_inspected_artifacts() -> None:
    artifact = {
        **FIXTURE,
        "n_features": 4,
        "feature_name_pattern": "first='cg_test_0007' last='cg_test_0009' n=4",
        "feature_names_head": ["cg_test_0007"],
        "n_nonzero_coefficients": 1,
    }
    text = audit.render_report([artifact], [])
    assert "4 features" in text
    assert "cg_test_0000" not in text
    assert "cg_test_0007" in text


def test_classify_feature_names() -> None:
    assert audit.classify_feature_names(["cg_test_0000", "cg_test_0001"]) == "FIXTURE"
    assert audit.classify_feature_names(["cg00000001", "cg00000002"]) == "CANDIDATE REAL"
    assert audit.classify_feature_names(["cg_test_0000", "cg00000001"]) == "MIXED"
    assert audit.classify_feature_names([]) == "OTHER (no feature_names_in_)"
