"""Offline tests for QTD→QTS HTML listing mapper (synthetic fixtures)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import ModuleType

import polars as pl
import typer
from typer.testing import CliRunner

_REPO = Path(__file__).resolve().parent.parent
_SCRIPT = _REPO / "scripts" / "eqtl" / "map_qts_to_qtd.py"
_SYNTHETIC_HTML = _REPO / "tests" / "fixtures" / "eqtl_catalogue_synthetic" / "html"
_DATASETS = _REPO / "analysis" / "eqtl" / "eqtl_datasets_used.tsv"


def _load() -> ModuleType:
    name = "map_qts_to_qtd"
    sys.modules.pop(name, None)
    spec = importlib.util.spec_from_file_location(name, _SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def _app(mod: ModuleType) -> typer.Typer:
    application = typer.Typer()
    application.command()(mod.main)
    return application


def test_parse_directory_hrefs_extracts_qts_and_files() -> None:
    mod = _load()
    html = (_SYNTHETIC_HTML / "sumstats_index.html").read_text(encoding="utf-8")
    names = mod.parse_directory_hrefs(html)
    assert "QTS000999" in names
    assert "QTS000998" in names


def test_dry_run_maps_14_of_14(tmp_path: Path) -> None:
    mod = _load()
    out = tmp_path / "eqtl_datasets_used_with_qts.tsv"
    provenance = tmp_path / "prov.json"
    result = CliRunner().invoke(
        _app(mod),
        [
            "--dry-run",
            "--datasets",
            str(_DATASETS),
            "--out",
            str(out),
            "--provenance-out",
            str(provenance),
        ],
    )
    assert result.exit_code == 0, result.output
    table = pl.read_csv(out, separator="\t")
    assert table.height == 14
    assert {"dataset_id", "qts_id", "file_type"}.issubset(set(table.columns))
    assert table.filter(pl.col("qts_id").is_null()).height == 0
    assert set(table["qts_id"].to_list()) == {"QTS000999"}
    assert set(table["file_type"].to_list()) == {".all.tsv.gz"}
    assert provenance.is_file()


def test_from_saved_html_mode(tmp_path: Path) -> None:
    mod = _load()
    out = tmp_path / "with_qts.tsv"
    provenance = tmp_path / "prov.json"
    result = CliRunner().invoke(
        _app(mod),
        [
            "--from-saved-html",
            str(_SYNTHETIC_HTML),
            "--datasets",
            str(_DATASETS),
            "--out",
            str(out),
            "--provenance-out",
            str(provenance),
        ],
    )
    assert result.exit_code == 0, result.output
    table = pl.read_csv(out, separator="\t")
    assert table.height == 14
