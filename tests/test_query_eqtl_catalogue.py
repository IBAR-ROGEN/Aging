"""Offline tests for eQTL Catalogue tabix query (synthetic dry-run)."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path
from types import ModuleType

import polars as pl
import typer
from typer.testing import CliRunner

_REPO = Path(__file__).resolve().parent.parent
_SCRIPT = _REPO / "scripts" / "eqtl" / "query_eqtl_catalogue.py"
_PANEL = _REPO / "data" / "processed" / "variants_47_input_v2.csv"


def _load() -> ModuleType:
    name = "query_eqtl_catalogue"
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


def test_expected_pairs_constant() -> None:
    mod = _load()
    assert mod.EXPECTED_PAIRS == 47 * 14 == 658


def test_parse_tabix_rows_matches_chrom_pos_alleles_not_rsid() -> None:
    mod = _load()
    variant = mod.PanelVariant(rsid="rsTEST", chrom="1", pos=100, ref="A", alt="G")
    dataset = mod.DatasetRow(dataset_id="QTD000146", qts_id="QTS000999", file_type=".all.tsv.gz")
    header = mod.SYNTHETIC_HEADER
    # Wrong rsid in file, correct chrom/pos/alleles → kept
    body = (
        "ENSG1\t1\t100\tA\tG\tchr1_100_A_G\t10\t0.1\t1e-4\t0.25\t0.05\tSNP\t"
        "20\t100\t1.0\tENSG1\t5.0\trsOTHER\n"
        # Wrong alt → dropped
        "ENSG1\t1\t100\tA\tT\tchr1_100_A_T\t10\t0.1\t1e-4\t0.25\t0.05\tSNP\t"
        "20\t100\t1.0\tENSG1\t5.0\trsOTHER\n"
    )
    rows = mod.parse_tabix_rows(header=header, body=body, variant=variant, dataset=dataset)
    assert len(rows) == 1
    assert rows[0]["canonical_rsid"] == "rsTEST"
    assert rows[0]["eqtl_alt_allele"] == "G"
    assert rows[0]["alt_matches_panel_alt"] == "Y"
    assert rows[0]["beta"] == "0.25"


def test_account_pairs_requires_658() -> None:
    mod = _load()
    statuses = [
        {"status": "returned"},
        {"status": "no data"},
        {"status": "error"},
    ] * 200
    statuses = statuses[:657]
    try:
        mod.account_pairs(statuses)
        raised = False
    except AssertionError:
        raised = True
    assert raised


def test_dry_run_accounts_658_and_uses_cache(tmp_path: Path) -> None:
    mod = _load()
    out = tmp_path / "eqtl_associations_full_rebuilt.csv"
    diff = tmp_path / "diff.csv"
    provenance = tmp_path / "prov.json"
    pair_status = tmp_path / "pairs.tsv"
    doc = tmp_path / "EQTL_CATALOGUE_REBUILD.md"
    cache_dir = tmp_path / "cache"
    result = CliRunner().invoke(
        _app(mod),
        [
            "--dry-run",
            "--panel",
            str(_PANEL),
            "--cache-dir",
            str(cache_dir),
            "--out",
            str(out),
            "--diff-out",
            str(diff),
            "--provenance-out",
            str(provenance),
            "--pair-status-out",
            str(pair_status),
            "--doc-out",
            str(doc),
            "--skip-integrated",
        ],
    )
    assert result.exit_code == 0, result.output
    pairs = pl.read_csv(pair_status, separator="\t")
    assert pairs.height == 658
    assert set(pairs["status"].unique().to_list()).issubset({"returned", "no data", "error"})
    assert out.is_file()
    assert provenance.is_file()
    payload = json.loads(provenance.read_text(encoding="utf-8"))
    assert payload["n_pairs"] == 658
    assert sum(payload["pair_counts"].values()) == 658
    assert payload["mode"] == "dry-run-synthetic"
    assert payload["catalogue_release_or_version"] == "not provided by source"

    # Second run should log cache hits (reuse).
    result2 = CliRunner().invoke(
        _app(mod),
        [
            "--dry-run",
            "--panel",
            str(_PANEL),
            "--cache-dir",
            str(cache_dir),
            "--out",
            str(out),
            "--diff-out",
            str(diff),
            "--provenance-out",
            str(provenance),
            "--pair-status-out",
            str(pair_status),
            "--doc-out",
            str(doc),
            "--skip-integrated",
        ],
    )
    assert result2.exit_code == 0, result2.output
    payload2 = json.loads(provenance.read_text(encoding="utf-8"))
    assert payload2["n_cache_hits"] == 658
