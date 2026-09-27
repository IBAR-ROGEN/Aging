"""Snapshot guards for the allele × dataset eQTL coverage matrix."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import ModuleType

import polars as pl
import pytest

_REPO_ROOT = Path(__file__).resolve().parent.parent
_SCRIPT = _REPO_ROOT / "scripts" / "build_eqtl_coverage.py"


def _load_script() -> ModuleType:
    name = "build_eqtl_coverage"
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, _SCRIPT)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load {_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


cov = _load_script()


def test_coverage_errors_survive_optimized_bytecode() -> None:
    assert not issubclass(cov.CoverageCheckError, AssertionError)


def test_missing_associations_file_is_not_an_assertion(tmp_path: Path) -> None:
    missing = tmp_path / "eqtl_associations_full.csv"
    with pytest.raises(FileNotFoundError, match="Associations file not found"):
        cov.load_associations(missing)


def test_thirteen_datasets_are_refused(tmp_path: Path) -> None:
    n_rows = cov.EXPECTED_ASSOC_ROWS
    n_rsids = cov.EXPECTED_ASSOC_RSIDS
    n_datasets = 13
    frame = pl.DataFrame(
        {
            "canonical_rsid": [f"rs{i % n_rsids}" for i in range(n_rows)],
            "eqtl_alt_allele": ["A"] * n_rows,
            "tissue": ["cerebellum"] * n_rows,
            "dataset_id": [f"DS{i % n_datasets}" for i in range(n_rows)],
            "n_samples": [100] * n_rows,
            "fdr_bh": [0.2] * n_rows,
        }
    )
    path = tmp_path / "eqtl_associations_full.csv"
    frame.write_csv(path)
    with pytest.raises(cov.CoverageCheckError, match="expected 14"):
        cov.load_associations(path)


def test_shared_tissue_stays_two_dataset_columns() -> None:
    associations = pl.DataFrame(
        {
            "canonical_rsid": ["rs1"],
            "eqtl_alt_allele": ["A"],
            "tissue": ["cerebellum"],
            "dataset_id": ["QTD000161"],
            "n_samples": [175],
            "fdr_bh": [0.01],
        }
    )
    extra = pl.DataFrame(
        {
            "canonical_rsid": ["rs9"],
            "eqtl_alt_allele": ["G"],
            "tissue": ["cerebellum"],
            "dataset_id": ["QTD000166"],
            "n_samples": [209],
            "fdr_bh": [0.4],
        }
    )
    associations = pl.concat([associations, extra])
    annotation = pl.DataFrame(
        {
            "gene": ["GENE", "OTHER"],
            "canonical_rsid": ["rs1", "rs2"],
            "alt_allele": ["A", "T"],
            "eqtl_query_status": [cov.STATUS_TESTED, cov.STATUS_VARIANT_ABSENT],
        }
    )
    long_cells = cov.classify_cells(annotation, associations, enforce_catalogue=False)
    datasets = cov.dataset_columns(associations, enforce_catalogue=False)
    wide = cov.pivot_matrix(long_cells, datasets)
    dataset_cols = [column for column in wide.columns if column not in cov.INDEX_COLUMNS]
    assert wide.height == 2
    assert len(dataset_cols) == 2
    assert all("cerebellum" in column for column in dataset_cols)
    assert all(wide[column].null_count() == 0 for column in dataset_cols)

    tested = wide.filter(pl.col("canonical_rsid") == "rs1")
    present = [column for column in dataset_cols if "QTD000161" in column][0]
    missing = [column for column in dataset_cols if "QTD000166" in column][0]
    assert tested[present].item() == cov.CELL_ASSOC
    assert tested[missing].item() == cov.CELL_VARIANT_ABSENT

    absent = wide.filter(pl.col("canonical_rsid") == "rs2")
    assert absent[present].item() == cov.CELL_VARIANT_ABSENT
    assert absent[missing].item() == cov.CELL_VARIANT_ABSENT
