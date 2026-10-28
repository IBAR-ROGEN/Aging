#!/usr/bin/env python3
"""Export Hannum coefficients, retrain without chrX/chrY, and compare sex models.

Example:
    uv run python scripts/clock/clock_sex_check.py
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any

import joblib
import numpy as np
import polars as pl
import typer
from scipy.stats import pearsonr
from sklearn.metrics import mean_absolute_error

from rogen_aging.clock.benchmark import (
    age_acceleration,
    load_gse40279_heldout_cohort,
    load_gse87571_cohort,
    predict_rogen_hannum,
    sex_acceleration_model,
)
from rogen_aging.clock.hannum import (
    gse87571_cpg_ids_from_geo,
    load_probe_chromosomes,
    select_probes_from_parquet,
    write_coefficient_table,
    write_provenance,
    write_reproducible_coefficient_table,
    write_training_matrix,
)
from rogen_aging.clock.train import train_clock
from rogen_aging.config import find_repo_root, get_config

app = typer.Typer(add_completion=False, no_args_is_help=False)
logger = logging.getLogger(__name__)

REPO_ROOT = find_repo_root()
DEFAULT_MODEL = REPO_ROOT / "models" / "gse40279_hannum450k_elasticnet.joblib"
DEFAULT_PROVENANCE = REPO_ROOT / "models" / "gse40279_hannum450k_elasticnet.provenance.json"
DEFAULT_COEFF = REPO_ROOT / "models" / "gse40279_hannum450k_elasticnet_coefficients.csv"
DEFAULT_NOXY_MODEL = REPO_ROOT / "models" / "gse40279_hannum450k_elasticnet_noXY.joblib"
DEFAULT_NOXY_PROVENANCE = (
    REPO_ROOT / "models" / "gse40279_hannum450k_elasticnet_noXY.provenance.json"
)
DEFAULT_NOXY_COEFF = REPO_ROOT / "models" / "gse40279_hannum450k_elasticnet_noXY_coefficients.csv"
DEFAULT_NOXY_CPGS = REPO_ROOT / "models" / "gse40279_hannum450k_elasticnet_noXY_cpgs.csv"
DEFAULT_NOXY_METRICS = REPO_ROOT / "outputs" / "gse40279_hannum450k_noXY_train_metrics.json"
DEFAULT_NOXY_MATRIX = REPO_ROOT / "data" / "methylation" / "GSE40279_training_columns_noXY.parquet"
DEFAULT_ZHOU = REPO_ROOT / "data" / "dryrun" / "reference" / "HM450.hg38.manifest.tsv.gz"
DEFAULT_DOC = REPO_ROOT / "docs" / "CLOCK_SEX_CHECK.md"
DEFAULT_TABLE = REPO_ROOT / "analysis" / "clock_sex_check" / "metrics_table.csv"
EXCLUDE_SEX_CHROMS = frozenset({"chrX", "chrY"})


def chromosome_share_report(model_path: Path, zhou_manifest: Path) -> dict[str, Any]:
    """Count non-zero coefficients on chrX/chrY and their absolute-weight share."""
    model = joblib.load(model_path)
    names = [str(name) for name in model.feature_names_in_]
    coef = np.ravel(model.named_steps["elasticnet"].coef_).astype(np.float64)
    nonzero_idx = np.flatnonzero(coef)
    chrom_map = load_probe_chromosomes(zhou_manifest)
    n_x = 0
    n_y = 0
    abs_xy = 0.0
    abs_all = float(np.abs(coef[nonzero_idx]).sum()) if nonzero_idx.size else 0.0
    unmapped = 0
    for index in nonzero_idx:
        chrom = chrom_map.get(names[index])
        if chrom is None:
            unmapped += 1
            continue
        weight = abs(float(coef[index]))
        if chrom == "chrX":
            n_x += 1
            abs_xy += weight
        elif chrom == "chrY":
            n_y += 1
            abs_xy += weight
    share = (abs_xy / abs_all) if abs_all else 0.0
    return {
        "n_nonzero": int(nonzero_idx.size),
        "n_chrX": n_x,
        "n_chrY": n_y,
        "n_unmapped": unmapped,
        "abs_coef_chrXY": abs_xy,
        "abs_coef_all_nonzero": abs_all,
        "abs_coef_chrXY_share": share,
    }


def evaluate_model_on_cohorts(
    model_path: Path,
    *,
    label: str,
    gse40279_parquet: Path,
    gse40279_matrix: Path,
    geo_cache_dir: Path,
    gse87571_matrix: Path | None,
) -> list[dict[str, Any]]:
    """MAE, r, mean signed error, and sex model for GSE40279 held-out and GSE87571."""
    model = joblib.load(model_path)
    feats = [str(name) for name in model.feature_names_in_]
    gse40279 = load_gse40279_heldout_cohort(
        parquet_path=gse40279_parquet,
        series_matrix=gse40279_matrix,
        restrict_cpgs=feats,
    )
    gse87571 = load_gse87571_cohort(
        geo_cache_dir=geo_cache_dir,
        series_matrix=gse87571_matrix if gse87571_matrix and gse87571_matrix.is_file() else None,
        restrict_cpgs=feats,
    )
    rows: list[dict[str, Any]] = []
    for cohort in (gse40279, gse87571):
        predicted = predict_rogen_hannum(cohort.frame, model_path).to_numpy(dtype=float)
        age = cohort.frame["chronological_age"].to_numpy(dtype=float)
        female = cohort.frame["female"].to_numpy(dtype=float)
        signed = predicted - age
        mae = float(mean_absolute_error(age, predicted))
        r_value = float(pearsonr(predicted, age).statistic)
        accel = age_acceleration(predicted, age)
        sex = sex_acceleration_model(accel=accel, female=female, age=age)
        rows.append(
            {
                "model": label,
                "cohort": cohort.name,
                "n": int(len(predicted)),
                "mae": mae,
                "pearson_r": r_value,
                "mean_signed_error": float(np.mean(signed)),
                "sex_coef_female_minus_male": sex["sex_coef_female_minus_male"],
                "sex_ci_low": sex["sex_ci_low"],
                "sex_ci_high": sex["sex_ci_high"],
                "sex_p": sex["sex_p"],
            }
        )
    return rows


def write_doc(
    path: Path,
    *,
    chrom_report: dict[str, Any],
    rows: list[dict[str, Any]],
) -> None:
    """Write CLOCK_SEX_CHECK.md with the metrics table and one short paragraph."""
    by_key = {(row["model"], row["cohort"]): row for row in rows}
    main_87571 = by_key[("main", "GSE87571")]
    noxy_87571 = by_key[("noXY", "GSE87571")]
    main_held = by_key[("main", "GSE40279_heldout")]
    noxy_held = by_key[("noXY", "GSE40279_heldout")]

    def sex_sig(row: dict[str, Any], alpha: float = 0.05) -> bool:
        return float(row["sex_p"]) < alpha

    main_has = sex_sig(main_87571) or sex_sig(main_held)
    noxy_has = sex_sig(noxy_87571) or sex_sig(noxy_held)
    if noxy_has and not main_has:
        verdict = (
            "The noXY model shows a sex difference (accel ~ female + age, p < 0.05 on at "
            "least one cohort) and the main model does not."
        )
    elif not noxy_has and not main_has:
        verdict = (
            "Neither the main model nor the noXY model shows a sex difference "
            "(accel ~ female + age, p < 0.05) on GSE40279 held-out or GSE87571."
        )
    elif main_has and not noxy_has:
        verdict = (
            "The main model shows a sex difference (accel ~ female + age, p < 0.05 on at "
            "least one cohort) and the noXY model does not."
        )
    else:
        verdict = (
            "Both the main model and the noXY model show a sex difference "
            "(accel ~ female + age, p < 0.05) on at least one cohort."
        )

    lines = [
        "# Clock sex check (chrX / chrY sensitivity)",
        "",
        "Non-zero coefficients in `models/gse40279_hannum450k_elasticnet.joblib` "
        "mapped with the Zhou HM450 hg38 manifest (`data/dryrun/reference/HM450.hg38.manifest.tsv.gz`).",
        "",
        "| Item | Value |",
        "|------|-------|",
        f"| Non-zero coefficients | {chrom_report['n_nonzero']} |",
        f"| Non-zero on chrX | {chrom_report['n_chrX']} |",
        f"| Non-zero on chrY | {chrom_report['n_chrY']} |",
        f"| Unmapped non-zero | {chrom_report['n_unmapped']} |",
        f"| Sum \\|coef\\| on chrX+chrY | {chrom_report['abs_coef_chrXY']:.6g} |",
        f"| Sum \\|coef\\| all non-zero | {chrom_report['abs_coef_all_nonzero']:.6g} |",
        f"| chrX+chrY share of sum \\|coef\\| | {chrom_report['abs_coef_chrXY_share']:.6g} |",
        "",
        "The noXY model uses the same GSE40279 data, split (`test_size=0.2`, `random_state=42`), "
        "missingness and age-correlation filters, probe cap (8000), and ElasticNetCV settings as "
        "the main provenance, except chrX and chrY probes are removed after the missingness "
        "filter and before the age-correlation filter. Artifact: "
        "`models/gse40279_hannum450k_elasticnet_noXY.joblib`.",
        "",
        "## Accuracy and sex model",
        "",
        "| Model | Cohort | n | MAE | r | Mean signed error | Sex coef (F−M) | Sex p |",
        "|-------|--------|---|-----|---|-------------------|----------------|-------|",
    ]
    for row in rows:
        lines.append(
            f"| {row['model']} | {row['cohort']} | {row['n']} | {row['mae']:.6g} | "
            f"{row['pearson_r']:.6g} | {row['mean_signed_error']:.6g} | "
            f"{row['sex_coef_female_minus_male']:.6g} | {row['sex_p']:.6g} |"
        )
    lines.extend(["", verdict, ""])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


@app.command()
def main(
    skip_retrain: bool = typer.Option(
        False,
        "--skip-retrain",
        help="Reuse an existing noXY joblib instead of retraining.",
    ),
) -> None:
    """Export coefficients, optional noXY retrain, evaluate, and write the doc."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    if not DEFAULT_MODEL.is_file():
        raise typer.Exit(f"Missing main model: {DEFAULT_MODEL}")
    if not DEFAULT_ZHOU.is_file():
        raise typer.Exit(f"Missing Zhou manifest: {DEFAULT_ZHOU}")

    typer.echo(f"Writing reproducible coefficients → {DEFAULT_COEFF}")
    nonzero = write_reproducible_coefficient_table(DEFAULT_MODEL, DEFAULT_COEFF)
    typer.echo(f"non-zero coefficients: {nonzero}")

    chrom_report = chromosome_share_report(DEFAULT_MODEL, DEFAULT_ZHOU)
    typer.echo(json.dumps(chrom_report, indent=2))

    if not skip_retrain or not DEFAULT_NOXY_MODEL.is_file():
        typer.echo(
            "Retraining noXY model (exclude chrX/chrY before age filter; "
            "GSE87571 probe set from GEO matrices, not the 8000-column parquet)…"
        )
        chrom_map = load_probe_chromosomes(DEFAULT_ZHOU)
        cfg = get_config().clock
        input_parquet = REPO_ROOT / "data" / "methylation" / "GSE40279_processed.parquet"
        geo_ids = set(gse87571_cpg_ids_from_geo(REPO_ROOT / "data" / "geo"))
        report = select_probes_from_parquet(
            input_parquet,
            geo_ids,
            test_size=float(cfg.test_size),
            random_state=int(cfg.random_state),
            probe_chromosomes=chrom_map,
            exclude_chromosomes=EXCLUDE_SEX_CHROMS,
        )
        typer.echo(
            f"noXY filter: trained={report.n_probes_trained} "
            f"after_age={report.n_probes_after_age_correlation} cv={report.cv}"
        )
        write_training_matrix(input_parquet, DEFAULT_NOXY_MATRIX, report.columns)
        metrics = train_clock(
            DEFAULT_NOXY_MATRIX,
            DEFAULT_NOXY_MODEL,
            DEFAULT_NOXY_METRICS,
            test_size=float(cfg.test_size),
            random_state=int(cfg.random_state),
            cv=report.cv,
        )
        write_coefficient_table(DEFAULT_NOXY_MODEL, DEFAULT_NOXY_CPGS)
        write_reproducible_coefficient_table(DEFAULT_NOXY_MODEL, DEFAULT_NOXY_COEFF)
        write_provenance(
            DEFAULT_NOXY_PROVENANCE,
            input_parquet=input_parquet,
            report=report,
            metrics=metrics,
            coefficient_csv=DEFAULT_NOXY_CPGS,
            model_path=DEFAULT_NOXY_MODEL,
            repo_root=REPO_ROOT,
            reproducible_coefficient_csv=DEFAULT_NOXY_COEFF,
            exclude_chromosomes=EXCLUDE_SEX_CHROMS,
        )
    else:
        typer.echo(f"Reusing existing noXY model: {DEFAULT_NOXY_MODEL}")

    gse40279_matrix = REPO_ROOT / "data" / "geo" / "GSE40279_series_matrix.txt.gz"
    gse87571_matrix = REPO_ROOT / "data" / "geo" / "GSE87571_series_matrix.txt.gz"
    geo_cache = REPO_ROOT / "data" / "geo"
    parquet = REPO_ROOT / "data" / "methylation" / "GSE40279_processed.parquet"

    rows: list[dict[str, Any]] = []
    rows.extend(
        evaluate_model_on_cohorts(
            DEFAULT_MODEL,
            label="main",
            gse40279_parquet=parquet,
            gse40279_matrix=gse40279_matrix,
            geo_cache_dir=geo_cache,
            gse87571_matrix=gse87571_matrix,
        )
    )
    rows.extend(
        evaluate_model_on_cohorts(
            DEFAULT_NOXY_MODEL,
            label="noXY",
            gse40279_parquet=parquet,
            gse40279_matrix=gse40279_matrix,
            geo_cache_dir=geo_cache,
            gse87571_matrix=gse87571_matrix,
        )
    )
    DEFAULT_TABLE.parent.mkdir(parents=True, exist_ok=True)
    pl.DataFrame(rows).write_csv(DEFAULT_TABLE)
    write_doc(DEFAULT_DOC, chrom_report=chrom_report, rows=rows)
    typer.echo(f"Wrote {DEFAULT_TABLE}")
    typer.echo(f"Wrote {DEFAULT_DOC}")


if __name__ == "__main__":
    app()
