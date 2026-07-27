"""Offline development fixtures for ROGEN annotation / integrative / clock pipelines.

Generates minimal synthetic (or locally derived) inputs so pipeline CLIs can be
exercised without live API access or large external datasets. Safe to re-run;
existing files are overwritten by default.
"""

from __future__ import annotations

import json
import pickle
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import polars as pl
from sklearn.linear_model import ElasticNet, ElasticNetCV

REPO_ROOT = Path(__file__).resolve().parent.parent.parent

JULY_VARIANTS = REPO_ROOT / "data" / "processed" / "prioritized_variants.csv"
JULY_ALPHAGENOME = REPO_ROOT / "data" / "scores" / "alphagenome_raw.parquet"
JULY_ALPHAMISSENSE = REPO_ROOT / "data" / "scores" / "alphamissense_raw.parquet"
JULY_LOCAL_VEP = REPO_ROOT / "data" / "processed" / "vep_local.jsonl"
JULY_CACHE_DIR = REPO_ROOT / "data" / "cache" / "july_annotation"

INTEGRATIVE_DIR = REPO_ROOT / "analysis" / "integrative" / "fixtures"
INTEGRATIVE_VARIANTS = INTEGRATIVE_DIR / "annotated_variants.parquet"
INTEGRATIVE_EQTLS = INTEGRATIVE_DIR / "eqtls.parquet"
INTEGRATIVE_PROBES = INTEGRATIVE_DIR / "probe_annotation.parquet"
INTEGRATIVE_SAMPLES = INTEGRATIVE_DIR / "sample_genotypes.parquet"

CLOCK_MODEL_PKL = REPO_ROOT / "models" / "ro_clock_elasticnet_gse40279.pkl"
CLOCK_MODEL_JOBLIB = REPO_ROOT / "models" / "methylation_clock_v1.joblib"
CLOCK_METH = REPO_ROOT / "data" / "methylation" / "GSE87571_processed.parquet"
CLOCK_META = REPO_ROOT / "data" / "methylation" / "GSE87571_meta.csv"

VEP_SOURCE = REPO_ROOT / "analysis" / "vep_annotation" / "la_snp_vep_annotations.csv"
AG_SOURCE = REPO_ROOT / "analysis" / "alphagenome" / "alphagenome_impact_analysis.csv"
GTEX_SOURCE = REPO_ROOT / "analysis" / "gtex_annotation" / "la_snp_gtex_eqtls.csv"

_DEMO_CPGS = ("cg00000001", "cg00000002", "cg00000003")


@dataclass(frozen=True)
class FixturePaths:
    """Resolved paths written by :func:`write_all_pipeline_fixtures`."""

    july_variants: Path
    july_alphagenome: Path
    july_alphamissense: Path
    july_local_vep: Path
    july_cache_dir: Path
    integrative_variants: Path
    integrative_eqtls: Path
    integrative_probes: Path
    integrative_samples: Path
    clock_model: Path
    clock_methylation: Path
    clock_meta: Path


def _ensure_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def _synthetic_variants(n: int = 8) -> pl.DataFrame:
    """Return a tiny prioritized-variant table for offline July/integrative demos."""
    rows = [
        {
            "chrom": "2",
            "pos": 206129125,
            "ref": "A",
            "alt": "G",
            "rsid": "rs6435324",
            "gene_symbol": "NDUFS1",
        },
        {
            "chrom": "2",
            "pos": 206141952,
            "ref": "T",
            "alt": "A",
            "rsid": "rs1801318",
            "gene_symbol": "NDUFS1",
        },
        {
            "chrom": "17",
            "pos": 33231561,
            "ref": "G",
            "alt": "A",
            "rsid": "rs12952455",
            "gene_symbol": "ASIC2",
        },
        {
            "chrom": "19",
            "pos": 44908684,
            "ref": "T",
            "alt": "C",
            "rsid": "rs429358",
            "gene_symbol": "APOE",
        },
        {
            "chrom": "6",
            "pos": 31575254,
            "ref": "G",
            "alt": "A",
            "rsid": "rs1046089",
            "gene_symbol": "HSPA1L",
        },
        {
            "chrom": "16",
            "pos": 56995236,
            "ref": "C",
            "alt": "T",
            "rsid": "rs5882",
            "gene_symbol": "CETP",
        },
        {
            "chrom": "11",
            "pos": 108093211,
            "ref": "G",
            "alt": "A",
            "rsid": "rs1800057",
            "gene_symbol": "ATM",
        },
        {
            "chrom": "1",
            "pos": 155204800,
            "ref": "C",
            "alt": "T",
            "rsid": "rs2230288",
            "gene_symbol": "GBA",
        },
    ]
    return pl.DataFrame(rows[:n])


def _variants_from_vep_csv(path: Path, limit: int = 12) -> pl.DataFrame | None:
    """Build prioritized variants from a local VEP annotation CSV when present."""
    if not path.is_file():
        return None
    frame = pl.read_csv(path, infer_schema_length=5_000)
    rename = {
        "rsID": "rsid",
        "chromosome": "chrom",
        "position_GRCh38": "pos",
        "ref_allele": "ref",
        "alt_allele": "alt",
        "gene_symbols": "gene_symbol",
    }
    present = {src: dst for src, dst in rename.items() if src in frame.columns}
    if not present:
        return None
    frame = frame.rename(present)
    needed = {"chrom", "pos", "ref", "alt", "rsid", "gene_symbol"}
    if not needed.issubset(frame.columns):
        return None
    frame = frame.select(list(needed)).with_columns(
        pl.col("alt").cast(pl.Utf8).str.split(",").list.first().alias("alt"),
        pl.col("gene_symbol")
        .cast(pl.Utf8)
        .str.split(";")
        .list.first()
        .alias("gene_symbol"),
    )
    return frame.drop_nulls(["chrom", "pos", "ref", "alt", "rsid"]).head(limit)


def write_july_fixtures(*, repo_root: Path = REPO_ROOT, limit: int = 12) -> dict[str, Path]:
    """Write July annotation inputs (variants, scores, local VEP, warm GTEx cache).

    Args:
        repo_root: Repository root used to resolve default paths.
        limit: Maximum number of demo variants.

    Returns:
        Mapping of logical name → written path.
    """
    variants_path = repo_root / "data" / "processed" / "prioritized_variants.csv"
    ag_path = repo_root / "data" / "scores" / "alphagenome_raw.parquet"
    am_path = repo_root / "data" / "scores" / "alphamissense_raw.parquet"
    vep_path = repo_root / "data" / "processed" / "vep_local.jsonl"
    cache_dir = repo_root / "data" / "cache" / "july_annotation"

    variants = _variants_from_vep_csv(repo_root / "analysis" / "vep_annotation" / "la_snp_vep_annotations.csv", limit)
    if variants is None:
        variants = _synthetic_variants(min(limit, 8))

    _ensure_parent(variants_path)
    variants.write_csv(variants_path)

    ag_src = repo_root / "analysis" / "alphagenome" / "alphagenome_impact_analysis.csv"
    if ag_src.is_file():
        ag = pl.read_csv(ag_src)
        if "snp" in ag.columns:
            keep = variants["rsid"].to_list()
            ag = ag.filter(pl.col("snp").is_in(keep))
        if ag.height == 0:
            ag = pl.DataFrame(
                {
                    "gene": variants["gene_symbol"],
                    "snp": variants["rsid"],
                    "ref_score": [0.1] * variants.height,
                    "alt_score": [0.12] * variants.height,
                    "diff": [0.02] * variants.height,
                    "perc_change": [20.0] * variants.height,
                    "abs_perc_change": [20.0] * variants.height,
                }
            )
    else:
        ag = pl.DataFrame(
            {
                "gene": variants["gene_symbol"],
                "snp": variants["rsid"],
                "ref_score": [0.1] * variants.height,
                "alt_score": [0.12] * variants.height,
                "diff": [0.02] * variants.height,
                "perc_change": [20.0] * variants.height,
                "abs_perc_change": [20.0] * variants.height,
            }
        )
    _ensure_parent(ag_path)
    ag.write_parquet(ag_path)

    rng = np.random.default_rng(42)
    am_scores = rng.uniform(0.05, 0.95, size=variants.height)
    am = pl.DataFrame(
        {
            "snp": variants["rsid"],
            "am_score": am_scores,
            "am_class": [
                "likely_pathogenic" if s > 0.5 else "likely_benign" for s in am_scores
            ],
        }
    )
    _ensure_parent(am_path)
    am.write_parquet(am_path)

    _write_local_vep_jsonl(variants, vep_path, repo_root=repo_root)
    _seed_gtex_cache(variants, cache_dir, repo_root=repo_root)

    return {
        "variants": variants_path,
        "alphagenome": ag_path,
        "alphamissense": am_path,
        "local_vep": vep_path,
        "cache_dir": cache_dir,
    }


def _write_local_vep_jsonl(
    variants: pl.DataFrame,
    path: Path,
    *,
    repo_root: Path,
) -> None:
    """Write synthetic Ensembl-like VEP JSONL keyed by rsID / variant_key."""
    vep_src = repo_root / "analysis" / "vep_annotation" / "la_snp_vep_annotations.csv"
    by_rsid: dict[str, dict[str, Any]] = {}
    if vep_src.is_file():
        src = pl.read_csv(vep_src)
        for row in src.to_dicts():
            rsid = str(row.get("rsID") or row.get("rsid") or "")
            if rsid:
                by_rsid[rsid] = row

    impact_for = {
        "missense_variant": "MODERATE",
        "stop_gained": "HIGH",
        "frameshift_variant": "HIGH",
        "intron_variant": "MODIFIER",
        "synonymous_variant": "LOW",
    }

    _ensure_parent(path)
    with path.open("w", encoding="utf-8") as handle:
        for row in variants.to_dicts():
            rsid = str(row["rsid"])
            chrom = str(row["chrom"]).removeprefix("chr")
            pos = int(row["pos"])
            ref = str(row["ref"]).upper()
            alt = str(row["alt"]).upper().split(",")[0]
            variant_key = f"{chrom}:{pos}:{ref}:{alt}"
            src = by_rsid.get(rsid, {})
            consequence = str(src.get("most_severe_consequence") or "intron_variant")
            impact = impact_for.get(consequence, "MODIFIER")
            gene = str(row["gene_symbol"])
            sift = src.get("SIFT") or src.get("sift") or ""
            poly = src.get("PolyPhen") or src.get("polyphen") or ""
            payload = [
                {
                    "id": rsid,
                    "most_severe_consequence": consequence,
                    "transcript_consequences": [
                        {
                            "impact": impact,
                            "consequence_terms": [consequence],
                            "transcript_id": "ENST_DEMO",
                            "gene_symbol": gene,
                            "canonical": 1,
                            "sift_prediction": sift or None,
                            "polyphen_prediction": poly or None,
                            "hgvsg": f"{chrom}:g.{pos}{ref}>{alt}",
                            "hgvsc": "c.1A>G",
                            "hgvsp": "p.Met1Val" if impact in {"HIGH", "MODERATE"} else None,
                        }
                    ],
                }
            ]
            record = {"rsid": rsid, "variant_key": variant_key, "payload": payload}
            handle.write(json.dumps(record) + "\n")


def _cache_path(cache_dir: Path, kind: str, key: str) -> Path:
    """Mirror july pipeline cache naming for warm demo caches."""
    import hashlib
    import re

    safe = re.sub(r"[^\w.\-]+", "_", key)
    if len(safe) > 120:
        digest = hashlib.sha256(key.encode("utf-8")).hexdigest()[:24]
        safe = f"{safe[:80]}__{digest}"
    return cache_dir / f"{kind}_{safe}.json"


def _seed_gtex_cache(
    variants: pl.DataFrame,
    cache_dir: Path,
    *,
    repo_root: Path,
) -> None:
    """Warm GTEx variant + eQTL caches so ``--cache-only`` demos stay offline."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    gtex_src = repo_root / "analysis" / "gtex_annotation" / "la_snp_gtex_eqtls.csv"
    eqtls_by_rsid: dict[str, list[dict[str, Any]]] = {}
    if gtex_src.is_file():
        src = pl.read_csv(gtex_src)
        rename = {"rsID": "rsid"}
        if "rsID" in src.columns:
            src = src.rename(rename)
        for row in src.to_dicts():
            rsid = str(row.get("rsid") or "")
            if not rsid:
                continue
            eqtls_by_rsid.setdefault(rsid, []).append(row)

    # Tissue list must match scripts/ukb/run_july_annotation_pipeline.TARGET_TISSUES
    # for aggregate cache key compatibility — import lazily to avoid cycles.
    target_tissues = (
        "Brain_Amygdala",
        "Brain_Anterior_cingulate_cortex_BA24",
        "Brain_Caudate_basal_ganglia",
        "Brain_Cerebellar_Hemisphere",
        "Brain_Cerebellum",
        "Brain_Cortex",
        "Brain_Frontal_Cortex_BA9",
        "Brain_Hippocampus",
        "Brain_Hypothalamus",
        "Brain_Nucleus_accumbens_basal_ganglia",
        "Brain_Putamen_basal_ganglia",
        "Brain_Spinal_cord_cervical_c-1",
        "Brain_Substantia_nigra",
        "Whole_Blood",
    )
    tissue_key = "_".join(target_tissues)

    for row in variants.to_dicts():
        rsid = str(row["rsid"])
        chrom = str(row["chrom"]).removeprefix("chr")
        pos = int(row["pos"])
        ref = str(row["ref"]).upper()
        alt = str(row["alt"]).upper().split(",")[0]
        gtex_variant_id = f"chr{chrom}_{pos}_{ref}_{alt}_b38"
        variant_payload = {
            "data": [
                {
                    "variantId": gtex_variant_id,
                    "snpId": rsid,
                    "chromosome": f"chr{chrom}",
                    "pos": pos,
                    "ref": ref,
                    "alt": alt,
                }
            ]
        }
        with _cache_path(cache_dir, "gtex_variant", rsid.lower()).open("w", encoding="utf-8") as handle:
            json.dump(variant_payload, handle)

        hits: list[dict[str, Any]] = []
        for eq in eqtls_by_rsid.get(rsid, []):
            tissue = str(eq.get("tissue") or "Whole_Blood")
            if tissue not in target_tissues:
                continue
            hits.append(
                {
                    "tissueSiteDetailId": tissue,
                    "nes": float(eq.get("nes") or 0.1),
                    "pValue": float(eq.get("p_value") or 1e-5),
                    "geneSymbol": str(eq.get("gene_symbol") or row["gene_symbol"]),
                }
            )
        if not hits:
            hits = [
                {
                    "tissueSiteDetailId": "Whole_Blood",
                    "nes": 0.15,
                    "pValue": 1e-6,
                    "geneSymbol": str(row["gene_symbol"]),
                },
                {
                    "tissueSiteDetailId": "Brain_Cortex",
                    "nes": -0.08,
                    "pValue": 2e-4,
                    "geneSymbol": str(row["gene_symbol"]),
                },
            ]
        aggregate = {"variantId": gtex_variant_id, "hits": hits}
        with _cache_path(cache_dir, "gtex_eqtl", f"{gtex_variant_id}__{tissue_key}").open(
            "w", encoding="utf-8"
        ) as handle:
            json.dump(aggregate, handle)


def write_integrative_fixtures(*, repo_root: Path = REPO_ROOT) -> dict[str, Path]:
    """Write annotated-variant / eQTL / probe / sample fixtures for integrative CLI.

    Args:
        repo_root: Repository root used to resolve default paths.

    Returns:
        Mapping of logical name → written path.
    """
    out_dir = repo_root / "analysis" / "integrative" / "fixtures"
    out_dir.mkdir(parents=True, exist_ok=True)

    variants = pl.DataFrame(
        {
            "chrom": ["chr17", "2", "19", "16"],
            "pos": [100, 200, 300, 400],
            "ref": ["A", "G", "C", "T"],
            "alt": ["G", "T", "T", "C"],
            "rsid": ["rs1", "rs2", "rs3", "rs4"],
            "gene_symbol": ["FOXO3", "APOE", "CETP", "HSPA1A"],
            "vep_impact": ["HIGH", "MODERATE", "LOW", "MODIFIER"],
            "alphagenome_abs_perc_change": [25.0, 10.0, 0.0, 5.0],
            "alphamissense_score": [0.8, 0.2, None, 0.55],
            "age_acceleration": [5.0, -2.0, 0.0, 1.5],
        }
    )
    eqtls = pl.DataFrame(
        {
            "rsid": ["rs1", "rs1", "rs2", "rs2", "rs4"],
            "tissue": [
                "Brain_Cortex",
                "Whole_Blood",
                "Brain_Hippocampus",
                "Lung",
                "Whole_Blood",
            ],
            "nes": [0.5, -0.2, 0.1, 0.9, 0.3],
            "p_value": [1e-8, 1e-4, 1e-3, 1e-10, 1e-5],
            "gene_symbol": ["FOXO3", "FOXO3", "APOE", "APOE", "HSPA1A"],
            "gtex_variant_id": ["v1", "v1", "v2", "v2", "v4"],
        }
    )
    probes = pl.DataFrame(
        {
            "IlmnID": ["cg0001", "cg0002", "cg0003"],
            "UCSC_RefGene_Name": ["FOXO3", "APOE;TOMM40", "HSPA1A"],
        }
    )
    samples = pl.DataFrame(
        {
            "sample_id": ["S1", "S1", "S2", "S2", "S3", "S3"],
            "rsid": ["rs1", "rs2", "rs1", "rs4", "rs2", "rs3"],
            "alt_dosage": [1.0, 0.0, 2.0, 1.0, 1.0, 0.0],
        }
    )

    variants_path = out_dir / "annotated_variants.parquet"
    eqtls_path = out_dir / "eqtls.parquet"
    probes_path = out_dir / "probe_annotation.parquet"
    samples_path = out_dir / "sample_genotypes.parquet"
    variants.write_parquet(variants_path)
    eqtls.write_parquet(eqtls_path)
    probes.write_parquet(probes_path)
    samples.write_parquet(samples_path)
    return {
        "variants": variants_path,
        "eqtls": eqtls_path,
        "probes": probes_path,
        "samples": samples_path,
    }


def _fit_demo_elasticnet() -> ElasticNet:
    """Fit a tiny bare ElasticNet with ``feature_names_in_`` for clock demos."""
    rng = np.random.default_rng(0)
    x = pl.DataFrame(
        {
            _DEMO_CPGS[0]: rng.uniform(0.1, 0.9, size=40),
            _DEMO_CPGS[1]: rng.uniform(0.1, 0.9, size=40),
            _DEMO_CPGS[2]: rng.uniform(0.1, 0.9, size=40),
        }
    ).to_pandas()
    y = (
        25.0
        + 30.0 * x[_DEMO_CPGS[0]]
        + 10.0 * x[_DEMO_CPGS[1]]
        + rng.normal(0.0, 0.5, size=40)
    )
    model = ElasticNet(alpha=0.01, l1_ratio=0.5, max_iter=10_000, random_state=0)
    model.fit(x, y)
    return model


def _load_clock_artifact(path: Path) -> Any | None:
    """Load a joblib/pickle clock if it looks usable by ``load_elasticnet_clock``."""
    if not path.is_file():
        return None
    suffix = path.suffix.lower()
    if suffix == ".joblib":
        import joblib

        obj = joblib.load(path)
    else:
        with path.open("rb") as handle:
            obj = pickle.load(handle)

    if type(obj) is ElasticNet:
        return obj
    if hasattr(obj, "named_steps"):
        final = list(obj.named_steps.values())[-1]
        if isinstance(final, (ElasticNet, ElasticNetCV)):
            return obj
    return None


def write_clock_fixtures(
    *,
    repo_root: Path = REPO_ROOT,
    force_synthetic: bool = False,
) -> dict[str, Path]:
    """Ensure clock model (+ optional synthetic methylation) exist for eval demos.

    Prefers re-serializing ``methylation_clock_v1.joblib`` (Pipeline or bare
    ElasticNet) to ``ro_clock_elasticnet_gse40279.pkl``. When methylation inputs
    are missing, writes a tiny synthetic cohort under ``data/methylation/``.

    Args:
        repo_root: Repository root used to resolve default paths.
        force_synthetic: Always write the synthetic cohort even if GSE files exist.

    Returns:
        Mapping of logical name → written path.
    """
    model_pkl = repo_root / "models" / "ro_clock_elasticnet_gse40279.pkl"
    model_joblib = repo_root / "models" / "methylation_clock_v1.joblib"
    meth_path = repo_root / "data" / "methylation" / "GSE87571_processed.parquet"
    meta_path = repo_root / "data" / "methylation" / "GSE87571_meta.csv"

    model = _load_clock_artifact(model_joblib)
    if model is None:
        model = _load_clock_artifact(model_pkl)
    if model is None:
        model = _fit_demo_elasticnet()

    _ensure_parent(model_pkl)
    with model_pkl.open("wb") as handle:
        pickle.dump(model, handle)

    if force_synthetic or not meth_path.is_file() or not meta_path.is_file():
        rng = np.random.default_rng(7)
        n = 24
        ids = [f"DEMO{i:04d}" for i in range(n)]
        feature_names = list(getattr(model, "feature_names_in_", ()))
        if not feature_names and hasattr(model, "named_steps"):
            feature_names = list(getattr(model, "feature_names_in_", ()))
        if not feature_names:
            final = (
                list(model.named_steps.values())[-1]
                if hasattr(model, "named_steps")
                else model
            )
            n_feat = int(getattr(final, "n_features_in_", len(_DEMO_CPGS)))
            feature_names = [f"cg{i:08d}" for i in range(1, min(n_feat, 50) + 1)]
        feature_names = [str(c) for c in feature_names[:50]]
        meth = pl.DataFrame(
            {
                "sample_id": ids,
                **{
                    name: rng.uniform(0.05, 0.95, size=n).tolist()
                    for name in feature_names
                },
            }
        )
        ages = np.concatenate(
            [
                rng.uniform(18.0, 29.0, size=8),
                rng.uniform(30.0, 60.0, size=8),
                rng.uniform(61.0, 85.0, size=8),
            ]
        )
        meta = pl.DataFrame({"sample_id": ids, "chronological_age": ages.tolist()})
        _ensure_parent(meth_path)
        meth.write_parquet(meth_path)
        _ensure_parent(meta_path)
        meta.write_csv(meta_path)

    return {
        "model": model_pkl,
        "methylation": meth_path,
        "meta": meta_path,
    }


def write_all_pipeline_fixtures(
    *,
    repo_root: Path = REPO_ROOT,
    july_limit: int = 12,
    force_synthetic_clock: bool = False,
) -> FixturePaths:
    """Write July, integrative, and clock development fixtures.

    Args:
        repo_root: Repository root.
        july_limit: Max variants for the July demo table.
        force_synthetic_clock: Force synthetic methylation cohort for clock eval.

    Returns:
        :class:`FixturePaths` pointing at all generated artefacts.
    """
    july = write_july_fixtures(repo_root=repo_root, limit=july_limit)
    integrative = write_integrative_fixtures(repo_root=repo_root)
    clock = write_clock_fixtures(
        repo_root=repo_root, force_synthetic=force_synthetic_clock
    )
    return FixturePaths(
        july_variants=july["variants"],
        july_alphagenome=july["alphagenome"],
        july_alphamissense=july["alphamissense"],
        july_local_vep=july["local_vep"],
        july_cache_dir=july["cache_dir"],
        integrative_variants=integrative["variants"],
        integrative_eqtls=integrative["eqtls"],
        integrative_probes=integrative["probes"],
        integrative_samples=integrative["samples"],
        clock_model=clock["model"],
        clock_methylation=clock["methylation"],
        clock_meta=clock["meta"],
    )
