"""Runtime configuration for the rogen_aging pipeline.

Loads layered YAML from ``config/default.yaml`` (and optional overrides) via
OmegaConf. Paths are resolved relative to the repository root.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any

from omegaconf import DictConfig, OmegaConf

_CONFIG: DictConfig | None = None
_REPO_ROOT: Path | None = None


def find_repo_root(start: Path | None = None) -> Path:
    """Locate the repository root (directory containing ``pyproject.toml``).

    Search order:
      1. Cached value from a previous successful lookup
      2. ``ROGEN_REPO_ROOT`` environment variable (when set)
      3. Walk parents of ``start`` (defaults to this module)
      4. Walk parents of the current working directory

    Args:
        start: Optional starting path. Defaults to this module's location.

    Returns:
        Absolute path to the repository root.

    Raises:
        FileNotFoundError: If no ``pyproject.toml`` is found.
    """
    global _REPO_ROOT
    if start is None and _REPO_ROOT is not None:
        return _REPO_ROOT

    env_root = os.environ.get("ROGEN_REPO_ROOT")
    if env_root:
        candidate = Path(env_root).expanduser().resolve()
        if (candidate / "pyproject.toml").is_file():
            _REPO_ROOT = candidate
            return candidate

    search_starts: list[Path] = []
    if start is not None:
        search_starts.append(Path(start).resolve())
    else:
        search_starts.append(Path(__file__).resolve())
    search_starts.append(Path.cwd().resolve())

    seen: set[Path] = set()
    for origin in search_starts:
        cur = origin.parent if origin.is_file() else origin
        for candidate in (cur, *cur.parents):
            if candidate in seen:
                continue
            seen.add(candidate)
            if (candidate / "pyproject.toml").is_file():
                _REPO_ROOT = candidate
                return candidate

    raise FileNotFoundError(
        "Could not locate repository root (pyproject.toml). "
        "Set ROGEN_REPO_ROOT or run from the repository."
    )


def default_config_dir(repo_root: Path | None = None) -> Path:
    """Return the repository ``config/`` directory."""
    return (repo_root or find_repo_root()) / "config"


def default_config_path(repo_root: Path | None = None) -> Path:
    """Return the path to ``config/default.yaml``."""
    return default_config_dir(repo_root) / "default.yaml"


def production_config_path(repo_root: Path | None = None) -> Path:
    """Return the path to ``config/production.yaml``."""
    return default_config_dir(repo_root) / "production.yaml"


def _as_dict(cfg: Any) -> dict[str, Any]:
    container = OmegaConf.to_container(cfg, resolve=True)
    if not isinstance(container, dict):
        raise TypeError(f"Expected mapping config node, got {type(container)!r}")
    return dict(container)


def resolve_repo_path(path: str | Path, *, repo_root: Path | None = None) -> Path:
    """Resolve a config path against the repository root.

    Absolute paths are returned unchanged. Relative paths are joined to
    ``repo_root``.
    """
    root = repo_root or find_repo_root()
    p = Path(path)
    return p if p.is_absolute() else (root / p)


def cfg_path(cfg: DictConfig, *keys: str, repo_root: Path | None = None) -> Path:
    """Read a nested string path from ``cfg`` and resolve it against the repo root.

    Args:
        cfg: Loaded OmegaConf config.
        *keys: Nested key path (e.g. ``\"paths\", \"models\", \"clock_elasticnet\"``).
        repo_root: Optional explicit repository root.

    Returns:
        Resolved filesystem path.
    """
    node: Any = cfg
    for key in keys:
        node = node[key]
    root = repo_root
    if root is None:
        repo_node = OmegaConf.select(cfg, "repo.root")
        root = Path(repo_node) if repo_node else find_repo_root()
    return resolve_repo_path(str(node), repo_root=root)


def load_config(
    config_path: str | Path | None = None,
    *,
    profile: str | None = None,
    overrides: dict[str, Any] | list[str] | None = None,
    repo_root: Path | None = None,
    set_active: bool = True,
) -> DictConfig:
    """Load pipeline configuration from YAML.

    Merge order:
      1. ``config/default.yaml``
      2. ``config/production.yaml`` when ``profile=\"production\"``
      3. Explicit ``config_path`` (if provided)
      4. ``overrides`` (dict deep-merge or OmegaConf dotlist)

    Args:
        config_path: Optional YAML file merged on top of defaults.
        profile: Optional named profile (``\"production\"``).
        overrides: Extra dict or OmegaConf dotlist overrides.
        repo_root: Repository root for relative paths; detected when omitted.
        set_active: When True (default), store result as the process-wide config.

    Returns:
        Merged ``DictConfig`` with ``repo.root`` set.
    """
    root = (repo_root or find_repo_root()).resolve()
    base_path = default_config_path(root)
    if not base_path.is_file():
        raise FileNotFoundError(f"Default config not found: {base_path}")

    cfg = OmegaConf.load(base_path)
    if profile == "production":
        prod = production_config_path(root)
        if not prod.is_file():
            raise FileNotFoundError(f"Production config not found: {prod}")
        cfg = OmegaConf.merge(cfg, OmegaConf.load(prod))
    if config_path is not None:
        path = Path(config_path)
        if not path.is_file():
            raise FileNotFoundError(f"Config file not found: {path}")
        cfg = OmegaConf.merge(cfg, OmegaConf.load(path))
    if overrides:
        if isinstance(overrides, dict):
            cfg = OmegaConf.merge(cfg, OmegaConf.create(overrides))
        else:
            cfg = OmegaConf.merge(cfg, OmegaConf.from_dotlist(list(overrides)))

    cfg = OmegaConf.create(OmegaConf.to_container(cfg, resolve=False))
    OmegaConf.set_struct(cfg, False)
    if not isinstance(cfg, DictConfig):
        raise TypeError(f"Expected DictConfig after merge, got {type(cfg)!r}")
    repo_cfg = cfg.get("repo")
    if not isinstance(repo_cfg, DictConfig):
        cfg["repo"] = {"root": str(root)}
    else:
        repo_cfg["root"] = str(root)

    if set_active:
        set_config(cfg)
    return cfg


def get_config() -> DictConfig:
    """Return the active config, loading ``config/default.yaml`` if needed."""
    if _CONFIG is None:
        return load_config()
    return _CONFIG


def set_config(cfg: DictConfig) -> DictConfig:
    """Set the process-wide active configuration and sync package defaults."""
    global _CONFIG, _REPO_ROOT
    _CONFIG = cfg
    repo_node = OmegaConf.select(cfg, "repo.root")
    if repo_node:
        _REPO_ROOT = Path(str(repo_node))
    _sync_package_defaults(cfg)
    return cfg


def reset_config() -> None:
    """Clear the active config (next ``get_config`` reloads defaults)."""
    global _CONFIG
    _CONFIG = None


def risk_weights(cfg: DictConfig | None = None) -> dict[str, float]:
    """Return integrative composite-risk channel weights."""
    c = cfg or get_config()
    return {str(k): float(v) for k, v in _as_dict(c.integrative.risk_weights).items()}


def vep_impact_scores(cfg: DictConfig | None = None) -> dict[str, float]:
    """Return VEP impact label → severity score mapping."""
    c = cfg or get_config()
    return {str(k): float(v) for k, v in _as_dict(c.integrative.vep_impact_scores).items()}


def target_tissues(cfg: DictConfig | None = None) -> tuple[str, ...]:
    """Return default GTEx target tissues for integrative mapping."""
    c = cfg or get_config()
    return tuple(str(t) for t in c.integrative.target_tissues)


def alphamissense_high_threshold(cfg: DictConfig | None = None) -> float:
    """Return the AlphaMissense high-pathogenicity threshold."""
    c = cfg or get_config()
    return float(c.integrative.alphamissense_high_threshold)


def _sync_package_defaults(cfg: DictConfig) -> None:
    """Refresh module-level DEFAULT_* aliases that consumers import directly.

    Only updates modules that are already imported to avoid circular imports
    during first load of :mod:`rogen_aging.config`.
    """
    import sys

    root = Path(str(cfg.repo.root))

    phenotype_mod = sys.modules.get("rogen_aging.integrative.phenotype_integrator")
    if phenotype_mod is not None:
        weights = risk_weights(cfg)
        scores = vep_impact_scores(cfg)
        existing_weights = getattr(phenotype_mod, "DEFAULT_WEIGHTS", None)
        if isinstance(existing_weights, dict):
            existing_weights.clear()
            existing_weights.update(weights)
        else:
            setattr(phenotype_mod, "DEFAULT_WEIGHTS", weights)
        existing_scores = getattr(phenotype_mod, "VEP_IMPACT_SCORES", None)
        if isinstance(existing_scores, dict):
            existing_scores.clear()
            existing_scores.update(scores)
        else:
            setattr(phenotype_mod, "VEP_IMPACT_SCORES", scores)
        setattr(
            phenotype_mod,
            "ALPHAMISSENSE_HIGH_THRESHOLD",
            alphamissense_high_threshold(cfg),
        )

    mapper_mod = sys.modules.get("rogen_aging.integrative.variant_tissue_mapper")
    if mapper_mod is not None:
        setattr(mapper_mod, "DEFAULT_TARGET_TISSUES", target_tissues(cfg))

    io_mod = sys.modules.get("rogen_aging.integrative.io")
    if io_mod is not None:
        setattr(io_mod, "REPO_ROOT", root)
        setattr(
            io_mod,
            "DEFAULT_JULY_XLSX",
            cfg_path(cfg, "paths", "integrative", "july_xlsx", repo_root=root),
        )
        setattr(
            io_mod,
            "DEFAULT_VARIANTS_PARQUET",
            cfg_path(cfg, "paths", "integrative", "variants_parquet", repo_root=root),
        )
        setattr(
            io_mod,
            "DEFAULT_EQTLS_PARQUET",
            cfg_path(cfg, "paths", "integrative", "eqtls_parquet", repo_root=root),
        )
        setattr(
            io_mod,
            "DEFAULT_EQTLS_CSV",
            cfg_path(cfg, "paths", "integrative", "eqtls_csv", repo_root=root),
        )
        setattr(
            io_mod,
            "DEFAULT_ALPHAGENOME",
            cfg_path(cfg, "paths", "integrative", "alphagenome", repo_root=root),
        )
        setattr(
            io_mod,
            "DEFAULT_OUTPUT_DIR",
            cfg_path(cfg, "paths", "integrative", "output_dir", repo_root=root),
        )


__all__ = [
    "alphamissense_high_threshold",
    "cfg_path",
    "default_config_dir",
    "default_config_path",
    "find_repo_root",
    "get_config",
    "load_config",
    "production_config_path",
    "reset_config",
    "resolve_repo_path",
    "risk_weights",
    "set_config",
    "target_tissues",
    "vep_impact_scores",
]
