#!/usr/bin/env python3
"""Read-only audit of clock serializations vs claimed GSE40279 provenance.

Does not rename, move, retrain, regenerate, or otherwise modify any model file.
Does not download data. Does not reconstruct a missing clock from metrics.

Example:
    uv run python scripts/audit_clock_artifacts.py
"""

from __future__ import annotations

import json
import pickle
import re
import subprocess
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
import typer

from rogen_aging.config import find_repo_root

REPO_ROOT = find_repo_root()
DEFAULT_REPORT = REPO_ROOT / "docs" / "CLOCK_ARTIFACT_AUDIT.md"
AUDIT_REPORT_NAME = "CLOCK_ARTIFACT_AUDIT.md"

SKIP_DIR_NAMES = frozenset(
    {
        ".git",
        ".venv",
        "venv",
        "node_modules",
        "__pycache__",
        ".ruff_cache",
        ".mypy_cache",
        ".pytest_cache",
        ".mplconfig",
        "tests",
    }
)
ARTIFACT_SUFFIXES = {".pkl", ".joblib", ".sav"}
FIXTURE_NAME_RE = re.compile(r"^cg_test_\d+$")
ILLUMINA_NAME_RE = re.compile(r"^cg\d{8}$")
SKLEARN_VERSION_RE = re.compile(rb"_sklearn_version.{1,20}?([0-9]+\.[0-9]+\.[0-9]+)")
MAE_JSON_RE = re.compile(
    r'"(?P<key>test_mae|mae|mae_overall|mae_overall_years)"\s*:\s*(?P<value>-?\d+(?:\.\d+)?)'
)
MAE_PROSE_RE = re.compile(
    r"(?P<label>MAE\s*(?:~|=|:)|test MAE|baseline MAE|MEAN ABSOLUTE ERROR \(MAE\))"
    r".{0,80}?(?P<value>\d+\.\d+)",
    re.IGNORECASE,
)
TEXT_SUFFIXES = {
    ".md",
    ".json",
    ".py",
    ".ipynb",
    ".txt",
    ".yml",
    ".yaml",
    ".csv",
    ".rst",
    ".toml",
}


def discover_artifacts(root: Path) -> list[Path]:
    """Return every ``.pkl`` / ``.joblib`` / ``.sav`` under ``root`` (not venv/git)."""
    found: list[Path] = []
    for path in root.rglob("*"):
        if not path.is_file():
            continue
        if path.suffix.lower() not in ARTIFACT_SUFFIXES:
            continue
        if any(part in SKIP_DIR_NAMES for part in path.parts):
            continue
        found.append(path)
    return sorted(found)


def _sklearn_version_from_bytes(raw: bytes) -> str | None:
    match = SKLEARN_VERSION_RE.search(raw)
    if match is None:
        return None
    return match.group(1).decode("ascii")


def _try_pickle(path: Path) -> tuple[Any | None, str]:
    try:
        with path.open("rb") as handle:
            return pickle.load(handle), "pickle"
    except Exception as exc:
        return None, f"pickle={type(exc).__name__}: {exc}"


def _try_joblib(path: Path) -> tuple[Any | None, str]:
    try:
        import joblib

        return joblib.load(path), "joblib"
    except Exception as exc:
        return None, f"joblib={type(exc).__name__}: {exc}"


def _load_estimator(path: Path) -> tuple[Any | None, str]:
    """Load a serialization without writing. Prefer native loader, then fallback."""
    if path.suffix.lower() == ".joblib":
        model, note = _try_joblib(path)
        if model is None:
            return None, f"unreadable: {note}"
        return model, note
    model, pickle_note = _try_pickle(path)
    if model is not None:
        return model, "pickle"
    model, joblib_note = _try_joblib(path)
    if model is not None:
        return model, f"joblib (pickle failed: {pickle_note})"
    return None, f"unreadable: {pickle_note}; {joblib_note}"


def _elasticnet_step(model: Any) -> Any | None:
    if hasattr(model, "named_steps"):
        steps = model.named_steps
        if "elasticnet" in steps:
            return steps["elasticnet"]
        return list(steps.values())[-1]
    if hasattr(model, "coef_"):
        return model
    return None


def _feature_names(model: Any) -> list[str]:
    names = getattr(model, "feature_names_in_", None)
    if names is not None:
        return [str(x) for x in names]
    if hasattr(model, "named_steps"):
        for step in reversed(list(model.named_steps.values())):
            raw = getattr(step, "feature_names_in_", None)
            if raw is not None:
                return [str(x) for x in raw]
    return []


def classify_feature_names(names: list[str]) -> str:
    """Return FIXTURE, CANDIDATE REAL, MIXED, or OTHER."""
    if not names:
        return "OTHER (no feature_names_in_)"
    fixture = all(FIXTURE_NAME_RE.match(n) is not None for n in names)
    illumina = all(ILLUMINA_NAME_RE.match(n) is not None for n in names)
    if fixture:
        return "FIXTURE"
    if illumina:
        return "CANDIDATE REAL"
    n_fix = sum(FIXTURE_NAME_RE.match(n) is not None for n in names)
    n_ilm = sum(ILLUMINA_NAME_RE.match(n) is not None for n in names)
    if n_fix and n_ilm:
        return "MIXED"
    if n_ilm:
        return "PARTIAL ILLUMINA (not all names match ^cg\\d{8}$)"
    return "OTHER"


def inspect_artifact(path: Path) -> dict[str, Any]:
    """Inspect one serialization. Never writes the file."""
    stat = path.stat()
    mtime = datetime.fromtimestamp(stat.st_mtime).isoformat(timespec="seconds")
    raw = path.read_bytes()
    sklearn_ver = _sklearn_version_from_bytes(raw)
    model, loader = _load_estimator(path)
    record: dict[str, Any] = {
        "path": str(path.relative_to(REPO_ROOT)),
        "mtime": mtime,
        "size_bytes": int(stat.st_size),
        "loader": loader,
        "sklearn_version_recorded": sklearn_ver or "not recorded on object",
        "type": None,
        "pipeline_steps": None,
        "n_features": None,
        "feature_name_pattern": None,
        "feature_flag": "OTHER",
        "n_nonzero_coefficients": None,
        "intercept": None,
        "feature_names_head": [],
        "error": None,
    }
    if model is None:
        record["error"] = loader
        return record

    record["type"] = f"{type(model).__module__}.{type(model).__name__}"
    if hasattr(model, "named_steps"):
        record["pipeline_steps"] = " → ".join(
            f"{name} ({type(step).__name__})" for name, step in model.named_steps.items()
        )
    else:
        record["pipeline_steps"] = "(bare estimator, no Pipeline)"

    names = _feature_names(model)
    record["n_features"] = len(names) if names else getattr(model, "n_features_in_", None)
    record["feature_names_head"] = names[:8]
    record["feature_flag"] = classify_feature_names(names)
    if names:
        if record["feature_flag"] == "FIXTURE":
            record["feature_name_pattern"] = f"^cg_test_\\d+$ (all {len(names)} names)"
        elif record["feature_flag"] == "CANDIDATE REAL":
            record["feature_name_pattern"] = f"^cg\\d{{8}}$ (all {len(names)} names)"
        else:
            record["feature_name_pattern"] = f"first={names[0]!r} last={names[-1]!r} n={len(names)}"
    else:
        record["feature_name_pattern"] = "unavailable"

    enet = _elasticnet_step(model)
    if enet is not None and hasattr(enet, "coef_"):
        coef = np.ravel(enet.coef_)
        record["n_nonzero_coefficients"] = int(np.count_nonzero(coef))
        intercept = getattr(enet, "intercept_", None)
        if intercept is not None:
            record["intercept"] = float(np.ravel(intercept)[0])
        step_ver = getattr(enet, "_sklearn_version", None)
        if step_ver and record["sklearn_version_recorded"] == "not recorded on object":
            record["sklearn_version_recorded"] = str(step_ver)
    return record


def _mentions_gse40279(path: Path, text: str) -> bool:
    name = path.as_posix().lower()
    if "gse40279" in name:
        return True
    return "gse40279" in text.lower()


def _skip_mae_file(path: Path, exclude: set[Path]) -> bool:
    """Skip the audit report so a later run does not cite its own table."""
    if path.name == AUDIT_REPORT_NAME:
        return True
    return path.resolve() in exclude


def find_mae_in_working_tree(
    root: Path,
    *,
    exclude: set[Path] | None = None,
) -> list[dict[str, str]]:
    """Find numeric MAE values in GSE40279 / clock-metrics files on disk."""
    excluded = {path.resolve() for path in (exclude or set())}
    hits: list[dict[str, str]] = []
    for path in root.rglob("*"):
        if not path.is_file() or path.suffix.lower() not in TEXT_SUFFIXES:
            continue
        if any(part in SKIP_DIR_NAMES for part in path.parts):
            continue
        if _skip_mae_file(path, excluded):
            continue
        try:
            text = path.read_text(encoding="utf-8")
        except (UnicodeDecodeError, OSError):
            continue
        rel = str(path.relative_to(root))
        gse = _mentions_gse40279(path, text)
        clock_metrics = path.name.endswith("clock_metrics.json") or path.name.endswith(
            "gse40279_train_metrics.json"
        )
        if not gse and not clock_metrics:
            # Activity 2.1.10 figure scripts state a target MAE for the clock workstream.
            if not re.search(r"MAE\s*~\s*\d", text) and "3.42 years" not in text:
                continue
        for match in MAE_JSON_RE.finditer(text):
            line_no = text[: match.start()].count("\n") + 1
            hits.append(
                {
                    "file": rel,
                    "line": str(line_no),
                    "value": match.group("value"),
                    "label": match.group("key"),
                    "source": "working-tree file",
                    "gse40279_context": "yes" if gse or clock_metrics else "clock-metrics file",
                    "snippet": text.splitlines()[line_no - 1].strip()[:200],
                }
            )
        for match in MAE_PROSE_RE.finditer(text):
            line_no = text[: match.start()].count("\n") + 1
            hits.append(
                {
                    "file": rel,
                    "line": str(line_no),
                    "value": match.group("value"),
                    "label": match.group("label").strip(),
                    "source": "working-tree file",
                    "gse40279_context": "yes" if gse else "prose MAE in clock-related file",
                    "snippet": text.splitlines()[line_no - 1].strip()[:200],
                }
            )
    # Deduplicate file+line+value
    uniq: dict[tuple[str, str, str], dict[str, str]] = {}
    for hit in hits:
        uniq[(hit["file"], hit["line"], hit["value"])] = hit
    return [uniq[k] for k in sorted(uniq)]


def find_mae_in_git(root: Path) -> list[dict[str, str]]:
    """Search commit messages and git-grep of the HEAD tree for GSE40279 MAE numbers."""
    hits: list[dict[str, str]] = []
    log = subprocess.run(
        ["git", "log", "--all", "--format=%h\t%s", "-i", "--grep=MAE"],
        cwd=root,
        check=False,
        capture_output=True,
        text=True,
    )
    for line in log.stdout.splitlines():
        if not line.strip():
            continue
        sha, _, subject = line.partition("\t")
        nums = re.findall(r"\d+\.\d+", subject)
        hits.append(
            {
                "file": f"git commit {sha}",
                "line": "subject",
                "value": nums[0] if nums else "(no numeric MAE in subject)",
                "label": "commit subject",
                "source": "git log --grep=MAE",
                "gse40279_context": "commit message",
                "snippet": subject.strip()[:200],
            }
        )
    gse_log = subprocess.run(
        ["git", "log", "--all", "--format=%h\t%s", "-i", "--grep=GSE40279"],
        cwd=root,
        check=False,
        capture_output=True,
        text=True,
    )
    for line in gse_log.stdout.splitlines():
        if not line.strip():
            continue
        sha, _, subject = line.partition("\t")
        nums = re.findall(r"MAE[^\d]{0,12}(\d+\.\d+)", subject, flags=re.I)
        hits.append(
            {
                "file": f"git commit {sha}",
                "line": "subject",
                "value": nums[0] if nums else "(no numeric MAE in subject)",
                "label": "GSE40279 commit",
                "source": "git log --grep=GSE40279",
                "gse40279_context": "yes",
                "snippet": subject.strip()[:200],
            }
        )
    grep = subprocess.run(
        [
            "git",
            "grep",
            "-n",
            "-I",
            "-E",
            r'"test_mae"|"mae":|MAE ~|MAE =',
            "--",
            "*.md",
            "*.json",
            "*.py",
            "*.ipynb",
        ],
        cwd=root,
        check=False,
        capture_output=True,
        text=True,
    )
    for line in grep.stdout.splitlines():
        if not line.strip():
            continue
        file_part, _, rest = line.partition(":")
        file_path = Path(file_part)
        if file_path.name == AUDIT_REPORT_NAME or "tests" in file_path.parts:
            continue
        line_no, _, snippet = rest.partition(":")
        val_m = re.search(r"(-?\d+\.\d+)", snippet)
        hits.append(
            {
                "file": file_part,
                "line": line_no,
                "value": val_m.group(1) if val_m else "(no numeric value on line)",
                "label": "git grep",
                "source": "git grep HEAD",
                "gse40279_context": "tracked file",
                "snippet": snippet.strip()[:200],
            }
        )
    return hits


def _json_object(path: Path) -> dict[str, Any]:
    if path.suffix.lower() != ".json" or not path.is_file():
        return {}
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return {}
    if isinstance(payload, dict):
        return payload
    return {}


def _recorded_feature_count(payload: dict[str, Any]) -> int | None:
    for key in ("n_cpgs_features", "n_features_used", "n_clock_cpgs"):
        value = payload.get(key)
        if isinstance(value, int):
            return value
    return None


def _probe_tokens(payload: dict[str, Any]) -> list[str]:
    tokens: list[str] = []
    for key in ("selected_cpgs", "missing_cpg_ids"):
        raw = payload.get(key)
        if isinstance(raw, list):
            tokens.extend(str(item) for item in raw)
    return tokens


def _points_at(payload: dict[str, Any], artifact_path: str) -> bool:
    model_path = payload.get("model_path")
    if not isinstance(model_path, str) or not artifact_path:
        return False
    return artifact_path in model_path.replace("\\", "/")


def _count_matches(payload: dict[str, Any], artifact: dict[str, Any]) -> bool:
    expected = artifact.get("n_features")
    if not isinstance(expected, int):
        return False
    return _recorded_feature_count(payload) == expected


def _matches_fixture(
    hit: dict[str, str], payload: dict[str, Any], artifact: dict[str, Any]
) -> bool:
    """True when this hit records the fixture's feature count and probe identity."""
    if not _count_matches(payload, artifact):
        return False
    tokens = _probe_tokens(payload)
    fixture_names = bool(tokens) and all(FIXTURE_NAME_RE.match(tok) for tok in tokens)
    if fixture_names or "cg_test_" in hit["snippet"]:
        return True
    return _points_at(payload, str(artifact["path"]))


def _matches_real(payload: dict[str, Any], artifact: dict[str, Any]) -> bool:
    """True when this hit records the candidate clock's feature count and Illumina ids."""
    if not _count_matches(payload, artifact):
        return False
    tokens = _probe_tokens(payload)
    if tokens and all(ILLUMINA_NAME_RE.match(tok) for tok in tokens):
        return True
    return _points_at(payload, str(artifact["path"]))


def _hardcoded_untraceable(hit: dict[str, str]) -> str | None:
    value = hit["value"]
    snippet = hit["snippet"].lower()
    path = hit["file"].lower()
    if value == "2.1" and ("simulat" in snippet or "proxy" in snippet or "mae ~" in snippet):
        return "hard-coded simulated target, not produced by any artifact"
    if value == "3.42" and ("dashboard" in path or "mockup" in path or "842" in snippet):
        return "dashboard mockup (N=842); no artifact has real Illumina probes or that N"
    return None


def verdict_for_mae(
    hit: dict[str, str],
    artifacts: list[dict[str, Any]],
    *,
    root: Path = REPO_ROOT,
) -> tuple[str, str]:
    """TRACEABLE only when the hit's own feature count and probe names match an artifact."""
    hardcoded = _hardcoded_untraceable(hit)
    if hardcoded is not None:
        return "UNTRACEABLE", hardcoded

    payload = _json_object(root / hit["file"])
    reals = [
        a for a in artifacts if a["feature_flag"] == "CANDIDATE REAL" and _matches_real(payload, a)
    ]
    if reals:
        return "TRACEABLE", ", ".join(str(a["path"]) for a in reals)

    fixtures = [
        a for a in artifacts if a["feature_flag"] == "FIXTURE" and _matches_fixture(hit, payload, a)
    ]
    if fixtures:
        names = ", ".join(str(a["path"]) for a in fixtures)
        return (
            "TRACEABLE",
            names + " (FIXTURE — not a GSE40279 clock; probe names are cg_test_*)",
        )
    return (
        "UNTRACEABLE",
        "no artifact shares this occurrence's feature count and probe names or model path",
    )


def _md_escape(text: str) -> str:
    return text.replace("|", "\\|").replace("\n", " ")


def _same_fit(left: dict[str, Any], right: dict[str, Any]) -> bool:
    keys = (
        "intercept",
        "n_nonzero_coefficients",
        "n_features",
        "feature_flag",
        "feature_names_head",
    )
    return all(left.get(key) == right.get(key) for key in keys)


def _headline(artifacts: list[dict[str, Any]]) -> str:
    reals = [a for a in artifacts if a["feature_flag"] == "CANDIDATE REAL"]
    fixtures = [a for a in artifacts if a["feature_flag"] == "FIXTURE"]
    if reals:
        return (
            "A candidate-real GSE40279 clock **may** exist: "
            + ", ".join(str(a["path"]) for a in reals)
            + ". Confirm against training data before treating it as the Hannum cohort model."
        )
    if not fixtures:
        return (
            "**No loadable clock with Illumina or fixture probe names was found in this checkout.**"
        )
    patterns = ", ".join(sorted({str(a["feature_name_pattern"]) for a in fixtures}))
    named = [a for a in fixtures if "gse40279" in str(a["path"]).lower()]
    named_bit = ""
    if named:
        named_bit = (
            " Paths containing `gse40279` that are flagged FIXTURE are demo clocks: "
            + ", ".join(f"`{a['path']}`" for a in named)
            + "."
        )
    return (
        "**No real GSE40279-trained clock exists in this checkout.** "
        "No artifact has Illumina probe IDs matching `^cg\\d{8}$`. "
        f"Loadable clocks flagged FIXTURE use {patterns}." + named_bit
    )


def _name_vs_contents(artifacts: list[dict[str, Any]]) -> list[str]:
    lines = ["### Name vs contents", ""]
    named = [a for a in artifacts if "gse40279" in str(a["path"]).lower()]
    if not named:
        lines.extend(["No artifact path contains `gse40279`.", ""])
    for artifact in named:
        lines.append(
            f"`{artifact['path']}` is named as a GSE40279-trained model. "
            f"Inspection recorded {artifact['n_features']} features "
            f"({artifact['feature_name_pattern']}), flag {artifact['feature_flag']}, "
            f"{artifact['n_nonzero_coefficients']} non-zero coefficients, "
            f"intercept {artifact['intercept']}. Loader: {artifact['loader']}."
        )
        if artifact["feature_flag"] == "FIXTURE":
            lines.append(
                "Those feature names are the demo probe set (`cg_test_*`, "
                "`test_data/mock_clock_wide.csv`), not Hannum 2013 450K IDs."
            )
        head = artifact.get("feature_names_head") or []
        if head:
            lines.append("First feature names: " + ", ".join(str(name) for name in head) + ".")
        lines.append("")
    seen: set[frozenset[str]] = set()
    for artifact in artifacts:
        paths = frozenset(str(other["path"]) for other in artifacts if _same_fit(artifact, other))
        if len(paths) < 2 or paths in seen:
            continue
        seen.add(paths)
        joined = ", ".join(f"`{path}`" for path in sorted(paths))
        lines.append(
            f"{joined} match on feature count, probe-name head, non-zero coefficient count, "
            "and intercept."
        )
        for path in sorted(paths):
            match = next(item for item in artifacts if item["path"] == path)
            if "pickle failed" in str(match["loader"]):
                lines.append(f"`{path}` did not load with pickle; joblib did ({match['loader']}).")
        lines.append("")
    return lines


def _is_gitignored(root: Path, relative: str) -> bool:
    result = subprocess.run(
        ["git", "check-ignore", "-q", "--", relative],
        cwd=root,
        check=False,
        capture_output=True,
    )
    return result.returncode == 0


def _mae_search_note(mae_hits: list[dict[str, str]], *, root: Path = REPO_ROOT) -> str:
    git_numeric = [hit for hit in mae_hits if hit["source"].startswith("git log")]
    if git_numeric:
        subjects = "; ".join(hit["snippet"] for hit in git_numeric)
        git_sentence = f"Commit subjects with a numeric MAE: {subjects}."
    else:
        git_sentence = "No git commit subject that mentions MAE or GSE40279 contains a numeric MAE."
    ignored = "analysis/gse40279_train_metrics.json"
    ignore_sentence = ""
    if _is_gitignored(root, ignored):
        ignore_sentence = f"`{ignored}` is gitignored. "
    return (
        "Searched working-tree text files, `git grep` on HEAD, and commit subjects matching "
        "MAE or GSE40279. `CLOCK_ARTIFACT_AUDIT.md` is excluded so a later run does not cite "
        "its own table. Format strings and commit subjects with no number are omitted. "
        + git_sentence
        + " "
        + ignore_sentence
        + "TRACEABLE means this occurrence records the same feature count as an artifact and "
        "either names that artifact or uses its probe-id pattern. A FIXTURE trace is not a "
        "Hannum/GSE40279 clock."
    )


def render_report(
    artifacts: list[dict[str, Any]],
    mae_hits: list[dict[str, str]],
    *,
    root: Path = REPO_ROOT,
) -> str:
    reals = [a for a in artifacts if a["feature_flag"] == "CANDIDATE REAL"]
    fixtures = [a for a in artifacts if a["feature_flag"] == "FIXTURE"]
    headline = _headline(artifacts)

    lines = [
        "# Clock artifact audit",
        "",
        headline,
        "",
        "Generated by `scripts/audit_clock_artifacts.py`. Models were opened read-only. "
        "Nothing was renamed, moved, retrained, or regenerated.",
        "",
        "## Artifacts (every `.pkl` / `.joblib` / `.sav` in this checkout)",
        "",
        "| Path | n features | Feature name pattern | Flag | n non-zero coefs | Intercept | Pipeline steps | sklearn version in object | mtime |",
        "|------|------------|----------------------|------|------------------|-----------|----------------|---------------------------|-------|",
    ]
    for a in artifacts:
        lines.append(
            "| "
            + " | ".join(
                _md_escape(str(x))
                for x in (
                    a["path"],
                    a["n_features"],
                    a["feature_name_pattern"],
                    a["feature_flag"],
                    a["n_nonzero_coefficients"],
                    a["intercept"],
                    a["pipeline_steps"],
                    a["sklearn_version_recorded"],
                    a["mtime"],
                )
            )
            + " |"
        )
    lines.append("")
    lines.extend(_name_vs_contents(artifacts))
    lines.extend(
        [
            "## GSE40279 baseline MAE occurrences",
            "",
            _mae_search_note(mae_hits, root=root),
            "",
            "| File | Line | Value | Verdict | Could an artifact here have produced it? |",
            "|------|------|-------|---------|------------------------------------------|",
        ]
    )
    for hit in mae_hits:
        verdict, reason = verdict_for_mae(hit, artifacts, root=root)
        lines.append(
            "| "
            + " | ".join(
                _md_escape(str(x))
                for x in (
                    hit["file"],
                    hit["line"],
                    hit["value"],
                    verdict,
                    reason,
                )
            )
            + " |"
        )
    n_fix = len(fixtures)
    n_real = len(reals)
    lines.extend(
        [
            "",
            "## Counts",
            "",
            f"- Artifacts inspected: {len(artifacts)}",
            f"- Flagged FIXTURE: {n_fix}",
            f"- Flagged CANDIDATE REAL: {n_real}",
            f"- MAE occurrence rows: {len(mae_hits)}",
            "",
            "If CANDIDATE REAL is zero, there is no GSE40279 clock to evaluate. "
            "Do not substitute the fixture.",
            "",
        ]
    )
    return "\n".join(lines) + "\n"


def _is_numeric_mae(value: str) -> bool:
    return bool(re.fullmatch(r"-?\d+\.\d+", value)) and value not in {"0.0", "1.0"}


def _filter_mae_hits(hits: list[dict[str, str]]) -> list[dict[str, str]]:
    """Keep numeric GSE40279 / clock MAE values; drop format strings and empty git subjects."""
    kept: list[dict[str, str]] = []
    seen: set[tuple[str, str, str]] = set()
    for hit in hits:
        path = hit["file"].replace("\\342\\200\\224", "—").strip('"')
        hit = {**hit, "file": path}
        if not _is_numeric_mae(hit["value"]):
            continue
        key = (hit["file"], hit["line"], hit["value"])
        if key in seen:
            continue
        snippet = hit["snippet"]
        gse = "gse40279" in path.lower() or "gse40279" in snippet.lower()
        clock_json = "clock_metrics.json" in path or "gse40279_train_metrics" in path
        simulated = "MAE ~" in snippet or "Simulate" in snippet
        dashboard = "3.42" in hit["value"]
        git_hist = hit["source"].startswith("git ")
        if not (gse or clock_json or simulated or dashboard or git_hist):
            continue
        seen.add(key)
        kept.append(hit)
    return kept


def main(
    report: Path = typer.Option(
        DEFAULT_REPORT,
        "--report",
        help="Markdown audit path. Default: docs/CLOCK_ARTIFACT_AUDIT.md",
    ),
) -> None:
    """Inspect clock serializations and write the GSE40279 provenance audit."""
    artifacts = [inspect_artifact(p) for p in discover_artifacts(REPO_ROOT)]
    mae_hits = _filter_mae_hits(
        find_mae_in_working_tree(REPO_ROOT, exclude={report.resolve()}) + find_mae_in_git(REPO_ROOT)
    )
    text = render_report(artifacts, mae_hits)
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text(text, encoding="utf-8")
    n_real = sum(a["feature_flag"] == "CANDIDATE REAL" for a in artifacts)
    typer.echo(text)
    typer.echo(f"Wrote {report}")
    if n_real == 0:
        typer.echo("CONCLUSION: no artifact with real Illumina probe IDs (^cg\\d{8}$).")


if __name__ == "__main__":
    typer.run(main)
