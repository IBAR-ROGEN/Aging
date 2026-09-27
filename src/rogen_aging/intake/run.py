"""Run checksum verification, pseudonymization, and the leak check.

Outputs are built in a temporary directory and moved into ``out_dir`` only
after the leak check passes. This module does not make network calls and
does not log the pseudonym key.
"""

from __future__ import annotations

import shutil
import tempfile
from pathlib import Path

from rogen_aging.config import find_repo_root
from rogen_aging.intake.errors import IntakeError
from rogen_aging.intake.key import load_pseudonym_key
from rogen_aging.intake.leakcheck import leak_check
from rogen_aging.intake.manifest import (
    ManifestRow,
    identifiers_to_scrub,
    load_manifest,
    verify_checksums,
)
from rogen_aging.intake.paths import is_within, resolve_lexical
from rogen_aging.intake.pseudonym import pseudonym_for
from rogen_aging.intake.report import (
    assert_report_excludes_original_ids,
    build_report,
    render_linkage,
    report_text,
)
from rogen_aging.intake.rewrite import rewrite_delivery


def run_intake(
    manifest: Path,
    delivery_dir: Path,
    out_dir: Path,
    custodian_dir: Path,
    repo_root: Path | None = None,
) -> Path:
    """Pseudonymize a delivery and write the report plus the linkage table.

    Args:
        manifest: TSV manifest path.
        delivery_dir: Directory containing the named delivery files.
        out_dir: Processing output. Receives pseudonymized files and
            ``intake_report.json`` only after the leak check passes.
        custodian_dir: Directory for ``linkage.tsv``. Must be outside
            ``out_dir`` and outside the repository.
        repo_root: Repository root. Discovered when omitted.

    Returns:
        ``out_dir`` after a successful move.

    Raises:
        IntakeError: On manifest, checksum, key, location, rewrite, or leak
            failure. Those failures leave ``out_dir`` unwritten.
    """
    root = repo_root if repo_root is not None else find_repo_root()
    rows = load_manifest(manifest)
    verify_checksums(delivery_dir, rows)
    key = load_pseudonym_key(root)
    out_resolved = _destination(out_dir)
    custodian_resolved = _destination(custodian_dir)
    _refuse_locations(out_resolved, custodian_resolved, root)
    if out_resolved.exists() and any(out_resolved.iterdir()):
        raise IntakeError("output directory is not empty")
    linkage_path = custodian_resolved / "linkage.tsv"
    if linkage_path.exists():
        raise IntakeError("linkage table already exists")

    pseudonyms = _pseudonyms(rows, key)
    created_out = not out_resolved.exists()
    moved: list[Path] = []
    with tempfile.TemporaryDirectory(prefix="rogen-intake-") as temporary:
        work = Path(temporary)
        produced = rewrite_delivery(work, delivery_dir, rows, pseudonyms)
        scrub_ids = identifiers_to_scrub(rows)
        leak_check(produced, scrub_ids)
        report = build_report(produced, root)
        text = report_text(report)
        assert_report_excludes_original_ids(text, scrub_ids)
        (work / "intake_report.json").write_text(text, encoding="utf-8")
        allowed = {item.output_path.name for item in produced}
        allowed.add("intake_report.json")
        unexpected = [child.name for child in work.iterdir() if child.name not in allowed]
        if unexpected:
            raise IntakeError("intake workspace has an unexpected file; refusing to publish")
        try:
            out_resolved.mkdir(parents=True, exist_ok=True)
            for child in list(work.iterdir()):
                dest = out_resolved / child.name
                shutil.move(str(child), dest)
                moved.append(dest)
            custodian_resolved.mkdir(parents=True, exist_ok=True)
            linkage_path.write_text(render_linkage(rows, pseudonyms), encoding="utf-8")
        except IntakeError:
            _rollback(moved, linkage_path, out_resolved, created_out)
            raise
        except OSError as exc:
            _rollback(moved, linkage_path, out_resolved, created_out)
            raise IntakeError(f"could not publish intake output: {exc}") from exc
    return out_resolved


def _pseudonyms(rows: list[ManifestRow], key: bytes) -> dict[str, str]:
    mapping: dict[str, str] = {}
    owners: dict[str, str] = {}
    for row in rows:
        pseudonym = pseudonym_for(key, row.original_sample_id)
        previous = owners.get(pseudonym)
        if previous is not None and previous != row.original_sample_id:
            raise IntakeError("pseudonym collision; refusing to write output")
        owners[pseudonym] = row.original_sample_id
        mapping[row.original_sample_id] = pseudonym
    return mapping


def _destination(path: Path) -> Path:
    if path.exists():
        return path.resolve()
    return resolve_lexical(path)


def _refuse_locations(out_dir: Path, custodian_dir: Path, repo_root: Path) -> None:
    if is_within(custodian_dir, repo_root):
        raise IntakeError("custodian directory must be outside the repository")
    if custodian_dir == out_dir or is_within(custodian_dir, out_dir):
        raise IntakeError("custodian directory must be outside the output directory")


def _rollback(
    moved: list[Path],
    linkage_path: Path,
    out_dir: Path,
    created_out: bool,
) -> None:
    for path in moved:
        if path.is_dir():
            shutil.rmtree(path)
        elif path.exists() or path.is_symlink():
            path.unlink()
    if linkage_path.exists():
        linkage_path.unlink()
    if created_out and out_dir.exists():
        shutil.rmtree(out_dir)
