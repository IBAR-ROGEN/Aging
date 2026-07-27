"""Regression tests for the UKB CI compliance audit script."""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parent.parent
_AUDIT = _REPO_ROOT / "scripts" / "dev" / "ukbb_ci_compliance_audit.sh"


@pytest.mark.skipif(not _AUDIT.is_file(), reason="audit script missing")
def test_ukbb_ci_compliance_audit_passes_on_repo() -> None:
    """Full-tree audit must exit 0 on the current working tree."""
    result = subprocess.run(
        [str(_AUDIT)],
        cwd=_REPO_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stdout + "\n" + result.stderr
    assert "No blocking violations" in result.stdout


@pytest.mark.skipif(not _AUDIT.is_file(), reason="audit script missing")
def test_ukbb_ci_compliance_audit_ignores_label_wrap_slash_fragments(
    tmp_path: Path,
) -> None:
    """Label-wrapping slash strings must not raise a false absolute-path FAIL.

    The historical CI failure matched ``"/\\n"``-style fragments inside
    ``str.replace`` helpers used for figure labels.
    """
    sample = tmp_path / "figure_labels.py"
    sample.write_text(
        'cluster = "A/B"\n'
        'cluster = cluster.replace("/", "/" + "\\n")\n'
        'other = "/\\n"\n'
        'allowed = "/tmp/scratch"\n',
        encoding="utf-8",
    )
    env = os.environ.copy()
    env["REPO_ROOT"] = str(tmp_path)
    result = subprocess.run(
        [str(_AUDIT)],
        cwd=_REPO_ROOT,
        capture_output=True,
        text=True,
        check=False,
        env=env,
    )
    assert result.returncode == 0, result.stdout + "\n" + result.stderr
    assert "disallowed absolute path" not in result.stdout
    assert "No blocking violations" in result.stdout
