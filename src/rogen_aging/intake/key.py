"""Pseudonym key loading and creation.

The key is read from ``ROGEN_PSEUDO_KEY_FILE``. It is never logged.
"""

from __future__ import annotations

import os
import secrets
import stat
import sys
from pathlib import Path

from loguru import logger

from rogen_aging.config import find_repo_root
from rogen_aging.intake.errors import IntakeError
from rogen_aging.intake.paths import is_within
from rogen_aging.intake.pseudonym import MIN_KEY_BYTES

KEY_ENV_VAR = "ROGEN_PSEUDO_KEY_FILE"
_GROUP_OR_OTHER_READ = stat.S_IRGRP | stat.S_IROTH


def load_pseudonym_key(repo_root: Path | None = None) -> bytes:
    """Read the HMAC key from ``ROGEN_PSEUDO_KEY_FILE``.

    Args:
        repo_root: Repository root. Discovered when omitted.

    Returns:
        Key bytes, length at least 32.

    Raises:
        IntakeError: If the variable is unset, the file is inside the
            repository, the POSIX mode allows group or other reads, or the
            key is shorter than 32 bytes.
    """
    root = repo_root if repo_root is not None else find_repo_root()
    raw = os.environ.get(KEY_ENV_VAR)
    if raw is None or not raw.strip():
        raise IntakeError(f"{KEY_ENV_VAR} is unset")
    path = Path(raw).expanduser()
    if is_within(path, root):
        raise IntakeError("refusing to read a key file inside the repository")
    if not path.is_file():
        raise IntakeError(f"{KEY_ENV_VAR} does not point to a file")
    if sys.platform == "win32":
        logger.warning(
            "POSIX mode check skipped on Windows; key file permissions were not verified"
        )
    else:
        mode = path.stat().st_mode
        if mode & _GROUP_OR_OTHER_READ:
            raise IntakeError("key file is readable by group or others; expected mode 600")
    data = path.read_bytes()
    if len(data) < MIN_KEY_BYTES:
        raise IntakeError("key is shorter than 32 bytes")
    return data


def make_key(out: Path, repo_root: Path | None = None) -> None:
    """Write a random 32-byte key with mode 600.

    Args:
        out: Destination path. Must be outside the repository.
        repo_root: Repository root. Discovered when omitted.

    Raises:
        IntakeError: If ``out`` is inside the repository or already exists.
    """
    root = repo_root if repo_root is not None else find_repo_root()
    if is_within(out, root):
        raise IntakeError("refusing to write a key inside the repository")
    target = out.expanduser()
    if not target.is_absolute():
        target = Path.cwd() / target
    parent = target.parent
    if is_within(parent, root) or is_within(target, root):
        raise IntakeError("refusing to write a key inside the repository")
    parent.mkdir(parents=True, exist_ok=True)
    if is_within(target, root):
        raise IntakeError("refusing to write a key inside the repository")
    key = secrets.token_bytes(MIN_KEY_BYTES)
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
    try:
        descriptor = os.open(target, flags, 0o600)
    except FileExistsError as exc:
        raise IntakeError(f"key file already exists: {target}") from exc
    try:
        os.write(descriptor, key)
    finally:
        os.close(descriptor)
    if sys.platform == "win32":
        logger.warning(
            "POSIX mode 600 cannot be enforced on Windows; key file permissions were not verified"
        )
    else:
        os.chmod(target, 0o600)
    logger.info("Wrote a 32-byte key with mode 600")
