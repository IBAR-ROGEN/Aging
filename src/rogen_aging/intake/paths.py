"""Path checks that keep keys and linkage outside the repository."""

from __future__ import annotations

import os
from pathlib import Path


def resolve_lexical(path: Path) -> Path:
    """Absolute path with ``..`` collapsed, without following the final symlink."""
    absolute = path.expanduser()
    if not absolute.is_absolute():
        absolute = Path.cwd() / absolute
    return Path(os.path.abspath(absolute))


def is_within(path: Path, root: Path) -> bool:
    """Return True when ``path`` is ``root`` or a path inside it.

    Both the lexical path and, when the path exists, the symlink target are
    checked. A key path that merely points through the repository is refused.
    """
    root_resolved = root.resolve()
    candidates = [resolve_lexical(path)]
    if path.exists() or path.is_symlink():
        candidates.append(path.resolve())
    else:
        candidates.append(path.expanduser().resolve())
    for candidate in candidates:
        if candidate == root_resolved:
            return True
        try:
            candidate.relative_to(root_resolved)
        except ValueError:
            continue
        return True
    return False
