"""Checks that run before an aggregate table leaves the processing zone."""

from __future__ import annotations

from rogen_aging.privacy.release_guard import (
    ReleaseGuardError,
    ReleaseGuardResult,
    run_release_guard,
)

__all__ = [
    "ReleaseGuardError",
    "ReleaseGuardResult",
    "run_release_guard",
]
