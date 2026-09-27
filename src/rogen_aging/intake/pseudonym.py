"""HMAC-SHA256 pseudonyms for original sample identifiers."""

from __future__ import annotations

import hashlib
import hmac

from rogen_aging.intake.errors import IntakeError

PSEUDONYM_PREFIX = "RO-"
PSEUDONYM_HEX_LEN = 12
MIN_KEY_BYTES = 32


def pseudonym_for(key: bytes, original_sample_id: str) -> str:
    """Return ``RO-`` plus the first 12 hex characters of HMAC-SHA256.

    Args:
        key: HMAC key. Must be at least 32 bytes.
        original_sample_id: Identifier from the delivery manifest.

    Returns:
        Pseudonym such as ``RO-0123456789ab``.

    Raises:
        IntakeError: If the key is shorter than 32 bytes.
    """
    if len(key) < MIN_KEY_BYTES:
        raise IntakeError("key is shorter than 32 bytes")
    digest = hmac.new(key, original_sample_id.encode("utf-8"), hashlib.sha256).hexdigest()
    return f"{PSEUDONYM_PREFIX}{digest[:PSEUDONYM_HEX_LEN]}"
