"""Download GSE40279 (Hannum 2013) and write a wide beta matrix.

There is no GSE40279 loader elsewhere in this package. GSE87571 helpers in
``external_data`` download and parse GEO series matrices; this module reuses
those helpers and does not invent betas when the download fails.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pyarrow as pa  # type: ignore[import-untyped]
import pyarrow.parquet as pq  # type: ignore[import-untyped]

from rogen_aging.clock.external_data import (
    _TABLE_BEGIN_MARKERS,
    _TABLE_END_MARKERS,
    _download_url,
    _open_text,
    _parse_chronological_ages,
    _split_soft_line,
    _strip_geo_field,
)

logger = logging.getLogger(__name__)

GSE40279_ID = "GSE40279"
EXPECTED_N_SAMPLES = 656
SERIES_MATRIX_FILENAME = f"{GSE40279_ID}_series_matrix.txt.gz"
SERIES_MATRIX_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE40nnn/"
    f"{GSE40279_ID}/matrix/{SERIES_MATRIX_FILENAME}"
)
_MISSING_TOKENS = frozenset({"", "NA", "NaN", "nan", "null", "NULL"})
_CHUNK_PROBES = 20_000


class Gse40279SampleCountError(RuntimeError):
    """Raised when the GEO series does not contain the published sample count."""

    def __init__(self, found: int, expected: int = EXPECTED_N_SAMPLES) -> None:
        self.found = found
        self.expected = expected
        difference = found - expected
        super().__init__(
            f"GSE40279 sample count is {found}, expected {expected} "
            f"(difference {difference:+d}). No samples were dropped or added."
        )


def _parse_beta_row(cells: list[str], n_samples: int) -> np.ndarray:
    if len(cells) != n_samples:
        raise ValueError(f"Beta row has {len(cells)} values, expected {n_samples}.")
    try:
        return np.asarray(cells, dtype=np.float32)
    except ValueError:
        out = np.empty(n_samples, dtype=np.float32)
        for index, cell in enumerate(cells):
            if cell in _MISSING_TOKENS:
                out[index] = np.nan
            else:
                out[index] = np.float32(cell)
        return out


def _read_series_metadata(path: Path) -> tuple[dict[str, object], list[str]]:
    """Read ``!Sample_*`` rows and the matrix header. Does not store probe rows."""
    sample_rows: dict[str, object] = {}
    header: list[str] | None = None
    in_table = False
    with _open_text(path) as handle:
        for line in handle:
            stripped = line.strip()
            lower = stripped.lower()
            if not in_table:
                if any(lower == marker.lower() for marker in _TABLE_BEGIN_MARKERS):
                    in_table = True
                    continue
                if stripped.startswith("!Sample_"):
                    parts = _split_soft_line(stripped)
                    if not parts:
                        continue
                    key = parts[0]
                    values = [_strip_geo_field(item) for item in parts[1:]]
                    if key == "!Sample_characteristics_ch1":
                        characteristics = sample_rows.setdefault(key, [])
                        if not isinstance(characteristics, list):
                            raise TypeError("characteristics bucket must be a list")
                        characteristics.append(values)
                    else:
                        sample_rows[key] = values
                continue
            if any(lower == marker.lower() for marker in _TABLE_END_MARKERS):
                break
            if not stripped or stripped.startswith("!"):
                continue
            parts = _split_soft_line(stripped)
            if not parts:
                continue
            header = [_strip_geo_field(item) for item in parts]
            break
    if header is None:
        raise ValueError(f"No series matrix table header found in {path}.")
    return sample_rows, header


def _align_samples(
    gsms: list[str],
    ages: list[float | None],
    header: list[str],
) -> tuple[list[str], list[float]]:
    if header[0].upper() != "ID_REF":
        raise ValueError(f"Expected first matrix column 'ID_REF', got {header[0]!r}.")
    matrix_samples = header[1:]
    if len(matrix_samples) != len(gsms):
        raise Gse40279SampleCountError(len(matrix_samples), expected=len(gsms))
    if set(matrix_samples) != set(gsms):
        raise ValueError(
            "Series matrix sample columns do not match !Sample_geo_accession. "
            f"only_metadata={sorted(set(gsms) - set(matrix_samples))[:8]} "
            f"only_matrix={sorted(set(matrix_samples) - set(gsms))[:8]}. "
            "No samples were dropped or added."
        )
    age_by_gsm = {gsm: age for gsm, age in zip(gsms, ages, strict=True)}
    ordered_ages: list[float] = []
    missing_age: list[str] = []
    for gsm in matrix_samples:
        age = age_by_gsm[gsm]
        if age is None:
            missing_age.append(gsm)
        else:
            ordered_ages.append(float(age))
    if missing_age:
        raise ValueError(
            f"{len(missing_age)} of {len(matrix_samples)} GSE40279 samples have no parsed "
            f"chronological age (first: {missing_age[:8]}). Refusing to drop them."
        )
    return matrix_samples, ordered_ages


def _stream_cg_matrix(path: Path, n_samples: int) -> tuple[np.ndarray, list[str]]:
    """Second pass: CpG rows only, float32, shape ``(n_probes, n_samples)``."""
    chunks: list[np.ndarray] = []
    probe_ids: list[str] = []
    buffer = np.empty((_CHUNK_PROBES, n_samples), dtype=np.float32)
    buffer_names: list[str] = []
    filled = 0
    seen: set[str] = set()
    in_table = False
    skipped_header = False
    with _open_text(path) as handle:
        for line in handle:
            stripped = line.strip()
            lower = stripped.lower()
            if not in_table:
                if any(lower == marker.lower() for marker in _TABLE_BEGIN_MARKERS):
                    in_table = True
                continue
            if any(lower == marker.lower() for marker in _TABLE_END_MARKERS):
                break
            if not stripped or stripped.startswith("!"):
                continue
            parts = [_strip_geo_field(item) for item in _split_soft_line(stripped)]
            if not skipped_header:
                skipped_header = True
                continue
            probe_id = parts[0]
            if not probe_id.startswith("cg"):
                continue
            if probe_id in seen:
                raise ValueError(f"Duplicate CpG id {probe_id} in {path}.")
            seen.add(probe_id)
            buffer[filled] = _parse_beta_row(parts[1:], n_samples)
            buffer_names.append(probe_id)
            filled += 1
            if filled == _CHUNK_PROBES:
                chunks.append(buffer.copy())
                probe_ids.extend(buffer_names)
                logger.info("Read %s cg probes", len(probe_ids))
                buffer_names = []
                filled = 0
    if not skipped_header:
        raise ValueError(f"Series matrix table in {path} has no header row.")
    if filled:
        chunks.append(buffer[:filled].copy())
        probe_ids.extend(buffer_names)
    if not chunks:
        raise ValueError(f"No cg* beta rows found in {path}.")
    matrix = np.concatenate(chunks, axis=0)
    return matrix, probe_ids


def _write_wide_parquet(
    matrix: np.ndarray,
    probe_ids: list[str],
    sample_ids: list[str],
    ages: list[float],
    output_path: Path,
) -> None:
    """Write samples × probes. ``matrix`` is ``(n_probes, n_samples)`` float32."""
    n_probes, n_samples = matrix.shape
    if n_probes != len(probe_ids) or n_samples != len(sample_ids) or n_samples != len(ages):
        raise ValueError("Beta matrix shape does not match probe or sample labels.")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        pa.field("sample_id", pa.string()),
        pa.field("chronological_age", pa.float32()),
    ]
    fields.extend(pa.field(name, pa.float32()) for name in probe_ids)
    schema = pa.schema(fields)
    ages_arr = np.asarray(ages, dtype=np.float32)
    writer = pq.ParquetWriter(output_path, schema, compression="zstd")
    batch_size = 64
    for start in range(0, n_samples, batch_size):
        stop = min(start + batch_size, n_samples)
        block = np.ascontiguousarray(matrix[:, start:stop].T)
        arrays: list[pa.Array] = [
            pa.array(sample_ids[start:stop]),
            pa.array(ages_arr[start:stop]),
        ]
        arrays.extend(pa.array(block[:, probe], type=pa.float32()) for probe in range(n_probes))
        writer.write_batch(pa.record_batch(arrays, schema=schema))
    writer.close()


def ensure_gse40279_series_matrix(geo_cache_dir: Path) -> Path:
    """Download the series matrix into ``geo_cache_dir`` if it is not already there."""
    dest = Path(geo_cache_dir) / SERIES_MATRIX_FILENAME
    if dest.is_file() and dest.stat().st_size > 0:
        return dest
    logger.info("Downloading %s", SERIES_MATRIX_URL)
    _download_url(SERIES_MATRIX_URL, dest, timeout_s=3600)
    return dest


def prepare_gse40279(
    output_path: Path,
    *,
    geo_cache_dir: Path,
    series_matrix_path: Path | None = None,
    expected_n_samples: int = EXPECTED_N_SAMPLES,
) -> dict[str, int]:
    """Download GSE40279 if needed and write the wide float32 beta table.

    Args:
        output_path: Destination parquet (``sample_id``, ``chronological_age``, ``cg*``).
        geo_cache_dir: Directory for the GEO series matrix download.
        series_matrix_path: Local series matrix. Downloaded when omitted.
        expected_n_samples: Published cohort size. A mismatch aborts.

    Returns:
        ``n_samples`` and ``n_probes`` written to ``output_path``.

    Raises:
        Gse40279SampleCountError: If the series sample count is not ``expected_n_samples``.
        RuntimeError: If the GEO download fails. The error includes the URL.
    """
    if series_matrix_path is None:
        matrix_path = ensure_gse40279_series_matrix(geo_cache_dir)
    else:
        matrix_path = Path(series_matrix_path)
        if not matrix_path.is_file():
            raise FileNotFoundError(matrix_path)

    sample_rows, header = _read_series_metadata(matrix_path)
    gsms_raw = sample_rows.get("!Sample_geo_accession", [])
    if not isinstance(gsms_raw, list):
        raise ValueError("!Sample_geo_accession was not parsed.")
    gsms = [str(item) for item in gsms_raw]
    if len(gsms) != expected_n_samples:
        raise Gse40279SampleCountError(len(gsms), expected=expected_n_samples)
    ages = _parse_chronological_ages(sample_rows, len(gsms))
    sample_ids, ordered_ages = _align_samples(gsms, ages, header)
    logger.info("Streaming %s cg beta rows for %s samples", GSE40279_ID, len(sample_ids))
    betas, probe_ids = _stream_cg_matrix(matrix_path, len(sample_ids))
    _write_wide_parquet(betas, probe_ids, sample_ids, ordered_ages, output_path)
    return {"n_samples": len(sample_ids), "n_probes": len(probe_ids)}
