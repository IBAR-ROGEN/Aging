# Note: GSE87571 is a public dataset suitable for external validation of methylation clocks trained on GSE40279.
"""Load GSE87571 (Illumina 450K whole blood) for external clock validation.

The GEO series matrix for GSE87571 contains sample metadata but no probe-level
beta matrix; betas are published in supplementary ``matrix1of2`` and
``matrix2of2`` files, which this module downloads and merges when needed.
"""

from __future__ import annotations

import gzip
import io
import re
import shutil
import urllib.error
import urllib.request
from pathlib import Path
from typing import IO, Annotated, Any, Sequence

import pandas as pd
import requests
import typer

GSE_ID = "GSE87571"
SERIES_MATRIX_FILENAME = f"{GSE_ID}_series_matrix.txt.gz"
MATRIX1_FILENAME = f"{GSE_ID}_matrix1of2.txt.gz"
MATRIX2_FILENAME = f"{GSE_ID}_matrix2of2.txt.gz"

SERIES_MATRIX_URL = (
    f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE87nnn/{GSE_ID}/matrix/{SERIES_MATRIX_FILENAME}"
)
SUPP_MATRIX1_URL = (
    f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE87nnn/{GSE_ID}/suppl/{MATRIX1_FILENAME}"
)
SUPP_MATRIX2_URL = (
    f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE87nnn/{GSE_ID}/suppl/{MATRIX2_FILENAME}"
)
GEO_QUERY_URL = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={GSE_ID}"

_TABLE_BEGIN_MARKERS = ("!series_matrix_table_begin", "!Series_matrix_table_begin")
_TABLE_END_MARKERS = ("!series_matrix_table_end", "!Series_matrix_table_end")

_TITLE_LABEL_RE = re.compile(r'^"?(X)(\d+)\b', re.IGNORECASE)
_AGE_VALUE_RE = re.compile(
    r"(?i)\bage\s*(?:\([^)]*\))?\s*[:=]\s*([\d.]+)",
)


def save_as_parquet(df: pd.DataFrame, output_path: str | Path) -> None:
    """Write ``df`` to Parquet using the pyarrow engine (requires pyarrow)."""
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(path, index=True, engine="pyarrow")


def _strip_geo_field(value: str) -> str:
    s = value.strip()
    if len(s) >= 2 and s[0] == '"' and s[-1] == '"':
        return s[1:-1]
    return s


def _open_text(path: Path) -> IO[str]:
    if path.suffix.lower() == ".gz":
        raw = gzip.open(path, "rb")
        return io.TextIOWrapper(raw, encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def _download_url(url: str, dest: Path, timeout_s: int = 600) -> None:
    """Stream ``url`` to ``dest`` using ``requests`` first, then stdlib ``urllib``."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".partial")
    attempts: list[str] = []

    def cleanup_tmp() -> None:
        if tmp.exists():
            tmp.unlink(missing_ok=True)

    try:
        with requests.get(url, stream=True, timeout=timeout_s) as resp:
            resp.raise_for_status()
            with tmp.open("wb") as out:
                for block in resp.iter_content(chunk_size=1024 * 1024):
                    if block:
                        out.write(block)
        tmp.replace(dest)
        return
    except (requests.RequestException, TimeoutError, OSError) as exc:
        cleanup_tmp()
        attempts.append(f"requests: {type(exc).__name__}: {exc}")

    try:
        with urllib.request.urlopen(url, timeout=timeout_s) as resp:
            with tmp.open("wb") as out:
                shutil.copyfileobj(resp, out, length=1024 * 1024)
        tmp.replace(dest)
        return
    except urllib.error.URLError as exc:
        cleanup_tmp()
        attempts.append(f"urllib: {type(exc).__name__}: {exc}")
    except TimeoutError as exc:
        cleanup_tmp()
        attempts.append(f"urllib: TimeoutError: {exc}")

    hint = (
        " If you see SSL or certificate errors, configure system trust / REQUESTS_CA_BUNDLE, "
        "or download the files manually from the GEO accession page."
    )
    raise RuntimeError(
        f"Failed to download {url!r} to {dest!r}.{hint} "
        f"Tried: {'; '.join(attempts)}. "
        f"Manual page: {GEO_QUERY_URL} (series matrix ``{SERIES_MATRIX_FILENAME}``; "
        f"supplementary ``{MATRIX1_FILENAME}`` and ``{MATRIX2_FILENAME}`` when the "
        f"embedded matrix table is empty)."
    )


def _try_geoparse_download_series_matrix(geo_cache_dir: Path) -> Path | None:
    # GEOparse is sometimes flaky with large series matrices; if the fetch fails,
    # recommend manual download from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE87571
    # and using local_path.
    geo_cache_dir.mkdir(parents=True, exist_ok=True)
    dest = geo_cache_dir / SERIES_MATRIX_FILENAME
    if dest.is_file():
        return dest
    try:
        import GEOparse  # noqa: PLC0415 — optional heavy import at runtime

        GEOparse.get_GEO(geo=GSE_ID, destdir=str(geo_cache_dir), silent=True)
    except ImportError:
        return None
    except Exception:
        return None
    if dest.is_file():
        return dest
    candidates = list(geo_cache_dir.glob(f"{GSE_ID}*_series_matrix.txt.gz"))
    if len(candidates) == 1:
        return candidates[0]
    if candidates:
        candidates.sort(key=lambda p: p.stat().st_mtime, reverse=True)
        return candidates[0]
    return None


def _ensure_series_matrix(geo_cache_dir: Path) -> Path:
    geo_cache_dir = Path(geo_cache_dir)
    dest = geo_cache_dir / SERIES_MATRIX_FILENAME
    if dest.is_file():
        return dest
    geoparse_path = _try_geoparse_download_series_matrix(geo_cache_dir)
    if geoparse_path is not None and geoparse_path.is_file():
        if geoparse_path.resolve() != dest.resolve():
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(geoparse_path, dest)
        return dest
    _download_url(SERIES_MATRIX_URL, dest)
    return dest


def _ensure_supplementary_matrix(url: str, filename: str, geo_cache_dir: Path) -> Path:
    geo_cache_dir = Path(geo_cache_dir)
    dest = geo_cache_dir / filename
    if dest.is_file():
        return dest
    _download_url(url, dest)
    return dest


def _split_soft_line(line: str) -> list[str]:
    line = line.rstrip("\n\r")
    if not line:
        return []
    return line.split("\t")


def _parse_title_to_gsm(sample_rows: dict[str, Any]) -> dict[str, str]:
    titles = sample_rows.get("!Sample_title", [])
    gsms = sample_rows.get("!Sample_geo_accession", [])
    if not titles or not gsms:
        raise ValueError(
            "Series matrix is missing !Sample_title or !Sample_geo_accession; cannot map columns."
        )
    if len(titles) != len(gsms):
        raise ValueError(
            f"!Sample_title length ({len(titles)}) != !Sample_geo_accession length ({len(gsms)})."
        )
    label_to_gsm: dict[str, str] = {}
    for title, gsm in zip(titles, gsms, strict=True):
        m = _TITLE_LABEL_RE.match(title.strip())
        if not m:
            continue
        label = f"X{m.group(2)}"
        label_to_gsm[label] = gsm
    if not label_to_gsm:
        raise ValueError("Could not parse any X-labels from !Sample_title for supplementary matrix.")
    return label_to_gsm


def _parse_chronological_ages(sample_rows: dict[str, Any], n_samples: int) -> list[float | None]:
    """GEO repeats ``!Sample_characteristics_ch1`` on multiple lines (gender, age, tissue, …)."""
    ages: list[float | None] = [None] * n_samples
    raw = sample_rows.get("!Sample_characteristics_ch1")
    rows: list[list[str]]
    if raw is None:
        return ages
    if isinstance(raw, list) and raw and isinstance(raw[0], list):
        rows = raw
    else:
        rows = [raw]
    for values in rows:
        if len(values) != n_samples:
            continue
        for i, cell in enumerate(values):
            m = _AGE_VALUE_RE.search(cell)
            if m:
                ages[i] = float(m.group(1))
    return ages


def _parse_series_matrix_file_full(path: Path) -> tuple[dict[str, list[str]], list[str], list[list[str]]]:
    """Parse file in one pass: metadata before table, then matrix block."""
    sample_rows: dict[str, list[str]] = {}
    header: list[str] | None = None
    body: list[list[str]] = []
    in_table = False
    with _open_text(path) as handle:
        for line in handle:
            stripped = line.strip()
            lower = stripped.lower()
            if not in_table:
                if any(lower == m.lower() for m in _TABLE_BEGIN_MARKERS):
                    in_table = True
                    continue
                if stripped.startswith("!Sample_"):
                    parts = _split_soft_line(stripped)
                    if parts:
                        key = parts[0]
                        vals = [_strip_geo_field(x) for x in parts[1:]]
                        if key == "!Sample_characteristics_ch1":
                            sample_rows.setdefault(key, []).append(vals)
                        else:
                            sample_rows[key] = vals
                continue
            if any(lower == m.lower() for m in _TABLE_END_MARKERS):
                break
            if not stripped or stripped.startswith("!"):
                continue
            parts = _split_soft_line(stripped)
            if not parts:
                continue
            if header is None:
                header = [_strip_geo_field(p) for p in parts]
                continue
            body.append([_strip_geo_field(p) for p in parts])
    if header is None:
        raise ValueError(f"No series matrix table begin marker found in {path}.")
    return sample_rows, header, body


def _is_beta_column_name(name: str) -> bool:
    return bool(re.match(r"^X\d+$", name))


def _load_supplementary_merged(
    geo_cache_dir: Path,
    restrict_to_cpgs: Sequence[str] | None,
) -> pd.DataFrame:
    p1 = _ensure_supplementary_matrix(SUPP_MATRIX1_URL, MATRIX1_FILENAME, geo_cache_dir)
    p2 = _ensure_supplementary_matrix(SUPP_MATRIX2_URL, MATRIX2_FILENAME, geo_cache_dir)
    read_kw: dict[str, Any] = {"sep": "\t", "compression": "gzip", "low_memory": False}
    if restrict_to_cpgs is not None:
        head_kw = {"sep": "\t", "compression": "gzip", "nrows": 0}
        head1 = pd.read_csv(p1, **head_kw)
        head2 = pd.read_csv(p2, **head_kw)
        cpg_set = {str(x) for x in restrict_to_cpgs if str(x).startswith("cg")}
        if not cpg_set:
            raise ValueError(
                "restrict_to_cpgs contained no probe IDs starting with 'cg'; nothing to load."
            )
        cols1 = ["ID_REF"] + [c for c in head1.columns if c != "ID_REF" and _is_beta_column_name(c)]
        cols2 = [c for c in head2.columns if c != "ID_REF" and _is_beta_column_name(c)]

        def read_filtered(path: Path, use_columns: list[str]) -> pd.DataFrame:
            chunks: list[pd.DataFrame] = []
            for chunk in pd.read_csv(path, usecols=use_columns, chunksize=50_000, **read_kw):
                sub = chunk.loc[chunk["ID_REF"].isin(cpg_set)]
                if not sub.empty:
                    chunks.append(sub)
            if not chunks:
                return pd.DataFrame(columns=use_columns)
            return pd.concat(chunks, ignore_index=True)

        left = read_filtered(p1, cols1)
        right = read_filtered(p2, ["ID_REF"] + cols2)
        merged = left.merge(right, on="ID_REF", how="outer")
        if merged.empty:
            raise ValueError(
                "No matching CpG rows were found in supplementary matrices for restrict_to_cpgs."
            )
    else:
        left = pd.read_csv(p1, **read_kw)
        right = pd.read_csv(p2, **read_kw)
        id2 = right.drop(columns=["ID_REF"])
        merged = pd.concat([left, id2], axis=1)
    merged = merged.set_index("ID_REF")
    beta_cols = [c for c in merged.columns if _is_beta_column_name(str(c))]
    return merged[beta_cols].apply(pd.to_numeric, errors="coerce")


def _probe_table_to_beta_wide(header: list[str], body: list[list[str]]) -> pd.DataFrame | None:
    if not body:
        return None
    if header[0].upper() != "ID_REF":
        raise ValueError(f"Expected first matrix column 'ID_REF', got {header[0]!r}.")
    data = pd.DataFrame(body, columns=header)
    data = data.set_index("ID_REF")
    numeric = data.apply(pd.to_numeric, errors="coerce")
    return numeric


def _build_output_frame(
    beta_wide_probes_samples: pd.DataFrame,
    label_to_gsm: dict[str, str],
    ages_by_gsm_order: list[float | None],
    gsms_ordered: list[str],
    restrict_to_cpgs: Sequence[str] | None,
) -> pd.DataFrame:
    rename_map: dict[str, str] = {}
    missing_labels: list[str] = []
    for col in beta_wide_probes_samples.columns:
        col_s = str(col)
        if not _is_beta_column_name(col_s):
            continue
        gsm = label_to_gsm.get(col_s)
        if gsm is None:
            missing_labels.append(col_s)
            continue
        rename_map[col_s] = gsm
    if missing_labels:
        raise ValueError(
            f"{len(missing_labels)} matrix columns have no GSM mapping from !Sample_title "
            f"(first few: {missing_labels[:8]})."
        )
    mat = beta_wide_probes_samples.rename(columns=rename_map)
    gsm_to_age: dict[str, float] = {}
    for gsm, age in zip(gsms_ordered, ages_by_gsm_order, strict=True):
        if age is not None:
            gsm_to_age[str(gsm)] = float(age)
    samples = [s for s in mat.columns if s in gsm_to_age]
    if not samples:
        raise ValueError("No samples with a parsed chronological age; check metadata.")
    out = mat[samples].T
    out.insert(0, "chronological_age", [gsm_to_age[s] for s in samples])
    if restrict_to_cpgs is not None:
        wanted = [c for c in restrict_to_cpgs if str(c).startswith("cg")]
        present = [c for c in wanted if c in out.columns]
        missing = sorted(set(wanted) - set(present))
        if missing:
            raise ValueError(
                f"{len(missing)} CpGs from restrict_to_cpgs are absent in GSE87571 after merge "
                f"(showing up to 12): {missing[:12]}"
            )
        out = out[["chronological_age", *present]]
    return out


def load_gse87571(
    local_path: str | Path | None = None,
    geo_cache_dir: str | Path = "./data/geo",
    restrict_to_cpgs: Sequence[str] | None = None,
) -> pd.DataFrame:
    """Load GSE87571 methylation beta values and chronological ages.

    Returns a DataFrame with rows = samples (GEO GSM accessions as the index),
    columns = CpG probe IDs (``cg*``) plus ``chronological_age``, suitable for
    ``scripts/validate_clock.py``.

    If ``local_path`` points to an existing series matrix (``.txt`` or
    ``.txt.gz``), it is used for metadata. When the embedded matrix table has
    no probe rows (true for GSE87571), supplementary ``matrix1of2`` /
    ``matrix2of2`` files are downloaded from NCBI FTP into ``geo_cache_dir``
    and merged.

    Parameters
    ----------
    local_path:
        Optional path to a local ``GSE87571_series_matrix.txt.gz`` (or plain
        ``.txt``) file.
    geo_cache_dir:
        Directory for GEO downloads and caches.
    restrict_to_cpgs:
        If set, keep only these CpG IDs (intersection must be non-empty).

    Raises
    ------
    FileNotFoundError
        If ``local_path`` is set but the file does not exist.
    RuntimeError
        On repeated network / GEOparse failures (see error text for manual steps).
    """
    geo_cache_dir = Path(geo_cache_dir)
    series_path: Path
    if local_path is not None:
        series_path = Path(local_path).expanduser().resolve()
        if not series_path.is_file():
            raise FileNotFoundError(
                f"local_path was set but the file does not exist: {series_path}"
            )
    else:
        series_path = _ensure_series_matrix(geo_cache_dir)

    sample_rows, table_header, table_body = _parse_series_matrix_file_full(series_path)
    gsms = sample_rows.get("!Sample_geo_accession", [])
    n = len(gsms)
    if n == 0:
        raise ValueError("No samples (!Sample_geo_accession) found in series matrix.")
    ages = _parse_chronological_ages(sample_rows, n)
    label_to_gsm = _parse_title_to_gsm(sample_rows)

    beta_from_series = _probe_table_to_beta_wide(table_header, table_body)
    if beta_from_series is not None and beta_from_series.shape[0] > 0:
        beta_wide = beta_from_series
        if restrict_to_cpgs is not None:
            cpg_set = {str(x) for x in restrict_to_cpgs if str(x).startswith("cg")}
            beta_wide = beta_wide.loc[beta_wide.index.astype(str).isin(cpg_set)]
            if beta_wide.shape[0] == 0:
                raise ValueError("restrict_to_cpgs did not match any probes in the series matrix table.")
    else:
        beta_wide = _load_supplementary_merged(geo_cache_dir, restrict_to_cpgs)

    out = _build_output_frame(
        beta_wide,
        label_to_gsm,
        ages,
        gsms,
        restrict_to_cpgs,
    )
    out.index.name = "sample_id"
    return out


def _parse_restrict_cpgs_file(path: Path) -> list[str]:
    lines = path.read_text(encoding="utf-8").splitlines()
    return [ln.strip() for ln in lines if ln.strip() and not ln.strip().startswith("#")]


def main(
    output: Annotated[
        Path,
        typer.Option(
            "--output",
            "-o",
            help="Output Parquet path (parent directories are created).",
        ),
    ],
    local_path: Annotated[
        Path | None,
        typer.Option(
            "--local-path",
            help=f"Local series matrix ({SERIES_MATRIX_FILENAME} or .txt) instead of downloading.",
        ),
    ] = None,
    geo_cache_dir: Annotated[
        Path,
        typer.Option("--geo-cache-dir", help="Cache directory for GEO / HTTPS downloads."),
    ] = Path("./data/geo"),
    restrict_cpgs_file: Annotated[
        Path | None,
        typer.Option(
            "--restrict-cpgs-file",
            help="Optional text file: one CpG ID per line (lines starting with # ignored).",
        ),
    ] = None,
) -> None:
    """CLI entrypoint: load GSE87571 and write ``--output`` Parquet."""
    restrict: list[str] | None = None
    if restrict_cpgs_file is not None:
        if not restrict_cpgs_file.is_file():
            raise FileNotFoundError(f"restrict_cpgs_file not found: {restrict_cpgs_file}")
        restrict = _parse_restrict_cpgs_file(restrict_cpgs_file)
    df = load_gse87571(
        local_path=local_path,
        geo_cache_dir=geo_cache_dir,
        restrict_to_cpgs=restrict,
    )
    save_as_parquet(df, output)
    typer.echo(f"Wrote {len(df)} samples x {df.shape[1]} columns to {output.resolve()}")


if __name__ == "__main__":
    typer.run(main)
