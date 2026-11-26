#!/usr/bin/env python3
"""Pool haplotype-split modkit bedMethyl files into one pileup (low memory).

wf-human-variation writes SAMPLE.wf_mods.1 / .2 / .ungrouped bedMethyl files.
Reads are partitioned between them, so per-site counts add up. This script
streams the sorted inputs, sums the count columns per (chrom, start, code),
recomputes the percentage from the summed counts, and writes one gzip file in
the same 18-column modkit layout. Percentages are never averaged.

Only codes m and h are kept (the clock scripts use only these). Input must be
strand-combined (strand column '.'). All inputs must list contigs in plain
text order (chr1, chr10, ..., chr9, chrM, chrX, chrY); any other order
stops the run. Contigs other than chr1-22, X, Y, M are counted
and skipped.

Usage:
  python scripts/clock/pool_wf_mods_bedmethyl.py \
      --out data/dryrun/ont_vs_bisulfite/SAMPLE.wf_mods.pooled.bedmethyl.gz \
      A.bedmethyl.gz B.bedmethyl.gz C.bedmethyl.gz
"""

from __future__ import annotations

import argparse
import datetime as dt
import gzip
import hashlib
import heapq
import json
import subprocess
import sys
from pathlib import Path

ALLOWED = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY", "chrM"}
# modkit / wf-human-variation write contigs in plain text order:
# chr1, chr10, ..., chr19, chr2, chr20, chr21, chr22, chr3, ..., chr9, chrM, chrX, chrY.
# Not every file contains every contig (the haplotype files have no chrM), so the
# order is fixed here, not learned from the files. The sort check below stops the
# run if any file uses a different order.
RANK = {c: i for i, c in enumerate(sorted(ALLOWED))}
KEEP = {"m", "h"}
# 0-based indexes of the summed count columns in the modkit layout
SUM_COLS = {
    9: "n_valid",
    11: "n_mod",
    12: "n_canonical",
    13: "n_other_mod",
    14: "n_delete",
    15: "n_fail",
    16: "n_diff",
    17: "n_nocall",
}


def _open(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_rows(path: Path, stats: dict):
    last = None
    with _open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            assert len(f) >= 18, f"STOP: expected 18 columns, got {len(f)} in {path}"
            if f[3] not in KEEP:
                stats["skipped_other_code"] += 1
                continue
            if f[0] not in ALLOWED:
                stats["skipped_other_contig"] += 1
                continue
            rank = RANK[f[0]]
            assert f[5] == ".", f"STOP: strand is {f[5]!r}, need strand-combined input: {path}"
            key = (rank, int(f[1]), f[3])
            assert last is None or key >= last, f"STOP: {path} not sorted at {key} after {last}"
            last = key
            stats["rows_in"] += 1
            yield key, f


def pool(paths: list[Path], out: Path) -> dict:
    stats = {"rows_in": 0, "rows_out": 0, "skipped_other_code": 0, "skipped_other_contig": 0}
    streams = [read_rows(p, stats) for p in paths]
    merged = heapq.merge(*streams, key=lambda kv: kv[0])
    cur_key, cur = None, None

    def flush(fh):
        key, f, sums = cur_key, cur_f, cur
        n_valid, n_mod = sums[9], sums[11]
        assert n_valid > 0, f"STOP: zero valid coverage at {key}"
        out_f = [
            f[0],
            f[1],
            f[2],
            f[3],
            str(n_valid),
            ".",
            f[1],
            f[2],
            "255,0,0",
            str(n_valid),
            f"{100.0 * n_mod / n_valid:.2f}",
            str(n_mod),
            str(sums[12]),
            str(sums[13]),
            str(sums[14]),
            str(sums[15]),
            str(sums[16]),
            str(sums[17]),
        ]
        fh.write("\t".join(out_f) + "\n")
        stats["rows_out"] += 1

    with gzip.open(out, "wt") as fh:
        cur_f = None
        for key, f in merged:
            if key != cur_key:
                if cur_key is not None:
                    flush(fh)
                cur_key, cur_f = key, f
                cur = {i: 0 for i in SUM_COLS}
            for i in SUM_COLS:
                cur[i] += int(f[i])
        if cur_key is not None:
            flush(fh)
    return stats


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("inputs", nargs="+", type=Path)
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args()
    for p in a.inputs:
        assert p.is_file(), f"STOP: missing input {p}"
    stats = pool(a.inputs, a.out)
    prov = {
        "script": "scripts/clock/pool_wf_mods_bedmethyl.py",
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "python": sys.version.split()[0],
        "git_commit": subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"], capture_output=True, text=True
        ).stdout.strip()
        or "unknown",
        "inputs": {str(p): _sha256(p) for p in a.inputs},
        "output": {str(a.out): _sha256(a.out)},
        "kept_codes": sorted(KEEP),
        "contig_order_expected": sorted(RANK, key=RANK.get),
        "method": "per-site sum of count columns across inputs; percent recomputed from summed counts",
        **stats,
    }
    Path(str(a.out) + ".provenance.json").write_text(json.dumps(prov, indent=2) + "\n")
    print(json.dumps(stats))


if __name__ == "__main__":
    main()
