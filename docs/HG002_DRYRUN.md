# HG002 short-read dry run

Public GIAB sample **HG002** (= NA24385 = GM24385). Pipeline:
download → `rogen-intake` → align → bcftools call → `rtg vcfeval` vs NISTv4.2.1.

bcftools is a baseline caller. Numbers below are **not** tuned for F1.

## Window

| Item | Value |
|------|-------|
| Region | `chr20:10000000-20000000` |
| Downsample | samtools `-s 0.1` (seed 0, ~30x target) |
| Aligner | `minimap2` |
| Mean depth | 28.761 |
| % mapped | 99.390 |
| % duplicates | 0.255 |

Reads were taken from the chr20 window of a pre-aligned GIAB BAM and realigned to chr20 only. The mapping rate is therefore close to 100% by construction and is not a whole-genome mapping rate. Precision and recall apply to this window, not to the whole genome.

## Benchmark (GIAB v4.2.1, BED ∩ region)

### SNVs

| TP | FP | FN | Precision | Recall | F1 |
|----|----|----|-----------|--------|----|
| 13214 | 71 | 93 | 0.995 | 0.993 | 0.994 |

### Indels

| TP | FP | FN | Precision | Recall | F1 |
|----|----|----|-----------|--------|----|
| 1898 | 182 | 172 | 0.916 | 0.917 | 0.916 |

## Tool versions

```text
{
  "samtools": "samtools 1.24",
  "bcftools": "bcftools 1.24",
  "minimap2": "2.31-r1302",
  "rtg": "Product: RTG Tools 3.13"
}
```

## Wall time (seconds)

```text
{
  "fetch": 3474.107,
  "intake": 3.169,
  "align": 25.96,
  "call": 30.933,
  "benchmark": 5.522
}
```

## Reproduce

```bash
uv run python scripts/dryrun/hg002_shortread_dryrun.py check-urls
uv run python scripts/dryrun/hg002_shortread_dryrun.py fetch
uv run python scripts/dryrun/hg002_shortread_dryrun.py intake
uv run python scripts/dryrun/hg002_shortread_dryrun.py align
uv run python scripts/dryrun/hg002_shortread_dryrun.py call
uv run python scripts/dryrun/hg002_shortread_dryrun.py benchmark
uv run python scripts/dryrun/hg002_shortread_dryrun.py report
```

Metrics JSON: [`analysis/dryrun/hg002_shortread_metrics.json`](../analysis/dryrun/hg002_shortread_metrics.json).
