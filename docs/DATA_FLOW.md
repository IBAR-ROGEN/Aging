Design document. Not yet reviewed by the institutional data protection officer.

# Data flow for sequencing intake

No Romanian sequencing data has arrived. Intake is exercised on small synthetic files that the tests create. This note says where each kind of file is supposed to live, who can see it, and what is allowed to move on.

The pseudonym key is not stored in this repository. `rogen-intake make-key` will not write it here, and `rogen-intake run` will not read a key file that sits inside the checkout.

## Custodian zone

The custodian is the person who receives the raw delivery and holds the link between a participant and a pseudonym. Analysts do not have this zone.

Files that live here:

- the raw delivery, with the original file names and the original sample identifiers (FASTQ, BAM, CRAM, VCF, POD5)
- the delivery manifest (original id, file name, file type, checksum, and an optional `aliases` column)
- the pseudonym key file, mode 600, path given only by `ROGEN_PSEUDO_KEY_FILE`
- `linkage.tsv`, original sample id next to the `RO-` pseudonym

What may leave: nothing that still has an original sample identifier, and not the key or the linkage table. The custodian runs intake on this machine. Intake checks checksums first and stops, without writing the processing directory, if a file is missing or the checksum does not match.

## Processing zone

People who work on pseudonymized sequences have this zone. They do not have the key or the linkage table.

Files that live here:

- renamed files, `<pseudonym>.<original extension>`
- `intake_report.json`: pseudonyms, file types, checksums that were verified, leak-check result, tool versions, timestamp, and git commit

The report is checked before it is published. If an original sample id appears in it, intake stops and the processing directory is not kept.

What may leave: those pseudonymized files and the intake report, and only after the leak check has passed. The check looks at output file names, BAM/CRAM/VCF headers, and the first 10,000 FASTQ read headers. It searches the manifest `original_sample_id` and every comma-separated alias on that row (for example `HG002` with aliases `NA24385` and `GM24385`). Header rewrite replaces those strings, including an `@RG SM` value that is an alias rather than the primary id. A hit deletes the temporary output and writes nothing. POD5 files are renamed after the checksum check. If the `pod5` package is not installed, the report says `pod5 metadata not scanned` and does not describe that metadata as clean.

## Release zone

The release zone is for sharing results outside the analysis group.

Files that live here: aggregate tables only. They are allowed in only after `rogen-release-check` passes. That guard is not part of intake. See [Release guard](#release-guard).

What may leave: those aggregate tables. Individual-level sequence files, the linkage table, and the key do not enter this zone.

## Release guard

`rogen-release-check` runs on every aggregate table before it leaves the processing zone. The minimum cell count is `release_guard.min_cell_count` in `config/default.yaml`. The project default is 5, to be confirmed with the data custodian.

```bash
uv run rogen-release-check \
  --table results/allele_counts.csv \
  --table-kind allele_counts \
  --n-samples 200 \
  --report results/release_report.json \
  --write-sanitized results/allele_counts_sanitized.csv
```

The command writes a JSON report with PASS or FAIL and a reason for each check.

1. Individual-level data. Fails if a column is named like `sample_id`, `eid`, `participant`, `pseudonym`, `iid`, or `fid`; if any cell matches a pseudonym (`RO-` plus 12 hex characters); or if the row count equals `--n-samples`.
2. Direct identifiers. Fails if a column looks like a date of birth, a full date, a postcode, or a name.
3. Allele-count tables (`AC` and `AN`, or `allele_count` and `allele_number`). Drops a row when `AC < k` or `AN - AC < k`, and rounds allele frequencies (`AF` or `allele_frequency`) to 3 decimals.
4. Strata tables. Suppresses phenotype-group counts below `k`, then suppresses one more count in that row when the row total would recover the hidden value.

A sanitized copy is written only when `--write-sanitized` is set and checks 1 and 2 pass. Those two checks are never rewritten. If they fail, release stops.

### Limitations

This guard does not protect against differencing attacks across several released tables, or against membership inference from summary statistics. These are open limitations.
