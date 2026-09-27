Design document. Not yet reviewed by the institutional data protection officer.

# Data flow for sequencing intake

No Romanian sequencing data has arrived. Intake is exercised on small synthetic files that the tests create. This note says where each kind of file is supposed to live, who can see it, and what is allowed to move on.

The pseudonym key is not stored in this repository. `rogen-intake make-key` will not write it here, and `rogen-intake run` will not read a key file that sits inside the checkout.

## Custodian zone

The custodian is the person who receives the raw delivery and holds the link between a participant and a pseudonym. Analysts do not have this zone.

Files that live here:

- the raw delivery, with the original file names and the original sample identifiers (FASTQ, BAM, CRAM, VCF, POD5)
- the delivery manifest (original id, file name, file type, checksum)
- the pseudonym key file, mode 600, path given only by `ROGEN_PSEUDO_KEY_FILE`
- `linkage.tsv`, original sample id next to the `RO-` pseudonym

What may leave: nothing that still has an original sample identifier, and not the key or the linkage table. The custodian runs intake on this machine. Intake checks checksums first and stops, without writing the processing directory, if a file is missing or the checksum does not match.

## Processing zone

People who work on pseudonymized sequences have this zone. They do not have the key or the linkage table.

Files that live here:

- renamed files, `<pseudonym>.<original extension>`
- `intake_report.json`: pseudonyms, file types, checksums that were verified, leak-check result, tool versions, timestamp, and git commit

The report is checked before it is published. If an original sample id appears in it, intake stops and the processing directory is not kept.

What may leave: those pseudonymized files and the intake report, and only after the leak check has passed. The check looks at output file names, BAM/CRAM/VCF headers, and the first 10,000 FASTQ read headers. A hit deletes the temporary output and writes nothing. POD5 files are renamed after the checksum check. If the `pod5` package is not installed, the report says `pod5 metadata not scanned` and does not describe that metadata as clean.

## Release zone

The release zone is for sharing results outside the analysis group.

Files that live here: aggregate tables only. They are allowed in only after the release guard that a later step will add. That guard is not part of intake.

What may leave: those aggregate tables. Individual-level sequence files, the linkage table, and the key do not enter this zone.
