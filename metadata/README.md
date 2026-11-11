# metadata

Small text files that describe the data. They are committed to Git.

- `samples.tsv`: one row per sample (ID, source accession, condition, file names).
- `sources.tsv`: where each raw file is downloaded from. Columns: file path, NASA OSDR download URL, size in bytes, and download date.
- `checksums.sha256`: checksums of the raw files you keep, so anyone can check. SHA-256 of each downloaded file. Written by the first run of `pixi run download`
  they have the same data (`shasum -a 256 -c metadata/checksums.sha256`).
- `README.md` (this file): a data dictionary that explains the files above.

Samples are labeled from the sample names in the counts table: '_GC_' means ground control (10 samples), and '_FLT_' means spaceflight (9 samples).

The KEGG gene sets are not checksummed because KEGG updates ;; their download date is saved in the file and shown in the report.
