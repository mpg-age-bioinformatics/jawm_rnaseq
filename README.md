# jawm_rnaseq

This is a jawm rnaseq module.

Installing jawm:
```
pip install git+ssh://git@github.com/mpg-age-bioinformatics/jawm.git
```
For more information on jawm please visit jawm's repo on [GitHub.com](https://github.com/mpg-age-bioinformatics/jawm/tree/main).

Example usage:
```
# clone this module
git clone git@github.com:mpg-age-bioinformatics/jawm_rnaseq.git

cd jawm_rnaseq

# download test data
jawm-test -r download

# docker
jawm rnaseq.py rnaseq -p ./yaml/docker.yaml

# slurm & apptainer with multiple yaml files
jawm rnaseq.py rnaseq -p ./yaml/vars.yaml ./yaml/hpc.yaml
```

Additional jawm workflows are available [here (GitHub.com)](https://github.com/mpg-age-bioinformatics?q=jawm_&type=all&language=&sort=).

## Raw-file checksums

The `raw_md5sums` Process uses `python:3.12-slim` to hash a list of files and
write one `<file>\t<md5sum>` row per file, without a header. It hashes the exact
file bytes, including compression for `.gz` files.

```python
checksums = _rnaseq.configure_raw_md5sums(files, "/path/to/checksums.tsv")
checksums.execute()
```

For local inputs, the RNA-seq workflow hashes the original filenames in the
first two columns of `groups`, resolved against `source_folder`. Each lane is
recorded separately, using its original path rather than a renamed or
concatenated read file. The job runs after raw-file listing, before FastQC, and writes
`metadata/raw-inputs.md5.tsv`. The file is included under `raw_inputs` in the
report file list.

For existing local inputs, set
`--global.var.md5sums_file /path/to/recorded-checksums.txt` to verify those
files against existing checksums by **basename**, ignoring recorded directories.
The reference may use standard `md5sum` format (`<md5sum>  <file>`) or the TSV
format above. Missing entries, mismatches, or conflicting entries for the same
basename fail the job before FastQC. Entries are compared against the original
raw filenames, including individual lanes before concatenation.

For service-provider downloads, `download_link` must also download the provider's
checksum file specified by `md5sums_file` (relative to `source_folder` or absolute).
The download step extracts its FASTQ entries to `fastq.md5sums.txt` (or
`fastq.md5sums.1.txt` if that file already exists), matches the flattened file
basenames, and checks them before preprocessing. Successful verification is
recorded in `fastq.md5sums.check.log`. The later checksum job records
the original raw files selected in `groups`.

Without an external file, existing checksum reports are skipped on resume.
With an external file, checksums are always verified. Missing inputs fail the
job; the report is published only after all files have been hashed and verified.
SRA downloads skip this checksum step because the SRA tools already verify them.
