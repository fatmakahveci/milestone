#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
workdir="${1:-/tmp/milestone-full-smoke}"

rm -rf "${workdir}"
mkdir -p "${workdir}"

schema_dir="${repo_root}/tests/fixtures/e2e_schema"
read1="${repo_root}/tests/fixtures/e2e_reads/sample_1.fastq"
read2="${repo_root}/tests/fixtures/e2e_reads/sample_2.fastq"
schema_output="${workdir}/schema"
allele_output="${workdir}/alleles"

python "${repo_root}/workflow/milestone.py" schema_creation \
  --reference smoke_ref \
  --schema_name "${schema_dir}" \
  --output "${schema_output}" \
  --threads 1

python "${repo_root}/workflow/milestone.py" allele_calling \
  --reference smoke_ref \
  --schema_name "${schema_dir}" \
  --output "${allele_output}" \
  --read1 "${read1}" \
  --read2 "${read2}" \
  --aligner vg \
  --threads 1 \
  --min-locus-coverage 10

test -f "${schema_output}/smoke_ref.fasta"
test -f "${schema_output}/smoke_ref.vcf.gz"
test -f "${allele_output}/vg/sample_wgmlst.tsv"
