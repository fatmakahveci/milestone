#!/usr/bin/env bash
set -euo pipefail

CHEWIE_BIN="${CHEWIE_BIN:-/tmp/milestone-bench-venv311/bin/chewBBACA.py}"
BLAST_BIN_DIR="${BLAST_BIN_DIR:-/tmp/milestone-blast/bin}"
WORKDIR="${1:-/tmp/chewie_real_smoke}"
export WORKDIR

rm -rf "${WORKDIR}"
mkdir -p "${WORKDIR}/schema_raw" "${WORKDIR}/input_cds"
cp tests/fixtures/e2e_schema/abc.fasta "${WORKDIR}/schema_raw/abc.fasta"

python3 - <<'PY'
from pathlib import Path
import os

workdir = Path(os.environ["WORKDIR"])
seq1 = "ATGAAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTTAAATAA"
seq2 = "ATGAAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTTAAACCCGCGTTTAAATAA"
(workdir / "input_cds" / "sample_a.fasta").write_text(f">abc\n{seq1}\n", encoding="ascii")
(workdir / "input_cds" / "sample_b.fasta").write_text(f">abc\n{seq2}\n", encoding="ascii")
PY

PATH="${BLAST_BIN_DIR}:$PATH" "${CHEWIE_BIN}" PrepExternalSchema \
  --schema-directory "${WORKDIR}/schema_raw" \
  --output-directory "${WORKDIR}/schema_chew" \
  --blast-path "${BLAST_BIN_DIR}" \
  --cpu 1

PATH="${BLAST_BIN_DIR}:$PATH" "${CHEWIE_BIN}" AlleleCall \
  --input-files "${WORKDIR}/input_cds" \
  --schema-directory "${WORKDIR}/schema_chew" \
  --output-directory "${WORKDIR}/allele_call" \
  --blast-path "${BLAST_BIN_DIR}" \
  --cpu 1 \
  --cds-input
