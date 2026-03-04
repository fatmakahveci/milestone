#!/usr/bin/env bash
set -euo pipefail

TMP_ROOT="${TMPDIR:-/tmp}/milestone-live-public-smoke"
rm -rf "$TMP_ROOT"
mkdir -p "$TMP_ROOT"

SCHEME_DB="pubmlst_neisseria_seqdef"
ISOLATE_DB="pubmlst_neisseria_isolates"

SCHEME_ID="$(
python workflow/scripts/import_pubmlst_scheme.py \
  --database "$SCHEME_DB" \
  --list-schemes \
  --scheme-type wgmlst | head -n 1 | cut -f1
)"

test -n "$SCHEME_ID"

python workflow/scripts/import_pubmlst_scheme.py \
  --database "$SCHEME_DB" \
  --scheme-id "$SCHEME_ID" \
  --output-dir "$TMP_ROOT/scheme"

test -f "$TMP_ROOT/scheme/scheme_metadata.json"
test -f "$TMP_ROOT/scheme/scheme_manifest.json"
find "$TMP_ROOT/scheme" -name '*.fasta' | head -n 1 | grep -q .

python workflow/scripts/import_pubmlst_benchmark_pack.py \
  --scheme-database "$SCHEME_DB" \
  --isolate-database "$ISOLATE_DB" \
  --scheme-id "$SCHEME_ID" \
  --discover-isolates 1 \
  --species "Neisseria spp." \
  --output-dir "$TMP_ROOT/benchmark_pack"

test -f "$TMP_ROOT/benchmark_pack/benchmark_manifest.json"
find "$TMP_ROOT/benchmark_pack/truth" -name '*_wgmlst.tsv' | head -n 1 | grep -q .

python workflow/scripts/import_enterobase_scheme.py \
  --database senterica \
  --list-schemes > "$TMP_ROOT/enterobase_schemes.tsv"

test -s "$TMP_ROOT/enterobase_schemes.tsv"

ENTEROBASE_SCHEME="$(
python - <<'PY' "$TMP_ROOT/enterobase_schemes.tsv"
from pathlib import Path
import sys
lines = Path(sys.argv[1]).read_text(encoding="utf-8").splitlines()
for line in lines:
    parts = line.split("\t")
    if len(parts) >= 3 and parts[0] and parts[2]:
        print(parts[0])
        break
PY
)"

test -n "$ENTEROBASE_SCHEME"

python workflow/scripts/import_enterobase_scheme.py \
  --database senterica \
  --scheme-name "$ENTEROBASE_SCHEME" \
  --output-dir "$TMP_ROOT/enterobase"

test -f "$TMP_ROOT/enterobase/enterobase_scheme_metadata.json"
find "$TMP_ROOT/enterobase" -name '*_profiles.tar.gz' | head -n 1 | grep -q .

echo "Live public import smoke completed successfully."
