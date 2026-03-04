#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
output_dir="${1:-/tmp/milestone-external-smoke}"

mkdir -p "${output_dir}"

if [[ -n "${CHEWIE_BIN:-}" && -n "${BLAST_BIN_DIR:-}" ]]; then
  echo "Running chewBBACA fixture smoke..."
  bash "${repo_root}/benchmarking/run_chewbbaca_fixture_smoke.sh" "${output_dir}/chewbbaca_fixture" || true
else
  echo "Skipping chewBBACA fixture smoke because CHEWIE_BIN or BLAST_BIN_DIR is unset."
fi

if command -v docker >/dev/null 2>&1; then
  echo "Running MentaLiST Docker smoke..."
  bash "${repo_root}/benchmarking/run_mentalist_docker_smoke.sh" || true
else
  echo "Skipping MentaLiST Docker smoke because docker is unavailable."
fi

python3 "${repo_root}/workflow/scripts/benchmark_typing_tools.py" \
  --config "${repo_root}/benchmarking/external_tool_smoke_suite.json" \
  --output-dir "${output_dir}/reports"
