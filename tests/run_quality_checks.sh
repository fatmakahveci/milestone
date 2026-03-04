#!/usr/bin/env bash
set -euo pipefail

python -m py_compile \
  workflow/milestone.py \
  workflow/scripts/compare_wgmlst_profiles.py \
  workflow/scripts/compare_wgmlst_matrix.py \
  workflow/scripts/compare_wgmlst_batch.py \
  workflow/scripts/benchmark_wgmlst_profiles.py \
  workflow/scripts/benchmark_typing_tools.py \
  workflow/scripts/build_validation_bundle.py \
  workflow/scripts/build_publication_package.py \
  workflow/scripts/freeze_public_benchmark_pack.py \
  workflow/scripts/import_pubmlst_benchmark_pack.py \
  workflow/scripts/import_enterobase_scheme.py \
  workflow/scripts/publication_readiness.py \
  workflow/scripts/compare_results.py \
  workflow/scripts/comp_milestone_chewie.py \
  workflow/scripts/create_reference.py \
  workflow/scripts/create_sample_info_txt.py \
  workflow/scripts/import_pubmlst_scheme.py \
  workflow/scripts/result_metrics.py \
  workflow/scripts/schema_manifest.py \
  workflow/scripts/schema_qc.py \
  workflow/scripts/wgmlst_utils.py \
  webapp/app.py \
  webapp/public_schemas.py \
  webapp/runner.py \
  webapp/job_worker.py \
  webapp/uploads.py

python -m ruff check workflow webapp tests
python -m mypy webapp workflow/scripts/result_metrics.py workflow/scripts/schema_qc.py workflow/scripts/wgmlst_utils.py workflow/milestone.py
