<div align="left">
  <h1>
    <img src="images/milestone.png" alt="milestone logo" width="96" />
    MILESTONE
  </h1>
</div>

[![Tests](https://github.com/fatmakahveci/milestone/actions/workflows/tests.yml/badge.svg)](https://github.com/fatmakahveci/milestone/actions/workflows/tests.yml)

Milestone is a local, reproducible workflow for bacterial `wgMLST` analysis:

- schema creation or public schema import
- allele calling
- profile comparison (`pairwise`, `batch`, `matrix`)
- schema QC and benchmark workflows
- publication-ready packaging and reporting
- Streamlit web UI for interactive runs

Live static site: https://fatmakahveci.github.io/milestone/

## Quick Start

### CLI setup

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r webapp/requirements.txt
```

### Run web app

```bash
streamlit run webapp/app.py
```

### Run tests

```bash
pytest tests
```

## Core Commands

### 1) Compare two profiles

```bash
python workflow/milestone.py profile_compare \
  --profile-a output/run_a/sample_a_wgmlst.tsv \
  --profile-b output/run_b/sample_b_wgmlst.tsv \
  --output-dir output/compare/sample_a_vs_b
```

Outputs:
- `wgmlst_profile_summary.json`
- `wgmlst_profile_comparison.tsv`

### 2) Build profile distance matrix

```bash
python workflow/milestone.py profile_matrix \
  --profiles output/run_a/sample_a_wgmlst.tsv output/run_b/sample_b_wgmlst.tsv output/run_c/sample_c_wgmlst.tsv \
  --output output/compare/matrix_abc
```

Outputs:
- `wgmlst_distance_matrix.tsv`
- `wgmlst_distance_summary.json`
- `wgmlst_distance_report.html`

### 3) Benchmark against truth profiles

```bash
python workflow/milestone.py profile_benchmark \
  --predicted-dir output/predicted_profiles \
  --truth-dir validation/truth_profiles \
  --output output/benchmark_run
```

Outputs:
- `wgmlst_benchmark_summary.json`
- `wgmlst_benchmark_per_sample.tsv`
- `wgmlst_benchmark_per_locus.tsv`
- `wgmlst_benchmark_report.html`

### 4) Schema QC

```bash
python workflow/scripts/schema_qc.py \
  --schema-dir tests/fixtures/smoke_schema \
  --output-dir output/schema_qc
```

Outputs:
- `schema_qc_summary.json`
- `schema_qc_issues.tsv`

## Public Schema / Benchmark Import

### PubMLST scheme import

```bash
python workflow/scripts/import_pubmlst_scheme.py \
  --database pubmlst_neisseria_seqdef \
  --scheme-id 1 \
  --output-dir schemas/neisseria_wgmlst \
  --include-profiles
```

### PubMLST benchmark pack import

```bash
python workflow/scripts/import_pubmlst_benchmark_pack.py \
  --scheme-database pubmlst_neisseria_seqdef \
  --isolate-database pubmlst_neisseria_isolates \
  --scheme-id 1 \
  --isolate-id 101 \
  --isolate-id 102 \
  --species "Neisseria meningitidis" \
  --output-dir benchmark_packs/neisseria_demo
```

### EnteroBase scheme import

```bash
python workflow/scripts/import_enterobase_scheme.py \
  --database senterica \
  --scheme-name wgMLST \
  --output-dir public_schemes/enterobase_senterica
```

If EnteroBase requires authentication, set one of:

```bash
export ENTEROBASE_TOKEN="..."
```

or

```bash
export ENTEROBASE_USERNAME="..."
export ENTEROBASE_PASSWORD="..."
```

## Publication / Validation Utilities

```bash
python workflow/milestone.py publication_readiness \
  --schema-qc-summary results/schema_qc/schema_qc_summary.json \
  --benchmark-summary results/benchmark/wgmlst_benchmark_summary.json \
  --schema-manifest results/schema/schema_manifest.json \
  --output results/publication_readiness
```

```bash
python workflow/milestone.py publication_package \
  --schema-qc-summary results/schema_qc/schema_qc_summary.json \
  --benchmark-summary results/benchmark/wgmlst_benchmark_summary.json \
  --schema-manifest results/schema/schema_manifest.json \
  --output results/milestone_publication_package.zip \
  --title "Milestone validation package"
```

## Docker

### Full stack

```bash
docker compose up --build
```

### UI only

```bash
docker compose up --build milestone-web
```

### Runtime image with full bioinformatics stack

```bash
docker compose up --build milestone-runtime
```

### External tool smoke benchmarking

```bash
docker compose run --rm milestone-benchmark
```

## Project Documentation

- [What Milestone Is](docs/scientific-positioning.md)
- [How To Use Milestone](docs/how-to-use-milestone.md)
- [Competitive Analysis](docs/competitive-analysis.md)
- [Public Benchmark Datasets](docs/public-benchmark-datasets.md)
- [Publication Checklist](docs/publication-checklist.md)
- [Publication Runbook](docs/publication-runbook.md)
- [Manuscript Sections](docs/manuscript-sections.md)
- [External Tool Benchmarking](docs/external-tool-benchmarking.md)
- [External Tool Install Notes](docs/external-tool-install.md)

## Notes

- Milestone outputs `*_wgmlst.tsv` profile files.
- Profile comparison is a schema-based discriminator; it is not a direct phylogeny engine.
- Species-specific validation is still required before strong biological claims.

## Citation Metadata

- [CITATION.cff](CITATION.cff)
- [codemeta.json](codemeta.json)
