# Milestone Publication Runbook

This runbook describes a defensible end-to-end workflow for taking a Milestone analysis from raw outputs to manuscript-ready supplementary material.

## 1. Freeze The Validation Corpus

If you already have a local benchmark-pack collection:

```bash
python workflow/milestone.py validation_corpus \
  --collection-dir tests/fixtures/benchmark/public_species_packs \
  --output results/validation_corpus
```

If you want to build a new public workspace from PubMLST-backed species definitions, create a JSON config and run:

```bash
python workflow/scripts/build_public_benchmark_workspace.py \
  --config configs/public_benchmark_workspace.example.json \
  --output-dir results/public_benchmark_workspace
```

This produces:

- frozen per-species benchmark packs
- a corpus bundle
- a workspace manifest documenting the imported species packs

## 2. Run Schema QC

```bash
python workflow/scripts/schema_qc.py \
  --schema-dir path/to/schema \
  --manifest path/to/schema_manifest.json \
  --output-dir results/schema_qc
```

Before submission, resolve all schema QC errors. Any remaining warnings should be discussed explicitly in the manuscript.

## 3. Run Profile Benchmarking

```bash
python workflow/milestone.py profile_benchmark \
  --benchmark-pack-dir results/public_benchmark_workspace/packs \
  --output results/benchmark
```

Key outputs:

- `wgmlst_benchmark_summary.json`
- `wgmlst_benchmark_per_sample.tsv`
- `wgmlst_benchmark_per_locus.tsv`
- `wgmlst_benchmark_report.html`

## 4. Run Cross-Tool Benchmark Harness

For local smoke benchmarking:

```bash
python workflow/scripts/benchmark_typing_tools.py \
  --config benchmarking/tool_profiles.json \
  --output-dir results/tool_benchmarks
```

Notes:

- `chewBBACA` and `MentaLiST` results are meaningful only when those tools are installed and configured with real commands and fixed schemas.
- Tool-version metadata and skip reasons are recorded automatically.

## 5. Evaluate Publication Readiness

```bash
python workflow/milestone.py publication_readiness \
  --schema-qc-summary results/schema_qc/schema_qc_summary.json \
  --benchmark-summary results/benchmark/wgmlst_benchmark_summary.json \
  --schema-manifest path/to/schema_manifest.json \
  --output results/publication_readiness
```

Review:

- `publication_readiness.json`
- `publication_readiness.md`

Do not ignore `fail` statuses. `warn` statuses should be documented in the manuscript.

## 6. Build A Publication Package

```bash
python workflow/milestone.py publication_package \
  --schema-qc-summary results/schema_qc/schema_qc_summary.json \
  --benchmark-summary results/benchmark/wgmlst_benchmark_summary.json \
  --schema-manifest path/to/schema_manifest.json \
  --output results/milestone_publication_package.zip \
  --title "Milestone validation package"
```

This ZIP includes:

- readiness report
- benchmark and QC summaries
- methods and limitations templates
- citation metadata
- environment metadata

## 7. Extract A Supplement Directory

```bash
python workflow/milestone.py manuscript_supplement \
  --publication-package results/milestone_publication_package.zip \
  --output results/manuscript_supplement
```

This creates a supplement-ready directory with `SUPPLEMENT_INDEX.md`.

## 8. What To Claim

Safe claims:

- Milestone performs schema-based wgMLST allele calling and profile discrimination.
- The workflow was validated against a frozen benchmark corpus.
- Unresolved loci were handled explicitly as non-comparable.

Unsafe claims:

- universal strain identification
- phylogenetic distance inference from wgMLST alone
- biological confirmation of novel alleles without orthogonal validation

## 9. Minimum Submission Package

Before submission, archive at least:

- schema manifest
- schema QC summary
- benchmark summary
- publication package ZIP
- manuscript supplement directory
- exact Milestone commit or release identifier
