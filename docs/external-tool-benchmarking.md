# External Tool Benchmarking Notes

This document records how Milestone was exercised against external tools and what was observed in the current development environment.

## Tools Considered

- `chewBBACA`
- `MentaLiST`

## chewBBACA

### Environment Used

- Python 3.11 virtual environment at `/tmp/milestone-bench-venv311`
- `chewBBACA 3.5.2`
- BLAST installed in `/tmp/milestone-blast/bin`

### Fixture Workflow

The repository now includes a reproducible smoke script:

```bash
bash benchmarking/run_chewbbaca_fixture_smoke.sh
```

What it does:

1. builds a minimal external schema from [abc.fasta](/Users/fatmakhv/Desktop/milestone/tests/fixtures/e2e_schema/abc.fasta)
2. adapts that schema with `PrepExternalSchema`
3. runs `AlleleCall` with `--cds-input` on two synthetic CDS FASTA files

### Observed Result

`PrepExternalSchema` completed successfully.

`AlleleCall` progressed through exact matching and result generation, but `chewBBACA 3.5.2` terminated with:

`TypeError: 'NoneType' object is not subscriptable`

in its internal `invalid_alleles` handling during wrap-up.

This is an external-tool failure, not a Milestone exception path.

### Scientific Interpretation

At this point, the repository supports:

- executable comparison harnesses
- real `chewBBACA` fixture invocation
- reproducible failure capture

It does **not** yet support a fully successful published `chewBBACA` allele-calling comparison in this environment.

## MentaLiST

### Installation Reality

`MentaLiST` was not available via `pip` or `bioconda` in the tested environment.

The repo now includes a Docker-based smoke script:

```bash
bash benchmarking/run_mentalist_docker_smoke.sh
```

This follows the official project direction more closely than the failed `pip` path.

### Observed Result In This Environment

Docker could not be used because the daemon was not running:

- error: failed to connect to Docker socket

So the repository currently supports:

- a documented `MentaLiST` smoke path
- benchmark harness integration points

but not a completed live `MentaLiST` benchmark in this environment.

## Honest Summary

What is now true:

- Milestone has real external benchmark scaffolding.
- `chewBBACA` was installed and exercised on a real minimal workload.
- `MentaLiST` has a documented reproducible smoke path.

What is not yet true:

- there is not yet a successful full external-tool comparison for both tools on the same frozen dataset in this environment
- `MentaLiST` has not been run successfully here
- `chewBBACA` finished schema adaptation but failed during allele-call wrap-up with its current version/environment combination

## Containerized Smoke Path

The repository now includes:

- [Dockerfile.benchmark](/Users/fatmakhv/Desktop/milestone/Dockerfile.benchmark)
- [run_external_tool_smoke_suite.sh](/Users/fatmakhv/Desktop/milestone/benchmarking/run_external_tool_smoke_suite.sh)
- [external_tool_smoke_suite.json](/Users/fatmakhv/Desktop/milestone/benchmarking/external_tool_smoke_suite.json)
- [external-benchmark.yml](/Users/fatmakhv/Desktop/milestone/.github/workflows/external-benchmark.yml)

This gives a repeatable Linux-oriented path for:

1. installing BLAST
2. installing `chewBBACA`
3. running the real `chewBBACA` fixture smoke
4. attempting the `MentaLiST` Docker smoke
5. generating a unified external benchmark report
