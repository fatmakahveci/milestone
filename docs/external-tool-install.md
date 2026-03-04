# External Tool Installation Notes

This note records practical installation paths for external comparison tools used with Milestone.

## chewBBACA

### Working Path Used In This Repository

A working local setup was achieved with:

- Python `3.11`
- `chewBBACA 3.5.2`
- BLAST installed separately

Example local setup:

```bash
/opt/homebrew/opt/python@3.11/bin/python3.11 -m venv /tmp/milestone-bench-venv311
/tmp/milestone-bench-venv311/bin/pip install setuptools wheel chewBBACA
/opt/homebrew/bin/micromamba create -y -p /tmp/milestone-blast -c conda-forge -c bioconda blast
```

Then:

```bash
CHEWIE_BIN=/tmp/milestone-bench-venv311/bin/chewBBACA.py \
BLAST_BIN_DIR=/tmp/milestone-blast/bin \
bash benchmarking/run_chewbbaca_fixture_smoke.sh
```

### Current Observed Limitation

`PrepExternalSchema` completed successfully, but `AlleleCall` failed during wrap-up with `chewBBACA 3.5.2` on the tested minimal fixture.

## MentaLiST

### What Did Not Work Here

- `pip install MentaLiST`
- `micromamba create ... mentalist`

Those routes did not provide a working install in the tested environment.

### Practical Recommended Path

Use Docker:

```bash
bash benchmarking/run_mentalist_docker_smoke.sh
```

This requires:

- a running Docker daemon
- network access to pull the image

For a reproducible Linux-style benchmark environment, the repository now also includes:

- [Dockerfile.benchmark](/Users/fatmakhv/Desktop/milestone/Dockerfile.benchmark)
- [run_external_tool_smoke_suite.sh](/Users/fatmakhv/Desktop/milestone/benchmarking/run_external_tool_smoke_suite.sh)

### Current Observed Limitation

In the tested environment, Docker was installed but the daemon/socket was not available, so the smoke run could not be completed.

## Reporting Recommendation

When publishing external-tool comparisons based on this repository, report:

- exact tool version
- exact installation route
- whether the comparison completed successfully or failed due to external-tool/runtime issues
- whether the comparison is a smoke benchmark or a real allele-calling workload
