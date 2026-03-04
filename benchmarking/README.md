# Benchmarking Profiles

This directory contains reproducible command profiles for running Milestone and comparison tools against fixed datasets.

## Included profiles

- `milestone_smoke.json`
- `chewbbaca_smoke.json`
- `mentalist_smoke.json`

These are not published benchmark results. They are executable harness definitions for environments where the external tools are installed.

## Example

```bash
python workflow/scripts/benchmark_typing_tools.py \
  --config benchmarking/milestone_smoke.json \
  --output-dir output/tool_benchmarks/milestone_smoke
```

Or compare multiple tools:

```bash
python workflow/scripts/benchmark_typing_tools.py \
  --config benchmarking/tool_profiles.json \
  --output-dir output/tool_benchmarks/all_tools
```

## Honest limit

The included `chewBBACA` and `MentaLiST` profiles are templates. They become meaningful only in a runtime where those tools and their input schemas are actually installed.

The benchmark harness now records:

- command availability and skip reasons
- tool version previews when available
- runtime summaries
- a comparison summary relative to the `Milestone` baseline when multiple tools complete successfully

Additional runnable helpers:

- [run_chewbbaca_fixture_smoke.sh](/Users/fatmakhv/Desktop/milestone/benchmarking/run_chewbbaca_fixture_smoke.sh): minimal real `chewBBACA` fixture workflow
- [run_mentalist_docker_smoke.sh](/Users/fatmakhv/Desktop/milestone/benchmarking/run_mentalist_docker_smoke.sh): Docker-based `MentaLiST` smoke path
- [chewbbaca_fixture_smoke.json](/Users/fatmakhv/Desktop/milestone/benchmarking/chewbbaca_fixture_smoke.json): runnable harness profile for the real `chewBBACA` fixture smoke
- [mentalist_docker_smoke.json](/Users/fatmakhv/Desktop/milestone/benchmarking/mentalist_docker_smoke.json): runnable harness profile for the Docker-based `MentaLiST` smoke

Observed external-tool behavior and current limits are documented in [external-tool-benchmarking.md](/Users/fatmakhv/Desktop/milestone/docs/external-tool-benchmarking.md).
Installation notes are documented in [external-tool-install.md](/Users/fatmakhv/Desktop/milestone/docs/external-tool-install.md).
