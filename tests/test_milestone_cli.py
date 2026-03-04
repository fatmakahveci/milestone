from __future__ import annotations

import importlib.util
import sys
from argparse import Namespace
from pathlib import Path


MODULE_PATH = Path(__file__).resolve().parent.parent / "workflow" / "milestone.py"
SPEC = importlib.util.spec_from_file_location("milestone_cli", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def test_build_snakefile_text_targets_wgmlst_output() -> None:
    args = Namespace(
        command="allele_calling",
        read1="samples/demo_1.fastq.gz",
        read2="samples/demo_2.fastq.gz",
        aligner="vg",
        threads=8,
        reference="demo_ref",
        schema_name="schema",
        output="output/demo",
        update_reference=False,
        min_locus_coverage=60.0,
        snakefile="Snakefile",
    )
    paths = MODULE.workflow_paths(args)

    snakefile_text = MODULE.build_snakefile_text(args, paths)

    assert "demo_wgmlst.tsv" in snakefile_text
    assert "allele_calling.smk" in snakefile_text


def test_build_profile_compare_command_invokes_compare_script(monkeypatch) -> None:
    args = Namespace(
        command="profile_compare",
        profile_a="out/a.tsv",
        profile_b="out/b.tsv",
        output="out/compare",
        label_a="A",
        label_b="B",
    )
    calls = []

    class FakeCompleted:
        returncode = 0

    def fake_run(command, check=False, cwd=None):
        calls.append({"command": command, "check": check, "cwd": cwd})
        return FakeCompleted()

    monkeypatch.setattr(MODULE.subprocess, "run", fake_run)

    returncode = MODULE.run_profile_compare(args)

    assert returncode == 0
    assert calls[0]["command"][1].endswith("compare_wgmlst_profiles.py")
    assert "--profile-a" in calls[0]["command"]
    assert calls[0]["cwd"] == MODULE.REPO_ROOT


def test_build_profile_matrix_command_invokes_matrix_script(monkeypatch) -> None:
    args = Namespace(
        command="profile_matrix",
        profiles=["out/a.tsv", "out/b.tsv", "out/c.tsv"],
        output="out/matrix",
        distance_threshold=0.1,
    )
    calls = []

    class FakeCompleted:
        returncode = 0

    def fake_run(command, check=False, cwd=None):
        calls.append({"command": command, "check": check, "cwd": cwd})
        return FakeCompleted()

    monkeypatch.setattr(MODULE.subprocess, "run", fake_run)

    returncode = MODULE.run_profile_matrix(args)

    assert returncode == 0
    assert calls[0]["command"][1].endswith("compare_wgmlst_matrix.py")
    assert "--profiles" in calls[0]["command"]
    assert "--distance-threshold" in calls[0]["command"]
    assert calls[0]["cwd"] == MODULE.REPO_ROOT


def test_build_profile_compare_batch_command_invokes_batch_script(monkeypatch) -> None:
    args = Namespace(
        command="profile_compare_batch",
        profiles=["out/a.tsv", "out/b.tsv"],
        output="out/batch",
    )
    calls = []

    class FakeCompleted:
        returncode = 0

    def fake_run(command, check=False, cwd=None):
        calls.append({"command": command, "check": check, "cwd": cwd})
        return FakeCompleted()

    monkeypatch.setattr(MODULE.subprocess, "run", fake_run)

    returncode = MODULE.run_profile_compare_batch(args)

    assert returncode == 0
    assert calls[0]["command"][1].endswith("compare_wgmlst_batch.py")


def test_build_profile_benchmark_command_invokes_benchmark_script(monkeypatch) -> None:
    args = Namespace(
        command="profile_benchmark",
        predicted_dir="predicted",
        truth_dir="truth",
        output="out/benchmark",
    )
    calls = []

    class FakeCompleted:
        returncode = 0

    def fake_run(command, check=False, cwd=None):
        calls.append({"command": command, "check": check, "cwd": cwd})
        return FakeCompleted()

    monkeypatch.setattr(MODULE.subprocess, "run", fake_run)

    returncode = MODULE.run_profile_benchmark(args)

    assert returncode == 0
    assert calls[0]["command"][1].endswith("benchmark_wgmlst_profiles.py")
    assert "--predicted-dir" in calls[0]["command"]
    assert "--truth-dir" in calls[0]["command"]


def test_build_publication_readiness_command_invokes_script(monkeypatch) -> None:
    args = Namespace(
        command="publication_readiness",
        schema_qc_summary="out/schema_qc_summary.json",
        benchmark_summary="out/wgmlst_benchmark_summary.json",
        schema_manifest="out/schema_manifest.json",
        output="out/readiness",
        max_non_comparable_rate=0.1,
    )
    calls = []

    class FakeCompleted:
        returncode = 0

    def fake_run(command, check=False, cwd=None):
        calls.append({"command": command, "check": check, "cwd": cwd})
        return FakeCompleted()

    monkeypatch.setattr(MODULE.subprocess, "run", fake_run)

    returncode = MODULE.run_publication_readiness(args)

    assert returncode == 0
    assert calls[0]["command"][1].endswith("publication_readiness.py")
    assert "--schema-qc-summary" in calls[0]["command"]
    assert "--schema-manifest" in calls[0]["command"]


def test_build_publication_package_command_invokes_script(monkeypatch) -> None:
    args = Namespace(
        command="publication_package",
        schema_qc_summary="out/schema_qc_summary.json",
        benchmark_summary="out/wgmlst_benchmark_summary.json",
        schema_manifest="out/schema_manifest.json",
        output="out/package.zip",
        title="Demo package",
        notes="demo",
        max_non_comparable_rate=0.05,
    )
    calls = []

    class FakeCompleted:
        returncode = 0

    def fake_run(command, check=False, cwd=None):
        calls.append({"command": command, "check": check, "cwd": cwd})
        return FakeCompleted()

    monkeypatch.setattr(MODULE.subprocess, "run", fake_run)

    returncode = MODULE.run_publication_package(args)

    assert returncode == 0
    assert calls[0]["command"][1].endswith("build_publication_package.py")
    assert "--title" in calls[0]["command"]
    assert "--notes" in calls[0]["command"]


def test_build_manuscript_supplement_command_invokes_script(monkeypatch) -> None:
    args = Namespace(
        command="manuscript_supplement",
        publication_package="out/package.zip",
        output="out/supplement",
    )
    calls = []

    class FakeCompleted:
        returncode = 0

    def fake_run(command, check=False, cwd=None):
        calls.append({"command": command, "check": check, "cwd": cwd})
        return FakeCompleted()

    monkeypatch.setattr(MODULE.subprocess, "run", fake_run)

    returncode = MODULE.run_manuscript_supplement(args)

    assert returncode == 0
    assert calls[0]["command"][1].endswith("build_manuscript_supplement.py")
    assert "--publication-package" in calls[0]["command"]


def test_build_validation_corpus_command_invokes_script(monkeypatch) -> None:
    args = Namespace(
        command="validation_corpus",
        collection_dir="tests/fixtures/benchmark/public_species_packs",
        output="out/corpus",
        zip_name="corpus.zip",
    )
    calls = []

    class FakeCompleted:
        returncode = 0

    def fake_run(command, check=False, cwd=None):
        calls.append({"command": command, "check": check, "cwd": cwd})
        return FakeCompleted()

    monkeypatch.setattr(MODULE.subprocess, "run", fake_run)

    returncode = MODULE.run_validation_corpus(args)

    assert returncode == 0
    assert calls[0]["command"][1].endswith("build_species_validation_corpus.py")
    assert "--collection-dir" in calls[0]["command"]


def test_build_snakemake_command_includes_selected_flags() -> None:
    args = Namespace(
        command="schema_creation",
        threads=4,
        forceall=True,
        dryrun=True,
        printshellcmds=True,
        ri=True,
        unlock=False,
        quiet=True,
        snakefile="Snakefile",
        reference="demo_ref",
        schema_name="schema",
        output="out",
        min_locus_coverage=72.5,
    )
    paths = MODULE.workflow_paths(args)

    command = MODULE.build_snakemake_command(args, paths)

    assert "--forceall" in command
    assert "--dryrun" in command
    assert "--printshellcmds" in command
    assert "--rerun-incomplete" in command
    assert "--quiet" in command


def test_build_config_text_records_min_locus_coverage() -> None:
    args = Namespace(
        command="schema_creation",
        threads=4,
        forceall=False,
        dryrun=False,
        printshellcmds=False,
        ri=False,
        unlock=False,
        quiet=False,
        snakefile="Snakefile",
        reference="demo_ref",
        schema_name="schema",
        output="out",
        min_locus_coverage=72.5,
        translation_table=11,
        allowed_start_codons="ATG,GTG",
    )
    paths = MODULE.workflow_paths(args)

    config_text = MODULE.build_config_text(args, paths)

    assert "min_locus_coverage: 72.5" in config_text
    assert 'allowed_start_codons: "ATG,GTG"' in config_text
