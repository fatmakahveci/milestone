from __future__ import annotations

import importlib.util
import io
import json
import sys
from pathlib import Path

RUNNER_PATH = Path(__file__).resolve().parent.parent / "webapp" / "runner.py"
RUNNER_SPEC = importlib.util.spec_from_file_location("runner", RUNNER_PATH)
RUNNER = importlib.util.module_from_spec(RUNNER_SPEC)
assert RUNNER_SPEC.loader is not None
sys.modules[RUNNER_SPEC.name] = RUNNER
RUNNER_SPEC.loader.exec_module(RUNNER)

UPLOADS_PATH = Path(__file__).resolve().parent.parent / "webapp" / "uploads.py"
UPLOADS_SPEC = importlib.util.spec_from_file_location("uploads", UPLOADS_PATH)
UPLOADS = importlib.util.module_from_spec(UPLOADS_SPEC)
assert UPLOADS_SPEC.loader is not None
sys.modules[UPLOADS_SPEC.name] = UPLOADS
UPLOADS_SPEC.loader.exec_module(UPLOADS)


class FakeUpload:
    def __init__(self, name: str, payload: bytes) -> None:
        self.name = name
        self._payload = payload

    def getbuffer(self) -> memoryview:
        return memoryview(self._payload)


class LargeFakeUpload(FakeUpload):
    pass


def test_start_background_job_writes_metadata_and_command_files(tmp_path: Path, monkeypatch) -> None:
    popen_calls = []

    class FakePopen:
        def __init__(self, command, cwd=None, stdout=None, stderr=None, start_new_session=None):
            popen_calls.append(
                {
                    "command": command,
                    "cwd": cwd,
                    "start_new_session": start_new_session,
                }
            )

    monkeypatch.setattr(RUNNER, "JOBS_ROOT", tmp_path / "jobs")
    monkeypatch.setattr(RUNNER.subprocess, "Popen", FakePopen)

    job_id = RUNNER.start_background_job("schema_creation", ["python", "workflow/milestone.py"], ["out.txt"])
    job_dir = RUNNER.JOBS_ROOT / job_id

    metadata = json.loads((job_dir / "metadata.json").read_text(encoding="utf-8"))
    command = json.loads((job_dir / "command.json").read_text(encoding="utf-8"))

    assert metadata["job_id"] == job_id
    assert metadata["label"] == "schema_creation"
    assert metadata["status"] == "queued"
    assert metadata["outputs"] == ["out.txt"]
    assert command["command"] == ["python", "workflow/milestone.py"]
    assert popen_calls[0]["cwd"] == RUNNER.REPO_ROOT


def test_get_job_snapshot_reads_existing_metadata(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(RUNNER, "JOBS_ROOT", tmp_path / "jobs")
    job_dir = RUNNER.JOBS_ROOT / "abc123"
    job_dir.mkdir(parents=True)
    (job_dir / "stdout.log").write_text("hello", encoding="utf-8")
    (job_dir / "stderr.log").write_text("", encoding="utf-8")
    (job_dir / "command.json").write_text(json.dumps({"command": ["python", "x.py"]}), encoding="utf-8")
    (job_dir / "metadata.json").write_text(
        json.dumps(
            {
                "job_id": "abc123",
                "label": "allele_calling",
                "status": "running",
                "command": ["python", "x.py"],
                "command_text": "python x.py",
                "created_at": "2026-03-04T00:00:00+00:00",
                "started_at": None,
                "finished_at": None,
                "pid": 1234,
                "returncode": None,
                "outputs": ["demo.txt"],
            }
        ),
        encoding="utf-8",
    )

    snapshot = RUNNER.get_job_snapshot("abc123")

    assert snapshot is not None
    assert snapshot.job_id == "abc123"
    assert snapshot.status == "running"
    assert snapshot.stdout_path.read_text(encoding="utf-8") == "hello"


def test_retry_job_requeues_same_command_and_outputs(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(RUNNER, "JOBS_ROOT", tmp_path / "jobs")
    queued = {}

    def fake_start_background_job(label, command, outputs):
        queued["label"] = label
        queued["command"] = command
        queued["outputs"] = outputs
        return "retry123"

    monkeypatch.setattr(RUNNER, "start_background_job", fake_start_background_job)

    snapshot = RUNNER.JobSnapshot(
        job_id="abc123",
        label="public_schema_import",
        status="failed",
        command=["python", "workflow/scripts/import_pubmlst_scheme.py"],
        command_text="python workflow/scripts/import_pubmlst_scheme.py",
        created_at="2026-03-04T00:00:00+00:00",
        started_at=None,
        finished_at=None,
        pid=None,
        returncode=1,
        stdout_path=tmp_path / "stdout.log",
        stderr_path=tmp_path / "stderr.log",
        outputs=["ci_schema"],
        job_dir=tmp_path / "jobdir",
    )

    new_job_id = RUNNER.retry_job(snapshot)

    assert new_job_id == "retry123"
    assert queued["label"] == "public_schema_import"
    assert queued["command"] == ["python", "workflow/scripts/import_pubmlst_scheme.py"]
    assert queued["outputs"] == ["ci_schema"]


def test_create_download_bundle_includes_existing_files(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(RUNNER, "REPO_ROOT", tmp_path)
    result_file = tmp_path / "demo.txt"
    result_file.write_text("payload", encoding="utf-8")

    bundle = RUNNER.create_download_bundle(["demo.txt", "missing.txt"])

    import zipfile

    with zipfile.ZipFile(io.BytesIO(bundle)) as archive:
        assert archive.namelist() == ["demo.txt"]
        assert archive.read("demo.txt") == b"payload"


def test_existing_outputs_and_path_resolution_use_repo_root(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(RUNNER, "REPO_ROOT", tmp_path)
    file_path = tmp_path / "nested" / "result.tsv"
    file_path.parent.mkdir(parents=True)
    file_path.write_text("ok", encoding="utf-8")

    assert RUNNER.path_exists("nested/result.tsv")
    assert RUNNER.existing_outputs(["nested/result.tsv", "missing.tsv"]) == ["nested/result.tsv"]


def test_upload_helpers_slugify_and_save_files(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(UPLOADS, "UPLOAD_ROOT", tmp_path / "uploads")
    upload_dir = UPLOADS.create_upload_dir("schema demo")
    saved = UPLOADS.save_uploaded_file(upload_dir, FakeUpload("sample 1.fastq.gz", b"reads"))

    assert upload_dir.exists()
    assert saved.name == "sample-1.fastq.gz"
    assert saved.read_bytes() == b"reads"


def test_upload_helpers_reject_parent_like_names(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(UPLOADS, "UPLOAD_ROOT", tmp_path / "uploads")
    upload_dir = UPLOADS.create_upload_dir("schema demo")
    saved = UPLOADS.save_uploaded_file(upload_dir, FakeUpload("..", b"reads"))

    assert saved.parent == upload_dir
    assert saved.name == "upload"


def test_upload_helpers_enforce_size_limit(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(UPLOADS, "UPLOAD_ROOT", tmp_path / "uploads")
    monkeypatch.setattr(UPLOADS, "MAX_UPLOAD_SIZE_MB", 0)
    upload_dir = UPLOADS.create_upload_dir("schema demo")

    try:
        UPLOADS.save_uploaded_file(upload_dir, LargeFakeUpload("sample.fastq.gz", b"reads"))
    except ValueError as exc:
        assert "exceeds" in str(exc)
    else:
        raise AssertionError("Expected upload size limit to raise ValueError")


def test_build_import_scheme_command_supports_chained_schema_creation() -> None:
    command = RUNNER.build_import_scheme_command(
        database="pubmlst_demo",
        scheme_id=7,
        output_dir="schemas/demo",
        include_profiles=True,
        scheme_type="wgmlst",
        run_schema_creation=True,
        reference="demo_ref",
        pipeline_output_dir="output/schema",
        threads=4,
    )

    assert "--database" in command
    assert "--scheme-type" in command
    assert "--run-schema-creation" in command
    assert "demo_ref" in command
    assert "output/schema" in command


def test_build_import_pubmlst_benchmark_pack_command_contains_expected_arguments() -> None:
    command = RUNNER.build_import_pubmlst_benchmark_pack_command(
        scheme_database="pubmlst_demo_seqdef",
        isolate_database="pubmlst_demo_isolates",
        scheme_id=7,
        isolate_ids=[101, 102],
        output_dir="benchmark_packs/demo",
        species="Neisseria meningitidis",
        download_contigs=True,
        discover_isolates=2,
    )

    assert "--scheme-database" in command
    assert "--isolate-database" in command
    assert "--scheme-id" in command
    assert command.count("--isolate-id") == 2
    assert "--discover-isolates" in command
    assert "--download-contigs" in command
    assert "benchmark_packs/demo" in command


def test_build_import_enterobase_scheme_command_contains_expected_arguments() -> None:
    command = RUNNER.build_import_enterobase_scheme_command(
        database="senterica",
        scheme_name="wgMLST",
        output_dir="public_schemes/enterobase_senterica",
    )

    assert "--database" in command
    assert "--scheme-name" in command
    assert "--output-dir" in command
    assert "wgMLST" in command


def test_build_import_enterobase_scheme_command_supports_list_mode() -> None:
    command = RUNNER.build_import_enterobase_scheme_command(
        database="senterica",
        scheme_name="wgMLST",
        output_dir="public_schemes/enterobase_senterica",
        list_schemes=True,
    )

    assert "--list-schemes" in command
    assert "--scheme-name" not in command


def test_build_profile_compare_command_contains_expected_arguments() -> None:
    command = RUNNER.build_profile_compare_command(
        profile_a="results/a.tsv",
        profile_b="results/b.tsv",
        output_dir="results/compare",
        label_a="sample_a",
        label_b="sample_b",
    )

    assert str(RUNNER.WORKFLOW_SCRIPT) in command
    assert "profile_compare" in command
    assert "--profile-a" in command
    assert "--profile-b" in command
    assert "results/compare" in command
    assert "sample_a" in command
    assert "sample_b" in command


def test_build_profile_matrix_and_batch_commands_contain_expected_arguments() -> None:
    matrix_command = RUNNER.build_profile_matrix_command(
        ["results/a.tsv", "results/b.tsv"],
        "results/matrix",
        distance_threshold=0.1,
    )
    batch_command = RUNNER.build_profile_compare_batch_command(
        ["results/a.tsv", "results/b.tsv"],
        "results/batch",
    )

    assert "profile_matrix" in matrix_command
    assert "--distance-threshold" in matrix_command
    assert "0.1" in matrix_command
    assert "profile_compare_batch" in batch_command


def test_build_profile_benchmark_command_contains_expected_arguments() -> None:
    command = RUNNER.build_profile_benchmark_command(
        predicted_dir="results/predicted",
        truth_dir="results/truth",
        output_dir="results/benchmark",
    )

    assert str(RUNNER.WORKFLOW_SCRIPT) in command
    assert "profile_benchmark" in command
    assert "--predicted-dir" in command
    assert "--truth-dir" in command
    assert "results/benchmark" in command


def test_build_profile_benchmark_command_supports_pack_dir() -> None:
    command = RUNNER.build_profile_benchmark_command(
        predicted_dir=None,
        truth_dir=None,
        output_dir="results/benchmark",
        benchmark_pack_dir="tests/fixtures/benchmark/public_species_packs",
    )

    assert "--benchmark-pack-dir" in command
    assert "tests/fixtures/benchmark/public_species_packs" in command


def test_build_schema_qc_command_contains_expected_arguments() -> None:
    command = RUNNER.build_schema_qc_command("schemas/demo", "output/qc")

    assert str(RUNNER.SCHEMA_QC_SCRIPT) in command
    assert "--schema-dir" in command
    assert "--output-dir" in command


def test_download_bundle_excludes_unsafe_absolute_paths(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(RUNNER, "REPO_ROOT", tmp_path / "repo")
    monkeypatch.setattr(RUNNER, "JOBS_ROOT", tmp_path / "repo" / "webapp_jobs")
    monkeypatch.setattr(RUNNER, "UPLOADS_ROOT", tmp_path / "repo" / "webapp_uploads")
    RUNNER.REPO_ROOT.mkdir(parents=True)

    unsafe_file = tmp_path / "outside.txt"
    unsafe_file.write_text("secret", encoding="utf-8")

    bundle = RUNNER.create_download_bundle([str(unsafe_file)])

    import zipfile

    with zipfile.ZipFile(io.BytesIO(bundle)) as archive:
        assert archive.namelist() == []


def test_is_safe_web_output_path_rejects_parent_traversal_and_outside_absolute(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(RUNNER, "REPO_ROOT", tmp_path / "repo")
    RUNNER.REPO_ROOT.mkdir(parents=True)

    assert RUNNER.is_safe_web_output_path("results/output")
    assert not RUNNER.is_safe_web_output_path("../outside")
    assert not RUNNER.is_safe_web_output_path(str(tmp_path / "outside"))


def test_existing_outputs_excludes_unsafe_absolute_paths(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(RUNNER, "REPO_ROOT", tmp_path / "repo")
    monkeypatch.setattr(RUNNER, "JOBS_ROOT", tmp_path / "repo" / "webapp_jobs")
    monkeypatch.setattr(RUNNER, "UPLOADS_ROOT", tmp_path / "repo" / "webapp_uploads")
    RUNNER.REPO_ROOT.mkdir(parents=True)
    safe = RUNNER.REPO_ROOT / "safe.txt"
    safe.write_text("ok", encoding="utf-8")
    outside = tmp_path / "outside.txt"
    outside.write_text("nope", encoding="utf-8")

    assert RUNNER.existing_outputs(["safe.txt", str(outside)]) == ["safe.txt"]


def test_load_result_metrics_reads_wgmlst_and_compare_outputs(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(RUNNER, "REPO_ROOT", tmp_path)
    monkeypatch.setattr(RUNNER, "JOBS_ROOT", tmp_path / "webapp_jobs")
    monkeypatch.setattr(RUNNER, "UPLOADS_ROOT", tmp_path / "webapp_uploads")
    (tmp_path / "demo_wgmlst.tsv").write_text("abcZ\t12\nadk\tLNF\n", encoding="utf-8")
    (tmp_path / "wgmlst_profile_summary.json").write_text(
        json.dumps(
            {
                "decision": "different",
                "comparable_loci": 1,
                "differing_loci": 1,
                "unresolved_loci": 1,
                "allele_distance": 1.0,
            }
        ),
        encoding="utf-8",
    )
    (tmp_path / "wgmlst_distance_summary.json").write_text(
        json.dumps(
            {
                "sample_count": 3,
                "pair_count": 3,
                "pairs": [
                    {"decision": "different"},
                    {"decision": "inconclusive"},
                    {"decision": "different"},
                ],
            }
        ),
        encoding="utf-8",
    )
    (tmp_path / "wgmlst_benchmark_summary.json").write_text(
        json.dumps(
            {
                "matched_samples": 2,
                "exact_allele_id_match_count": 1,
                "mean_concordance": 0.94,
            }
        ),
        encoding="utf-8",
    )
    (tmp_path / "sample_novel_alleles.tsv").write_text(
        "locus\tallele_id\tsequence_length\tstatus\nabc\t9\t9\tnovel_candidate\n",
        encoding="utf-8",
    )
    (tmp_path / "wgmlst_benchmark_pack_summary.json").write_text(
        json.dumps(
            {
                "species_count": 2,
                "mean_concordance": 0.97,
            }
        ),
        encoding="utf-8",
    )

    metrics = RUNNER.load_result_metrics(
        [
            "demo_wgmlst.tsv",
            "sample_novel_alleles.tsv",
            "wgmlst_profile_summary.json",
            "wgmlst_distance_summary.json",
            "wgmlst_benchmark_summary.json",
            "wgmlst_benchmark_pack_summary.json",
        ]
    )

    assert {metric.kind for metric in metrics} == {
        "wgmlst_profile",
        "novel_alleles",
        "profile_compare",
        "profile_matrix",
        "profile_benchmark",
        "profile_benchmark_collection",
    }


def test_prune_jobs_removes_old_completed_jobs(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(RUNNER, "JOBS_ROOT", tmp_path / "jobs")
    job_dir = RUNNER.JOBS_ROOT / "oldjob"
    job_dir.mkdir(parents=True)
    (job_dir / "stdout.log").write_text("", encoding="utf-8")
    (job_dir / "stderr.log").write_text("", encoding="utf-8")
    (job_dir / "command.json").write_text(json.dumps({"command": ["python", "x.py"]}), encoding="utf-8")
    (job_dir / "metadata.json").write_text(
        json.dumps(
            {
                "job_id": "oldjob",
                "label": "schema_qc",
                "status": "completed",
                "command": ["python", "x.py"],
                "command_text": "python x.py",
                "created_at": "2020-01-01T00:00:00+00:00",
                "started_at": None,
                "finished_at": None,
                "pid": None,
                "returncode": 0,
                "outputs": [],
            }
        ),
        encoding="utf-8",
    )

    removed = RUNNER.prune_jobs(older_than_hours=1)

    assert removed == 1
    assert not job_dir.exists()
