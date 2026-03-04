from __future__ import annotations

import io
import json
import os
import shlex
import subprocess
import sys
import zipfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from uuid import uuid4

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from workflow.scripts.result_metrics import (
    summarize_batch_compare_summary,
    summarize_distance_summary,
    summarize_benchmark_summary,
    summarize_benchmark_pack_summary,
    summarize_comparison_summary,
    summarize_novel_alleles,
    summarize_schema_qc,
    summarize_wgmlst_profile,
)

WORKFLOW_SCRIPT = REPO_ROOT / "workflow" / "milestone.py"
IMPORT_SCHEME_SCRIPT = REPO_ROOT / "workflow" / "scripts" / "import_pubmlst_scheme.py"
IMPORT_BENCHMARK_PACK_SCRIPT = REPO_ROOT / "workflow" / "scripts" / "import_pubmlst_benchmark_pack.py"
IMPORT_ENTEROBASE_SCHEME_SCRIPT = REPO_ROOT / "workflow" / "scripts" / "import_enterobase_scheme.py"
SCHEMA_QC_SCRIPT = REPO_ROOT / "workflow" / "scripts" / "schema_qc.py"
JOB_WORKER_SCRIPT = Path(__file__).resolve().parent / "job_worker.py"
JOBS_ROOT = REPO_ROOT / "webapp_jobs"
UPLOADS_ROOT = REPO_ROOT / "webapp_uploads"


@dataclass(frozen=True)
class JobSnapshot:
    job_id: str
    label: str
    status: str
    command: list[str]
    command_text: str
    created_at: str
    started_at: str | None
    finished_at: str | None
    pid: int | None
    returncode: int | None
    stdout_path: Path
    stderr_path: Path
    outputs: list[str]
    job_dir: Path


@dataclass(frozen=True)
class ResultMetric:
    kind: str
    values: dict[str, object]


def timestamp() -> str:
    return datetime.now(timezone.utc).isoformat()


def ensure_jobs_root() -> Path:
    JOBS_ROOT.mkdir(parents=True, exist_ok=True)
    return JOBS_ROOT


def build_schema_command(
    reference: str,
    schema_name: str,
    output: str,
    threads: int,
    dryrun: bool,
) -> list[str]:
    command = [
        sys.executable,
        str(WORKFLOW_SCRIPT),
        "schema_creation",
        "--reference",
        reference,
        "--schema_name",
        schema_name,
        "--output",
        output,
        "--threads",
        str(threads),
    ]
    if dryrun:
        command.append("--dryrun")
    return command


def command_to_text(command: list[str]) -> str:
    return shlex.join(command)


def build_allele_command(
    reference: str,
    schema_name: str,
    output: str,
    read1: str,
    read2: str,
    threads: int,
    aligner: str,
    update_reference: bool,
    dryrun: bool,
    translation_table: int = 11,
    allowed_start_codons: str | None = None,
) -> list[str]:
    command = [
        sys.executable,
        str(WORKFLOW_SCRIPT),
        "allele_calling",
        "--reference",
        reference,
        "--schema_name",
        schema_name,
        "--output",
        output,
        "--read1",
        read1,
        "--read2",
        read2,
        "--threads",
        str(threads),
        "--aligner",
        aligner,
        "--translation-table",
        str(translation_table),
    ]
    if allowed_start_codons:
        command.extend(["--allowed-start-codons", allowed_start_codons])
    if update_reference:
        command.append("--update_reference")
    if dryrun:
        command.append("--dryrun")
    return command


def build_import_scheme_command(
    database: str,
    scheme_id: int,
    output_dir: str,
    include_profiles: bool = False,
    scheme_type: str | None = None,
    run_schema_creation: bool = False,
    reference: str | None = None,
    pipeline_output_dir: str | None = None,
    threads: int = 1,
) -> list[str]:
    command = [
        sys.executable,
        str(IMPORT_SCHEME_SCRIPT),
        "--database",
        database,
        "--scheme-id",
        str(scheme_id),
        "--output-dir",
        output_dir,
    ]
    if include_profiles:
        command.append("--include-profiles")
    if scheme_type:
        command.extend(["--scheme-type", scheme_type])
    if run_schema_creation:
        command.append("--run-schema-creation")
        command.extend(["--reference", reference or "public_reference"])
        command.extend(["--pipeline-output-dir", pipeline_output_dir or "public_output"])
        command.extend(["--threads", str(threads)])
    return command


def build_import_pubmlst_benchmark_pack_command(
    scheme_database: str,
    isolate_database: str,
    scheme_id: int,
    isolate_ids: list[int],
    output_dir: str,
    species: str | None = None,
    download_contigs: bool = False,
    discover_isolates: int = 0,
) -> list[str]:
    command = [
        sys.executable,
        str(IMPORT_BENCHMARK_PACK_SCRIPT),
        "--scheme-database",
        scheme_database,
        "--isolate-database",
        isolate_database,
        "--scheme-id",
        str(scheme_id),
        "--output-dir",
        output_dir,
    ]
    for isolate_id in isolate_ids:
        command.extend(["--isolate-id", str(isolate_id)])
    if species:
        command.extend(["--species", species])
    if discover_isolates > 0:
        command.extend(["--discover-isolates", str(discover_isolates)])
    if download_contigs:
        command.append("--download-contigs")
    return command


def build_import_enterobase_scheme_command(
    database: str,
    scheme_name: str,
    output_dir: str,
    list_schemes: bool = False,
) -> list[str]:
    command = [
        sys.executable,
        str(IMPORT_ENTEROBASE_SCHEME_SCRIPT),
        "--database",
        database,
    ]
    if list_schemes:
        command.append("--list-schemes")
    else:
        command.extend(["--scheme-name", scheme_name, "--output-dir", output_dir])
    return command


def build_profile_compare_command(
    profile_a: str,
    profile_b: str,
    output_dir: str,
    label_a: str,
    label_b: str,
) -> list[str]:
    return [
        sys.executable,
        str(WORKFLOW_SCRIPT),
        "profile_compare",
        "--profile-a",
        profile_a,
        "--profile-b",
        profile_b,
        "--output-dir",
        output_dir,
        "--label-a",
        label_a,
        "--label-b",
        label_b,
    ]


def build_profile_benchmark_command(
    predicted_dir: str | None,
    truth_dir: str | None,
    output_dir: str,
    benchmark_pack_dir: str | None = None,
) -> list[str]:
    command = [
        sys.executable,
        str(WORKFLOW_SCRIPT),
        "profile_benchmark",
        "--output",
        output_dir,
    ]
    if benchmark_pack_dir:
        command.extend(["--benchmark-pack-dir", benchmark_pack_dir])
    else:
        command.extend(["--predicted-dir", predicted_dir or "", "--truth-dir", truth_dir or ""])
    return command


def build_profile_matrix_command(
    profiles: list[str],
    output_dir: str,
    distance_threshold: float = 0.0,
) -> list[str]:
    command = [
        sys.executable,
        str(WORKFLOW_SCRIPT),
        "profile_matrix",
        "--profiles",
        *profiles,
        "--output",
        output_dir,
        "--distance-threshold",
        str(distance_threshold),
    ]
    return command


def build_profile_compare_batch_command(
    profiles: list[str],
    output_dir: str,
) -> list[str]:
    return [
        sys.executable,
        str(WORKFLOW_SCRIPT),
        "profile_compare_batch",
        "--profiles",
        *profiles,
        "--output",
        output_dir,
    ]


def build_schema_qc_command(
    schema_dir: str,
    output_dir: str,
    translation_table: int = 11,
    allowed_start_codons: str | None = None,
    manifest_path: str | None = None,
) -> list[str]:
    command = [
        sys.executable,
        str(SCHEMA_QC_SCRIPT),
        "--schema-dir",
        schema_dir,
        "--output-dir",
        output_dir,
        "--translation-table",
        str(translation_table),
    ]
    if allowed_start_codons:
        command.extend(["--allowed-start-codons", allowed_start_codons])
    if manifest_path:
        command.extend(["--manifest", manifest_path])
    return command


def existing_outputs(paths: list[str]) -> list[str]:
    existing = []
    for path in paths:
        resolved = resolve_downloadable_path(path)
        if resolved is not None and resolved.exists():
            existing.append(path)
    return existing


def path_exists(path_text: str) -> bool:
    return resolve_repo_path(path_text).exists()


def resolve_repo_path(path_text: str) -> Path:
    path = Path(path_text)
    if path.is_absolute():
        return path
    return REPO_ROOT / path


def _is_relative_to(path: Path, root: Path) -> bool:
    try:
        path.resolve().relative_to(root.resolve())
    except ValueError:
        return False
    return True


def is_safe_download_path(path: Path) -> bool:
    return any(
        _is_relative_to(path, root)
        for root in [REPO_ROOT, JOBS_ROOT, UPLOADS_ROOT]
    )


def is_safe_web_output_path(path_text: str) -> bool:
    path = resolve_repo_path(path_text)
    if path.is_absolute() and not _is_relative_to(path, REPO_ROOT):
        return False
    normalized_parts = [part for part in Path(path_text).parts if part not in {"."}]
    return ".." not in normalized_parts


def resolve_downloadable_path(path_text: str) -> Path | None:
    path = resolve_repo_path(path_text)
    if not is_safe_download_path(path):
        return None
    return path


def read_text_tail(path: Path, max_chars: int = 12000) -> str:
    if not path.exists():
        return ""
    content = path.read_text(encoding="utf-8", errors="replace")
    if len(content) <= max_chars:
        return content
    return content[-max_chars:]


def load_result_metric(path_text: str) -> ResultMetric | None:
    path = resolve_downloadable_path(path_text)
    if path is None or not path.exists() or not path.is_file():
        return None
    if path.name.endswith("_wgmlst.tsv"):
        return ResultMetric("wgmlst_profile", summarize_wgmlst_profile(path))
    if path.name.endswith("_novel_alleles.tsv"):
        return ResultMetric("novel_alleles", summarize_novel_alleles(path))
    if path.name == "wgmlst_profile_summary.json":
        return ResultMetric("profile_compare", summarize_comparison_summary(path))
    if path.name == "wgmlst_distance_summary.json":
        return ResultMetric("profile_matrix", summarize_distance_summary(path))
    if path.name == "wgmlst_batch_compare_summary.json":
        return ResultMetric("profile_compare_batch", summarize_batch_compare_summary(path))
    if path.name == "wgmlst_benchmark_summary.json":
        return ResultMetric("profile_benchmark", summarize_benchmark_summary(path))
    if path.name == "wgmlst_benchmark_pack_summary.json":
        return ResultMetric("profile_benchmark_collection", summarize_benchmark_pack_summary(path))
    if path.name == "schema_qc_summary.json":
        return ResultMetric("schema_qc", summarize_schema_qc(path))
    return None


def load_result_metrics(paths: list[str]) -> list[ResultMetric]:
    metrics = []
    for path_text in paths:
        metric = load_result_metric(path_text)
        if metric is not None:
            metrics.append(metric)
    return metrics


def create_download_bundle(paths: list[str]) -> bytes:
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for path_text in paths:
            path = resolve_downloadable_path(path_text)
            if path is not None and path.exists() and path.is_file():
                archive.write(path, arcname=path.name)
    buffer.seek(0)
    return buffer.getvalue()


def _job_paths(job_dir: Path) -> dict[str, Path]:
    return {
        "metadata": job_dir / "metadata.json",
        "stdout": job_dir / "stdout.log",
        "stderr": job_dir / "stderr.log",
        "command": job_dir / "command.json",
    }


def start_background_job(label: str, command: list[str], outputs: list[str]) -> str:
    ensure_jobs_root()
    job_id = uuid4().hex[:12]
    job_dir = JOBS_ROOT / job_id
    job_dir.mkdir(parents=True, exist_ok=True)
    paths = _job_paths(job_dir)

    metadata = {
        "job_id": job_id,
        "label": label,
        "status": "queued",
        "command": command,
        "command_text": shlex.join(command),
        "created_at": timestamp(),
        "started_at": None,
        "finished_at": None,
        "pid": None,
        "returncode": None,
        "outputs": outputs,
    }
    paths["metadata"].write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    paths["command"].write_text(json.dumps({"command": command}, indent=2), encoding="utf-8")
    paths["stdout"].write_text("", encoding="utf-8")
    paths["stderr"].write_text("", encoding="utf-8")

    subprocess.Popen(
        [sys.executable, str(JOB_WORKER_SCRIPT), str(job_dir)],
        cwd=REPO_ROOT,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        start_new_session=True,
    )
    return job_id


def retry_job(snapshot: JobSnapshot) -> str:
    return start_background_job(snapshot.label, snapshot.command, snapshot.outputs)


def list_jobs(limit: int = 10) -> list[JobSnapshot]:
    ensure_jobs_root()
    snapshots = []
    for metadata_path in sorted(JOBS_ROOT.glob("*/metadata.json"), reverse=True):
        snapshot = get_job_snapshot(metadata_path.parent.name)
        if snapshot is not None:
            snapshots.append(snapshot)
    snapshots.sort(key=lambda snapshot: snapshot.created_at, reverse=True)
    return snapshots[:limit]


def prune_jobs(older_than_hours: int = 24, statuses: tuple[str, ...] = ("completed", "failed")) -> int:
    ensure_jobs_root()
    cutoff_seconds = older_than_hours * 3600
    removed = 0
    now = datetime.now(timezone.utc)
    for metadata_path in JOBS_ROOT.glob("*/metadata.json"):
        snapshot = get_job_snapshot(metadata_path.parent.name)
        if snapshot is None or snapshot.status not in statuses:
            continue
        created_at = datetime.fromisoformat(snapshot.created_at)
        age_seconds = (now - created_at).total_seconds()
        if age_seconds >= cutoff_seconds:
            for path in snapshot.job_dir.glob("*"):
                if path.is_file():
                    path.unlink(missing_ok=True)
            snapshot.job_dir.rmdir()
            removed += 1
    return removed


def get_job_snapshot(job_id: str) -> JobSnapshot | None:
    job_dir = JOBS_ROOT / job_id
    paths = _job_paths(job_dir)
    if not paths["metadata"].exists():
        return None
    metadata = json.loads(paths["metadata"].read_text(encoding="utf-8"))
    return JobSnapshot(
        job_id=metadata["job_id"],
        label=metadata["label"],
        status=metadata["status"],
        command=metadata["command"],
        command_text=metadata["command_text"],
        created_at=metadata["created_at"],
        started_at=metadata.get("started_at"),
        finished_at=metadata.get("finished_at"),
        pid=metadata.get("pid"),
        returncode=metadata.get("returncode"),
        stdout_path=paths["stdout"],
        stderr_path=paths["stderr"],
        outputs=metadata.get("outputs", []),
        job_dir=job_dir,
    )


def refresh_job_status(snapshot: JobSnapshot) -> JobSnapshot:
    latest = get_job_snapshot(snapshot.job_id)
    if latest is None:
        return snapshot
    if latest.status == "running" and latest.pid is not None and not _pid_is_alive(latest.pid):
        return get_job_snapshot(snapshot.job_id) or latest
    return latest


def _pid_is_alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except OSError:
        return False
    return True
