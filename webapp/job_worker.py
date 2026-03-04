from __future__ import annotations

import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

JOBS_ROOT = Path(__file__).resolve().parent.parent / "webapp_jobs"
REPO_ROOT = Path(__file__).resolve().parent.parent
ALLOWED_SCRIPTS = {
    REPO_ROOT / "workflow" / "milestone.py",
    REPO_ROOT / "workflow" / "scripts" / "import_pubmlst_scheme.py",
    REPO_ROOT / "workflow" / "scripts" / "import_pubmlst_benchmark_pack.py",
    REPO_ROOT / "workflow" / "scripts" / "import_enterobase_scheme.py",
    REPO_ROOT / "workflow" / "scripts" / "schema_qc.py",
}


def timestamp() -> str:
    return datetime.now(timezone.utc).isoformat()


def update_metadata(metadata_path: Path, **updates) -> dict:
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    metadata.update(updates)
    metadata_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    return metadata


def main() -> int:
    job_dir = Path(sys.argv[1]).resolve()
    try:
        job_dir.relative_to(JOBS_ROOT.resolve())
    except ValueError as exc:
        raise SystemExit(f"Refusing to run job outside {JOBS_ROOT}") from exc

    metadata_path = job_dir / "metadata.json"
    command_path = job_dir / "command.json"
    stdout_path = job_dir / "stdout.log"
    stderr_path = job_dir / "stderr.log"

    command = json.loads(command_path.read_text(encoding="utf-8"))["command"]
    if len(command) < 2:
        raise SystemExit("Refusing to run empty or malformed command.")
    script_path = Path(command[1]).resolve()
    if script_path not in ALLOWED_SCRIPTS:
        raise SystemExit(f"Refusing to run non-whitelisted command target: {script_path}")
    update_metadata(metadata_path, status="running", started_at=timestamp())

    with stdout_path.open("a", encoding="utf-8") as stdout_handle, stderr_path.open(
        "a", encoding="utf-8"
    ) as stderr_handle:
        process = subprocess.Popen(
            command,
            cwd=job_dir.parent.parent,
            text=True,
            stdout=stdout_handle,
            stderr=stderr_handle,
        )
        update_metadata(metadata_path, pid=process.pid)
        returncode = process.wait()

    status = "completed" if returncode == 0 else "failed"
    update_metadata(
        metadata_path,
        status=status,
        finished_at=timestamp(),
        returncode=returncode,
    )
    return returncode


if __name__ == "__main__":
    raise SystemExit(main())
