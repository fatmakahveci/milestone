#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import subprocess
import time
from pathlib import Path

DEFAULT_TOOL_CONFIG = {
    "milestone": {
        "label": "Milestone",
        "command": ["python", "workflow/milestone.py", "profile_compare", "--help"],
        "workdir": ".",
        "timeout_seconds": 60,
        "required_files": ["workflow/milestone.py"],
    },
    "chewbbaca": {
        "label": "chewBBACA",
        "command": ["chewBBACA.py", "--help"],
        "timeout_seconds": 60,
        "notes": "Install chewBBACA in PATH to collect real runtime measurements.",
    },
    "mentalist": {
        "label": "MentaLiST",
        "command": ["mentalist", "--help"],
        "timeout_seconds": 60,
        "notes": "Install MentaLiST in PATH to collect real runtime measurements.",
    },
}


def load_config(path: str | None) -> dict[str, dict[str, object]]:
    if not path:
        return DEFAULT_TOOL_CONFIG
    return json.loads(Path(path).read_text(encoding="utf-8"))


def _command_available(command: list[str]) -> bool:
    executable = command[0]
    if executable in {"python", "python3"}:
        return True
    return shutil.which(executable) is not None


def _required_files_exist(required_files: list[str], workdir: str) -> bool:
    for path_text in required_files:
        path = Path(workdir) / path_text
        if not path.exists():
            return False
    return True


def _run_command(
    command: list[str],
    *,
    workdir: str,
    env: dict[str, str],
    timeout_seconds: int,
) -> tuple[str, int | None, str, str, float]:
    start = time.perf_counter()
    try:
        completed = subprocess.run(
            command,
            capture_output=True,
            text=True,
            check=False,
            cwd=workdir,
            env=env,
            timeout=timeout_seconds,
        )
        status = "completed" if completed.returncode == 0 else "failed"
        return (
            status,
            completed.returncode,
            completed.stdout[:200],
            completed.stderr[:200],
            time.perf_counter() - start,
        )
    except subprocess.TimeoutExpired as exc:
        return (
            "timed_out",
            None,
            (exc.stdout or "")[:200],
            (exc.stderr or "")[:200],
            time.perf_counter() - start,
        )


def _resolve_version_command(entry: dict[str, object], command: list[str]) -> list[str] | None:
    version_command = entry.get("version_command")
    if version_command:
        return [str(part) for part in version_command]
    if command and command[0] in {"python", "python3"}:
        return None
    if command:
        return [command[0], "--version"]
    return None


def benchmark_tools(config: dict[str, dict[str, object]]) -> dict[str, object]:
    results = []
    for tool_id, entry in config.items():
        command = [str(part) for part in entry["command"]]
        workdir = str(entry.get("workdir", "."))
        timeout_seconds = int(entry.get("timeout_seconds", 60))
        env = {**os.environ, **{str(key): str(value) for key, value in entry.get("env", {}).items()}}
        required_files = [str(path) for path in entry.get("required_files", [])]

        if not _command_available(command):
            results.append(
                {
                    "tool": tool_id,
                    "label": entry.get("label", tool_id),
                    "status": "skipped",
                    "reason": "command_not_found",
                    "notes": entry.get("notes", ""),
                }
            )
            continue
        if required_files and not _required_files_exist(required_files, workdir):
            results.append(
                {
                    "tool": tool_id,
                    "label": entry.get("label", tool_id),
                    "status": "skipped",
                    "reason": "required_files_missing",
                    "notes": entry.get("notes", ""),
                }
            )
            continue

        version = None
        version_command = _resolve_version_command(entry, command)
        if version_command and _command_available(version_command):
            version_status, _, version_stdout, version_stderr, _ = _run_command(
                version_command,
                workdir=workdir,
                env=env,
                timeout_seconds=min(timeout_seconds, 30),
            )
            if version_status == "completed":
                version = (version_stdout or version_stderr).strip() or None

        status, returncode, stdout_preview, stderr_preview, elapsed = _run_command(
            command,
            workdir=workdir,
            env=env,
            timeout_seconds=timeout_seconds,
        )

        results.append(
            {
                "tool": tool_id,
                "label": entry.get("label", tool_id),
                "status": status,
                "returncode": returncode,
                "runtime_seconds": round(elapsed, 4),
                "workdir": workdir,
                "command": command,
                "stdout_preview": stdout_preview,
                "stderr_preview": stderr_preview,
                "notes": entry.get("notes", ""),
                "version": version,
            }
        )
    return {
        "kind": "tool_benchmark_suite",
        "tool_count": len(results),
        "tools": results,
        "comparison_summary": summarize_comparison(results),
    }


def summarize_comparison(results: list[dict[str, object]]) -> dict[str, object]:
    milestone = next((result for result in results if result.get("tool") == "milestone"), None)
    baseline_runtime = milestone.get("runtime_seconds") if milestone and milestone.get("status") == "completed" else None
    completed_tools = [result for result in results if result.get("status") == "completed"]
    slower_than_milestone = []
    faster_than_milestone = []
    if baseline_runtime is not None:
        for result in completed_tools:
            if result.get("tool") == "milestone":
                continue
            runtime = result.get("runtime_seconds")
            if runtime is None:
                continue
            if runtime > baseline_runtime:
                slower_than_milestone.append(result["label"])
            elif runtime < baseline_runtime:
                faster_than_milestone.append(result["label"])
    return {
        "completed_tool_count": len(completed_tools),
        "skipped_tool_count": sum(1 for result in results if result.get("status") == "skipped"),
        "failed_tool_count": sum(1 for result in results if result.get("status") == "failed"),
        "baseline_tool": milestone.get("label") if milestone else None,
        "baseline_runtime_seconds": baseline_runtime,
        "faster_than_milestone": faster_than_milestone,
        "slower_than_milestone": slower_than_milestone,
    }


def write_report(report: dict[str, object], output_dir: str | Path) -> tuple[Path, Path, Path]:
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    json_path = output_path / "tool_benchmark_report.json"
    json_path.write_text(json.dumps(report, indent=2), encoding="utf-8")

    tsv_path = output_path / "tool_benchmark_report.tsv"
    with tsv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["tool", "status", "runtime_seconds", "returncode", "reason", "version", "notes"])
        for tool in report["tools"]:
            writer.writerow(
                [
                    tool["label"],
                    tool["status"],
                    tool.get("runtime_seconds"),
                    tool.get("returncode"),
                    tool.get("reason", ""),
                    tool.get("version", ""),
                    tool.get("notes", ""),
                ]
            )

    html_path = output_path / "tool_benchmark_report.html"
    rows = "".join(
        f"<tr><td>{tool['label']}</td><td>{tool['status']}</td><td>{tool.get('runtime_seconds')}</td>"
        f"<td>{tool.get('returncode', '')}</td><td>{tool.get('reason', '')}</td><td>{tool.get('version', '')}</td><td>{tool.get('notes', '')}</td></tr>"
        for tool in report["tools"]
    )
    html_path.write_text(
        "<html><head><meta charset='utf-8'><title>Typing tool benchmark suite</title>"
        "<meta name='description' content='Milestone benchmark harness for wgMLST typing tools including chewBBACA and MentaLiST.'>"
        "<style>body{font-family:Arial,sans-serif;margin:2rem;}table{border-collapse:collapse;width:100%;}"
        "th,td{border:1px solid #ccc;padding:0.4rem;text-align:left;}th{background:#f5f5f5;}</style></head>"
        "<body><h1>Typing tool benchmark suite</h1>"
        f"<p>Completed tools: {report.get('comparison_summary', {}).get('completed_tool_count')}</p>"
        "<table><tr><th>Tool</th><th>Status</th><th>Runtime (s)</th><th>Return code</th><th>Reason</th><th>Version</th><th>Notes</th></tr>"
        f"{rows}</table></body></html>",
        encoding="utf-8",
    )
    return json_path, tsv_path, html_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run a reproducible benchmark harness for wgMLST typing tools.")
    parser.add_argument("--config", help="Optional JSON config overriding the default tool commands.")
    parser.add_argument("--output-dir", required=True, help="Directory for benchmark outputs.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    report = benchmark_tools(load_config(args.config))
    json_path, tsv_path, html_path = write_report(report, args.output_dir)
    print(f"Tool benchmark report written to {json_path}")
    print(f"TSV report written to {tsv_path}")
    print(f"HTML report written to {html_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
