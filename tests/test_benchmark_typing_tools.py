from __future__ import annotations

import json
import sys
from importlib import util
from pathlib import Path

MODULE_PATH = Path(__file__).resolve().parent.parent / "workflow" / "scripts" / "benchmark_typing_tools.py"
SPEC = util.spec_from_file_location("benchmark_typing_tools", MODULE_PATH)
MODULE = util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def test_benchmark_tools_marks_missing_commands_as_skipped() -> None:
    report = MODULE.benchmark_tools({"missing": {"label": "Missing", "command": ["definitely-not-installed"]}})
    assert report["tools"][0]["status"] == "skipped"
    assert report["tools"][0]["reason"] == "command_not_found"


def test_benchmark_tools_marks_missing_required_files_as_skipped(tmp_path: Path) -> None:
    report = MODULE.benchmark_tools(
        {
            "milestone": {
                "label": "Milestone",
                "command": ["python", "-c", "print('ok')"],
                "workdir": str(tmp_path),
                "required_files": ["missing.txt"],
            }
        }
    )
    assert report["tools"][0]["status"] == "skipped"
    assert report["tools"][0]["reason"] == "required_files_missing"


def test_write_report_creates_json_tsv_and_html(tmp_path: Path) -> None:
    report = {
        "kind": "tool_benchmark_suite",
        "tool_count": 1,
        "comparison_summary": {"completed_tool_count": 1},
        "tools": [{"tool": "milestone", "label": "Milestone", "status": "completed", "runtime_seconds": 0.1}],
    }
    json_path, tsv_path, html_path = MODULE.write_report(report, tmp_path)
    assert json.loads(json_path.read_text(encoding="utf-8"))["tool_count"] == 1
    assert "Milestone" in tsv_path.read_text(encoding="utf-8")
    assert "Milestone" in html_path.read_text(encoding="utf-8")


def test_load_config_reads_benchmarking_profile(tmp_path: Path) -> None:
    config_path = tmp_path / "profile.json"
    config_path.write_text(json.dumps({"demo": {"label": "Demo", "command": ["python", "--help"]}}), encoding="utf-8")
    loaded = MODULE.load_config(str(config_path))
    assert loaded["demo"]["label"] == "Demo"


def test_benchmark_tools_collects_comparison_summary() -> None:
    report = MODULE.benchmark_tools(
        {
            "milestone": {"label": "Milestone", "command": ["python", "-c", "print('ok')"]},
            "demo": {"label": "Demo", "command": ["python", "-c", "print('ok')"]},
        }
    )
    assert report["comparison_summary"]["completed_tool_count"] == 2
