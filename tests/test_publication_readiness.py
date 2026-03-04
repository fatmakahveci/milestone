from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path


MODULE_PATH = Path(__file__).resolve().parent.parent / "workflow" / "scripts" / "publication_readiness.py"
SPEC = importlib.util.spec_from_file_location("publication_readiness", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def test_evaluate_readiness_reports_ready_with_complete_inputs() -> None:
    report = MODULE.evaluate_readiness(
        {
            "error_count": 0,
            "warning_count": 0,
        },
        {
            "matched_samples": 3,
            "mean_concordance": 0.98,
            "mean_non_comparable_rate": 0.02,
            "exact_allele_id_match_count": 2,
        },
        {
            "schema_name": "demo",
            "schema_type": "wgmlst",
            "schema_version": "1",
            "species": "Demo species",
            "source": "fixture",
            "created_at": "2026-03-04T00:00:00+00:00",
        },
    )

    assert report["overall_status"] in {"ready", "ready_with_warnings"}
    assert report["fail_count"] == 0


def test_write_report_creates_json_and_markdown(tmp_path: Path) -> None:
    report = {
        "overall_status": "ready_with_warnings",
        "fail_count": 0,
        "warning_count": 1,
        "checks": [{"name": "schema_qc_warnings", "status": "warn", "detail": "warning"}],
        "schema_qc_summary": {},
        "benchmark_summary": {},
        "schema_manifest": None,
        "max_non_comparable_rate": 0.1,
    }

    json_path, markdown_path = MODULE.write_report(report, tmp_path)

    assert json.loads(json_path.read_text(encoding="utf-8"))["overall_status"] == "ready_with_warnings"
    assert "Publication Readiness" in markdown_path.read_text(encoding="utf-8")
