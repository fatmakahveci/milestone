from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path


MODULE_PATH = Path(__file__).resolve().parent.parent / "workflow" / "scripts" / "compare_wgmlst_matrix.py"
SPEC = importlib.util.spec_from_file_location("compare_wgmlst_matrix", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def test_build_distance_matrix_reports_pairwise_decisions() -> None:
    profiles = {
        "a": {"abc": "1", "adk": "2"},
        "b": {"abc": "1", "adk": "2"},
        "c": {"abc": "1", "adk": "7"},
    }

    labels, rows, comparisons = MODULE.build_distance_matrix(profiles)

    assert labels == ["a", "b", "c"]
    assert rows[0][1] == 0.0
    assert comparisons[0]["decision"] == "indistinguishable"
    assert any(pair["decision"] == "different" for pair in comparisons)


def test_build_clustering_summary_groups_identical_profiles() -> None:
    clustering = MODULE.build_clustering_summary(
        ["a", "b", "c"],
        [
            {"sample_a": "a", "sample_b": "b", "decision": "indistinguishable", "allele_distance": 0.0},
            {"sample_a": "a", "sample_b": "c", "decision": "different", "allele_distance": 0.5},
            {"sample_a": "b", "sample_b": "c", "decision": "different", "allele_distance": 0.5},
        ],
    )

    assert clustering["cluster_count"] == 2
    assert clustering["largest_cluster_size"] == 2


def test_write_outputs_creates_matrix_and_summary(tmp_path: Path) -> None:
    labels = ["a", "b"]
    rows = [
        ["a", 0.0, 0.5],
        ["b", 0.5, 0.0],
    ]
    comparisons = [
        {
            "sample_a": "a",
            "sample_b": "b",
            "decision": "different",
            "comparable_loci": 2,
            "differing_loci": 1,
            "unresolved_loci": 0,
            "allele_distance": 0.5,
        }
    ]

    matrix_path, summary_path = MODULE.write_outputs(labels, rows, comparisons, tmp_path)

    assert "sample\ta\tb" in matrix_path.read_text(encoding="utf-8")
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    assert summary["sample_count"] == 2
    assert summary["pairs"][0]["decision"] == "different"
    assert "clustering" in summary
    assert (tmp_path / "wgmlst_distance_report.html").exists()
