from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

MODULE_PATH = Path(__file__).resolve().parent.parent / "workflow" / "scripts" / "compare_wgmlst_batch.py"
SPEC = importlib.util.spec_from_file_location("compare_wgmlst_batch", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def test_compare_batch_summarizes_decisions(tmp_path: Path) -> None:
    profile_a = tmp_path / "a_wgmlst.tsv"
    profile_b = tmp_path / "b_wgmlst.tsv"
    profile_a.write_text("abc\t1\nadk\t2\n", encoding="utf-8")
    profile_b.write_text("abc\t1\nadk\t7\n", encoding="utf-8")

    report = MODULE.compare_batch([str(profile_a), str(profile_b)])

    assert report["pair_count"] == 1
    assert report["decision_counts"]["different"] == 1


def test_write_outputs_creates_json_and_tsv(tmp_path: Path) -> None:
    report = {
        "kind": "profile_compare_batch",
        "sample_count": 2,
        "pair_count": 1,
        "decision_counts": {"different": 1, "indistinguishable": 0, "inconclusive": 0},
        "pairs": [
            {
                "sample_a": "a",
                "sample_b": "b",
                "decision": "different",
                "comparable_loci": 2,
                "differing_loci": 1,
                "unresolved_loci": 0,
                "allele_distance": 0.5,
            }
        ],
    }

    json_path, tsv_path = MODULE.write_outputs(report, tmp_path)

    assert json.loads(json_path.read_text(encoding="utf-8"))["pair_count"] == 1
    assert "a\tb\tdifferent" in tsv_path.read_text(encoding="utf-8")
