from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path


MODULE_PATH = Path(__file__).resolve().parent.parent / "workflow" / "scripts" / "compare_wgmlst_profiles.py"
MODULE_SPEC = importlib.util.spec_from_file_location("compare_wgmlst_profiles", MODULE_PATH)
MODULE = importlib.util.module_from_spec(MODULE_SPEC)
assert MODULE_SPEC.loader is not None
sys.modules[MODULE_SPEC.name] = MODULE
MODULE_SPEC.loader.exec_module(MODULE)


def test_compare_profiles_marks_different_when_called_locus_differs() -> None:
    summary = MODULE.compare_profiles(
        {"abcZ": "12", "adk": "4"},
        {"abcZ": "12", "adk": "7"},
    )

    assert summary["decision"] == "different"
    assert summary["comparable_loci"] == 2
    assert summary["differing_loci"] == 1


def test_compare_profiles_marks_indistinguishable_when_all_calls_match() -> None:
    summary = MODULE.compare_profiles(
        {"abcZ": "12", "adk": "4"},
        {"abcZ": "12", "adk": "4"},
    )

    assert summary["decision"] == "indistinguishable"
    assert summary["unresolved_loci"] == 0
    assert summary["allele_distance"] == 0


def test_compare_profiles_marks_inconclusive_when_missing_calls_exist() -> None:
    summary = MODULE.compare_profiles(
        {"abcZ": "12", "adk": "LNF"},
        {"abcZ": "12", "adk": "4"},
    )

    assert summary["decision"] == "inconclusive"
    assert summary["identical_loci"] == 1
    assert summary["sample_b_only_calls"] == 1


def test_parse_and_write_summary_round_trip(tmp_path: Path) -> None:
    profile_a = tmp_path / "a.tsv"
    profile_b = tmp_path / "b.tsv"
    profile_a.write_text("abcZ\t12\nadk\t4\n", encoding="utf-8")
    profile_b.write_text("abcZ\t12\nadk\t7\n", encoding="utf-8")

    summary = MODULE.compare_profiles(MODULE.parse_profile(profile_a), MODULE.parse_profile(profile_b))
    summary_path, details_path = MODULE.write_summary(summary, tmp_path / "out", "A", "B")

    loaded = json.loads(summary_path.read_text(encoding="utf-8"))
    details = details_path.read_text(encoding="utf-8")

    assert loaded["decision"] == "different"
    assert loaded["differing_loci"] == 1
    assert "locus\tA\tB\tstatus" in details
    assert "adk\t4\t7\tdifferent" in details


def test_parse_profile_rejects_duplicate_loci(tmp_path: Path) -> None:
    profile = tmp_path / "dup.tsv"
    profile.write_text("abcZ\t12\nabcZ\t7\n", encoding="utf-8")

    try:
        MODULE.parse_profile(profile)
    except ValueError as exc:
        assert "Duplicate locus" in str(exc)
    else:
        raise AssertionError("Expected duplicate loci to raise ValueError")
