from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path


MODULE_PATH = Path(__file__).resolve().parent.parent / "workflow" / "scripts" / "schema_qc.py"
SPEC = importlib.util.spec_from_file_location("schema_qc", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def test_inspect_schema_dir_reports_mismatches_and_invalid_orfs(tmp_path: Path) -> None:
    schema_dir = tmp_path / "schema"
    schema_dir.mkdir()
    (schema_dir / "abc.fasta").write_text(
        ">abc_1\nATGAAATAA\n>xyz_2\nATGAAA\n",
        encoding="ascii",
    )

    report = MODULE.inspect_schema_dir(schema_dir)

    assert report["locus_count"] == 1
    assert report["issue_count"] >= 1
    issue_types = {issue["issue"] for issue in report["issues"]}
    assert "locus_name_mismatch" in issue_types or "orf_invalid" in issue_types


def test_inspect_schema_dir_reports_duplicate_sequences_and_invalid_ids(tmp_path: Path) -> None:
    schema_dir = tmp_path / "schema"
    schema_dir.mkdir()
    (schema_dir / "abc.fasta").write_text(
        ">abc_1\nATGAAATAA\n>abc_bad id\nATGAAATAA\n",
        encoding="ascii",
    )

    report = MODULE.inspect_schema_dir(schema_dir)

    issue_types = {issue["issue"] for issue in report["issues"]}
    assert "duplicate_sequence" in issue_types
    assert "invalid_allele_id_format" in issue_types


def test_inspect_schema_dir_flags_invalid_reference_orf(tmp_path: Path) -> None:
    schema_dir = tmp_path / "schema"
    schema_dir.mkdir()
    (schema_dir / "abc.fasta").write_text(
        ">abc_1\nATGTAACCCTAA\n>abc_2\nATGAAACCCTAA\n",
        encoding="ascii",
    )

    report = MODULE.inspect_schema_dir(schema_dir)

    assert any(
        issue["issue"] == "orf_invalid" and "abc_1" in issue["detail"]
        for issue in report["issues"]
    )


def test_inspect_schema_dir_validates_manifest_and_translation_table(tmp_path: Path) -> None:
    schema_dir = tmp_path / "schema"
    schema_dir.mkdir()
    (schema_dir / "abc.fasta").write_text(">abc_1\nATGAAATAA\n", encoding="ascii")
    manifest = tmp_path / "schema_manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "manifest_version": "1",
                "schema_name": "demo",
                "schema_type": "wgmlst",
                "schema_version": "1",
                "species": "Demo",
                "source": "fixture",
                "source_url": "local",
                "created_at": "2026-03-04T00:00:00+00:00",
                "locus_count": 1,
                "loci": ["abc"],
                "translation_table": 4,
            }
        ),
        encoding="utf-8",
    )

    report = MODULE.inspect_schema_dir(schema_dir, translation_table=11, manifest_path=manifest)

    issue_types = {issue["issue"] for issue in report["issues"]}
    assert "translation_table_mismatch" in issue_types
    assert report["manifest_valid"] is True


def test_write_report_outputs_json_and_tsv(tmp_path: Path) -> None:
    report = {
        "schema_dir": "schema",
        "locus_count": 1,
        "issue_count": 1,
        "error_count": 1,
        "warning_count": 0,
        "loci": [{"locus": "abc", "allele_count": 2, "reference_allele": "abc_1", "mode_length": 9}],
        "issues": [{"locus": "abc", "issue": "duplicate_allele_id", "detail": "abc_1", "severity": "error"}],
    }

    json_path, tsv_path = MODULE.write_report(report, tmp_path / "out")

    assert json.loads(json_path.read_text(encoding="utf-8"))["issue_count"] == 1
    assert "duplicate_allele_id" in tsv_path.read_text(encoding="utf-8")
