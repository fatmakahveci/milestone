from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

MODULE_PATH = Path(__file__).resolve().parent.parent / "workflow" / "scripts" / "schema_manifest.py"
SPEC = importlib.util.spec_from_file_location("schema_manifest", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def test_build_schema_manifest_tracks_loci_and_metadata(tmp_path: Path) -> None:
    (tmp_path / "abc.fasta").write_text(">abc_1\nATG\n", encoding="utf-8")
    (tmp_path / "adk.fasta").write_text(">adk_1\nATG\n", encoding="utf-8")

    manifest = MODULE.build_schema_manifest(
        tmp_path,
        schema_name="demo_schema",
        schema_type="wgmlst",
        species="Escherichia coli",
        source="custom",
        source_url="local",
        schema_version="2026.1",
    )

    assert manifest["locus_count"] == 2
    assert manifest["loci"] == ["abc", "adk"]
    assert manifest["schema_name"] == "demo_schema"


def test_validate_schema_manifest_rejects_directory_mismatch(tmp_path: Path) -> None:
    (tmp_path / "abc.fasta").write_text(">abc_1\nATG\n", encoding="utf-8")
    manifest_path = tmp_path / "scheme_manifest.json"
    manifest_path.write_text(
        json.dumps(
            {
                "manifest_version": "1",
                "schema_name": "demo",
                "schema_type": "wgmlst",
                "species": "Test species",
                "source": "custom",
                "source_url": "local",
                "created_at": "2026-03-04T00:00:00+00:00",
                "locus_count": 2,
                "loci": ["abc", "adk"],
            }
        ),
        encoding="utf-8",
    )

    try:
        MODULE.validate_schema_manifest(manifest_path, tmp_path)
    except ValueError as exc:
        assert "do not match" in str(exc)
    else:
        raise AssertionError("Expected manifest validation mismatch to raise ValueError")
