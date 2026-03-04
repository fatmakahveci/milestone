from __future__ import annotations

import importlib.util
import json
import sys
import zipfile
from pathlib import Path

READY_MODULE_PATH = Path(__file__).resolve().parent.parent / "workflow" / "scripts" / "publication_readiness.py"
READY_SPEC = importlib.util.spec_from_file_location("publication_readiness", READY_MODULE_PATH)
READY_MODULE = importlib.util.module_from_spec(READY_SPEC)
assert READY_SPEC.loader is not None
sys.modules[READY_SPEC.name] = READY_MODULE
READY_SPEC.loader.exec_module(READY_MODULE)

MODULE_PATH = Path(__file__).resolve().parent.parent / "workflow" / "scripts" / "build_publication_package.py"
SPEC = importlib.util.spec_from_file_location("build_publication_package", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def test_create_publication_package_writes_zip_with_manifest_and_summary(tmp_path: Path) -> None:
    schema_qc = tmp_path / "schema_qc_summary.json"
    benchmark = tmp_path / "wgmlst_benchmark_summary.json"
    manifest = tmp_path / "schema_manifest.json"
    schema_qc.write_text(json.dumps({"error_count": 0, "warning_count": 0}), encoding="utf-8")
    benchmark.write_text(
        json.dumps(
            {
                "matched_samples": 2,
                "mean_concordance": 0.95,
                "mean_non_comparable_rate": 0.01,
                "exact_allele_id_match_count": 1,
            }
        ),
        encoding="utf-8",
    )
    manifest.write_text(
        json.dumps(
            {
                "schema_name": "demo",
                "schema_type": "wgmlst",
                "schema_version": "1",
                "species": "Demo species",
                "source": "fixture",
                "created_at": "2026-03-04T00:00:00+00:00",
            }
        ),
        encoding="utf-8",
    )

    bundle = MODULE.create_publication_package(schema_qc, benchmark, tmp_path / "publication.zip", manifest)

    with zipfile.ZipFile(bundle) as archive:
        package_manifest = json.loads(archive.read("publication_package_manifest.json").decode("utf-8"))
        readiness = json.loads(archive.read("publication_readiness.json").decode("utf-8"))
        summary = archive.read("publication_summary.md").decode("utf-8")

        assert package_manifest["readiness_overall_status"] in {"ready", "ready_with_warnings"}
        assert readiness["fail_count"] == 0
        assert "Publication Package Summary" in summary
