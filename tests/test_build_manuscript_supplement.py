from __future__ import annotations

import importlib.util
import json
import sys
import zipfile
from pathlib import Path


READY_PATH = Path(__file__).resolve().parent.parent / "workflow" / "scripts" / "publication_readiness.py"
READY_SPEC = importlib.util.spec_from_file_location("publication_readiness", READY_PATH)
READY_MODULE = importlib.util.module_from_spec(READY_SPEC)
assert READY_SPEC.loader is not None
sys.modules[READY_SPEC.name] = READY_MODULE
READY_SPEC.loader.exec_module(READY_MODULE)

PACKAGE_PATH = Path(__file__).resolve().parent.parent / "workflow" / "scripts" / "build_publication_package.py"
PACKAGE_SPEC = importlib.util.spec_from_file_location("build_publication_package", PACKAGE_PATH)
PACKAGE_MODULE = importlib.util.module_from_spec(PACKAGE_SPEC)
assert PACKAGE_SPEC.loader is not None
sys.modules[PACKAGE_SPEC.name] = PACKAGE_MODULE
PACKAGE_SPEC.loader.exec_module(PACKAGE_MODULE)

MODULE_PATH = Path(__file__).resolve().parent.parent / "workflow" / "scripts" / "build_manuscript_supplement.py"
SPEC = importlib.util.spec_from_file_location("build_manuscript_supplement", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def test_create_supplement_extracts_package_and_index(tmp_path: Path) -> None:
    schema_qc = tmp_path / "schema_qc_summary.json"
    benchmark = tmp_path / "wgmlst_benchmark_summary.json"
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
    package = PACKAGE_MODULE.create_publication_package(schema_qc, benchmark, tmp_path / "publication.zip")

    result = MODULE.create_supplement(package, tmp_path / "supplement")

    index = Path(result["index_path"]).read_text(encoding="utf-8")
    assert "Supplement Index" in index
    assert (tmp_path / "supplement" / "publication_summary.md").exists()
