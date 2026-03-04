from __future__ import annotations

import importlib.util
import io
import json
import sys
import zipfile
from pathlib import Path


MODULE_PATH = Path(__file__).resolve().parent.parent / "workflow" / "scripts" / "freeze_public_benchmark_pack.py"
SPEC = importlib.util.spec_from_file_location("freeze_public_benchmark_pack", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class FakeResponse(io.BytesIO):
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        self.close()
        return False


def test_freeze_pubmlst_benchmark_pack_writes_bundle_and_manifest(tmp_path: Path, monkeypatch) -> None:
    import_module = MODULE.sys.modules["workflow.scripts.import_pubmlst_benchmark_pack"]

    responses = {
        "https://rest.pubmlst.org/db/demo_seqdef/schemes/7": {"description": "Demo wgMLST"},
        "https://rest.pubmlst.org/db/demo_isolates/isolates/101": {"id": "iso-101"},
        "https://rest.pubmlst.org/db/demo_isolates/isolates/101/schemes/7/allele_ids": {"abc": 1, "adk": 2},
    }

    def fake_urlopen(url: str, timeout=None):
        return FakeResponse(json.dumps(responses[url]).encode("utf-8"))

    monkeypatch.setattr(import_module, "urlopen", fake_urlopen)

    result = MODULE.freeze_pubmlst_benchmark_pack(
        scheme_database="demo_seqdef",
        isolate_database="demo_isolates",
        scheme_id=7,
        isolate_ids=[101],
        output_dir=tmp_path,
        zip_name="bundle.zip",
    )

    manifest = json.loads(Path(result["manifest_path"]).read_text(encoding="utf-8"))
    assert manifest["frozen"] is True
    assert manifest["frozen_bundle"] == "bundle.zip"
    with zipfile.ZipFile(result["bundle_path"]) as archive:
        assert "benchmark_manifest.json" in archive.namelist()
        assert "iso-101_wgmlst.tsv" in archive.namelist()
