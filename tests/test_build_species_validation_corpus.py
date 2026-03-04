from __future__ import annotations

import importlib.util
import json
import sys
import zipfile
from pathlib import Path


MODULE_PATH = Path(__file__).resolve().parent.parent / "workflow" / "scripts" / "build_species_validation_corpus.py"
SPEC = importlib.util.spec_from_file_location("build_species_validation_corpus", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def test_create_corpus_package_writes_manifest_and_bundle(tmp_path: Path) -> None:
    collection = tmp_path / "packs"
    species = collection / "species_a"
    (species / "predicted").mkdir(parents=True)
    (species / "truth").mkdir(parents=True)
    (species / "benchmark_manifest.json").write_text(
        json.dumps(
            {
                "species": "Species A",
                "schema_name": "demo_schema",
                "schema_type": "wgmlst",
                "source_name": "fixture",
                "source_snapshot": "fixture:1",
            }
        ),
        encoding="utf-8",
    )
    (species / "predicted" / "a.tsv").write_text("abc\t1\n", encoding="utf-8")
    (species / "truth" / "a.tsv").write_text("abc\t1\n", encoding="utf-8")

    result = MODULE.create_corpus_package(collection, tmp_path / "out", "corpus.zip")

    summary = json.loads(Path(result["summary_path"]).read_text(encoding="utf-8"))
    assert summary["species_count"] == 1
    with zipfile.ZipFile(result["bundle_path"]) as archive:
        assert "validation_corpus_manifest.json" in archive.namelist()
        assert "frozen_packs/species_a.zip" in archive.namelist()
