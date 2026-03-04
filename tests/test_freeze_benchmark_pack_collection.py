from __future__ import annotations

import importlib.util
import json
import sys
import zipfile
from pathlib import Path


MODULE_PATH = Path(__file__).resolve().parent.parent / "workflow" / "scripts" / "freeze_benchmark_pack_collection.py"
SPEC = importlib.util.spec_from_file_location("freeze_benchmark_pack_collection", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def test_freeze_collection_creates_archives_and_manifest(tmp_path: Path) -> None:
    collection = tmp_path / "packs"
    pack_a = collection / "species_a"
    pack_b = collection / "species_b"
    (pack_a / "predicted").mkdir(parents=True)
    (pack_b / "truth").mkdir(parents=True)
    (pack_a / "benchmark_manifest.json").write_text('{"species":"A"}', encoding="utf-8")
    (pack_b / "benchmark_manifest.json").write_text('{"species":"B"}', encoding="utf-8")
    (pack_a / "predicted" / "a.tsv").write_text("abc\t1\n", encoding="utf-8")
    (pack_b / "truth" / "b.tsv").write_text("abc\t2\n", encoding="utf-8")

    result = MODULE.freeze_collection(collection, tmp_path / "out")

    assert result["pack_count"] == 2
    manifest = json.loads(Path(result["manifest_path"]).read_text(encoding="utf-8"))
    assert manifest["packs"] == ["species_a", "species_b"]
    with zipfile.ZipFile(tmp_path / "out" / "species_a.zip") as archive:
        assert "species_a/predicted/a.tsv" in archive.namelist()
