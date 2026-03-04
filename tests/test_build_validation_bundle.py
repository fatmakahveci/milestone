from __future__ import annotations

import importlib.util
import json
import sys
import zipfile
from pathlib import Path


MODULE_PATH = Path(__file__).resolve().parent.parent / "workflow" / "scripts" / "build_validation_bundle.py"
SPEC = importlib.util.spec_from_file_location("build_validation_bundle", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def test_create_bundle_writes_manifest_and_files(tmp_path: Path) -> None:
    input_file = tmp_path / "result.tsv"
    input_file.write_text("abc\t1\n", encoding="utf-8")
    bundle_path = MODULE.create_bundle([str(input_file)], tmp_path / "bundle.zip", notes="demo")

    with zipfile.ZipFile(bundle_path) as archive:
        manifest = json.loads(archive.read("validation_bundle_manifest.json").decode("utf-8"))
        assert manifest["notes"] == "demo"
        assert archive.read("result.tsv") == b"abc\t1\n"
