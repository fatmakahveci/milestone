from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path


MODULE_PATH = Path(__file__).resolve().parent.parent / "workflow" / "scripts" / "build_public_benchmark_workspace.py"
SPEC = importlib.util.spec_from_file_location("build_public_benchmark_workspace", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def test_build_workspace_creates_manifest_and_corpus(tmp_path: Path, monkeypatch) -> None:
    config = {
        "base_url": "https://rest.pubmlst.org",
        "species": [
            {
                "slug": "species_a",
                "species": "Species A",
                "scheme_database": "demo_seqdef",
                "isolate_database": "demo_isolates",
                "scheme_id": 7,
                "isolate_ids": [101],
            }
        ],
    }

    def fake_freeze_pubmlst_benchmark_pack(**kwargs):
        output_dir = kwargs["output_dir"]
        output_dir.mkdir(parents=True, exist_ok=True)
        manifest = output_dir / "benchmark_manifest.json"
        manifest.write_text(
            json.dumps(
                {
                    "species": kwargs["species"],
                    "schema_name": "demo",
                    "schema_type": "wgmlst",
                    "source_name": kwargs["source_name"],
                    "source_snapshot": "fixture:1",
                }
            ),
            encoding="utf-8",
        )
        return {"manifest_path": str(manifest), "bundle_path": str(output_dir / "bundle.zip")}

    monkeypatch.setattr(MODULE, "freeze_pubmlst_benchmark_pack", fake_freeze_pubmlst_benchmark_pack)

    result = MODULE.build_workspace(config, tmp_path / "workspace")

    assert result["species_count"] == 1
    assert Path(tmp_path / "workspace" / "public_benchmark_workspace.json").exists()
    assert Path(result["corpus"]["bundle_path"]).exists()
