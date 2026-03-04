#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import sys
import zipfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from workflow.scripts.build_validation_bundle import build_manifest


def discover_pack_dirs(collection_dir: str | Path) -> list[Path]:
    root = Path(collection_dir)
    return sorted(path for path in root.iterdir() if path.is_dir() and (path / "benchmark_manifest.json").exists())


def freeze_pack(pack_dir: Path, output_dir: Path) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    archive_path = output_dir / f"{pack_dir.name}.zip"
    with zipfile.ZipFile(archive_path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(pack_dir.rglob("*")):
            if path.is_file():
                archive.write(path, arcname=f"{pack_dir.name}/{path.relative_to(pack_dir)}")
    return archive_path


def freeze_collection(collection_dir: str | Path, output_dir: str | Path) -> dict[str, object]:
    collection_path = Path(collection_dir)
    output_path = Path(output_dir)
    pack_dirs = discover_pack_dirs(collection_path)
    archives = [freeze_pack(pack_dir, output_path) for pack_dir in pack_dirs]
    manifest = build_manifest([str(path) for path in archives], notes="Frozen benchmark pack collection")
    manifest["pack_count"] = len(pack_dirs)
    manifest["collection_dir"] = str(collection_dir)
    manifest["packs"] = [pack_dir.name for pack_dir in pack_dirs]
    manifest_path = output_path / "frozen_collection_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return {
        "pack_count": len(pack_dirs),
        "archives": [str(path) for path in archives],
        "manifest_path": str(manifest_path),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Freeze a directory of benchmark packs into individual ZIP archives.")
    parser.add_argument("--collection-dir", required=True, help="Directory containing species benchmark pack subdirectories.")
    parser.add_argument("--output-dir", required=True, help="Directory for frozen ZIP archives.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    result = freeze_collection(args.collection_dir, args.output_dir)
    print(f"Frozen packs: {result['pack_count']}")
    print(f"Collection manifest written to {result['manifest_path']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
