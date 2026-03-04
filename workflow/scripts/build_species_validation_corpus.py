#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import zipfile
from collections import Counter
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from workflow.scripts.freeze_benchmark_pack_collection import discover_pack_dirs, freeze_collection


def load_manifest(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def summarize_collection(collection_dir: str | Path) -> dict[str, object]:
    root = Path(collection_dir)
    pack_dirs = discover_pack_dirs(root)
    manifests = [load_manifest(pack_dir / "benchmark_manifest.json") for pack_dir in pack_dirs]
    schema_types = Counter(str(manifest.get("schema_type", "unknown")) for manifest in manifests)
    return {
        "kind": "species_validation_corpus",
        "collection_dir": str(collection_dir),
        "species_count": len(pack_dirs),
        "species": [
            manifest.get("species", pack_dir.name)
            for pack_dir, manifest in zip(pack_dirs, manifests, strict=False)
        ],
        "schema_types": dict(schema_types),
        "packs": manifests,
    }


def write_outputs(summary: dict[str, object], output_dir: str | Path) -> tuple[Path, Path]:
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    json_path = output_path / "validation_corpus_manifest.json"
    json_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    markdown_path = output_path / "validation_corpus_summary.md"
    lines = [
        "# Species Validation Corpus",
        "",
        f"- Species count: `{summary['species_count']}`",
        f"- Schema types: `{summary['schema_types']}`",
        "",
        "## Included Packs",
        "",
    ]
    for pack in summary["packs"]:
        lines.append(
            f"- `{pack.get('species', 'unknown')}` | schema `{pack.get('schema_name', 'unknown')}` | "
            f"source `{pack.get('source_name', 'unknown')}` | snapshot `{pack.get('source_snapshot', 'unknown')}`"
        )
    lines.append("")
    markdown_path.write_text("\n".join(lines), encoding="utf-8")
    return json_path, markdown_path


def create_corpus_package(collection_dir: str | Path, output_dir: str | Path, zip_name: str) -> dict[str, str]:
    summary = summarize_collection(collection_dir)
    json_path, markdown_path = write_outputs(summary, output_dir)
    frozen = freeze_collection(collection_dir, Path(output_dir) / "frozen_packs")
    zip_path = Path(output_dir) / zip_name
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        archive.write(json_path, arcname=json_path.name)
        archive.write(markdown_path, arcname=markdown_path.name)
        archive.write(Path(frozen["manifest_path"]), arcname="frozen_collection_manifest.json")
        for archive_path in frozen["archives"]:
            path = Path(archive_path)
            archive.write(path, arcname=f"frozen_packs/{path.name}")
    return {
        "summary_path": str(json_path),
        "markdown_path": str(markdown_path),
        "frozen_manifest_path": str(frozen["manifest_path"]),
        "bundle_path": str(zip_path),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build a species-specific validation corpus bundle from benchmark packs.")
    parser.add_argument("--collection-dir", required=True, help="Directory containing benchmark pack subdirectories.")
    parser.add_argument("--output-dir", required=True, help="Directory for corpus outputs.")
    parser.add_argument("--zip-name", default="species_validation_corpus.zip", help="ZIP filename for the corpus package.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    result = create_corpus_package(args.collection_dir, args.output_dir, args.zip_name)
    print(f"Summary written to {result['summary_path']}")
    print(f"Corpus bundle written to {result['bundle_path']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
