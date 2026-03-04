#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from workflow.scripts.build_validation_bundle import create_bundle
from workflow.scripts.import_pubmlst_benchmark_pack import import_benchmark_pack
from workflow.scripts.import_pubmlst_scheme import DEFAULT_BASE_URL, validate_base_url


def freeze_pubmlst_benchmark_pack(
    *,
    scheme_database: str,
    isolate_database: str,
    scheme_id: int,
    isolate_ids: list[int],
    output_dir: Path,
    base_url: str = DEFAULT_BASE_URL,
    species: str | None = None,
    source_name: str = "Frozen PubMLST benchmark pack",
    selection_criteria: str | None = None,
    download_contigs: bool = False,
    zip_name: str = "frozen_benchmark_pack.zip",
) -> dict[str, str]:
    downloaded = import_benchmark_pack(
        scheme_database=scheme_database,
        isolate_database=isolate_database,
        scheme_id=scheme_id,
        isolate_ids=isolate_ids,
        output_dir=output_dir,
        base_url=base_url,
        species=species,
        source_name=source_name,
        selection_criteria=selection_criteria,
        download_contigs=download_contigs,
    )
    manifest_path = output_dir / "benchmark_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["frozen_bundle"] = zip_name
    manifest["frozen"] = True
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    bundle_inputs = [str(path) for path in downloaded if path.exists() and path.is_file()]
    bundle_path = create_bundle(
        bundle_inputs,
        output_dir / zip_name,
        notes=f"Frozen benchmark pack for {manifest['species']} from {manifest['source_snapshot']}",
    )
    return {"manifest_path": str(manifest_path), "bundle_path": str(bundle_path)}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Download and freeze a PubMLST benchmark pack into a ZIP bundle.")
    parser.add_argument("--scheme-database", required=True, help="PubMLST scheme database.")
    parser.add_argument("--isolate-database", required=True, help="PubMLST isolate database.")
    parser.add_argument("--scheme-id", type=int, required=True, help="Scheme identifier.")
    parser.add_argument("--isolate-id", dest="isolate_ids", type=int, action="append", required=True, help="Isolate ID to freeze. Repeat per isolate.")
    parser.add_argument("--output-dir", type=Path, required=True, help="Output directory for the frozen pack.")
    parser.add_argument("--base-url", default=DEFAULT_BASE_URL, help=f"REST API base URL. Default: {DEFAULT_BASE_URL}")
    parser.add_argument("--species", help="Optional species label.")
    parser.add_argument("--source-name", default="Frozen PubMLST benchmark pack", help="Human-readable source label.")
    parser.add_argument("--selection-criteria", default=None, help="Free-text isolate selection criteria.")
    parser.add_argument("--download-contigs", action="store_true", help="Include contigs FASTA for isolates.")
    parser.add_argument("--zip-name", default="frozen_benchmark_pack.zip", help="Output ZIP filename.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    result = freeze_pubmlst_benchmark_pack(
        scheme_database=args.scheme_database,
        isolate_database=args.isolate_database,
        scheme_id=args.scheme_id,
        isolate_ids=list(args.isolate_ids),
        output_dir=args.output_dir,
        base_url=validate_base_url(args.base_url),
        species=args.species,
        source_name=args.source_name,
        selection_criteria=args.selection_criteria,
        download_contigs=args.download_contigs,
        zip_name=args.zip_name,
    )
    print(f"Frozen manifest written to {result['manifest_path']}")
    print(f"Frozen bundle written to {result['bundle_path']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
