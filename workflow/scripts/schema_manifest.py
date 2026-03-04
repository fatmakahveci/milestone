#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


REQUIRED_FIELDS = {
    "manifest_version",
    "schema_name",
    "schema_type",
    "species",
    "source",
    "source_url",
    "created_at",
    "locus_count",
    "loci",
}


def utc_timestamp() -> str:
    return datetime.now(timezone.utc).isoformat()


def discover_locus_files(schema_dir: str | Path) -> list[Path]:
    schema_path = Path(schema_dir)
    return sorted(
        path
        for path in schema_path.glob("*.fasta")
        if path.is_file()
    )


def locus_names(schema_dir: str | Path) -> list[str]:
    return [path.stem for path in discover_locus_files(schema_dir)]


def build_schema_manifest(
    schema_dir: str | Path,
    *,
    schema_name: str,
    schema_type: str,
    species: str,
    source: str,
    source_url: str,
    schema_version: str = "1",
    notes: str | None = None,
) -> dict[str, object]:
    loci = locus_names(schema_dir)
    manifest = {
        "manifest_version": "1",
        "schema_name": schema_name,
        "schema_type": schema_type,
        "schema_version": schema_version,
        "species": species,
        "source": source,
        "source_url": source_url,
        "created_at": utc_timestamp(),
        "locus_count": len(loci),
        "loci": loci,
    }
    if notes:
        manifest["notes"] = notes
    return manifest


def write_schema_manifest(manifest: dict[str, object], output_path: str | Path) -> Path:
    path = Path(output_path)
    path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return path


def validate_schema_manifest(manifest_path: str | Path, schema_dir: str | Path | None = None) -> dict[str, object]:
    path = Path(manifest_path)
    payload = json.loads(path.read_text(encoding="utf-8"))
    missing = sorted(REQUIRED_FIELDS - payload.keys())
    if missing:
        raise ValueError(f"Manifest is missing required fields: {', '.join(missing)}")

    loci = payload["loci"]
    if not isinstance(loci, list) or not all(isinstance(item, str) and item for item in loci):
        raise ValueError("Manifest field 'loci' must be a non-empty list of locus names.")
    if len(set(loci)) != len(loci):
        raise ValueError("Manifest locus names must be unique.")
    if payload["locus_count"] != len(loci):
        raise ValueError("Manifest locus_count does not match locus list length.")

    if schema_dir is not None:
        discovered = locus_names(schema_dir)
        if sorted(discovered) != sorted(loci):
            raise ValueError("Manifest loci do not match the schema directory FASTA files.")

    return payload


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create or validate a Milestone schema manifest.")
    parser.add_argument("--schema-dir", required=True, help="Directory containing one FASTA per locus.")
    parser.add_argument("--manifest-path", required=True, help="Path to write or validate the manifest JSON.")
    parser.add_argument("--validate", action="store_true", help="Validate an existing manifest instead of creating one.")
    parser.add_argument("--schema-name", help="Human-readable schema name.")
    parser.add_argument("--schema-type", default="wgmlst", help="Schema type, e.g. wgmlst or cgmlst.")
    parser.add_argument("--species", default="unspecified", help="Species label for the schema.")
    parser.add_argument("--source", default="custom", help="Schema provenance source.")
    parser.add_argument("--source-url", default="local", help="Source URL or identifier.")
    parser.add_argument("--schema-version", default="1", help="Schema version string.")
    parser.add_argument("--notes", help="Optional free-text notes.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.validate:
        validate_schema_manifest(args.manifest_path, args.schema_dir)
        print(f"Validated schema manifest: {args.manifest_path}")
        return 0

    if not args.schema_name:
        raise SystemExit("--schema-name is required unless --validate is used.")

    manifest = build_schema_manifest(
        args.schema_dir,
        schema_name=args.schema_name,
        schema_type=args.schema_type,
        species=args.species,
        source=args.source,
        source_url=args.source_url,
        schema_version=args.schema_version,
        notes=args.notes,
    )
    output_path = write_schema_manifest(manifest, args.manifest_path)
    print(f"Wrote schema manifest to {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
