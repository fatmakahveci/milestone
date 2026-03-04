#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import platform
import zipfile
from datetime import datetime, timezone
from pathlib import Path


def utc_timestamp() -> str:
    return datetime.now(timezone.utc).isoformat()


def build_manifest(inputs: list[str], notes: str | None = None) -> dict[str, object]:
    return {
        "bundle_version": "1",
        "created_at": utc_timestamp(),
        "platform": platform.platform(),
        "python_version": platform.python_version(),
        "inputs": inputs,
        "notes": notes or "",
    }


def create_bundle(inputs: list[str], output_path: str | Path, notes: str | None = None) -> Path:
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    manifest = build_manifest(inputs, notes)
    with zipfile.ZipFile(output, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        archive.writestr("validation_bundle_manifest.json", json.dumps(manifest, indent=2))
        for input_path in inputs:
            path = Path(input_path)
            if path.exists() and path.is_file():
                archive.write(path, arcname=path.name)
    return output


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create a reproducible validation bundle ZIP.")
    parser.add_argument("--inputs", nargs="+", required=True, help="Input files to include in the bundle.")
    parser.add_argument("--output", required=True, help="Path for the output ZIP file.")
    parser.add_argument("--notes", help="Optional notes stored in the bundle manifest.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    bundle_path = create_bundle(args.inputs, args.output, args.notes)
    print(f"Validation bundle written to {bundle_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
