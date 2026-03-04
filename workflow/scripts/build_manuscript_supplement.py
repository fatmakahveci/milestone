#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import zipfile
from pathlib import Path


def ensure_dir(path: str | Path) -> Path:
    output = Path(path)
    output.mkdir(parents=True, exist_ok=True)
    return output


def extract_publication_package(package_path: str | Path, output_dir: str | Path) -> Path:
    output = ensure_dir(output_dir)
    with zipfile.ZipFile(package_path) as archive:
        archive.extractall(output)
    return output


def build_supplement_index(extracted_dir: Path) -> str:
    readiness = extracted_dir / "publication_readiness.json"
    readiness_payload = json.loads(readiness.read_text(encoding="utf-8")) if readiness.exists() else {}
    lines = [
        "# Supplement Index",
        "",
        "This directory was generated from a Milestone publication package.",
        "",
        f"- Overall readiness: `{readiness_payload.get('overall_status', 'unknown')}`",
        f"- Fail count: `{readiness_payload.get('fail_count', 'unknown')}`",
        f"- Warning count: `{readiness_payload.get('warning_count', 'unknown')}`",
        "",
        "## Included Core Files",
        "",
    ]
    for name in [
        "publication_summary.md",
        "publication_readiness.json",
        "README.md",
        "docs/methods-template.md",
        "docs/limitations.md",
        "docs/publication-checklist.md",
        "CITATION.cff",
        "codemeta.json",
    ]:
        path = extracted_dir / name
        if path.exists():
            lines.append(f"- `{name}`")
    lines.append("")
    return "\n".join(lines)


def create_supplement(package_path: str | Path, output_dir: str | Path) -> dict[str, str]:
    extracted_dir = extract_publication_package(package_path, output_dir)
    index_path = extracted_dir / "SUPPLEMENT_INDEX.md"
    index_path.write_text(build_supplement_index(extracted_dir), encoding="utf-8")
    return {
        "supplement_dir": str(extracted_dir),
        "index_path": str(index_path),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Extract a Milestone publication package into a manuscript supplement directory.")
    parser.add_argument("--publication-package", required=True, help="Path to the publication package ZIP.")
    parser.add_argument("--output-dir", required=True, help="Directory for the extracted supplement.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    result = create_supplement(args.publication_package, args.output_dir)
    print(f"Supplement directory created at {result['supplement_dir']}")
    print(f"Supplement index written to {result['index_path']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
