#!/usr/bin/env python3

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import zipfile
from datetime import datetime, timezone
from pathlib import Path

from publication_readiness import evaluate_readiness, load_json


REPO_ROOT = Path(__file__).resolve().parents[2]


def utc_timestamp() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256sum(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def collect_supporting_files() -> list[Path]:
    candidates = [
        REPO_ROOT / "README.md",
        REPO_ROOT / "LICENSE.md",
        REPO_ROOT / "CITATION.cff",
        REPO_ROOT / "codemeta.json",
        REPO_ROOT / "pyproject.toml",
        REPO_ROOT / "workflow" / "environment.yml",
        REPO_ROOT / "docs" / "methods-template.md",
        REPO_ROOT / "docs" / "limitations.md",
        REPO_ROOT / "docs" / "publication-checklist.md",
    ]
    return [path for path in candidates if path.exists()]


def build_package_manifest(
    included_paths: list[Path],
    readiness_report: dict[str, object],
    title: str | None = None,
    notes: str | None = None,
) -> dict[str, object]:
    return {
        "package_version": "1",
        "title": title or "Milestone publication package",
        "created_at": utc_timestamp(),
        "platform": platform.platform(),
        "python_version": platform.python_version(),
        "notes": notes or "",
        "readiness_overall_status": readiness_report["overall_status"],
        "artifacts": [
            {
                "path": str(path.relative_to(REPO_ROOT)) if path.is_relative_to(REPO_ROOT) else path.name,
                "sha256": sha256sum(path),
                "size_bytes": path.stat().st_size,
            }
            for path in included_paths
        ],
    }


def build_publication_summary(readiness_report: dict[str, object]) -> str:
    lines = [
        "# Publication Package Summary",
        "",
        f"- Overall readiness: `{readiness_report['overall_status']}`",
        f"- Schema QC errors: `{readiness_report['schema_qc_summary'].get('error_count', 0)}`",
        f"- Schema QC warnings: `{readiness_report['schema_qc_summary'].get('warning_count', 0)}`",
        f"- Matched benchmark samples: `{readiness_report['benchmark_summary'].get('matched_samples', 0)}`",
        f"- Mean concordance: `{readiness_report['benchmark_summary'].get('mean_concordance')}`",
        f"- Mean non-comparable rate: `{readiness_report['benchmark_summary'].get('mean_non_comparable_rate')}`",
        "",
        "## Mandatory Reporting Notes",
        "",
        "- Report this work as schema-based wgMLST profiling and discrimination, not as universal strain identification.",
        "- State the schema provenance, version, and validation dataset.",
        "- Describe unresolved loci handling explicitly.",
        "- Do not interpret allele-ID equality across different schemas or tools as guaranteed biological equivalence.",
        "",
        "## Readiness Checks",
        "",
    ]
    for check in readiness_report["checks"]:
        lines.append(f"- `{check['status']}` `{check['name']}`: {check['detail']}")
    lines.append("")
    return "\n".join(lines)


def create_publication_package(
    schema_qc_summary_path: str | Path,
    benchmark_summary_path: str | Path,
    output_path: str | Path,
    schema_manifest_path: str | Path | None = None,
    title: str | None = None,
    notes: str | None = None,
    max_non_comparable_rate: float = 0.1,
) -> Path:
    schema_qc_summary = load_json(schema_qc_summary_path)
    benchmark_summary = load_json(benchmark_summary_path)
    schema_manifest = load_json(schema_manifest_path) if schema_manifest_path else None
    readiness_report = evaluate_readiness(
        schema_qc_summary,
        benchmark_summary,
        schema_manifest=schema_manifest,
        max_non_comparable_rate=max_non_comparable_rate,
    )

    included_paths = [
        Path(schema_qc_summary_path),
        Path(benchmark_summary_path),
    ]
    if schema_manifest_path is not None:
        included_paths.append(Path(schema_manifest_path))
    included_paths.extend(collect_supporting_files())

    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    package_manifest = build_package_manifest(included_paths, readiness_report, title=title, notes=notes)
    summary_markdown = build_publication_summary(readiness_report)

    with zipfile.ZipFile(output, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        archive.writestr("publication_package_manifest.json", json.dumps(package_manifest, indent=2))
        archive.writestr("publication_summary.md", summary_markdown)
        archive.writestr("publication_readiness.json", json.dumps(readiness_report, indent=2))
        for path in included_paths:
            arcname = str(path.relative_to(REPO_ROOT)) if path.is_relative_to(REPO_ROOT) else path.name
            archive.write(path, arcname=arcname)
    return output


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Bundle Milestone publication-support artifacts into a ZIP package.")
    parser.add_argument("--schema-qc-summary", required=True, help="Path to schema_qc_summary.json.")
    parser.add_argument("--benchmark-summary", required=True, help="Path to wgmlst_benchmark_summary.json.")
    parser.add_argument("--schema-manifest", help="Optional schema_manifest.json path.")
    parser.add_argument("--output", required=True, help="Output ZIP path.")
    parser.add_argument("--title", help="Optional package title.")
    parser.add_argument("--notes", help="Optional notes to store in the package manifest.")
    parser.add_argument(
        "--max-non-comparable-rate",
        type=float,
        default=0.1,
        help="Warn in readiness report when mean_non_comparable_rate exceeds this value. Default: 0.1",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    bundle = create_publication_package(
        args.schema_qc_summary,
        args.benchmark_summary,
        args.output,
        schema_manifest_path=args.schema_manifest,
        title=args.title,
        notes=args.notes,
        max_non_comparable_rate=args.max_non_comparable_rate,
    )
    print(f"Publication package written to {bundle}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
