#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]


@dataclass(frozen=True)
class ReadinessCheck:
    name: str
    status: str
    detail: str


def load_json(path: str | Path) -> dict[str, object]:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def evaluate_readiness(
    schema_qc_summary: dict[str, object],
    benchmark_summary: dict[str, object],
    schema_manifest: dict[str, object] | None = None,
    max_non_comparable_rate: float = 0.1,
) -> dict[str, object]:
    checks: list[ReadinessCheck] = []

    error_count = int(schema_qc_summary.get("error_count", 0))
    if error_count == 0:
        checks.append(ReadinessCheck("schema_qc_errors", "pass", "Schema QC reported no errors."))
    else:
        checks.append(
            ReadinessCheck(
                "schema_qc_errors",
                "fail",
                f"Schema QC reported {error_count} error(s).",
            )
        )

    warning_count = int(schema_qc_summary.get("warning_count", 0))
    if warning_count == 0:
        checks.append(ReadinessCheck("schema_qc_warnings", "pass", "Schema QC reported no warnings."))
    else:
        checks.append(
            ReadinessCheck(
                "schema_qc_warnings",
                "warn",
                f"Schema QC reported {warning_count} warning(s); these should be discussed in the manuscript.",
            )
        )

    matched_samples = int(benchmark_summary.get("matched_samples", 0))
    if matched_samples > 0:
        checks.append(
            ReadinessCheck(
                "benchmark_sample_count",
                "pass",
                f"Benchmark summary contains {matched_samples} matched sample(s).",
            )
        )
    else:
        checks.append(ReadinessCheck("benchmark_sample_count", "fail", "Benchmark summary contains no matched samples."))

    mean_concordance = benchmark_summary.get("mean_concordance")
    if mean_concordance is None:
        checks.append(
            ReadinessCheck(
                "benchmark_mean_concordance",
                "warn",
                "Mean concordance is unavailable; benchmark interpretation will be limited.",
            )
        )
    else:
        checks.append(
            ReadinessCheck(
                "benchmark_mean_concordance",
                "pass",
                f"Mean concordance is {float(mean_concordance):.4f}.",
            )
        )

    non_comparable_rate = benchmark_summary.get("mean_non_comparable_rate")
    if non_comparable_rate is None:
        checks.append(
            ReadinessCheck(
                "benchmark_non_comparable_rate",
                "warn",
                "Mean non-comparable rate is unavailable.",
            )
        )
    elif float(non_comparable_rate) <= max_non_comparable_rate:
        checks.append(
            ReadinessCheck(
                "benchmark_non_comparable_rate",
                "pass",
                f"Mean non-comparable rate ({float(non_comparable_rate):.4f}) is within the configured threshold.",
            )
        )
    else:
        checks.append(
            ReadinessCheck(
                "benchmark_non_comparable_rate",
                "warn",
                (
                    f"Mean non-comparable rate ({float(non_comparable_rate):.4f}) exceeds the configured "
                    f"threshold ({max_non_comparable_rate:.4f})."
                ),
            )
        )

    exact_match_count = int(benchmark_summary.get("exact_allele_id_match_count", 0))
    checks.append(
        ReadinessCheck(
            "benchmark_exact_allele_id_matches",
            "pass",
            (
                f"Exact allele-ID matches: {exact_match_count}. This metric only reflects shared nomenclature, "
                "not biological equivalence across different schemas or tools."
            ),
        )
    )

    if schema_manifest is None:
        checks.append(
            ReadinessCheck(
                "schema_manifest",
                "warn",
                "No schema manifest was supplied; provenance reporting will be incomplete.",
            )
        )
    else:
        required_fields = {"schema_name", "schema_type", "schema_version", "species", "source", "created_at"}
        missing = sorted(field for field in required_fields if not schema_manifest.get(field))
        if missing:
            checks.append(
                ReadinessCheck(
                    "schema_manifest",
                    "warn",
                    f"Schema manifest is present but missing fields: {', '.join(missing)}.",
                )
            )
        else:
            checks.append(ReadinessCheck("schema_manifest", "pass", "Schema manifest contains core provenance fields."))

    doc_paths = {
        "methods_template": REPO_ROOT / "docs" / "methods-template.md",
        "limitations_doc": REPO_ROOT / "docs" / "limitations.md",
        "citation_file": REPO_ROOT / "CITATION.cff",
        "codemeta_file": REPO_ROOT / "codemeta.json",
    }
    for name, path in doc_paths.items():
        status = "pass" if path.exists() else "warn"
        detail = f"Found {path.name}." if path.exists() else f"Missing required publication-support file: {path.name}."
        checks.append(ReadinessCheck(name, status, detail))

    fail_count = sum(1 for check in checks if check.status == "fail")
    warn_count = sum(1 for check in checks if check.status == "warn")
    overall_status = "ready" if fail_count == 0 and warn_count == 0 else "ready_with_warnings" if fail_count == 0 else "not_ready"
    return {
        "kind": "publication_readiness",
        "overall_status": overall_status,
        "fail_count": fail_count,
        "warning_count": warn_count,
        "checks": [asdict(check) for check in checks],
        "schema_qc_summary": schema_qc_summary,
        "benchmark_summary": benchmark_summary,
        "schema_manifest": schema_manifest,
        "max_non_comparable_rate": max_non_comparable_rate,
    }


def write_report(report: dict[str, object], output_dir: str | Path) -> tuple[Path, Path]:
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    json_path = output_path / "publication_readiness.json"
    json_path.write_text(json.dumps(report, indent=2), encoding="utf-8")

    markdown_path = output_path / "publication_readiness.md"
    lines = [
        "# Publication Readiness",
        "",
        f"- Overall status: `{report['overall_status']}`",
        f"- Fail count: `{report['fail_count']}`",
        f"- Warning count: `{report['warning_count']}`",
        "",
        "## Checks",
        "",
    ]
    for check in report["checks"]:
        lines.append(f"- `{check['status']}` `{check['name']}`: {check['detail']}")
    lines.append("")
    markdown_path.write_text("\n".join(lines), encoding="utf-8")
    return json_path, markdown_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Evaluate whether a Milestone analysis is ready to support publication.")
    parser.add_argument("--schema-qc-summary", required=True, help="Path to schema_qc_summary.json.")
    parser.add_argument("--benchmark-summary", required=True, help="Path to wgmlst_benchmark_summary.json.")
    parser.add_argument("--schema-manifest", help="Optional path to schema_manifest.json.")
    parser.add_argument("--output-dir", required=True, help="Directory for readiness outputs.")
    parser.add_argument(
        "--max-non-comparable-rate",
        type=float,
        default=0.1,
        help="Warn when benchmark mean_non_comparable_rate exceeds this value. Default: 0.1",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    schema_manifest = load_json(args.schema_manifest) if args.schema_manifest else None
    report = evaluate_readiness(
        load_json(args.schema_qc_summary),
        load_json(args.benchmark_summary),
        schema_manifest=schema_manifest,
        max_non_comparable_rate=args.max_non_comparable_rate,
    )
    json_path, markdown_path = write_report(report, args.output_dir)
    print(f"Overall status: {report['overall_status']}")
    print(f"JSON written to {json_path}")
    print(f"Markdown written to {markdown_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
