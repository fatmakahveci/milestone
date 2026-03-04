#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
from collections import Counter
from dataclasses import asdict, dataclass
from pathlib import Path
import sys

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from wgmlst_utils import (
    describe_quality_check,
    get_allele_id_from_allele_name,
    get_cds_name_from_allele_name,
    is_valid_allele_identifier,
    select_reference_record,
)
from schema_manifest import validate_schema_manifest


@dataclass(frozen=True)
class SchemaIssue:
    locus: str
    issue: str
    detail: str
    severity: str


@dataclass(frozen=True)
class SchemaRecord:
    id: str
    sequence: str


def parse_fasta(path: Path) -> list[SchemaRecord]:
    records: list[SchemaRecord] = []
    current_id: str | None = None
    current_seq: list[str] = []
    with path.open("r", encoding="ascii") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    records.append(SchemaRecord(current_id, "".join(current_seq)))
                current_id = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
    if current_id is not None:
        records.append(SchemaRecord(current_id, "".join(current_seq)))
    return records


def inspect_locus(
    locus_fasta: Path,
    translation_table: int = 11,
    allowed_start_codons: set[str] | None = None,
) -> tuple[dict[str, object], list[SchemaIssue]]:
    records = parse_fasta(locus_fasta)
    locus_name = locus_fasta.stem
    issues: list[SchemaIssue] = []
    if not records:
        issues.append(SchemaIssue(locus_name, "empty_locus", "No FASTA records found.", "error"))
        return {"locus": locus_name, "allele_count": 0, "reference_allele": None}, issues

    record_ids = [record.id for record in records]
    duplicate_ids = [record_id for record_id, count in Counter(record_ids).items() if count > 1]
    for record_id in duplicate_ids:
        issues.append(SchemaIssue(locus_name, "duplicate_allele_id", record_id, "error"))

    duplicate_sequences = [
        sequence
        for sequence, count in Counter(record.sequence for record in records).items()
        if count > 1
    ]
    for sequence in duplicate_sequences:
        duplicate_members = [record.id for record in records if record.sequence == sequence]
        issues.append(
            SchemaIssue(
                locus_name,
                "duplicate_sequence",
                f"Alleles share the same sequence: {', '.join(sorted(duplicate_members))}.",
                "warning",
            )
        )

    for record in records:
        if get_cds_name_from_allele_name(record.id) != locus_name:
            issues.append(
                SchemaIssue(
                    locus_name,
                    "locus_name_mismatch",
                    f"Record '{record.id}' does not match locus file name '{locus_name}'.",
                    "error",
                )
            )
        allele_id = get_allele_id_from_allele_name(record.id)
        if not is_valid_allele_identifier(allele_id):
            issues.append(
                SchemaIssue(
                    locus_name,
                    "invalid_allele_id_format",
                    f"Allele identifier '{allele_id}' contains unsupported characters.",
                    "error",
                )
            )

    reference_record = select_reference_record(records)
    lengths = [len(record.sequence) for record in records]
    mode_length = Counter(lengths).most_common(1)[0][0]

    for record in records:
        qc_result = describe_quality_check(
            record.sequence,
            reference_record.sequence,
            mode_length,
            translation_table=translation_table,
            allowed_start_codons=allowed_start_codons,
        )
        if qc_result["status"] == "LNF":
            issues.append(
                SchemaIssue(
                    locus_name,
                    "orf_invalid",
                    f"Allele '{record.id}' fails coding-sequence validation: {qc_result['reason']}.",
                    "warning",
                )
            )

    return {
        "locus": locus_name,
        "allele_count": len(records),
        "reference_allele": reference_record.id,
        "mode_length": mode_length,
    }, issues


def inspect_schema_dir(
    schema_dir: Path,
    translation_table: int = 11,
    allowed_start_codons: set[str] | None = None,
    manifest_path: Path | None = None,
) -> dict[str, object]:
    locus_fastas = sorted(schema_dir.glob("*.fasta"))
    summaries = []
    issues: list[SchemaIssue] = []
    for locus_fasta in locus_fastas:
        summary, locus_issues = inspect_locus(locus_fasta, translation_table, allowed_start_codons)
        summaries.append(summary)
        issues.extend(locus_issues)

    manifest_payload = None
    if manifest_path is not None:
        try:
            manifest_payload = validate_schema_manifest(manifest_path, schema_dir)
            manifest_translation_table = manifest_payload.get("translation_table")
            if manifest_translation_table is not None and int(manifest_translation_table) != translation_table:
                issues.append(
                    SchemaIssue(
                        locus="schema_manifest",
                        issue="translation_table_mismatch",
                        detail=(
                            f"Manifest translation_table={manifest_translation_table} does not match "
                            f"requested QC table={translation_table}."
                        ),
                        severity="warning",
                    )
                )
        except ValueError as exc:
            issues.append(SchemaIssue("schema_manifest", "manifest_invalid", str(exc), "error"))

    severity_counts = Counter(issue.severity for issue in issues)
    return {
        "schema_dir": str(schema_dir),
        "translation_table": translation_table,
        "locus_count": len(summaries),
        "issue_count": len(issues),
        "error_count": severity_counts.get("error", 0),
        "warning_count": severity_counts.get("warning", 0),
        "manifest_path": str(manifest_path) if manifest_path is not None else None,
        "manifest_valid": manifest_payload is not None,
        "loci": summaries,
        "issues": [asdict(issue) for issue in issues],
    }


def write_report(report: dict[str, object], output_dir: Path) -> tuple[Path, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = output_dir / "schema_qc_summary.json"
    json_path.write_text(json.dumps(report, indent=2), encoding="utf-8")

    tsv_path = output_dir / "schema_qc_issues.tsv"
    with tsv_path.open("w", encoding="utf-8") as handle:
        handle.write("locus\tseverity\tissue\tdetail\n")
        for issue in report["issues"]:
            handle.write(f"{issue['locus']}\t{issue['severity']}\t{issue['issue']}\t{issue['detail']}\n")
    return json_path, tsv_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run schema-level QC checks for Milestone loci.")
    parser.add_argument("--schema-dir", type=Path, required=True, help="Directory containing locus FASTA files.")
    parser.add_argument("--output-dir", type=Path, required=True, help="Directory for QC report outputs.")
    parser.add_argument("--translation-table", type=int, default=11, help="NCBI translation table for ORF validation.")
    parser.add_argument("--manifest", type=Path, default=None, help="Optional schema_manifest.json to validate.")
    parser.add_argument(
        "--allowed-start-codons",
        type=str,
        default=None,
        help="Optional comma-separated start codons overriding translation-table defaults.",
    )
    return parser.parse_args()


def parse_start_codons(raw_value: str | None) -> set[str] | None:
    if not raw_value:
        return None
    values = {item.strip().upper() for item in raw_value.split(",") if item.strip()}
    return values or None


def main() -> int:
    args = parse_args()
    report = inspect_schema_dir(
        args.schema_dir,
        translation_table=args.translation_table,
        allowed_start_codons=parse_start_codons(args.allowed_start_codons),
        manifest_path=args.manifest,
    )
    json_path, tsv_path = write_report(report, args.output_dir)
    print(f"Loci inspected: {report['locus_count']}")
    print(f"Issues found: {report['issue_count']}")
    print(f"Summary written to {json_path}")
    print(f"Issues written to {tsv_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
