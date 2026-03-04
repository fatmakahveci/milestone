#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path

NON_COMPARABLE_CALLS = {"LNF", "ASM", "ALM", "NIPH", "NIPHEM"}


@dataclass(frozen=True)
class ComparisonRow:
    locus: str
    allele_a: str
    allele_b: str
    status: str


def parse_profile(path: str | Path) -> dict[str, str]:
    profile = {}
    with Path(path).open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 2:
                raise ValueError(f"Invalid wgMLST row: {line}")
            if fields[0] in profile:
                raise ValueError(f"Duplicate locus in wgMLST profile: {fields[0]}")
            profile[fields[0]] = fields[1]
    return profile


def is_comparable_call(allele: str) -> bool:
    stripped = allele.strip()
    if not stripped:
        return False
    if stripped in NON_COMPARABLE_CALLS:
        return False
    if "PLOT" in stripped:
        return False
    return True


def compare_profiles(profile_a: dict[str, str], profile_b: dict[str, str]) -> dict[str, object]:
    rows: list[ComparisonRow] = []
    identical_loci = 0
    differing_loci = 0
    comparable_loci = 0
    unresolved_loci = 0
    sample_a_only_calls = 0
    sample_b_only_calls = 0

    for locus in sorted(set(profile_a) | set(profile_b)):
        allele_a = profile_a.get(locus, "MISSING")
        allele_b = profile_b.get(locus, "MISSING")
        call_a = is_comparable_call(allele_a)
        call_b = is_comparable_call(allele_b)

        if call_a and call_b:
            comparable_loci += 1
            if allele_a == allele_b:
                status = "identical"
                identical_loci += 1
            else:
                status = "different"
                differing_loci += 1
        else:
            status = "unresolved"
            unresolved_loci += 1
            if call_a and not call_b:
                sample_a_only_calls += 1
            elif call_b and not call_a:
                sample_b_only_calls += 1

        rows.append(ComparisonRow(locus=locus, allele_a=allele_a, allele_b=allele_b, status=status))

    decision = classify_comparison(
        comparable_loci=comparable_loci,
        differing_loci=differing_loci,
        unresolved_loci=unresolved_loci,
    )
    allele_distance = differing_loci / comparable_loci if comparable_loci else None
    similarity = identical_loci / comparable_loci if comparable_loci else None
    return {
        "decision": decision,
        "total_loci": len(rows),
        "comparable_loci": comparable_loci,
        "identical_loci": identical_loci,
        "differing_loci": differing_loci,
        "unresolved_loci": unresolved_loci,
        "sample_a_only_calls": sample_a_only_calls,
        "sample_b_only_calls": sample_b_only_calls,
        "allele_distance": allele_distance,
        "allele_similarity": similarity,
        "rows": rows,
    }


def classify_comparison(comparable_loci: int, differing_loci: int, unresolved_loci: int) -> str:
    if differing_loci > 0:
        return "different"
    if comparable_loci == 0:
        return "inconclusive"
    if unresolved_loci > 0:
        return "inconclusive"
    return "indistinguishable"


def write_summary(summary: dict[str, object], output_dir: str | Path, label_a: str, label_b: str) -> tuple[Path, Path]:
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    details_path = output_path / "wgmlst_profile_comparison.tsv"
    with details_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["locus", label_a, label_b, "status"])
        for row in summary["rows"]:
            row_obj = row if isinstance(row, ComparisonRow) else ComparisonRow(**row)
            writer.writerow([row_obj.locus, row_obj.allele_a, row_obj.allele_b, row_obj.status])

    serializable = {
        key: [asdict(row) for row in value] if key == "rows" else value
        for key, value in summary.items()
    }
    summary_path = output_path / "wgmlst_profile_summary.json"
    summary_path.write_text(json.dumps(serializable, indent=2), encoding="utf-8")
    return summary_path, details_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare two wgMLST profiles and report whether they are distinguishable."
    )
    parser.add_argument("--profile-a", required=True, help="Path to first wgMLST TSV profile.")
    parser.add_argument("--profile-b", required=True, help="Path to second wgMLST TSV profile.")
    parser.add_argument("--output-dir", required=True, help="Directory for comparison outputs.")
    parser.add_argument("--label-a", default="sample_a", help="Column label for first profile.")
    parser.add_argument("--label-b", default="sample_b", help="Column label for second profile.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    profile_a = parse_profile(args.profile_a)
    profile_b = parse_profile(args.profile_b)
    summary = compare_profiles(profile_a, profile_b)
    summary_path, details_path = write_summary(summary, args.output_dir, args.label_a, args.label_b)
    print(f"Decision: {summary['decision']}")
    print(f"Comparable loci: {summary['comparable_loci']}")
    print(f"Differing loci: {summary['differing_loci']}")
    print(f"Summary written to {summary_path}")
    print(f"Details written to {details_path}")


if __name__ == "__main__":
    main()
