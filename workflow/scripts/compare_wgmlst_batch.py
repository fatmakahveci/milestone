#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from workflow.scripts.compare_wgmlst_profiles import compare_profiles, parse_profile


def load_profiles(profile_paths: list[str]) -> dict[str, dict[str, str]]:
    payload = {}
    for path_text in profile_paths:
        path = Path(path_text)
        payload[path.stem] = parse_profile(path)
    return payload


def compare_batch(profile_paths: list[str]) -> dict[str, object]:
    profiles = load_profiles(profile_paths)
    labels = sorted(profiles)
    pairs = []
    decision_counts = {"different": 0, "indistinguishable": 0, "inconclusive": 0}
    for index, label_a in enumerate(labels):
        for label_b in labels[index + 1 :]:
            summary = compare_profiles(profiles[label_a], profiles[label_b])
            row = {
                "sample_a": label_a,
                "sample_b": label_b,
                "decision": summary["decision"],
                "comparable_loci": summary["comparable_loci"],
                "differing_loci": summary["differing_loci"],
                "unresolved_loci": summary["unresolved_loci"],
                "allele_distance": summary["allele_distance"],
            }
            pairs.append(row)
            decision_counts[summary["decision"]] += 1
    return {
        "kind": "profile_compare_batch",
        "sample_count": len(labels),
        "pair_count": len(pairs),
        "decision_counts": decision_counts,
        "pairs": pairs,
    }


def write_outputs(report: dict[str, object], output_dir: str | Path) -> tuple[Path, Path]:
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    json_path = output_path / "wgmlst_batch_compare_summary.json"
    json_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    tsv_path = output_path / "wgmlst_batch_compare.tsv"
    with tsv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            ["sample_a", "sample_b", "decision", "comparable_loci", "differing_loci", "unresolved_loci", "allele_distance"]
        )
        for pair in report["pairs"]:
            writer.writerow(
                [
                    pair["sample_a"],
                    pair["sample_b"],
                    pair["decision"],
                    pair["comparable_loci"],
                    pair["differing_loci"],
                    pair["unresolved_loci"],
                    pair["allele_distance"],
                ]
            )
    return json_path, tsv_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Batch-compare multiple wgMLST profiles.")
    parser.add_argument("--profiles", nargs="+", required=True, help="Two or more wgMLST profile TSV files.")
    parser.add_argument("--output-dir", required=True, help="Directory for batch comparison outputs.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if len(args.profiles) < 2:
        raise SystemExit("At least two profiles are required.")
    report = compare_batch(args.profiles)
    json_path, tsv_path = write_outputs(report, args.output_dir)
    print(f"Samples compared: {report['sample_count']}")
    print(f"Pairwise decisions: {report['pair_count']}")
    print(f"Summary written to {json_path}")
    print(f"Details written to {tsv_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
