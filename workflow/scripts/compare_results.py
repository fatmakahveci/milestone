#!/usr/bin/env python3

from __future__ import annotations

import argparse
import sys
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from compare_wgmlst_profiles import compare_profiles


def load_legacy_merged_results(path: str) -> tuple[dict[str, str], dict[str, str]]:
    chew_profile: dict[str, str] = {}
    milestone_profile: dict[str, str] = {}
    with open(path, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            fields = line.replace(" ", "\t").split("\t")
            if len(fields) < 3:
                raise ValueError(f"Invalid legacy comparison row: {line}")
            locus, chew_allele, milestone_allele = fields[:3]
            chew_profile[locus] = chew_allele
            milestone_profile[locus] = milestone_allele
    return chew_profile, milestone_profile


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Legacy chewBBACA-vs-Milestone comparison wrapper built on the new wgMLST comparator."
    )
    parser.add_argument(
        "-f",
        "--file",
        required=True,
        help="Tabular file with locus, chewBBACA allele, and Milestone allele columns.",
    )
    parser.add_argument(
        "-a",
        "--aligner",
        required=True,
        help="Aligner label to display in the report.",
    )
    args = parser.parse_args()

    chew_profile, milestone_profile = load_legacy_merged_results(args.file)
    summary = compare_profiles(chew_profile, milestone_profile)

    print("Legacy compare_results.py now uses the unified wgMLST comparison engine.")
    print(f"Decision: {summary['decision']}")
    print(f"Comparable loci: {summary['comparable_loci']}")
    print(f"Identical loci: {summary['identical_loci']}")
    print(f"Differing loci: {summary['differing_loci']}")
    print(f"Unresolved loci: {summary['unresolved_loci']}")
    print(f"Compared profiles: chewBBACA vs {args.aligner}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
