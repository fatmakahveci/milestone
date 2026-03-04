#!/usr/bin/env python3

from __future__ import annotations

import argparse
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from compare_wgmlst_profiles import compare_profiles, parse_profile, write_summary


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Deprecated helper retained as a wrapper around the unified wgMLST comparison engine."
    )
    parser.add_argument("--milestone-profile", required=True, help="Path to Milestone wgMLST TSV profile.")
    parser.add_argument("--chew-profile", required=True, help="Path to chewBBACA wgMLST TSV profile.")
    parser.add_argument("--output-dir", required=True, help="Directory for comparison outputs.")
    args = parser.parse_args()

    milestone_profile = parse_profile(args.milestone_profile)
    chew_profile = parse_profile(args.chew_profile)
    summary = compare_profiles(milestone_profile, chew_profile)
    summary_path, details_path = write_summary(summary, args.output_dir, "milestone", "chewbbaca")

    print("Deprecated comp_milestone_chewie.py redirected to the unified wgMLST comparison engine.")
    print(f"Decision: {summary['decision']}")
    print(f"Summary written to {summary_path}")
    print(f"Details written to {details_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
