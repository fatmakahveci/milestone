#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import html
import json
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from workflow.scripts.compare_wgmlst_profiles import compare_profiles, parse_profile


def load_profiles(profile_paths: list[str]) -> dict[str, dict[str, str]]:
    profiles: dict[str, dict[str, str]] = {}
    for path_text in profile_paths:
        path = Path(path_text)
        profiles[path.stem] = parse_profile(path)
    return profiles


def build_distance_matrix(
    profiles: dict[str, dict[str, str]],
) -> tuple[list[str], list[list[str | float]], list[dict[str, object]]]:
    labels = sorted(profiles)
    rows: list[list[str | float]] = []
    comparisons: list[dict[str, object]] = []
    for label_a in labels:
        matrix_row: list[str | float] = [label_a]
        for label_b in labels:
            if label_a == label_b:
                matrix_row.append(0.0)
                continue
            summary = compare_profiles(profiles[label_a], profiles[label_b])
            matrix_row.append(summary["allele_distance"] if summary["allele_distance"] is not None else "NA")
            if label_a < label_b:
                comparisons.append(
                    {
                        "sample_a": label_a,
                        "sample_b": label_b,
                        "decision": summary["decision"],
                        "comparable_loci": summary["comparable_loci"],
                        "differing_loci": summary["differing_loci"],
                        "unresolved_loci": summary["unresolved_loci"],
                        "allele_distance": summary["allele_distance"],
                    }
                )
        rows.append(matrix_row)
    return labels, rows, comparisons


def build_clustering_summary(
    labels: list[str],
    comparisons: list[dict[str, object]],
    distance_threshold: float = 0.0,
) -> dict[str, object]:
    adjacency = {label: set() for label in labels}
    for comparison in comparisons:
        distance = comparison["allele_distance"]
        if comparison["decision"] == "indistinguishable" or (
            distance is not None and distance <= distance_threshold and comparison["decision"] != "inconclusive"
        ):
            adjacency[comparison["sample_a"]].add(comparison["sample_b"])
            adjacency[comparison["sample_b"]].add(comparison["sample_a"])

    clusters = []
    visited: set[str] = set()
    for label in labels:
        if label in visited:
            continue
        stack = [label]
        component = []
        while stack:
            current = stack.pop()
            if current in visited:
                continue
            visited.add(current)
            component.append(current)
            stack.extend(sorted(adjacency[current] - visited))
        clusters.append(sorted(component))
    clusters.sort(key=lambda cluster: (-len(cluster), cluster))
    return {
        "cluster_count": len(clusters),
        "largest_cluster_size": max((len(cluster) for cluster in clusters), default=0),
        "singleton_count": sum(1 for cluster in clusters if len(cluster) == 1),
        "clusters": clusters,
    }


def write_outputs(
    labels: list[str],
    rows: list[list[str | float]],
    comparisons: list[dict[str, object]],
    output_dir: str | Path,
    distance_threshold: float = 0.0,
) -> tuple[Path, Path]:
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    matrix_path = output_path / "wgmlst_distance_matrix.tsv"
    with matrix_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", *labels])
        for row in rows:
            writer.writerow(row)

    summary = {
        "kind": "profile_matrix",
        "sample_count": len(labels),
        "pair_count": len(comparisons),
        "pairs": comparisons,
        "clustering": build_clustering_summary(labels, comparisons, distance_threshold=distance_threshold),
    }
    summary_path = output_path / "wgmlst_distance_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    html_path = output_path / "wgmlst_distance_report.html"
    pair_rows = "".join(
        f"<tr><td>{html.escape(str(pair['sample_a']))}</td><td>{html.escape(str(pair['sample_b']))}</td>"
        f"<td>{pair['decision']}</td><td>{pair['allele_distance']}</td></tr>"
        for pair in comparisons
    )
    html_path.write_text(
        (
            "<html><head><meta charset='utf-8'><title>Milestone wgMLST distance matrix</title>"
            "<style>body{font-family:Arial,sans-serif;margin:2rem;}table{border-collapse:collapse;width:100%;}"
            "th,td{border:1px solid #ccc;padding:0.4rem;text-align:left;}th{background:#f5f5f5;}</style></head>"
            "<body><h1>Milestone wgMLST distance matrix</h1>"
            f"<p>Samples: {len(labels)}</p>"
            f"<p>Clusters: {summary['clustering']['cluster_count']}</p>"
            "<table><tr><th>Sample A</th><th>Sample B</th><th>Decision</th><th>Allele distance</th></tr>"
            f"{pair_rows}</table></body></html>"
        ),
        encoding="utf-8",
    )
    return matrix_path, summary_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a pairwise wgMLST allele-distance matrix from multiple profiles."
    )
    parser.add_argument(
        "--profiles",
        nargs="+",
        required=True,
        help="Two or more wgMLST TSV profiles.",
    )
    parser.add_argument("--output-dir", required=True, help="Directory for matrix outputs.")
    parser.add_argument(
        "--distance-threshold",
        type=float,
        default=0.0,
        help="Threshold used to build connected clusters from pairwise allele distances. Default: 0.0",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if len(args.profiles) < 2:
        raise SystemExit("At least two profiles are required.")
    profiles = load_profiles(args.profiles)
    labels, rows, comparisons = build_distance_matrix(profiles)
    matrix_path, summary_path = write_outputs(
        labels,
        rows,
        comparisons,
        args.output_dir,
        distance_threshold=args.distance_threshold,
    )
    print(f"Samples compared: {len(labels)}")
    print(f"Pairwise comparisons: {len(comparisons)}")
    print(f"Matrix written to {matrix_path}")
    print(f"Summary written to {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
