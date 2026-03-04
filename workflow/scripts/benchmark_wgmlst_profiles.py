#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import html
import json
import sys
from datetime import date
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from workflow.scripts.compare_wgmlst_profiles import compare_profiles, parse_profile


def load_benchmark_manifest(manifest_path: str | Path) -> dict[str, object]:
    return json.loads(Path(manifest_path).read_text(encoding="utf-8"))


def normalize_profile_name(path: Path) -> str:
    stem = path.stem
    if stem.endswith("_wgmlst"):
        return stem[:-7]
    return stem


def load_profile_dir(profile_dir: str | Path) -> dict[str, Path]:
    return {
        normalize_profile_name(path): path
        for path in sorted(Path(profile_dir).glob("*.tsv"))
        if path.is_file()
    }


def summarize_locus_mismatches(predicted_dir: str | Path, truth_dir: str | Path, shared_samples: list[str]) -> list[dict[str, object]]:
    locus_stats: dict[str, dict[str, int]] = {}
    for sample in shared_samples:
        summary = compare_profiles(
            parse_profile(load_profile_dir(predicted_dir)[sample]),
            parse_profile(load_profile_dir(truth_dir)[sample]),
        )
        for row in summary["rows"]:
            locus = row.locus
            stats = locus_stats.setdefault(
                locus,
                {"different_count": 0, "unresolved_count": 0, "identical_count": 0, "sample_count": 0},
            )
            stats["sample_count"] += 1
            stats[f"{row.status}_count"] += 1
    ranked = []
    for locus, stats in locus_stats.items():
        comparable = stats["different_count"] + stats["identical_count"]
        discordance = stats["different_count"] / comparable if comparable else None
        non_comparable_rate = stats["unresolved_count"] / stats["sample_count"] if stats["sample_count"] else None
        ranked.append(
            {
                "locus": locus,
                **stats,
                "discordance_rate": discordance,
                "non_comparable_rate": non_comparable_rate,
            }
        )
    ranked.sort(
        key=lambda row: (
            row["discordance_rate"] if row["discordance_rate"] is not None else -1,
            row["non_comparable_rate"] if row["non_comparable_rate"] is not None else -1,
            row["locus"],
        ),
        reverse=True,
    )
    return ranked


def _safe_ratio(numerator: int, denominator: int) -> float | None:
    return numerator / denominator if denominator else None


def benchmark_directories(predicted_dir: str | Path, truth_dir: str | Path) -> dict[str, object]:
    predicted = load_profile_dir(predicted_dir)
    truth = load_profile_dir(truth_dir)
    shared = sorted(set(predicted) & set(truth))
    results: list[dict[str, object]] = []

    for sample in shared:
        summary = compare_profiles(parse_profile(predicted[sample]), parse_profile(truth[sample]))
        comparable = int(summary["comparable_loci"])
        unresolved = int(summary["unresolved_loci"])
        concordance = _safe_ratio(int(summary["identical_loci"]), comparable)
        results.append(
            {
                "sample": sample,
                "decision": summary["decision"],
                "comparable_loci": summary["comparable_loci"],
                "identical_loci": summary["identical_loci"],
                "differing_loci": summary["differing_loci"],
                "unresolved_loci": summary["unresolved_loci"],
                "allele_distance": summary["allele_distance"],
                "concordance": concordance,
                "non_comparable_rate": _safe_ratio(unresolved, int(summary["total_loci"])),
                "exact_allele_id_match": summary["decision"] == "indistinguishable",
            }
        )

    overall_concordance_values = [row["concordance"] for row in results if row["concordance"] is not None]
    overall_non_comparable = [row["non_comparable_rate"] for row in results if row["non_comparable_rate"] is not None]
    return {
        "kind": "profile_benchmark",
        "predicted_dir": str(predicted_dir),
        "truth_dir": str(truth_dir),
        "report_date": str(date.today()),
        "matched_samples": len(shared),
        "missing_predicted_samples": sorted(set(truth) - set(predicted)),
        "missing_truth_samples": sorted(set(predicted) - set(truth)),
        "exact_allele_id_match_count": sum(1 for row in results if row["exact_allele_id_match"]),
        "mean_concordance": (
            sum(overall_concordance_values) / len(overall_concordance_values)
            if overall_concordance_values
            else None
        ),
        "mean_non_comparable_rate": (
            sum(overall_non_comparable) / len(overall_non_comparable)
            if overall_non_comparable
            else None
        ),
        "locus_mismatch_summary": summarize_locus_mismatches(predicted_dir, truth_dir, shared),
        "samples": results,
    }


def benchmark_pack(pack_dir: str | Path) -> dict[str, object]:
    pack_path = Path(pack_dir)
    manifest = load_benchmark_manifest(pack_path / "benchmark_manifest.json")
    report = benchmark_directories(pack_path / "predicted", pack_path / "truth")
    report["benchmark_pack"] = manifest
    report["species"] = manifest.get("species")
    report["schema_name"] = manifest.get("schema_name")
    return report


def benchmark_pack_collection(collection_dir: str | Path) -> dict[str, object]:
    collection_path = Path(collection_dir)
    reports = []
    for pack_dir in sorted(path for path in collection_path.iterdir() if path.is_dir()):
        manifest_path = pack_dir / "benchmark_manifest.json"
        if manifest_path.exists():
            reports.append(benchmark_pack(pack_dir))

    concordances = [report["mean_concordance"] for report in reports if report["mean_concordance"] is not None]
    non_comparable_rates = [
        report["mean_non_comparable_rate"]
        for report in reports
        if report["mean_non_comparable_rate"] is not None
    ]
    return {
        "kind": "profile_benchmark_collection",
        "collection_dir": str(collection_dir),
        "species_count": len(reports),
        "mean_concordance": sum(concordances) / len(concordances) if concordances else None,
        "mean_non_comparable_rate": (
            sum(non_comparable_rates) / len(non_comparable_rates) if non_comparable_rates else None
        ),
        "packs": reports,
    }


def write_benchmark_report(report: dict[str, object], output_dir: str | Path) -> tuple[Path, Path]:
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    json_path = output_path / "wgmlst_benchmark_summary.json"
    json_path.write_text(json.dumps(report, indent=2), encoding="utf-8")

    tsv_path = output_path / "wgmlst_benchmark_per_sample.tsv"
    with tsv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "sample",
                "decision",
                "comparable_loci",
                "identical_loci",
                "differing_loci",
                "unresolved_loci",
                "allele_distance",
                "concordance",
                "non_comparable_rate",
                "exact_allele_id_match",
            ]
        )
        for row in report["samples"]:
            writer.writerow(
                [
                    row["sample"],
                    row["decision"],
                    row["comparable_loci"],
                    row["identical_loci"],
                    row["differing_loci"],
                    row["unresolved_loci"],
                    row["allele_distance"],
                    row["concordance"],
                    row["non_comparable_rate"],
                    row["exact_allele_id_match"],
                ]
            )
    _write_benchmark_locus_tsv(report, output_path / "wgmlst_benchmark_per_locus.tsv")
    _write_html_report(report, output_path / "wgmlst_benchmark_report.html", "Milestone wgMLST benchmark")
    return json_path, tsv_path


def write_benchmark_collection_report(report: dict[str, object], output_dir: str | Path) -> tuple[Path, Path]:
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    json_path = output_path / "wgmlst_benchmark_pack_summary.json"
    json_path.write_text(json.dumps(report, indent=2), encoding="utf-8")

    tsv_path = output_path / "wgmlst_benchmark_species.tsv"
    with tsv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "species",
                "schema_name",
                "matched_samples",
                "exact_allele_id_match_count",
                "mean_concordance",
                "mean_non_comparable_rate",
                "source_name",
                "source_url",
                "retrieval_date",
                "source_snapshot",
                "selection_criteria",
            ]
        )
        for pack in report["packs"]:
            manifest = pack.get("benchmark_pack", {})
            writer.writerow(
                [
                    pack.get("species"),
                    pack.get("schema_name"),
                    pack.get("matched_samples"),
                    pack.get("exact_allele_id_match_count"),
                    pack.get("mean_concordance"),
                    pack.get("mean_non_comparable_rate"),
                    manifest.get("source_name"),
                    manifest.get("source_url"),
                    manifest.get("retrieval_date"),
                    manifest.get("source_snapshot"),
                    manifest.get("selection_criteria"),
                ]
            )
    _write_html_report(
        report,
        output_path / "wgmlst_benchmark_collection_report.html",
        "Milestone wgMLST benchmark collection",
    )
    return json_path, tsv_path


def _write_benchmark_locus_tsv(report: dict[str, object], output_path: Path) -> None:
    rows = report.get("locus_mismatch_summary", [])
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "locus",
                "sample_count",
                "identical_count",
                "different_count",
                "unresolved_count",
                "discordance_rate",
                "non_comparable_rate",
            ]
        )
        for row in rows:
            writer.writerow(
                [
                    row["locus"],
                    row["sample_count"],
                    row["identical_count"],
                    row["different_count"],
                    row["unresolved_count"],
                    row["discordance_rate"],
                    row["non_comparable_rate"],
                ]
            )


def _write_html_report(report: dict[str, object], output_path: Path, title: str) -> None:
    if report["kind"] == "profile_benchmark":
        summary_rows = "".join(
            f"<tr><td>{html.escape(str(sample['sample']))}</td><td>{sample['decision']}</td>"
            f"<td>{sample['concordance']}</td><td>{sample['non_comparable_rate']}</td></tr>"
            for sample in report["samples"]
        )
        locus_rows = "".join(
            f"<tr><td>{html.escape(str(row['locus']))}</td><td>{row['different_count']}</td>"
            f"<td>{row['unresolved_count']}</td><td>{row['discordance_rate']}</td></tr>"
            for row in report.get("locus_mismatch_summary", [])[:20]
        )
        body = (
            f"<p>Matched samples: {report['matched_samples']}</p>"
            f"<p>Mean concordance: {report['mean_concordance']}</p>"
            f"<p>Mean non-comparable rate: {report['mean_non_comparable_rate']}</p>"
            "<h2>Per-sample summary</h2>"
            "<table><tr><th>Sample</th><th>Decision</th><th>Concordance</th><th>Non-comparable rate</th></tr>"
            f"{summary_rows}</table>"
            "<h2>Top discordant loci</h2>"
            "<table><tr><th>Locus</th><th>Different</th><th>Unresolved</th><th>Discordance</th></tr>"
            f"{locus_rows}</table>"
        )
    else:
        pack_rows = "".join(
            f"<tr><td>{html.escape(str(pack.get('species')))}</td><td>{pack.get('mean_concordance')}</td>"
            f"<td>{pack.get('mean_non_comparable_rate')}</td></tr>"
            for pack in report.get("packs", [])
        )
        body = (
            f"<p>Species packs: {report['species_count']}</p>"
            f"<p>Mean concordance: {report['mean_concordance']}</p>"
            f"<p>Mean non-comparable rate: {report['mean_non_comparable_rate']}</p>"
            "<table><tr><th>Species</th><th>Mean concordance</th><th>Mean non-comparable rate</th></tr>"
            f"{pack_rows}</table>"
        )
    output_path.write_text(
        (
            "<html><head><meta charset='utf-8'>"
            f"<title>{html.escape(title)}</title>"
            "<style>body{font-family:Arial,sans-serif;margin:2rem;}table{border-collapse:collapse;width:100%;}"
            "th,td{border:1px solid #ccc;padding:0.4rem;text-align:left;}th{background:#f5f5f5;}</style>"
            f"</head><body><h1>{html.escape(title)}</h1>{body}</body></html>"
        ),
        encoding="utf-8",
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Benchmark predicted wgMLST profiles against truth-set wgMLST profiles."
    )
    parser.add_argument("--predicted-dir", help="Directory containing predicted wgMLST TSV files.")
    parser.add_argument("--truth-dir", help="Directory containing truth-set wgMLST TSV files.")
    parser.add_argument(
        "--benchmark-pack-dir",
        help="Directory containing one or more species benchmark packs with manifest/predicted/truth layout.",
    )
    parser.add_argument("--output-dir", required=True, help="Directory for benchmark outputs.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.benchmark_pack_dir:
        report = benchmark_pack_collection(args.benchmark_pack_dir)
        json_path, tsv_path = write_benchmark_collection_report(report, args.output_dir)
        print(f"Species packs benchmarked: {report['species_count']}")
        print(f"Pack summary written to {json_path}")
        print(f"Per-species results written to {tsv_path}")
        return 0
    if not args.predicted_dir or not args.truth_dir:
        raise SystemExit("--predicted-dir and --truth-dir are required unless --benchmark-pack-dir is used.")
    report = benchmark_directories(args.predicted_dir, args.truth_dir)
    json_path, tsv_path = write_benchmark_report(report, args.output_dir)
    print(f"Matched samples: {report['matched_samples']}")
    print(f"Exact allele-ID matches: {report['exact_allele_id_match_count']}")
    print(f"Summary written to {json_path}")
    print(f"Per-sample results written to {tsv_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
