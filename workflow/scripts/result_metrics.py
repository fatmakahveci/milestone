from __future__ import annotations

import json
from pathlib import Path


NON_COMPARABLE = {"LNF", "ASM", "ALM", "NIPH", "NIPHEM"}


def summarize_wgmlst_profile(path: str | Path) -> dict[str, float | int | str]:
    rows = []
    with Path(path).open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            locus, allele = line.split("\t")[:2]
            rows.append((locus, allele))

    total = len(rows)
    unresolved = sum(1 for _, allele in rows if allele in NON_COMPARABLE or "PLOT" in allele)
    novel = sum(1 for _, allele in rows if allele.startswith("NEW_") or allele.isdigit() and int(allele) > 1)
    reference = sum(1 for _, allele in rows if allele == "1")
    assigned = total - unresolved
    return {
        "kind": "wgmlst_profile",
        "total_loci": total,
        "assigned_loci": assigned,
        "unresolved_loci": unresolved,
        "novel_calls": novel,
        "reference_calls": reference,
    }


def summarize_comparison_summary(path: str | Path) -> dict[str, float | int | str | None]:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    return {
        "kind": "profile_compare",
        "decision": payload["decision"],
        "comparable_loci": payload["comparable_loci"],
        "differing_loci": payload["differing_loci"],
        "unresolved_loci": payload["unresolved_loci"],
        "allele_distance": payload.get("allele_distance"),
    }


def summarize_schema_qc(path: str | Path) -> dict[str, float | int | str]:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    return {
        "kind": "schema_qc",
        "locus_count": payload["locus_count"],
        "issue_count": payload["issue_count"],
        "error_count": payload["error_count"],
        "warning_count": payload["warning_count"],
    }


def summarize_distance_summary(path: str | Path) -> dict[str, float | int | str]:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    different = sum(1 for pair in payload["pairs"] if pair["decision"] == "different")
    inconclusive = sum(1 for pair in payload["pairs"] if pair["decision"] == "inconclusive")
    return {
        "kind": "profile_matrix",
        "sample_count": payload["sample_count"],
        "pair_count": payload["pair_count"],
        "different_pairs": different,
        "inconclusive_pairs": inconclusive,
        "cluster_count": payload.get("clustering", {}).get("cluster_count"),
    }


def summarize_benchmark_summary(path: str | Path) -> dict[str, float | int | str | None]:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    return {
        "kind": "profile_benchmark",
        "matched_samples": payload["matched_samples"],
        "exact_allele_id_match_count": payload["exact_allele_id_match_count"],
        "mean_concordance": payload["mean_concordance"],
        "mean_non_comparable_rate": payload.get("mean_non_comparable_rate"),
    }


def summarize_benchmark_pack_summary(path: str | Path) -> dict[str, float | int | str | None]:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    return {
        "kind": "profile_benchmark_collection",
        "species_count": payload["species_count"],
        "mean_concordance": payload["mean_concordance"],
        "mean_non_comparable_rate": payload.get("mean_non_comparable_rate"),
    }


def summarize_batch_compare_summary(path: str | Path) -> dict[str, float | int | str | None]:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    decision_counts = payload.get("decision_counts", {})
    return {
        "kind": "profile_compare_batch",
        "sample_count": payload["sample_count"],
        "pair_count": payload["pair_count"],
        "different_pairs": decision_counts.get("different", 0),
        "indistinguishable_pairs": decision_counts.get("indistinguishable", 0),
        "inconclusive_pairs": decision_counts.get("inconclusive", 0),
    }


def summarize_novel_alleles(path: str | Path) -> dict[str, float | int | str]:
    rows = []
    with Path(path).open("r", encoding="utf-8") as handle:
        for index, raw_line in enumerate(handle):
            line = raw_line.strip()
            if not line:
                continue
            if index == 0 and line.startswith("locus\tallele_id"):
                continue
            rows.append(line.split("\t"))
    return {
        "kind": "novel_alleles",
        "novel_allele_count": len(rows),
        "unique_loci": len({row[0] for row in rows}),
        "orf_reason_count": len({row[6] for row in rows if len(row) > 6 and row[6]}),
    }
