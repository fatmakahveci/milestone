from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path


MODULE_PATH = Path(__file__).resolve().parent.parent / "workflow" / "scripts" / "benchmark_wgmlst_profiles.py"
SPEC = importlib.util.spec_from_file_location("benchmark_wgmlst_profiles", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def test_benchmark_directories_reports_exact_and_different_samples(tmp_path: Path) -> None:
    predicted = tmp_path / "predicted"
    truth = tmp_path / "truth"
    predicted.mkdir()
    truth.mkdir()
    (predicted / "sample_a_wgmlst.tsv").write_text("abc\t1\nadk\t2\n", encoding="utf-8")
    (truth / "sample_a_wgmlst.tsv").write_text("abc\t1\nadk\t2\n", encoding="utf-8")
    (predicted / "sample_b_wgmlst.tsv").write_text("abc\t1\nadk\t7\n", encoding="utf-8")
    (truth / "sample_b_wgmlst.tsv").write_text("abc\t1\nadk\t2\n", encoding="utf-8")

    report = MODULE.benchmark_directories(predicted, truth)

    assert report["matched_samples"] == 2
    assert report["exact_allele_id_match_count"] == 1
    assert report["mean_concordance"] == 0.75
    assert "report_date" in report
    assert report["mean_non_comparable_rate"] == 0.0
    assert report["locus_mismatch_summary"][0]["locus"] in {"abc", "adk"}


def test_write_benchmark_report_creates_json_and_tsv(tmp_path: Path) -> None:
    report = {
        "kind": "profile_benchmark",
        "predicted_dir": "pred",
        "truth_dir": "truth",
        "matched_samples": 1,
        "missing_predicted_samples": [],
        "missing_truth_samples": [],
        "exact_allele_id_match_count": 1,
        "mean_concordance": 1.0,
        "mean_non_comparable_rate": 0.0,
        "locus_mismatch_summary": [
            {
                "locus": "abc",
                "sample_count": 1,
                "identical_count": 1,
                "different_count": 0,
                "unresolved_count": 0,
                "discordance_rate": 0.0,
                "non_comparable_rate": 0.0,
            }
        ],
        "samples": [
            {
                "sample": "sample_a",
                "decision": "indistinguishable",
                "comparable_loci": 2,
                "identical_loci": 2,
                "differing_loci": 0,
                "unresolved_loci": 0,
                "allele_distance": 0.0,
                "concordance": 1.0,
                "non_comparable_rate": 0.0,
                "exact_allele_id_match": True,
            }
        ],
    }

    json_path, tsv_path = MODULE.write_benchmark_report(report, tmp_path)

    loaded = json.loads(json_path.read_text(encoding="utf-8"))
    assert loaded["exact_allele_id_match_count"] == 1
    assert "sample\tsample_a" not in tsv_path.read_text(encoding="utf-8")
    assert "sample_a\tindistinguishable" in tsv_path.read_text(encoding="utf-8")
    assert (tmp_path / "wgmlst_benchmark_per_locus.tsv").exists()
    assert (tmp_path / "wgmlst_benchmark_report.html").exists()


def test_benchmark_pack_collection_summarizes_species_packs(tmp_path: Path) -> None:
    pack_root = tmp_path / "packs"
    ecoli = pack_root / "escherichia_coli"
    ecoli_pred = ecoli / "predicted"
    ecoli_truth = ecoli / "truth"
    ecoli_pred.mkdir(parents=True)
    ecoli_truth.mkdir(parents=True)
    (ecoli / "benchmark_manifest.json").write_text(
        json.dumps(
            {
                "species": "Escherichia coli",
                "schema_name": "demo_ecoli_wgmlst",
                "source_name": "fixture",
                "source_url": "https://pubmlst.org/",
                "retrieval_date": "2026-03-04",
                "source_snapshot": "fixture:1",
                "selection_criteria": "regression",
            }
        ),
        encoding="utf-8",
    )
    (ecoli_pred / "ecoli_a_wgmlst.tsv").write_text("abc\t1\nadk\t2\n", encoding="utf-8")
    (ecoli_truth / "ecoli_a_wgmlst.tsv").write_text("abc\t1\nadk\t2\n", encoding="utf-8")

    report = MODULE.benchmark_pack_collection(pack_root)

    assert report["species_count"] == 1
    assert report["packs"][0]["species"] == "Escherichia coli"
    assert report["mean_non_comparable_rate"] == 0.0


def test_write_benchmark_collection_report_creates_summary_files(tmp_path: Path) -> None:
    report = {
        "kind": "profile_benchmark_collection",
        "collection_dir": "packs",
        "species_count": 1,
        "mean_concordance": 1.0,
        "mean_non_comparable_rate": 0.0,
        "packs": [
            {
                "species": "Escherichia coli",
                "schema_name": "demo_ecoli_wgmlst",
                "matched_samples": 1,
                "exact_allele_id_match_count": 1,
                "mean_concordance": 1.0,
                "mean_non_comparable_rate": 0.0,
                "benchmark_pack": {
                    "source_name": "fixture",
                    "source_url": "https://pubmlst.org/",
                    "retrieval_date": "2026-03-04",
                    "source_snapshot": "fixture:1",
                    "selection_criteria": "regression",
                },
            }
        ],
    }

    json_path, tsv_path = MODULE.write_benchmark_collection_report(report, tmp_path)

    assert json.loads(json_path.read_text(encoding="utf-8"))["species_count"] == 1
    assert "Escherichia coli\tdemo_ecoli_wgmlst\t1\t1\t1.0\t0.0\tfixture" in tsv_path.read_text(encoding="utf-8")
    assert (tmp_path / "wgmlst_benchmark_collection_report.html").exists()
