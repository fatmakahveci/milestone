from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class DemoRun:
    title: str
    summary: str
    key_metric: str
    value: str


SCHEMA_RUNS = [
    DemoRun(
        title="Reference bundle",
        summary="Reference FASTA, VCF, and info files are assembled from locus FASTA inputs.",
        key_metric="Loci indexed",
        value="1,742",
    ),
    DemoRun(
        title="Novel allele ready",
        summary="Reference artifacts are prepared for later allele discovery and updates.",
        key_metric="Info entries",
        value="18,406",
    ),
]


ALLELE_RUNS = [
    DemoRun(
        title="Clinical sample ACME-42",
        summary="Reads align against the graph reference and produce a wgMLST summary table.",
        key_metric="Assigned loci",
        value="1,703 / 1,742",
    ),
    DemoRun(
        title="Reference update candidate",
        summary="Novel allele candidates are flagged for optional reference enrichment.",
        key_metric="Novel alleles",
        value="14",
    ),
]

def schema_outputs(reference: str, output: str) -> list[str]:
    base = output.rstrip("/")
    return [
        f"{base}/{reference}.fasta",
        f"{base}/{reference}.vcf.gz",
        f"{base}/{reference}_info.txt",
        f"{base}/logs/schema_creation.log",
    ]


def allele_outputs(reference: str, output: str, aligner: str, sample_name: str) -> list[str]:
    base = output.rstrip("/")
    return [
        f"{base}/{aligner}/{sample_name}_wgmlst.tsv",
        f"{base}/{aligner}/{sample_name}_novel_alleles.fasta",
        f"{base}/{aligner}/{sample_name}_novel_alleles.tsv",
        f"{base}/{aligner}/{reference}.xg",
        f"{base}/logs/allele_calling.log",
    ]


def compare_outputs(output: str) -> list[str]:
    base = output.rstrip("/")
    return [
        f"{base}/wgmlst_profile_summary.json",
        f"{base}/wgmlst_profile_comparison.tsv",
    ]


def schema_qc_outputs(output: str) -> list[str]:
    base = output.rstrip("/")
    return [
        f"{base}/schema_qc_summary.json",
        f"{base}/schema_qc_issues.tsv",
    ]


def benchmark_outputs(output: str) -> list[str]:
    base = output.rstrip("/")
    return [
        f"{base}/wgmlst_benchmark_summary.json",
        f"{base}/wgmlst_benchmark_per_sample.tsv",
        f"{base}/wgmlst_benchmark_per_locus.tsv",
        f"{base}/wgmlst_benchmark_report.html",
        f"{base}/wgmlst_benchmark_pack_summary.json",
        f"{base}/wgmlst_benchmark_species.tsv",
        f"{base}/wgmlst_benchmark_collection_report.html",
    ]


def matrix_outputs(output: str) -> list[str]:
    base = output.rstrip("/")
    return [
        f"{base}/wgmlst_distance_matrix.tsv",
        f"{base}/wgmlst_distance_summary.json",
        f"{base}/wgmlst_distance_report.html",
    ]


def batch_compare_outputs(output: str) -> list[str]:
    base = output.rstrip("/")
    return [
        f"{base}/wgmlst_batch_compare_summary.json",
        f"{base}/wgmlst_batch_compare.tsv",
    ]


def species_benchmark_pack_dir() -> str:
    return "tests/fixtures/benchmark/public_species_packs"


def pubmlst_benchmark_pack_outputs(output: str) -> list[str]:
    base = output.rstrip("/")
    return [
        f"{base}/benchmark_manifest.json",
        f"{base}/isolate_metadata.json",
        f"{base}/predicted/README.md",
    ]


def enterobase_scheme_outputs(output: str, scheme_name: str) -> list[str]:
    base = output.rstrip("/")
    return [
        f"{base}/{scheme_name}_profiles.tar.gz",
        f"{base}/enterobase_scheme_metadata.json",
    ]


def sample_wgmlst_table() -> list[dict[str, str]]:
    return [
        {"locus": "abcZ", "allele_id": "12", "status": "match"},
        {"locus": "adk", "allele_id": "4", "status": "match"},
        {"locus": "aroE", "allele_id": "LNF", "status": "low coverage"},
        {"locus": "fumC", "allele_id": "NEW_14", "status": "novel allele"},
        {"locus": "gdh", "allele_id": "1", "status": "reference allele"},
    ]


def sample_comparison_table() -> list[dict[str, str]]:
    return [
        {"locus": "abcZ", "sample_a": "12", "sample_b": "12", "status": "identical"},
        {"locus": "adk", "sample_a": "4", "sample_b": "7", "status": "different"},
        {"locus": "aroE", "sample_a": "LNF", "sample_b": "3", "status": "unresolved"},
        {"locus": "fumC", "sample_a": "NEW_14", "sample_b": "NEW_14", "status": "identical"},
    ]


def sample_novel_allele_table() -> list[dict[str, str]]:
    return [
        {
            "locus": "fumC",
            "allele_id": "14",
            "sequence_length": "1503",
            "status": "novel_candidate",
            "reference_allele": "fumC_1",
            "variation_summary": "187:C>T,411:G>A",
            "orf_reason": "passes_orf_qc",
        },
        {
            "locus": "icd",
            "allele_id": "22",
            "sequence_length": "1134",
            "status": "novel_candidate",
            "reference_allele": "icd_1",
            "variation_summary": "92:A>G",
            "orf_reason": "passes_orf_qc",
        },
    ]


def sample_benchmark_table() -> list[dict[str, str]]:
    return [
        {
            "sample": "strain_a",
            "decision": "indistinguishable",
            "concordance": "1.00",
            "non_comparable_rate": "0.00",
            "exact_allele_id_match": "True",
        },
        {
            "sample": "strain_b",
            "decision": "different",
            "concordance": "0.88",
            "non_comparable_rate": "0.06",
            "exact_allele_id_match": "False",
        },
    ]


def available_species_packs() -> list[dict[str, str]]:
    return [
        {
            "species": "Escherichia coli",
            "path": "tests/fixtures/benchmark/public_species_packs/escherichia_coli",
        },
        {
            "species": "Neisseria meningitidis",
            "path": "tests/fixtures/benchmark/public_species_packs/neisseria_meningitidis",
        },
        {
            "species": "Acinetobacter baumannii",
            "path": "tests/fixtures/benchmark/public_species_packs/acinetobacter_baumannii",
        },
    ]
