#!/usr/bin/env python3

from __future__ import annotations

import argparse
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent.parent
WORKFLOW_DIR = Path(__file__).resolve().parent
RULES_DIR = WORKFLOW_DIR / "rules"
COMPARE_SCRIPT = WORKFLOW_DIR / "scripts" / "compare_wgmlst_profiles.py"
MATRIX_SCRIPT = WORKFLOW_DIR / "scripts" / "compare_wgmlst_matrix.py"
COMPARE_BATCH_SCRIPT = WORKFLOW_DIR / "scripts" / "compare_wgmlst_batch.py"
BENCHMARK_SCRIPT = WORKFLOW_DIR / "scripts" / "benchmark_wgmlst_profiles.py"
PUBLICATION_READINESS_SCRIPT = WORKFLOW_DIR / "scripts" / "publication_readiness.py"
PUBLICATION_PACKAGE_SCRIPT = WORKFLOW_DIR / "scripts" / "build_publication_package.py"
SUPPLEMENT_SCRIPT = WORKFLOW_DIR / "scripts" / "build_manuscript_supplement.py"
VALIDATION_CORPUS_SCRIPT = WORKFLOW_DIR / "scripts" / "build_species_validation_corpus.py"
CONFIG_PATH = WORKFLOW_DIR / "config.yaml"


@dataclass(frozen=True)
class WorkflowPaths:
    log_dir: Path
    schema_log: Path
    allele_log: Path
    snakefile: Path


def workflow_paths(args: argparse.Namespace) -> WorkflowPaths:
    log_dir = Path(args.output) / "logs"
    return WorkflowPaths(
        log_dir=log_dir,
        schema_log=log_dir / "schema_creation.log",
        allele_log=log_dir / "allele_calling.log",
        snakefile=WORKFLOW_DIR / args.snakefile,
    )


def sample_name_from_read(read1: str) -> str:
    return Path(read1).name.split("_1")[0]


def build_config_text(args: argparse.Namespace, paths: WorkflowPaths) -> str:
    lines = [
        f"# milestone.py created {CONFIG_PATH}.",
        "",
        f'logs: "{paths.log_dir}"',
        f'working_dir: "{WORKFLOW_DIR}"',
        "",
        "parameters: ",
        f" threads: {args.threads}",
        f" min_locus_coverage: {args.min_locus_coverage}",
        f" translation_table: {getattr(args, 'translation_table', 11)}",
        f' allowed_start_codons: "{getattr(args, "allowed_start_codons", "") or ""}"',
        f'reference: "{Path(args.output) / args.reference}"',
        f'reference_vcf: "{Path(args.output) / f"{args.reference}.vcf"}"',
        f'reference_vcf_gz: "{Path(args.output) / f"{args.reference}.vcf.gz"}"',
        f'reference_vcf_gz_tbi: "{Path(args.output) / f"{args.reference}.vcf.gz.tbi"}"',
        f'reference_info_txt: "{Path(args.output) / f"{args.reference}_info.txt"}"',
        f'reference_fasta: "{Path(args.output) / f"{args.reference}.fasta"}"',
    ]
    if args.command == "schema_creation":
        lines.extend(
            [
                f'schema_dir: "{str(args.schema_name).rstrip("/") + "/"}"',
                f'schema_creation_log_file: "{paths.schema_log}"',
                f'output_dir: "{str(args.output).rstrip("/")}"',
            ]
        )
    elif args.command == "allele_calling":
        lines.extend(
            [
                f'schema_dir: "{str(args.schema_name).rstrip("/") + "/"}"',
                f'sample: "{sample_name_from_read(args.read1)}"',
                "samples:",
                f' sample1: "{args.read1}"',
                f' sample2: "{args.read2}"',
                f'aligner: "{args.aligner}"',
                f'output_dir: "{str(args.output).rstrip("/")}"',
                f'aligner_reference: "{Path(args.output) / args.aligner / args.reference}"',
                f'update_reference: "{args.update_reference}"',
                f'allele_calling_log_file: "{paths.allele_log}"',
            ]
        )
    return "\n".join(lines) + "\n"


def write_config(args: argparse.Namespace, paths: WorkflowPaths) -> None:
    CONFIG_PATH.write_text(build_config_text(args, paths), encoding="utf-8")


def build_snakefile_text(args: argparse.Namespace, paths: WorkflowPaths) -> str:
    lines = [
        "import glob, os, sys",
        "",
        f'configfile: "{CONFIG_PATH}"',
        "",
    ]
    if args.command == "schema_creation":
        lines.extend(
            [
                f'include: "{RULES_DIR / "schema_creation.smk"}"',
                "",
                "rule all:",
                "\tinput:",
                f'\t\treference_vcf_gz = "{Path(args.output) / f"{args.reference}.vcf.gz"}"',
            ]
        )
    elif args.command == "allele_calling":
        sample_name = sample_name_from_read(args.read1)
        out_dir = str(args.output).rstrip("/")
        lines.extend(
            [
                f'include: "{RULES_DIR / "allele_calling.smk"}"',
                "",
                "rule all:",
                "\tinput:",
                f'\t\tsample_wgmlst = "{out_dir}/{args.aligner}/{sample_name}_wgmlst.tsv"',
            ]
        )
    return "\n".join(lines) + "\n"


def write_snakefile(args: argparse.Namespace, paths: WorkflowPaths) -> None:
    paths.snakefile.write_text(build_snakefile_text(args, paths), encoding="utf-8")
    if args.command == "schema_creation":
        paths.schema_log.touch(exist_ok=True)
    elif args.command == "allele_calling":
        paths.allele_log.touch(exist_ok=True)


def build_snakemake_command(args: argparse.Namespace, paths: WorkflowPaths) -> list[str]:
    command = [
        "snakemake",
        "-p",
        "--use-conda",
        "--configfile",
        str(CONFIG_PATH),
        "--cores",
        str(args.threads),
        "--snakefile",
        str(paths.snakefile),
    ]
    if args.forceall:
        command.append("--forceall")
    if args.dryrun:
        command.append("--dryrun")
    if args.printshellcmds:
        command.append("--printshellcmds")
    if args.ri:
        command.append("--rerun-incomplete")
    if args.unlock:
        command.append("--unlock")
    if args.quiet:
        command.append("--quiet")
    return command


def run_snakemake(args: argparse.Namespace, paths: WorkflowPaths) -> int:
    completed = subprocess.run(build_snakemake_command(args, paths), check=False, cwd=REPO_ROOT)
    return completed.returncode


def run_profile_compare(args: argparse.Namespace) -> int:
    command = [
        sys.executable,
        str(COMPARE_SCRIPT),
        "--profile-a",
        args.profile_a,
        "--profile-b",
        args.profile_b,
        "--output-dir",
        args.output,
        "--label-a",
        args.label_a,
        "--label-b",
        args.label_b,
    ]
    completed = subprocess.run(command, check=False, cwd=REPO_ROOT)
    return completed.returncode


def run_profile_matrix(args: argparse.Namespace) -> int:
    command = [
        sys.executable,
        str(MATRIX_SCRIPT),
        "--profiles",
        *args.profiles,
        "--output-dir",
        args.output,
    ]
    command.extend(["--distance-threshold", str(args.distance_threshold)])
    completed = subprocess.run(command, check=False, cwd=REPO_ROOT)
    return completed.returncode


def run_profile_compare_batch(args: argparse.Namespace) -> int:
    command = [
        sys.executable,
        str(COMPARE_BATCH_SCRIPT),
        "--profiles",
        *args.profiles,
        "--output-dir",
        args.output,
    ]
    completed = subprocess.run(command, check=False, cwd=REPO_ROOT)
    return completed.returncode


def run_profile_benchmark(args: argparse.Namespace) -> int:
    command = [
        sys.executable,
        str(BENCHMARK_SCRIPT),
        "--output-dir",
        args.output,
    ]
    benchmark_pack_dir = getattr(args, "benchmark_pack_dir", None)
    if benchmark_pack_dir:
        command.extend(["--benchmark-pack-dir", benchmark_pack_dir])
    else:
        command.extend(["--predicted-dir", args.predicted_dir, "--truth-dir", args.truth_dir])
    completed = subprocess.run(command, check=False, cwd=REPO_ROOT)
    return completed.returncode


def run_publication_readiness(args: argparse.Namespace) -> int:
    command = [
        sys.executable,
        str(PUBLICATION_READINESS_SCRIPT),
        "--schema-qc-summary",
        args.schema_qc_summary,
        "--benchmark-summary",
        args.benchmark_summary,
        "--output-dir",
        args.output,
        "--max-non-comparable-rate",
        str(args.max_non_comparable_rate),
    ]
    if args.schema_manifest:
        command.extend(["--schema-manifest", args.schema_manifest])
    completed = subprocess.run(command, check=False, cwd=REPO_ROOT)
    return completed.returncode


def run_publication_package(args: argparse.Namespace) -> int:
    command = [
        sys.executable,
        str(PUBLICATION_PACKAGE_SCRIPT),
        "--schema-qc-summary",
        args.schema_qc_summary,
        "--benchmark-summary",
        args.benchmark_summary,
        "--output",
        args.output,
        "--max-non-comparable-rate",
        str(args.max_non_comparable_rate),
    ]
    if args.schema_manifest:
        command.extend(["--schema-manifest", args.schema_manifest])
    if args.title:
        command.extend(["--title", args.title])
    if args.notes:
        command.extend(["--notes", args.notes])
    completed = subprocess.run(command, check=False, cwd=REPO_ROOT)
    return completed.returncode


def run_manuscript_supplement(args: argparse.Namespace) -> int:
    command = [
        sys.executable,
        str(SUPPLEMENT_SCRIPT),
        "--publication-package",
        args.publication_package,
        "--output-dir",
        args.output,
    ]
    completed = subprocess.run(command, check=False, cwd=REPO_ROOT)
    return completed.returncode


def run_validation_corpus(args: argparse.Namespace) -> int:
    command = [
        sys.executable,
        str(VALIDATION_CORPUS_SCRIPT),
        "--collection-dir",
        args.collection_dir,
        "--output-dir",
        args.output,
        "--zip-name",
        args.zip_name,
    ]
    completed = subprocess.run(command, check=False, cwd=REPO_ROOT)
    return completed.returncode


def build_parent_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument(
        "-n",
        "--dryrun",
        "--dry-run",
        default=False,
        action="store_true",
        required=False,
        help="Snakemake dry-run mode.",
    )
    parser.add_argument(
        "-p",
        "--printshellcmds",
        default=False,
        action="store_true",
        required=False,
        help="Print Snakemake shell commands before execution.",
    )
    parser.add_argument(
        "-s",
        "--snakefile",
        default="Snakefile",
        required=False,
        help="Snakefile name to generate inside workflow/.",
    )
    parser.add_argument(
        "-t",
        "--threads",
        "--set-threads",
        required=False,
        type=int,
        default=1,
        help="Thread count for Snakemake rules.",
    )
    parser.add_argument(
        "--min-locus-coverage",
        dest="min_locus_coverage",
        required=False,
        type=float,
        default=60.0,
        help="Minimum locus breadth coverage percentage before calling a locus as present. Default: 60.",
    )
    parser.add_argument(
        "--translation-table",
        dest="translation_table",
        required=False,
        type=int,
        default=11,
        help="NCBI translation table used for ORF validation. Default: 11.",
    )
    parser.add_argument(
        "--allowed-start-codons",
        dest="allowed_start_codons",
        required=False,
        default=None,
        help="Optional comma-separated start codons overriding translation-table defaults.",
    )
    parser.add_argument(
        "-F",
        "--forceall",
        default=False,
        action="store_true",
        required=False,
        help="Force all dependent rules to rerun.",
    )
    parser.add_argument(
        "--ri",
        "--rerun-incomplete",
        default=False,
        action="store_true",
        required=False,
        help="Rerun incomplete jobs.",
    )
    parser.add_argument(
        "--unlock",
        default=False,
        action="store_true",
        required=False,
        help="Unlock Snakemake working directory.",
    )
    parser.add_argument(
        "-q",
        "--quiet",
        default=False,
        action="store_true",
        required=False,
        help="Reduce Snakemake progress output.",
    )
    parser.add_argument(
        "-r",
        "--reference",
        required=True,
        help="Reference prefix without extension.",
    )
    return parser


def parse_arguments() -> argparse.Namespace:
    parent_parser = build_parent_parser()
    parser = argparse.ArgumentParser(add_help=True)
    subparsers = parser.add_subparsers(title="commands", dest="command", required=True)

    schema_creation_parser = subparsers.add_parser(
        "schema_creation",
        parents=[parent_parser],
        description="Run schema_creation workflow to create reference FASTA and VCF files.",
        help="schema_creation",
    )
    schema_creation_parser.add_argument(
        "-sn",
        "--schema_name",
        type=str,
        required=True,
        help="Directory containing user-provided coding sequences and their alleles.",
    )
    schema_creation_parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Directory for schema_creation outputs.",
    )

    allele_calling_parser = subparsers.add_parser(
        "allele_calling",
        parents=[parent_parser],
        description="Run graph-based allele calling and optional reference update.",
        help="allele_calling",
    )
    allele_calling_parser.add_argument(
        "--aligner",
        default="vg",
        required=False,
        help="Graph aligner to use. Default: vg.",
    )
    allele_calling_parser.add_argument(
        "-e",
        "--read1",
        type=str,
        required=True,
        help="First paired-end read file.",
    )
    allele_calling_parser.add_argument(
        "-E",
        "--read2",
        type=str,
        required=True,
        help="Second paired-end read file.",
    )
    allele_calling_parser.add_argument(
        "--ur",
        "--update_reference",
        dest="update_reference",
        default=False,
        action="store_true",
        required=False,
        help="Update reference VCF and info files after calling novel alleles.",
    )
    allele_calling_parser.add_argument(
        "-sn",
        "--schema_name",
        type=str,
        required=True,
        help="Directory containing schema sequences and alleles.",
    )
    allele_calling_parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Directory for allele_calling outputs.",
    )

    compare_parser = subparsers.add_parser(
        "profile_compare",
        description="Compare two wgMLST profiles and report whether they are distinguishable.",
        help="profile_compare",
    )
    compare_parser.add_argument("--profile-a", required=True, help="First wgMLST TSV profile.")
    compare_parser.add_argument("--profile-b", required=True, help="Second wgMLST TSV profile.")
    compare_parser.add_argument("-o", "--output", required=True, help="Output directory for comparison results.")
    compare_parser.add_argument("--label-a", default="sample_a", help="Display label for profile A.")
    compare_parser.add_argument("--label-b", default="sample_b", help="Display label for profile B.")

    matrix_parser = subparsers.add_parser(
        "profile_matrix",
        description="Build a pairwise wgMLST distance matrix from multiple profiles.",
        help="profile_matrix",
    )
    matrix_parser.add_argument(
        "--profiles",
        nargs="+",
        required=True,
        help="Two or more wgMLST TSV profiles.",
    )
    matrix_parser.add_argument(
        "--distance-threshold",
        type=float,
        default=0.0,
        help="Distance threshold used for clustering profiles in the matrix summary.",
    )
    matrix_parser.add_argument("-o", "--output", required=True, help="Output directory for matrix results.")

    compare_batch_parser = subparsers.add_parser(
        "profile_compare_batch",
        description="Compare multiple wgMLST profiles and emit pairwise decisions.",
        help="profile_compare_batch",
    )
    compare_batch_parser.add_argument("--profiles", nargs="+", required=True, help="Two or more wgMLST TSV profiles.")
    compare_batch_parser.add_argument("-o", "--output", required=True, help="Output directory for batch comparisons.")

    benchmark_parser = subparsers.add_parser(
        "profile_benchmark",
        description="Benchmark predicted wgMLST profiles against truth-set profiles.",
        help="profile_benchmark",
    )
    benchmark_parser.add_argument("--predicted-dir", help="Directory containing predicted wgMLST TSVs.")
    benchmark_parser.add_argument("--truth-dir", help="Directory containing truth wgMLST TSVs.")
    benchmark_parser.add_argument(
        "--benchmark-pack-dir",
        help="Directory containing species benchmark packs with manifest/predicted/truth layout.",
    )
    benchmark_parser.add_argument("-o", "--output", required=True, help="Output directory for benchmark results.")

    readiness_parser = subparsers.add_parser(
        "publication_readiness",
        description="Evaluate whether a Milestone analysis has the documentation and validation expected for publication.",
        help="publication_readiness",
    )
    readiness_parser.add_argument("--schema-qc-summary", required=True, help="Path to schema_qc_summary.json.")
    readiness_parser.add_argument("--benchmark-summary", required=True, help="Path to wgmlst_benchmark_summary.json.")
    readiness_parser.add_argument("--schema-manifest", help="Optional path to schema_manifest.json.")
    readiness_parser.add_argument(
        "--max-non-comparable-rate",
        type=float,
        default=0.1,
        help="Warn when benchmark mean_non_comparable_rate exceeds this value.",
    )
    readiness_parser.add_argument("-o", "--output", required=True, help="Output directory for readiness outputs.")

    package_parser = subparsers.add_parser(
        "publication_package",
        description="Create a ZIP package containing publication-support artifacts and metadata.",
        help="publication_package",
    )
    package_parser.add_argument("--schema-qc-summary", required=True, help="Path to schema_qc_summary.json.")
    package_parser.add_argument("--benchmark-summary", required=True, help="Path to wgmlst_benchmark_summary.json.")
    package_parser.add_argument("--schema-manifest", help="Optional path to schema_manifest.json.")
    package_parser.add_argument("--title", help="Optional package title.")
    package_parser.add_argument("--notes", help="Optional package notes.")
    package_parser.add_argument(
        "--max-non-comparable-rate",
        type=float,
        default=0.1,
        help="Warn in readiness report when mean_non_comparable_rate exceeds this value.",
    )
    package_parser.add_argument("-o", "--output", required=True, help="Output ZIP path.")

    supplement_parser = subparsers.add_parser(
        "manuscript_supplement",
        description="Extract a publication package into a supplement-ready directory.",
        help="manuscript_supplement",
    )
    supplement_parser.add_argument("--publication-package", required=True, help="Publication package ZIP.")
    supplement_parser.add_argument("-o", "--output", required=True, help="Output directory for supplement files.")

    corpus_parser = subparsers.add_parser(
        "validation_corpus",
        description="Build a species-specific validation corpus package from benchmark packs.",
        help="validation_corpus",
    )
    corpus_parser.add_argument("--collection-dir", required=True, help="Directory containing benchmark pack subdirectories.")
    corpus_parser.add_argument("--zip-name", default="species_validation_corpus.zip", help="Output ZIP filename.")
    corpus_parser.add_argument("-o", "--output", required=True, help="Output directory for corpus artifacts.")

    return parser.parse_args()


def main() -> int:
    args = parse_arguments()
    if args.command == "profile_compare":
        return run_profile_compare(args)
    if args.command == "profile_matrix":
        return run_profile_matrix(args)
    if args.command == "profile_compare_batch":
        return run_profile_compare_batch(args)
    if args.command == "profile_benchmark":
        return run_profile_benchmark(args)
    if args.command == "publication_readiness":
        return run_publication_readiness(args)
    if args.command == "publication_package":
        return run_publication_package(args)
    if args.command == "manuscript_supplement":
        return run_manuscript_supplement(args)
    if args.command == "validation_corpus":
        return run_validation_corpus(args)

    paths = workflow_paths(args)
    paths.log_dir.mkdir(parents=True, exist_ok=True)
    write_config(args, paths)
    write_snakefile(args, paths)

    if not CONFIG_PATH.exists():
        print(f"Path to configfile does not exist: {CONFIG_PATH}")
        return 1
    return run_snakemake(args, paths)


if __name__ == "__main__":
    raise SystemExit(main())
