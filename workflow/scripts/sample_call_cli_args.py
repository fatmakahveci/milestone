from __future__ import annotations

import argparse


def parse_cli_args(default_min_locus_coverage: float) -> argparse.Namespace:
    parser = argparse.ArgumentParser(add_help=True)

    parser.add_argument(
        "--reference_fasta",
        type=str,
        required=True,
        help="Reference's fasta file name with its directory.",
    )
    parser.add_argument(
        "--reference_info",
        type=str,
        required=True,
        help="Reference's info file name with its directory.",
    )
    parser.add_argument(
        "--reference_vcf",
        type=str,
        required=True,
        help="Reference's vcf file name with its directory.",
    )
    parser.add_argument(
        "--sample_depth",
        type=str,
        required=True,
        help="Sample's depth file name with its directory.",
    )
    parser.add_argument(
        "--sample_vcf",
        type=str,
        required=True,
        help="Sample's vcf file name with its directory.",
    )
    parser.add_argument(
        "--schema_dir",
        type=str,
        required=True,
        help="Directory of schema's to write novel alleles.",
    )
    parser.add_argument(
        "--sample_sam",
        type=str,
        required=True,
        help="Sample's sam file name with its directory.",
    )
    parser.add_argument(
        "--threads",
        type=str,
        required=True,
        help="Number of threads.",
    )
    parser.add_argument(
        "--update_reference",
        type=str,
        required=False,
        help="Update reference's vcf and info file for the further analysis. False if it is not given.",
    )
    parser.add_argument(
        "--min-locus-coverage",
        dest="min_locus_coverage",
        type=float,
        required=False,
        default=default_min_locus_coverage,
        help="Minimum locus breadth coverage percentage before assigning a non-LNF allele.",
    )
    parser.add_argument(
        "--translation-table",
        dest="translation_table",
        type=int,
        required=False,
        default=11,
        help="NCBI translation table used for ORF validation. Default: 11.",
    )
    parser.add_argument(
        "--allowed-start-codons",
        dest="allowed_start_codons",
        type=str,
        required=False,
        default=None,
        help="Optional comma-separated start codons overriding the translation-table defaults.",
    )

    return parser.parse_args()
