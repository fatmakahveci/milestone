from __future__ import annotations

from collections import Counter
from pathlib import Path

from script_utils import ParsedVariationInfo


def remove_common_suffixes(var1: str, var2: str) -> tuple[str, str]:
    check_len = min(len(var1), len(var2))
    var1 = var1[::-1]
    var2 = var2[::-1]

    i = 0
    while i < check_len and var1[i] == var2[i]:
        i += 1

    return var1[i:][::-1], var2[i:][::-1]


def remove_common_prefixes(pos: int, var1: str, var2: str) -> tuple[int, str, str]:
    check_len = min(len(var1), len(var2))

    i = 0
    while i < check_len and var1[i] == var2[i]:
        i += 1

    return pos + i, var1[i:], var2[i:]


def remove_common_mid(
    pos: int,
    var1: str,
    var2: str,
    qual: int,
) -> tuple[list[int], list[str], list[str], list[int], int]:
    split_idx_list, pos_list, ref_list, alt_list, qual_list = [], [], [], [], []

    i = 0
    while i < len(var1):
        if var1[i] != var2[i]:
            split_idx_list.append(i)
        i += 1

    for i in split_idx_list:
        pos_list.append(pos + i)
        ref_list.append(var1[i])
        alt_list.append(var2[i])
        qual_list.append(qual)

    return pos_list, ref_list, alt_list, qual_list, len(split_idx_list)


def remove_redundance(variations: ParsedVariationInfo) -> ParsedVariationInfo:
    variations_to_replace = []

    for i in range(len(variations.pos_list)):
        if len(variations.ref_list[i]) > 1 and len(variations.alt_list[i]) > 1:
            variations.ref_list[i], variations.alt_list[i] = remove_common_suffixes(
                var1=variations.ref_list[i],
                var2=variations.alt_list[i],
            )

        if len(variations.ref_list[i]) > 1 and len(variations.alt_list[i]) > 1:
            (
                variations.pos_list[i],
                variations.ref_list[i],
                variations.alt_list[i],
            ) = remove_common_prefixes(
                pos=variations.pos_list[i],
                var1=variations.ref_list[i],
                var2=variations.alt_list[i],
            )

        if len(variations.ref_list[i]) == len(variations.alt_list[i]) and len(variations.ref_list[i]) > 2:
            variations_to_replace.append(
                [
                    i,
                    remove_common_mid(
                        pos=variations.pos_list[i],
                        var1=variations.ref_list[i],
                        var2=variations.alt_list[i],
                        qual=variations.qual_list[i],
                    ),
                ]
            )

    if len(variations_to_replace) > 0:
        for variation in variations_to_replace[::-1]:
            del variations.pos_list[variation[0]]
            del variations.ref_list[variation[0]]
            del variations.alt_list[variation[0]]
            del variations.qual_list[variation[0]]

            idx_to_put = variation[0]
            var = variation[1]

            for i in range(len(var[0]))[::-1]:
                variations.pos_list.insert(idx_to_put, var[0][i])
                variations.ref_list.insert(idx_to_put, var[1][i])
                variations.alt_list.insert(idx_to_put, var[2][i])
                variations.qual_list.insert(idx_to_put, var[3][i])

    return variations


def merge_variations(variations: ParsedVariationInfo) -> ParsedVariationInfo:
    if len(variations.pos_list) == 1:
        return variations

    merged_variations_list = []

    for pos, ref, alt, qual in zip(
        variations.pos_list,
        variations.ref_list,
        variations.alt_list,
        variations.qual_list,
        strict=False,
    ):
        merged_variations_list.append(f"{pos}*{ref}>{alt}-{qual}")

    merged_list = []
    var_list = remove_redundance(ParsedVariationInfo(",".join(merged_variations_list)))

    for pos, ref, alt, qual in zip(
        var_list.pos_list,
        var_list.ref_list,
        var_list.alt_list,
        var_list.qual_list,
        strict=False,
    ):
        merged_list.append(f"{pos}*{ref}>{alt}-{qual}")

    return ParsedVariationInfo(",".join(merged_list))


def build_variation_signature(
    variation_info: ParsedVariationInfo,
) -> Counter[tuple[int, str, str]]:
    return Counter(
        zip(
            variation_info.pos_list,
            variation_info.ref_list,
            variation_info.alt_list,
            strict=False,
        )
    )


def write_novel_allele_exports(
    sample_vcf_path: str | Path,
    novel_alleles: list[dict[str, str]],
) -> tuple[Path, Path]:
    sample_vcf_path = Path(sample_vcf_path)
    fasta_path = sample_vcf_path.with_name(f"{sample_vcf_path.stem}_novel_alleles.fasta")
    tsv_path = sample_vcf_path.with_name(f"{sample_vcf_path.stem}_novel_alleles.tsv")

    with fasta_path.open("w", encoding="utf-8") as fasta_handle:
        for record in novel_alleles:
            fasta_handle.write(f">{record['locus']}_{record['allele_id']}\n{record['sequence']}\n")

    with tsv_path.open("w", encoding="utf-8") as tsv_handle:
        tsv_handle.write("locus\tallele_id\tsequence_length\tstatus\treference_allele\tvariation_summary\torf_reason\n")
        for record in novel_alleles:
            tsv_handle.write(
                f"{record['locus']}\t{record['allele_id']}\t{len(record['sequence'])}\t{record['status']}\t"
                f"{record.get('reference_allele', '')}\t{record.get('variation_summary', '')}\t"
                f"{record.get('orf_reason', '')}\n"
            )

    return fasta_path, tsv_path
