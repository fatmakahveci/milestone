from __future__ import annotations

import re


def get_var_type(info: str) -> str:
    start = info.index("TYPE=") + 5

    if not info.strip("/n").endswith(";") and ";" not in info[start:-1]:
        return info[start:]

    end = info.index(";", start)
    return info[start:end]


def get_cigar(info: str) -> str:
    start = info.index("CIGAR=") + 6
    end = info.index(";", start)
    return info[start:end]


def get_cigar_info(info: str) -> tuple[str, int]:
    cigar = get_cigar(info=info).split(",")[0]

    for sep in ["M", "D", "I", "S", "H", "=", "X"]:
        cigar = cigar.replace(sep, ".")

    return get_cigar(info=info), sum(list(map(int, cigar.rstrip(".").split("."))))


def resolve_cigar(
    vcf_line,
    cigar: str,
    base_ratio_lookup,
    reference_target_lookup,
) -> tuple[list, list, list, list]:
    pos_list, ref_list, alt_list, qual_list = [], [], [], []

    cigar_list = list(zip(list(map(int, re.findall("[0-9]+", cigar))), [s for s in cigar if not s.isdigit()]))

    i, j = 0, 0

    for k, (case_count, case) in enumerate(cigar_list):
        if case == "D":
            pos_list.append(vcf_line.pos + i - 1)
            ref_list.append(vcf_line.ref[i - 1 : i + case_count])
            alt_list.append(vcf_line.alt[j - 1])
            qual_list.append(vcf_line.qual)
            i += case_count
        elif case == "I":
            pos_list.append(vcf_line.pos + i - 1)
            ref_list.append(vcf_line.ref[i - 1])
            alt_list.append(vcf_line.alt[j - 1 : j + case_count])
            qual_list.append(vcf_line.qual)
            j += case_count
        elif case == "X":
            if k == len(cigar_list) - 1 or cigar_list[k + 1][1] == "M":
                for n in range(case_count):
                    if vcf_line.ref[i + n] != base_ratio_lookup(
                        variation_pos=vcf_line.pos + i + n,
                        cds_name=reference_target_lookup(vcf_line.chr),
                    ):
                        pos_list.append(vcf_line.pos + i + n)
                        ref_list.append(vcf_line.ref[i + n])
                        alt_list.append(vcf_line.alt[j + n])
                        qual_list.append(vcf_line.qual)
            elif cigar_list[k + 1][1] in {"I", "D"}:
                for n in range(case_count):
                    pos_list.append(vcf_line.pos + i - 1)
                    ref_list.append(vcf_line.ref[i + n - 1])
                    alt_list.append(vcf_line.alt[j + n - 1])
                    qual_list.append(vcf_line.qual)
            i += case_count
            j += case_count
        else:
            i += case_count
            j += case_count

    return pos_list, ref_list, alt_list, qual_list


def get_sample_sam_dict(sample_sam_path: str, sam_class) -> dict:
    sample_sam_dict = {}

    with open(sample_sam_path, "r") as file:
        for line in file.readlines():
            if line.startswith("@"):
                continue
            sam = sam_class(line)
            sample_sam_dict.setdefault(sam.target_sequence_name, []).append(
                [sam.position, sam.sequence, sam.cigar]
            )

    return sample_sam_dict


def create_sample_variation_dict(
    sample_vcf_path: str,
    vcf_class,
    variation_factory,
    merge_variations,
    base_ratio_lookup,
    reference_target_lookup,
) -> dict:
    sample_variation_dict = {}

    with open(sample_vcf_path, "r") as file:
        for line in file.readlines():
            if line.startswith("#"):
                continue

            vcf_line = vcf_class(line)
            var_type = get_var_type(info=vcf_line.info)

            variation_pos_list, variation_ref_list, variation_alt_list, variation_qual_list = [], [], [], []

            if var_type in ["complex", "mnp", "del", "ins"]:
                cigar, _ = get_cigar_info(info=vcf_line.info)
                (
                    variation_pos_list,
                    variation_ref_list,
                    variation_alt_list,
                    variation_qual_list,
                ) = resolve_cigar(
                    vcf_line=vcf_line,
                    cigar=cigar,
                    base_ratio_lookup=base_ratio_lookup,
                    reference_target_lookup=reference_target_lookup,
                )
            elif var_type == "snp":
                if len(vcf_line.alt) > 1 and len(vcf_line.alt) == len(vcf_line.ref):
                    for sb in range(len(vcf_line.alt)):
                        variation_pos_list.append(vcf_line.pos + sb)
                        variation_ref_list.append(vcf_line.ref[sb])
                        variation_alt_list.append(vcf_line.alt[sb])
                        variation_qual_list.append(vcf_line.qual)
                else:
                    variation_pos_list.append(vcf_line.pos)
                    variation_ref_list.append(vcf_line.ref)
                    alt_idx = 0 if len(vcf_line.alt) == 1 else int(vcf_line.sample[0])
                    variation_alt_list.append(vcf_line.alt.split(",")[alt_idx])
                    variation_qual_list.append(vcf_line.qual)

            if variation_pos_list:
                is_new_locus = vcf_line.chr not in sample_variation_dict
                if is_new_locus:
                    sample_variation_dict[vcf_line.chr] = variation_factory(
                        f"{variation_pos_list[0]}*{variation_ref_list[0]}>{variation_alt_list[0]}-{variation_qual_list[0]}"
                    )

                start_idx = 1 if is_new_locus else 0
                for i in range(start_idx, len(variation_pos_list)):
                    sample_variation_dict[vcf_line.chr].pos_list.append(variation_pos_list[i])
                    sample_variation_dict[vcf_line.chr].ref_list.append(variation_ref_list[i])
                    sample_variation_dict[vcf_line.chr].alt_list.append(variation_alt_list[i])
                    sample_variation_dict[vcf_line.chr].qual_list.append(variation_qual_list[i])

    for sample_cds, variations in sample_variation_dict.items():
        sample_variation_dict[sample_cds] = merge_variations(variations=variations)

    return sample_variation_dict
