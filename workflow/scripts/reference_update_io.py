from __future__ import annotations

import os


def write_allele_sequence_to_schema_seed(
    schema_dir: str,
    sample_cds: str,
    cds_allele_id: str,
    sequence_with_variation: str,
) -> None:
    with open(os.path.join(schema_dir, sample_cds + ".fasta"), "a") as file:
        file.write(f">{sample_cds}_{cds_allele_id}\n")
        file.write(f"{sequence_with_variation}\n")


def write_variations_to_reference_info_file(
    reference_info_path: str,
    cds: str,
    allele_id: str,
    cds_variation,
) -> None:
    with open(reference_info_path, "a") as file:
        line = []
        for pos, ref, alt, qual in zip(
            cds_variation.pos_list,
            cds_variation.ref_list,
            cds_variation.alt_list,
            cds_variation.qual_list,
            strict=False,
        ):
            if type(alt) is list:
                alt = ";".join(alt)
            line.append(f"{pos}*{ref}>{alt}-{qual}")

        file.write(f"{cds}_{allele_id}\t{','.join(line)}\n")


def write_variations_to_reference_vcf_file(
    sample_vcf_path: str,
    cds: str,
    reference_allele_name: str,
    temp_sample_vcf_dir: str,
    cds_variation,
    vcf_class,
    run_command,
) -> None:
    temp_sample_vcf_file_name = os.path.join(temp_sample_vcf_dir, f"{cds}.vcf")

    with open(temp_sample_vcf_file_name, "w") as temp_sample_vcf_file:
        with open(sample_vcf_path, "r") as file:
            for line in file.readlines():
                if line.startswith("#"):
                    if line.startswith("#CHROM"):
                        temp_sample_vcf_file.write("\t".join(line.rstrip("\n").split("\t")[:-1]) + "\tREFERENCE\n")
                    else:
                        temp_sample_vcf_file.write(line)
                else:
                    if line.startswith(f"{reference_allele_name}\t"):
                        vcf_line = vcf_class(line)
                        if vcf_line.pos in cds_variation.pos_list:
                            vcf_line.chr = reference_allele_name
                            vcf_line.info = "."
                            vcf_line.sample_format = "GT"
                            vcf_line.sample = "1"
                            temp_sample_vcf_file.write(str(vcf_line))

    run_command(["bgzip", "-f", temp_sample_vcf_file_name])
    run_command(["tabix", "-p", "vcf", f"{temp_sample_vcf_file_name}.gz"])
