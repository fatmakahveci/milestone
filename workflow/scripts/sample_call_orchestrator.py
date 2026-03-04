from __future__ import annotations

from typing import Callable


class SampleCallingResult:
    def __init__(self, output_path: str, allele_calls: dict[str, str]):
        self.output_path = output_path
        self.allele_calls = allele_calls


def build_wgmlst_output_path(sample_vcf_path: str) -> str:
    return f"{sample_vcf_path[:-4]}_wgmlst.tsv"


def run_sample_calling_workflow(
    sample_vcf_path: str,
    init_state: Callable[[], None],
    collect_allele_calls: Callable[[], dict[str, str]],
    get_cds_name_from_allele_name: Callable[[str], str],
) -> SampleCallingResult:
    init_state()
    allele_calls = collect_allele_calls()
    output_path = build_wgmlst_output_path(sample_vcf_path)

    with open(output_path, "w", encoding="utf-8") as file:
        for sample_cds, allele_id in allele_calls.items():
            file.write(f"{get_cds_name_from_allele_name(sample_cds)}\t{allele_id}\n")

    return SampleCallingResult(output_path=output_path, allele_calls=allele_calls)
