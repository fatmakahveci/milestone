from __future__ import annotations

import subprocess
from pathlib import Path


def call_variants_of_allele(
    reference: str,
    sample: str,
    threads: int,
    run_command,
) -> None:
    with open(f"{sample}.paf", "w") as paf_handle:
        run_command(
            [
                "minimap2",
                "-c",
                "--cs=long",
                "-t",
                str(threads),
                "-A",
                "1",
                "-B",
                "4",
                "-O",
                "6",
                "-E",
                "1",
                f"{reference}.fasta",
                f"{sample}.fasta",
            ],
            stdout=paf_handle,
            stderr=subprocess.DEVNULL,
        )

    with open(f"{sample}.paf", "r", encoding="utf-8") as paf_input, open(f"{sample}.vcf", "w", encoding="utf-8") as vcf_output:
        sort_process = subprocess.Popen(
            ["sort", "-k6,6", "-k8,8n"],
            stdin=paf_input,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
        )
        assert sort_process.stdout is not None
        call_process = subprocess.Popen(
            ["paftools.js", "call", "-L0", "-l0", "-f", f"{reference}.fasta", "-s", sample, "-"],
            stdin=sort_process.stdout,
            stdout=vcf_output,
            stderr=subprocess.DEVNULL,
            text=True,
        )
        sort_process.stdout.close()
        call_returncode = call_process.wait()
        sort_returncode = sort_process.wait()
        if sort_returncode != 0 or call_returncode != 0:
            raise subprocess.CalledProcessError(call_returncode or sort_returncode, ["sort", "|", "paftools.js", "call"])

    run_command(
        ["bcftools", "norm", "-a", "-d", "none", f"{sample}.vcf", "-Ov", "-o", f"{sample}.vcf."],
        stderr=subprocess.DEVNULL,
    )
    Path(f"{sample}.vcf.").replace(f"{sample}.vcf")


def sort_zip_and_index_vcf_files(vcf_file: str, run_command) -> None:
    run_command(["bcftools", "sort", vcf_file, "-Ov", "-o", f"{vcf_file}."], stderr=subprocess.DEVNULL)
    Path(f"{vcf_file}.").replace(vcf_file)
    run_command(["bgzip", "-f", vcf_file], stderr=subprocess.DEVNULL)
    run_command(["tabix", "-f", "-p", "vcf", f"{vcf_file}.gz"], stderr=subprocess.DEVNULL)
    run_command(
        ["bcftools", "concat", "-a", "--rm-dups", "none", f"{vcf_file}.gz", "-Ov", "-o", vcf_file],
        stderr=subprocess.DEVNULL,
    )
    run_command(["bgzip", "-f", vcf_file], stderr=subprocess.DEVNULL)
    run_command(["tabix", "-f", "-p", "vcf", f"{vcf_file}.gz"], stderr=subprocess.DEVNULL)
