from __future__ import annotations


def merge_reference_with_novel_vcfs(
    temp_vcfs: list[str],
    sample_vcf_path: str,
    reference_vcf_path: str,
    threads: str | int,
    run_command,
) -> None:
    if not temp_vcfs:
        return

    threads_str = str(threads)
    run_command(
        ["bcftools", "concat", *temp_vcfs, "--threads", threads_str, "-Oz", "-o", f"{sample_vcf_path}.gz"]
    )
    run_command(["tabix", "-f", "-p", "vcf", f"{sample_vcf_path}.gz"])
    run_command(
        [
            "bcftools",
            "concat",
            "-a",
            "--threads",
            threads_str,
            f"{reference_vcf_path}.gz",
            f"{sample_vcf_path}.gz",
            "-Ov",
            "-o",
            reference_vcf_path,
        ]
    )
    run_command(["bcftools", "sort", reference_vcf_path, "-Oz", "-o", f"{reference_vcf_path}.gz"])
    run_command(["bcftools", "norm", f"{reference_vcf_path}.gz", "-m", "+any", "-Ov", "-o", reference_vcf_path])
    run_command(["bgzip", "-f", reference_vcf_path])
    run_command(["tabix", "-f", "-p", "vcf", f"{reference_vcf_path}.gz"])
