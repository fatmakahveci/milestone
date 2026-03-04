from __future__ import annotations

import importlib.util
from argparse import Namespace
from pathlib import Path

MODULE_PATH = Path(__file__).resolve().parent.parent / "workflow" / "scripts" / "create_sample_info_txt.py"
SPEC = importlib.util.spec_from_file_location("create_sample_info_txt", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def test_get_reference_allele_name_prefers_allele_one() -> None:
    cds_seq_dict = {"abc_1": "ATGTAA", "abc_7": "ATGAAATAA"}
    assert MODULE.get_reference_allele_name("abc", cds_seq_dict) == "abc_1"


def test_get_reference_allele_name_uses_unique_reference_candidate() -> None:
    cds_seq_dict = {"abc_ref": "ATGTAA", "xyz_1": "ATGCCCTAA"}
    assert MODULE.get_reference_allele_name("abc", cds_seq_dict) == "abc_ref"


def test_get_reference_allele_id_uses_inferred_reference_candidate() -> None:
    cds_seq_dict = {"abc_ref": "ATGTAA", "xyz_1": "ATGCCCTAA"}
    assert MODULE.get_reference_allele_id("abc", cds_seq_dict) == "ref"


def test_get_reference_allele_name_rejects_ambiguous_reference_candidates() -> None:
    cds_seq_dict = {"abc_refA": "ATGTAA", "abc_refB": "ATGAAATAA"}
    try:
        MODULE.get_reference_allele_name("abc", cds_seq_dict)
    except ValueError as exc:
        assert "Unable to determine a unique reference allele" in str(exc)
    else:
        raise AssertionError("Expected ValueError for ambiguous reference candidates")


def test_quality_check_accepts_common_bacterial_start_codons() -> None:
    sequence = "GTGAAACCCGGGTAA"
    assert MODULE.quality_check(sequence, "ATGAAACCCGGGTAA", locus_mode_length=len(sequence)) == "Q"


def test_quality_check_classifies_short_and_long_alleles_against_mode() -> None:
    assert MODULE.quality_check("ATGAAATAA", "ATGCCCCCCCTAA", locus_mode_length=12) == "ASM"
    assert MODULE.quality_check("ATGAAAAAACCCTAA", "ATGCCCCCCCTAA", locus_mode_length=12) == "ALM"


def test_quality_check_rejects_internal_stop_codons() -> None:
    assert MODULE.quality_check("ATGTAACCCTAA", "ATGAAACCCTAA", locus_mode_length=12) == "LNF"


def test_quality_check_rejects_invalid_reference_sequence_even_when_equal() -> None:
    assert MODULE.quality_check("ATGTAACCCTAA", "ATGTAACCCTAA", locus_mode_length=12) == "LNF"


def test_quality_check_respects_translation_table_four_for_tga() -> None:
    assert MODULE.quality_check(
        "ATGTGACCCTAA",
        "ATGAAACCCTAA",
        locus_mode_length=12,
        translation_table=4,
    ) == "Q"


def test_quality_check_respects_translation_table_two_for_aga_stop() -> None:
    assert MODULE.quality_check(
        "ATGAAAAGATAA",
        "ATGAAACCCTAA",
        locus_mode_length=12,
        translation_table=2,
    ) == "LNF"


def test_quality_check_respects_translation_table_twentyfive_for_tga_sense() -> None:
    assert MODULE.quality_check(
        "ATGTGACCCTAA",
        "ATGAAACCCTAA",
        locus_mode_length=12,
        translation_table=25,
    ) == "Q"


def test_quality_check_allows_custom_start_codons_override() -> None:
    assert MODULE.quality_check(
        "ATAAAACCCTAA",
        "ATGAAACCCTAA",
        locus_mode_length=12,
        allowed_start_codons={"ATA"},
    ) == "Q"


def test_get_locus_length_mode_uses_schema_fasta_lengths(tmp_path: Path) -> None:
    schema_dir = tmp_path / "schema"
    schema_dir.mkdir()
    (schema_dir / "abc.fasta").write_text(
        ">abc_1\nATGAAATAA\n>abc_2\nATGCCCCCCCAA\n>abc_3\nATGTTTTTTTAA\n",
        encoding="ascii",
    )
    MODULE.args = Namespace(schema_dir=str(schema_dir))
    assert MODULE.get_locus_length_mode("abc", "ATGAAATAA") == 12


def test_get_locus_length_mode_falls_back_to_reference_length(tmp_path: Path) -> None:
    MODULE.args = Namespace(schema_dir=str(tmp_path / "missing"))
    assert MODULE.get_locus_length_mode("abc", "ATGAAATAA") == 9


def test_get_next_novel_allele_id_uses_schema_fasta_not_reference_info(tmp_path: Path) -> None:
    schema_dir = tmp_path / "schema"
    schema_dir.mkdir()
    (schema_dir / "abc.fasta").write_text(
        ">abc_ref\nATGAAATAA\n>abc_2\nATGCCCCCCCAA\n>abc_7\nATGTTTTTTTAA\n",
        encoding="ascii",
    )
    (tmp_path / "reference_info.txt").write_text("", encoding="ascii")
    MODULE.args = Namespace(schema_dir=str(schema_dir), reference_info=str(tmp_path / "reference_info.txt"))

    assert MODULE.get_next_novel_allele_id("abc") == 8
    assert MODULE.get_allele_ids_of_cds_in_reference_info_txt()["abc"] == 8


def test_get_reference_alignment_target_name_uses_unique_sam_target() -> None:
    MODULE.sample_sam_dict = {"abc_ref": [], "xyz_1": []}
    assert MODULE.get_reference_alignment_target_name("abc") == "abc_ref"


def test_locus_name_parsing_supports_internal_underscores() -> None:
    assert MODULE.get_cds_name_from_allele_name("abc_def_12") == "abc_def"
    assert MODULE.get_allele_id_from_allele_name("abc_def_12") == "12"


def test_write_variations_to_reference_vcf_file_respects_non_numeric_reference_name(
    tmp_path: Path, monkeypatch
) -> None:
    sample_vcf = tmp_path / "sample.vcf"
    sample_vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        "abc_ref\t3\t.\tG\tA\t60\tPASS\tTYPE=snp\tGT\t1\n",
        encoding="ascii",
    )
    temp_dir = tmp_path / "temp"
    temp_dir.mkdir()
    MODULE.args = Namespace(sample_vcf=str(sample_vcf))

    commands = []

    def fake_run_command(command):
        commands.append(command)

    monkeypatch.setattr(MODULE, "run_command", fake_run_command)

    MODULE.write_variations_to_reference_vcf_file(
        cds="abc",
        reference_allele_name="abc_ref",
        temp_sample_vcf_dir=str(temp_dir),
        cds_variation=MODULE.Info("3*G>A-60"),
    )

    created = (temp_dir / "abc.vcf").read_text(encoding="ascii")
    assert "abc_ref\t3\t.\tG\tA\t60.0\tPASS\t.\tGT\t1" in created
    assert commands == [
        ["bgzip", "-f", str(temp_dir / "abc.vcf")],
        ["tabix", "-p", "vcf", str(temp_dir / "abc.vcf.gz")],
    ]


def test_write_novel_allele_exports_creates_fasta_and_tsv(tmp_path: Path) -> None:
    MODULE.args = Namespace(sample_vcf=str(tmp_path / "sample.vcf"))

    fasta_path, tsv_path = MODULE.write_novel_allele_exports(
        [
            {
                "locus": "abc",
                "allele_id": "9",
                "sequence": "ATGAAATAA",
                "status": "novel_candidate",
                "reference_allele": "abc_1",
                "variation_summary": "3:G>A",
                "orf_reason": "passes_orf_qc",
            }
        ]
    )

    assert ">abc_9" in fasta_path.read_text(encoding="utf-8")
    assert "abc\t9\t9\tnovel_candidate\tabc_1\t3:G>A\tpasses_orf_qc" in tsv_path.read_text(encoding="utf-8")


def test_quality_check_details_returns_reason() -> None:
    details = MODULE.quality_check_details("ATGTAACCCTAA", "ATGAAACCCTAA", locus_mode_length=12)
    assert details["status"] == "LNF"
    assert details["reason"] == "internal_stop_codon"


def test_compare_ref_to_sample_variations_compares_linked_variant_tuples(tmp_path: Path) -> None:
    schema_dir = tmp_path / "schema"
    schema_dir.mkdir()
    (schema_dir / "abc.fasta").write_text(">abc_1\nATGAAACCCGGGTAA\n", encoding="ascii")
    MODULE.args = Namespace(schema_dir=str(schema_dir), allowed_start_codons=None, translation_table=11)
    MODULE.novel_allele_id_of_cds_dict = {"abc": 9}

    allele_id, is_novel = MODULE.compare_ref_to_sample_variations(
        "abc",
        {"abc_1": "ATGAAACCCGGGTAA"},
        {
            "2": MODULE.Info("4*A>C-60,7*C>G-60"),
            "3": MODULE.Info("4*A>G-60,7*C>A-60"),
        },
        MODULE.Info("4*A>C-60,7*C>A-60"),
    )

    assert allele_id == "9"
    assert is_novel is True
