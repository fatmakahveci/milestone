import importlib.util
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent.parent / "workflow" / "scripts"
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))


MODULE_PATH = Path(__file__).resolve().parent.parent / "workflow" / "scripts" / "script_utils.py"
SPEC = importlib.util.spec_from_file_location("script_utils", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)

VARIATION_MODULE_PATH = Path(__file__).resolve().parent.parent / "workflow" / "scripts" / "variation_utils.py"
VARIATION_SPEC = importlib.util.spec_from_file_location("variation_utils", VARIATION_MODULE_PATH)
VARIATION_MODULE = importlib.util.module_from_spec(VARIATION_SPEC)
assert VARIATION_SPEC.loader is not None
VARIATION_SPEC.loader.exec_module(VARIATION_MODULE)

REFERENCE_FORMATS_PATH = Path(__file__).resolve().parent.parent / "workflow" / "scripts" / "reference_formats.py"
REFERENCE_FORMATS_SPEC = importlib.util.spec_from_file_location("reference_formats", REFERENCE_FORMATS_PATH)
REFERENCE_FORMATS_MODULE = importlib.util.module_from_spec(REFERENCE_FORMATS_SPEC)
assert REFERENCE_FORMATS_SPEC.loader is not None
REFERENCE_FORMATS_SPEC.loader.exec_module(REFERENCE_FORMATS_MODULE)

SAMPLE_FORMATS_PATH = Path(__file__).resolve().parent.parent / "workflow" / "scripts" / "sample_formats.py"
SAMPLE_FORMATS_SPEC = importlib.util.spec_from_file_location("sample_formats", SAMPLE_FORMATS_PATH)
SAMPLE_FORMATS_MODULE = importlib.util.module_from_spec(SAMPLE_FORMATS_SPEC)
assert SAMPLE_FORMATS_SPEC.loader is not None
SAMPLE_FORMATS_SPEC.loader.exec_module(SAMPLE_FORMATS_MODULE)

REFERENCE_UPDATE_IO_PATH = Path(__file__).resolve().parent.parent / "workflow" / "scripts" / "reference_update_io.py"
REFERENCE_UPDATE_IO_SPEC = importlib.util.spec_from_file_location("reference_update_io", REFERENCE_UPDATE_IO_PATH)
REFERENCE_UPDATE_IO_MODULE = importlib.util.module_from_spec(REFERENCE_UPDATE_IO_SPEC)
assert REFERENCE_UPDATE_IO_SPEC.loader is not None
REFERENCE_UPDATE_IO_SPEC.loader.exec_module(REFERENCE_UPDATE_IO_MODULE)

REFERENCE_MERGE_PATH = Path(__file__).resolve().parent.parent / "workflow" / "scripts" / "reference_merge_pipeline.py"
REFERENCE_MERGE_SPEC = importlib.util.spec_from_file_location("reference_merge_pipeline", REFERENCE_MERGE_PATH)
REFERENCE_MERGE_MODULE = importlib.util.module_from_spec(REFERENCE_MERGE_SPEC)
assert REFERENCE_MERGE_SPEC.loader is not None
REFERENCE_MERGE_SPEC.loader.exec_module(REFERENCE_MERGE_MODULE)

ORCHESTRATOR_PATH = Path(__file__).resolve().parent.parent / "workflow" / "scripts" / "sample_call_orchestrator.py"
ORCHESTRATOR_SPEC = importlib.util.spec_from_file_location("sample_call_orchestrator", ORCHESTRATOR_PATH)
ORCHESTRATOR_MODULE = importlib.util.module_from_spec(ORCHESTRATOR_SPEC)
assert ORCHESTRATOR_SPEC.loader is not None
ORCHESTRATOR_SPEC.loader.exec_module(ORCHESTRATOR_MODULE)


def test_parsed_variation_info_parses_multiple_variants() -> None:
    info = MODULE.ParsedVariationInfo("3*G>A-60,7*C>T-45")

    assert info.pos_list == [3, 7]
    assert info.ref_list == ["G", "C"]
    assert info.alt_list == ["A", "T"]
    assert info.qual_list == ["60", "45"]


def test_variation_utils_builds_linked_signature() -> None:
    info = MODULE.ParsedVariationInfo("3*G>A-60,7*C>T-45")
    signature = VARIATION_MODULE.build_variation_signature(info)

    assert signature[(3, "G", "A")] == 1
    assert signature[(7, "C", "T")] == 1


def test_reference_formats_vcf_parses_locus_name() -> None:
    vcf = REFERENCE_FORMATS_MODULE.Vcf("abc_def_12\t4\t.\tA\tC\t60\tPASS\t.\tGT\t1")

    assert vcf.chr == "abc_def"
    assert vcf.alt == "C"


def test_sample_formats_coverage_and_vcf_parse() -> None:
    coverage = SAMPLE_FORMATS_MODULE.Coverage("abc\t.\t10\t7\t8\t80.0\t12.5\t36.1\t48.0")
    vcf = SAMPLE_FORMATS_MODULE.Vcf("abc_def_12\t4\t.\tA\tC,G\t60\tPASS\tTYPE=snp\tGT\t1")

    assert coverage.cds == "abc"
    assert coverage.coverage == 80.0
    assert vcf.chr == "abc_def"
    assert vcf.alt == "C"


def test_reference_update_io_writes_reference_info(tmp_path: Path) -> None:
    info_path = tmp_path / "reference_info.txt"
    variation = MODULE.ParsedVariationInfo("3*G>A-60,7*C>T-45")

    REFERENCE_UPDATE_IO_MODULE.write_variations_to_reference_info_file(
        reference_info_path=str(info_path),
        cds="abc",
        allele_id="9",
        cds_variation=variation,
    )

    content = info_path.read_text(encoding="utf-8")
    assert content.strip() == "abc_9\t3*G>A-60,7*C>T-45"


def test_reference_merge_pipeline_emits_expected_command_sequence() -> None:
    commands = []

    def fake_run_command(command):
        commands.append(command)

    REFERENCE_MERGE_MODULE.merge_reference_with_novel_vcfs(
        temp_vcfs=["tmp/a.vcf.gz", "tmp/b.vcf.gz"],
        sample_vcf_path="sample.vcf",
        reference_vcf_path="reference.vcf",
        threads=4,
        run_command=fake_run_command,
    )

    assert commands == [
        ["bcftools", "concat", "tmp/a.vcf.gz", "tmp/b.vcf.gz", "--threads", "4", "-Oz", "-o", "sample.vcf.gz"],
        ["tabix", "-f", "-p", "vcf", "sample.vcf.gz"],
        ["bcftools", "concat", "-a", "--threads", "4", "reference.vcf.gz", "sample.vcf.gz", "-Ov", "-o", "reference.vcf"],
        ["bcftools", "sort", "reference.vcf", "-Oz", "-o", "reference.vcf.gz"],
        ["bcftools", "norm", "reference.vcf.gz", "-m", "+any", "-Ov", "-o", "reference.vcf"],
        ["bgzip", "-f", "reference.vcf"],
        ["tabix", "-f", "-p", "vcf", "reference.vcf.gz"],
    ]


def test_sample_call_orchestrator_writes_output(tmp_path: Path) -> None:
    init_calls = []

    def init_state() -> None:
        init_calls.append("init")

    def collect_alleles() -> dict[str, str]:
        return {"abc_1": "1", "xyz_4": "Q"}

    result = ORCHESTRATOR_MODULE.run_sample_calling_workflow(
        sample_vcf_path=str(tmp_path / "sample.vcf"),
        init_state=init_state,
        collect_allele_calls=collect_alleles,
        get_cds_name_from_allele_name=lambda name: name.split("_")[0],
    )

    content = Path(result.output_path).read_text(encoding="utf-8")
    assert init_calls == ["init"]
    assert "abc\t1" in content
    assert "xyz\tQ" in content
