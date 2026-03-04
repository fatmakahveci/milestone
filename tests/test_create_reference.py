from __future__ import annotations

import importlib.util
import sys
import types
from argparse import Namespace
from dataclasses import dataclass
from pathlib import Path

MODULE_PATH = Path(__file__).resolve().parent.parent / "workflow" / "scripts" / "create_reference.py"


def fasta_parse(handle, _format: str):
    header = None
    sequence_chunks: list[str] = []
    for raw_line in handle:
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                yield FakeRecord(header, "".join(sequence_chunks))
            header = line[1:]
            sequence_chunks = []
            continue
        sequence_chunks.append(line)
    if header is not None:
        yield FakeRecord(header, "".join(sequence_chunks))


bio_module = types.ModuleType("Bio")
bio_module.SeqIO = types.SimpleNamespace(parse=fasta_parse)
sys.modules.setdefault("Bio", bio_module)

SPEC = importlib.util.spec_from_file_location("create_reference", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


@dataclass
class FakeRecord:
    id: str
    seq: str


def test_select_reference_record_prefers_allele_one() -> None:
    records = [FakeRecord("abc_7", "AAAA"), FakeRecord("abc_1", "TTTT")]
    selected = MODULE.select_reference_record(records)
    assert selected.id == "abc_1"


def test_select_reference_record_falls_back_deterministically() -> None:
    records = [FakeRecord("abc_7", "TTTT"), FakeRecord("abc_ref", "AAAA"), FakeRecord("abc_2", "CCCC")]
    selected = MODULE.select_reference_record(records)
    assert selected.id == "abc_ref"


def test_select_reference_record_supports_locus_names_with_internal_underscores() -> None:
    records = [FakeRecord("abc_def_7", "AAAA"), FakeRecord("abc_def_2", "TTTT")]
    selected = MODULE.select_reference_record(records)
    assert selected.id == "abc_def_2"


def test_create_cds_list_uses_fallback_reference_when_allele_one_is_missing(
    tmp_path: Path, monkeypatch
) -> None:
    schema_dir = tmp_path / "schema"
    schema_dir.mkdir()
    locus_fasta = schema_dir / "abc.fasta"
    locus_fasta.write_text(
        ">abc_ref\nATGAAATAA\n>abc_2\nATGCCCCCCCAA\n",
        encoding="ascii",
    )

    MODULE.args = Namespace(schema_dir=str(schema_dir), threads="1")
    MODULE.create_dirs_to_split_sequences_to_call_variants()

    captured_calls: list[tuple[str, str]] = []

    def fake_create_allele_dict_for_a_cds(
        write_dir: str,
        allele_name: str,
        cds_dir: str,
        cds_name: str,
        reference_allele_name: str,
    ) -> None:
        captured_calls.append((allele_name, reference_allele_name))

    monkeypatch.setattr(MODULE, "create_allele_dict_for_a_cds", fake_create_allele_dict_for_a_cds)

    cds_to_merge = MODULE.create_cds_list(str(schema_dir), "abc.fasta", [])

    reference_path = schema_dir / "references" / "abc_ref.fasta"
    allele_path = schema_dir / "alleles" / "abc" / "abc_2.fasta"

    assert cds_to_merge == []
    assert reference_path.exists()
    assert allele_path.exists()
    assert captured_calls == [("abc_2", "abc_ref")]
