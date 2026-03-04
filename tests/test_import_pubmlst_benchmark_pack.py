from __future__ import annotations

import importlib.util
import io
import json
import sys
from pathlib import Path

MODULE_PATH = Path(__file__).resolve().parent.parent / "workflow" / "scripts" / "import_pubmlst_benchmark_pack.py"
SPEC = importlib.util.spec_from_file_location("import_pubmlst_benchmark_pack", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class FakeResponse(io.BytesIO):
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        self.close()
        return False


def test_import_benchmark_pack_writes_truth_profiles_and_manifest(tmp_path: Path, monkeypatch) -> None:
    responses = {
        "https://rest.pubmlst.org/db/demo_seqdef/schemes/7": {"description": "Demo wgMLST"},
        "https://rest.pubmlst.org/db/demo_isolates/isolates/101": {"id": "iso-101"},
        "https://rest.pubmlst.org/db/demo_isolates/isolates/102": {"id": "iso-102"},
        "https://rest.pubmlst.org/db/demo_isolates/isolates/101/schemes/7/allele_ids": {"abc": 1, "adk": 2},
        "https://rest.pubmlst.org/db/demo_isolates/isolates/102/schemes/7/allele_ids": {"abc": 1, "adk": 7},
    }

    def fake_urlopen(url: str, timeout=None):
        payload = responses[url]
        return FakeResponse(json.dumps(payload).encode("utf-8"))

    monkeypatch.setattr(MODULE, "urlopen", fake_urlopen)

    downloaded = MODULE.import_benchmark_pack(
        scheme_database="demo_seqdef",
        isolate_database="demo_isolates",
        scheme_id=7,
        isolate_ids=[101, 102],
        output_dir=tmp_path,
        species="Neisseria meningitidis",
    )

    assert (tmp_path / "truth" / "iso-101_wgmlst.tsv").read_text(encoding="utf-8") == "abc\t1\nadk\t2\n"
    assert (tmp_path / "predicted" / "README.md").exists()
    manifest = json.loads((tmp_path / "benchmark_manifest.json").read_text(encoding="utf-8"))
    assert manifest["species"] == "Neisseria meningitidis"
    assert manifest["truth_profile_count"] == 2
    assert manifest["public_profile_count"] == 2
    assert manifest["selection_criteria"] == "explicit_isolate_ids"
    assert "retrieval_date" in manifest
    assert len(downloaded) >= 4


def test_import_benchmark_pack_can_download_contigs(tmp_path: Path, monkeypatch) -> None:
    responses = {
        "https://rest.pubmlst.org/db/demo_seqdef/schemes/7": {"description": "Demo wgMLST"},
        "https://rest.pubmlst.org/db/demo_isolates/isolates/101": {"id": "iso-101"},
        "https://rest.pubmlst.org/db/demo_isolates/isolates/101/schemes/7/allele_ids": {"abc": 1},
        "https://rest.pubmlst.org/db/demo_isolates/isolates/101/contigs_fasta": ">contig1\nATGC\n",
    }

    def fake_urlopen(url: str, timeout=None):
        payload = responses[url]
        if isinstance(payload, dict):
            return FakeResponse(json.dumps(payload).encode("utf-8"))
        return FakeResponse(payload.encode("utf-8"))

    monkeypatch.setattr(MODULE, "urlopen", fake_urlopen)

    MODULE.import_benchmark_pack(
        scheme_database="demo_seqdef",
        isolate_database="demo_isolates",
        scheme_id=7,
        isolate_ids=[101],
        output_dir=tmp_path,
        download_contigs=True,
    )

    assert (tmp_path / "assemblies" / "iso-101.fasta").read_text(encoding="utf-8") == ">contig1\nATGC\n"


def test_fetch_isolate_ids_reads_first_n_uris(monkeypatch) -> None:
    payload = {
        "isolates": [
            "https://rest.pubmlst.org/db/demo_isolates/isolates/101",
            "https://rest.pubmlst.org/db/demo_isolates/isolates/102",
        ]
    }

    def fake_urlopen(url: str, timeout=None):
        assert url.endswith("/db/demo_isolates/isolates?page_size=2")
        return FakeResponse(json.dumps(payload).encode("utf-8"))

    monkeypatch.setattr(MODULE, "urlopen", fake_urlopen)

    assert MODULE.fetch_isolate_ids("demo_isolates", "https://rest.pubmlst.org", 2) == [101, 102]


def test_main_supports_discover_isolates_without_explicit_ids(monkeypatch, tmp_path: Path) -> None:
    monkeypatch.setattr(
        MODULE,
        "fetch_isolate_ids",
        lambda database, base_url, limit: [101, 102],
    )
    captured = {}

    def fake_import_benchmark_pack(**kwargs):
        captured.update(kwargs)
        return []

    monkeypatch.setattr(MODULE, "import_benchmark_pack", fake_import_benchmark_pack)
    monkeypatch.setattr(
        MODULE,
        "parse_args",
        lambda: MODULE.argparse.Namespace(
            database=None,
            scheme_database="demo_seqdef",
            isolate_database="demo_isolates",
            scheme_id=7,
            isolate_ids=[],
            discover_isolates=2,
            output_dir=tmp_path,
            base_url="https://rest.pubmlst.org",
            species=None,
            source_name="fixture",
            selection_criteria=None,
            download_contigs=False,
        ),
    )

    assert MODULE.main() == 0
    assert captured["isolate_ids"] == [101, 102]
