from __future__ import annotations

import importlib.util
import io
import json
import sys
from pathlib import Path

MODULE_PATH = Path(__file__).resolve().parent.parent / "workflow" / "scripts" / "import_pubmlst_scheme.py"
SPEC = importlib.util.spec_from_file_location("import_pubmlst_scheme", MODULE_PATH)
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


def _request_url(value) -> str:
    return value.full_url if hasattr(value, "full_url") else value


def test_list_schemes_reads_pubmlst_payload(monkeypatch) -> None:
    payload = [{"scheme": "https://rest.pubmlst.org/db/demo/schemes/42", "description": "wgMLST"}]

    def fake_urlopen(url: str, timeout=None):
        resolved = _request_url(url)
        assert resolved.endswith("/db/demo/schemes")
        return FakeResponse(json.dumps(payload).encode("utf-8"))

    monkeypatch.setattr(MODULE, "urlopen", fake_urlopen)
    assert MODULE.list_schemes("demo") == payload


def test_filter_schemes_by_type_matches_wgmlst_and_cgmlst() -> None:
    schemes = [
        {"description": "Neisseria wgMLST scheme"},
        {"description": "Core genome cgMLST scheme"},
        {"description": "Seven gene MLST"},
    ]
    assert len(MODULE.filter_schemes_by_type(schemes, "wgmlst")) == 1
    assert len(MODULE.filter_schemes_by_type(schemes, "cgmlst")) == 1


def test_download_scheme_writes_locus_fastas_and_metadata(tmp_path: Path, monkeypatch) -> None:
    responses = {
        "https://rest.pubmlst.org/db/demo/schemes/7": {
            "description": "Demo wgMLST",
            "locus_count": 2,
            "profiles_csv": "https://rest.pubmlst.org/db/demo/schemes/7/profiles_csv",
        },
        "https://rest.pubmlst.org/db/demo/schemes/7/loci": {
            "loci": [
                "https://rest.pubmlst.org/db/demo/loci/abcZ",
                "https://rest.pubmlst.org/db/demo/loci/adk",
            ]
        },
        "https://rest.pubmlst.org/db/demo/loci/abcZ/alleles_fasta": ">abcZ_1\nATG\n",
        "https://rest.pubmlst.org/db/demo/loci/adk/alleles_fasta": ">adk_1\nATG\n",
        "https://rest.pubmlst.org/db/demo/schemes/7/profiles_csv": "ST\tabcZ\tadk\n1\t1\t1\n",
    }

    def fake_urlopen(url: str, timeout=None):
        payload = responses[_request_url(url)]
        if isinstance(payload, dict):
            return FakeResponse(json.dumps(payload).encode("utf-8"))
        return FakeResponse(payload.encode("utf-8"))

    monkeypatch.setattr(MODULE, "urlopen", fake_urlopen)

    downloaded = MODULE.download_scheme("demo", 7, tmp_path, include_profiles=True)

    assert (tmp_path / "abcZ.fasta").read_text(encoding="utf-8") == ">abcZ_1\nATG\n"
    assert (tmp_path / "adk.fasta").read_text(encoding="utf-8") == ">adk_1\nATG\n"
    metadata = json.loads((tmp_path / "scheme_metadata.json").read_text(encoding="utf-8"))
    assert metadata["description"] == "Demo wgMLST"
    assert metadata["scheme_id"] == 7
    manifest = json.loads((tmp_path / "scheme_manifest.json").read_text(encoding="utf-8"))
    assert manifest["schema_type"] == "wgmlst"
    assert manifest["locus_count"] == 2
    assert (tmp_path / "profiles.tsv").exists()
    assert len(downloaded) == 3


def test_download_scheme_rejects_mismatched_expected_scheme_type(tmp_path: Path, monkeypatch) -> None:
    responses = {
        "https://rest.pubmlst.org/db/demo/schemes/7": {
            "description": "Demo cgMLST",
            "locus_count": 1,
        },
        "https://rest.pubmlst.org/db/demo/schemes/7/loci": {
            "loci": ["https://rest.pubmlst.org/db/demo/loci/abcZ"]
        },
    }

    def fake_urlopen(url: str, timeout=None):
        payload = responses[_request_url(url)]
        return FakeResponse(json.dumps(payload).encode("utf-8"))

    monkeypatch.setattr(MODULE, "urlopen", fake_urlopen)

    try:
        MODULE.download_scheme("demo", 7, tmp_path, expected_scheme_type="wgmlst")
    except ValueError as exc:
        assert "not 'wgmlst'" in str(exc)
    else:
        raise AssertionError("Expected scheme type mismatch to raise ValueError")


def test_locus_name_from_uri_extracts_terminal_segment() -> None:
    assert MODULE.locus_name_from_uri("https://rest.pubmlst.org/db/demo/loci/abcZ") == "abcZ"


def test_validate_base_url_rejects_non_http_schemes() -> None:
    try:
        MODULE.validate_base_url("file:///tmp/demo")
    except ValueError as exc:
        assert "Unsupported base URL scheme" in str(exc)
    else:
        raise AssertionError("Expected invalid base URL scheme to raise ValueError")


def test_validate_database_name_rejects_invalid_characters() -> None:
    try:
        MODULE.validate_database_name("ci/reference")
    except ValueError as exc:
        assert "letters, numbers, underscores, and hyphens" in str(exc)
    else:
        raise AssertionError("Expected invalid database name to raise ValueError")


def test_fetch_json_surfaces_clear_not_found_error(monkeypatch) -> None:
    from urllib.error import HTTPError

    def fake_urlopen(url: str, timeout=None):
        raise HTTPError(url, 404, "Not Found", hdrs=None, fp=None)

    monkeypatch.setattr(MODULE, "urlopen", fake_urlopen)

    try:
        MODULE.fetch_json("https://rest.pubmlst.org/db/demo/schemes/7")
    except RuntimeError as exc:
        assert "may be invalid" in str(exc)
    else:
        raise AssertionError("Expected remote fetch failure to raise RuntimeError")


def test_fetch_url_retries_on_temporary_network_errors(monkeypatch) -> None:
    calls = {"count": 0}

    def fake_urlopen(request, timeout=None):
        calls["count"] += 1
        if calls["count"] < 3:
            raise MODULE.URLError("temporary failure")
        return FakeResponse(b"[]")

    monkeypatch.setattr(MODULE, "urlopen", fake_urlopen)
    monkeypatch.setattr(MODULE.time, "sleep", lambda _: None)

    payload = MODULE.fetch_json("https://rest.pubmlst.org/db/demo/schemes")
    assert payload == []
    assert calls["count"] == 3


def test_download_locus_alleles_sanitizes_output_filename(tmp_path: Path, monkeypatch) -> None:
    def fake_fetch_text(url: str) -> str:
        return ">abc_1\nATG\n"

    monkeypatch.setattr(MODULE, "fetch_text", fake_fetch_text)
    result = MODULE.download_locus_alleles("demo", "https://rest.pubmlst.org/db/demo/loci/abc:Z", tmp_path)

    assert result.name == "abc_Z.fasta"
    assert result.read_text(encoding="utf-8") == ">abc_1\nATG\n"


def test_build_schema_creation_command_points_to_milestone_workflow(tmp_path: Path) -> None:
    command = MODULE.build_schema_creation_command(
        schema_dir=tmp_path / "schema",
        reference="demo_ref",
        output_dir=tmp_path / "output",
        threads=8,
    )
    assert command[1].endswith("workflow/milestone.py")
    assert command[-1] == "8"


def test_run_schema_creation_returns_subprocess_exit_code(tmp_path: Path, monkeypatch) -> None:
    calls = []

    class FakeCompleted:
        returncode = 0

    def fake_run(command, check=False):
        calls.append({"command": command, "check": check})
        return FakeCompleted()

    monkeypatch.setattr(MODULE.subprocess, "run", fake_run)
    returncode = MODULE.run_schema_creation(tmp_path / "schema", "demo_ref", tmp_path / "out", 2)
    assert returncode == 0
    assert calls and calls[0]["check"] is False


def test_download_scheme_accepts_manifest_overrides(tmp_path: Path, monkeypatch) -> None:
    responses = {
        "https://rest.pubmlst.org/db/demo/schemes/7": {
            "description": "Demo wgMLST",
            "locus_count": 1,
        },
        "https://rest.pubmlst.org/db/demo/schemes/7/loci": {
            "loci": ["https://rest.pubmlst.org/db/demo/loci/abcZ"]
        },
        "https://rest.pubmlst.org/db/demo/loci/abcZ/alleles_fasta": ">abcZ_1\nATG\n",
    }

    def fake_urlopen(url: str, timeout=None):
        payload = responses[_request_url(url)]
        if isinstance(payload, dict):
            return FakeResponse(json.dumps(payload).encode("utf-8"))
        return FakeResponse(payload.encode("utf-8"))

    monkeypatch.setattr(MODULE, "urlopen", fake_urlopen)

    MODULE.download_scheme(
        "demo",
        7,
        tmp_path,
        schema_name="Custom Demo",
        species="Neisseria meningitidis",
        schema_version="2026.1",
    )

    manifest = json.loads((tmp_path / "scheme_manifest.json").read_text(encoding="utf-8"))
    assert manifest["schema_name"] == "Custom Demo"
    assert manifest["species"] == "Neisseria meningitidis"
    assert manifest["schema_version"] == "2026.1"
