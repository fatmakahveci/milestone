from __future__ import annotations

import importlib.util
import io
import json
import sys
from pathlib import Path


MODULE_PATH = Path(__file__).resolve().parent.parent / "workflow" / "scripts" / "import_enterobase_scheme.py"
SPEC = importlib.util.spec_from_file_location("import_enterobase_scheme", MODULE_PATH)
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


def test_list_schemes_reads_enterobase_api_payload(monkeypatch) -> None:
    payload = {
        "Schemes": [
            {"scheme_name": "wgMLST", "label": "EnteroBase wgMLST", "download_sts_link": "https://example.org/wg.tar.gz"}
        ]
    }

    def fake_urlopen(url: str, timeout=None):
        resolved = _request_url(url)
        assert resolved.endswith("/api/v2.0/senterica/schemes?limit=1000")
        return FakeResponse(json.dumps(payload).encode("utf-8"))

    monkeypatch.setattr(MODULE, "urlopen", fake_urlopen)
    schemes = MODULE.list_schemes("senterica")
    assert schemes[0]["scheme_name"] == "wgMLST"


def test_download_scheme_profiles_writes_archive_and_metadata(tmp_path: Path, monkeypatch) -> None:
    payload = {
        "Schemes": [
            {
                "scheme_name": "wgMLST",
                "label": "EnteroBase wgMLST",
                "download_sts_link": "https://example.org/wg.tar.gz",
            }
        ]
    }

    def fake_urlopen(url: str, timeout=None):
        resolved = _request_url(url)
        if resolved.endswith("/api/v2.0/senterica/schemes?limit=1000"):
            return FakeResponse(json.dumps(payload).encode("utf-8"))
        if resolved == "https://example.org/wg.tar.gz":
            return FakeResponse(b"archive-bytes")
        raise AssertionError(f"Unexpected URL: {resolved}")

    monkeypatch.setattr(MODULE, "urlopen", fake_urlopen)

    downloaded = MODULE.download_scheme_profiles("senterica", "wgMLST", tmp_path)

    assert (tmp_path / "wgMLST_profiles.tar.gz").read_bytes() == b"archive-bytes"
    metadata = json.loads((tmp_path / "enterobase_scheme_metadata.json").read_text(encoding="utf-8"))
    assert metadata["scheme_name"] == "wgMLST"
    assert len(downloaded) == 2


def test_validate_database_name_rejects_invalid_characters() -> None:
    try:
        MODULE.validate_database_name("senterica/api")
    except ValueError as exc:
        assert "letters, numbers, underscores, and hyphens" in str(exc)
    else:
        raise AssertionError("Expected invalid database name to raise ValueError")


def test_fetch_json_surfaces_clear_unauthorized_error(monkeypatch) -> None:
    from urllib.error import HTTPError

    def fake_urlopen(url: str, timeout=None):
        raise HTTPError(url, 401, "Unauthorized", hdrs=None, fp=None)

    monkeypatch.setattr(MODULE, "urlopen", fake_urlopen)

    try:
        MODULE.fetch_json("https://enterobase.warwick.ac.uk/api/v2.0/senterica/schemes?limit=1000")
    except RuntimeError as exc:
        assert "require authentication" in str(exc) or "allow anonymous access" in str(exc)
    else:
        raise AssertionError("Expected unauthorized fetch to raise RuntimeError")


def test_build_headers_supports_bearer_and_basic_auth(monkeypatch) -> None:
    monkeypatch.setenv("ENTEROBASE_TOKEN", "secret-token")
    headers = MODULE.build_headers()
    assert headers["Authorization"] == "Bearer secret-token"

    monkeypatch.delenv("ENTEROBASE_TOKEN")
    monkeypatch.setenv("ENTEROBASE_USERNAME", "demo")
    monkeypatch.setenv("ENTEROBASE_PASSWORD", "pass")
    headers = MODULE.build_headers()
    assert headers["Authorization"].startswith("Basic ")


def test_fetch_url_retries_on_temporary_network_errors(monkeypatch) -> None:
    calls = {"count": 0}

    def fake_urlopen(request, timeout=None):
        calls["count"] += 1
        if calls["count"] < 3:
            raise MODULE.URLError("temporary failure")
        return FakeResponse(b'{"Schemes": []}')

    monkeypatch.setattr(MODULE, "urlopen", fake_urlopen)
    monkeypatch.setattr(MODULE.time, "sleep", lambda _: None)
    payload = MODULE.fetch_json("https://enterobase.warwick.ac.uk/api/v2.0/senterica/schemes?limit=1000")
    assert payload == {"Schemes": []}
    assert calls["count"] == 3
