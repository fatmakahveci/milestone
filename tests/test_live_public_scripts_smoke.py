from __future__ import annotations

from pathlib import Path


def test_live_public_import_smoke_script_exists() -> None:
    script = Path(__file__).resolve().parent / "run_live_public_import_smoke.sh"
    assert script.exists()
    assert script.read_text(encoding="utf-8").startswith("#!/usr/bin/env bash")
