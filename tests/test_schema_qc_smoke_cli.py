from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPT = REPO_ROOT / "workflow" / "scripts" / "schema_qc.py"


def test_schema_qc_cli_smoke(tmp_path: Path) -> None:
    schema_dir = tmp_path / "schema"
    schema_dir.mkdir()
    (schema_dir / "abc.fasta").write_text(">abc_1\nATGAAATAA\n>abc_2\nATGCCCCCCCAA\n", encoding="ascii")
    output_dir = tmp_path / "out"

    completed = subprocess.run(
        [sys.executable, str(SCRIPT), "--schema-dir", str(schema_dir), "--output-dir", str(output_dir)],
        check=False,
        capture_output=True,
        text=True,
        cwd=REPO_ROOT,
    )

    assert completed.returncode == 0, completed.stderr
    summary = json.loads((output_dir / "schema_qc_summary.json").read_text(encoding="utf-8"))
    assert summary["locus_count"] == 1
