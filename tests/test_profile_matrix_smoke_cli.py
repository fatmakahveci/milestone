from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
WORKFLOW_SCRIPT = REPO_ROOT / "workflow" / "milestone.py"
FIXTURE_DIR = Path(__file__).resolve().parent / "fixtures" / "smoke"


def test_profile_matrix_cli_smoke(tmp_path: Path) -> None:
    output_dir = tmp_path / "matrix"
    command = [
        sys.executable,
        str(WORKFLOW_SCRIPT),
        "profile_matrix",
        "--profiles",
        str(FIXTURE_DIR / "strain_a_wgmlst.tsv"),
        str(FIXTURE_DIR / "strain_b_wgmlst.tsv"),
        "--output",
        str(output_dir),
    ]

    completed = subprocess.run(command, check=False, capture_output=True, text=True, cwd=REPO_ROOT)

    assert completed.returncode == 0, completed.stderr
    assert "Samples compared: 2" in completed.stdout
    summary = json.loads((output_dir / "wgmlst_distance_summary.json").read_text(encoding="utf-8"))
    matrix = (output_dir / "wgmlst_distance_matrix.tsv").read_text(encoding="utf-8")
    assert summary["pair_count"] == 1
    assert summary["pairs"][0]["decision"] == "different"
    assert "sample\tstrain_a_wgmlst\tstrain_b_wgmlst" in matrix
