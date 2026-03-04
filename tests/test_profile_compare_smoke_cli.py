from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent.parent
WORKFLOW_SCRIPT = REPO_ROOT / "workflow" / "milestone.py"
FIXTURE_DIR = Path(__file__).resolve().parent / "fixtures" / "smoke"


def test_profile_compare_cli_smoke(tmp_path: Path) -> None:
    output_dir = tmp_path / "compare"
    command = [
        sys.executable,
        str(WORKFLOW_SCRIPT),
        "profile_compare",
        "--profile-a",
        str(FIXTURE_DIR / "strain_a_wgmlst.tsv"),
        "--profile-b",
        str(FIXTURE_DIR / "strain_b_wgmlst.tsv"),
        "--label-a",
        "strain_a",
        "--label-b",
        "strain_b",
        "--output",
        str(output_dir),
    ]

    completed = subprocess.run(command, check=False, capture_output=True, text=True, cwd=REPO_ROOT)

    assert completed.returncode == 0, completed.stderr
    assert "Decision: different" in completed.stdout
    summary = json.loads((output_dir / "wgmlst_profile_summary.json").read_text(encoding="utf-8"))
    details = (output_dir / "wgmlst_profile_comparison.tsv").read_text(encoding="utf-8")
    assert summary["decision"] == "different"
    assert summary["differing_loci"] == 1
    assert "adk\t4\t7\tdifferent" in details
