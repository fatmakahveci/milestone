from __future__ import annotations

import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPT_DIR = REPO_ROOT / "workflow" / "scripts"
FIXTURE_DIR = Path(__file__).resolve().parent / "fixtures" / "smoke"


def test_compare_results_wrapper_runs_on_legacy_merged_table(tmp_path: Path) -> None:
    merged = tmp_path / "merged.tsv"
    merged.write_text(
        "abcZ\t12\t12\n"
        "adk\t4\t7\n"
        "aroE\tLNF\t3\n",
        encoding="utf-8",
    )

    completed = subprocess.run(
        [sys.executable, str(SCRIPT_DIR / "compare_results.py"), "--file", str(merged), "--aligner", "vg"],
        check=False,
        capture_output=True,
        text=True,
        cwd=REPO_ROOT,
    )

    assert completed.returncode == 0
    assert "Decision: different" in completed.stdout


def test_comp_milestone_chewie_wrapper_runs(tmp_path: Path) -> None:
    output_dir = tmp_path / "legacy"
    completed = subprocess.run(
        [
            sys.executable,
            str(SCRIPT_DIR / "comp_milestone_chewie.py"),
            "--milestone-profile",
            str(FIXTURE_DIR / "strain_a_wgmlst.tsv"),
            "--chew-profile",
            str(FIXTURE_DIR / "strain_b_wgmlst.tsv"),
            "--output-dir",
            str(output_dir),
        ],
        check=False,
        capture_output=True,
        text=True,
        cwd=REPO_ROOT,
    )

    assert completed.returncode == 0
    assert "Decision: different" in completed.stdout
    assert (output_dir / "wgmlst_profile_summary.json").exists()
