from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
WORKFLOW_SCRIPT = REPO_ROOT / "workflow" / "milestone.py"
FIXTURE_DIR = Path(__file__).resolve().parent / "fixtures" / "benchmark" / "demo_species"
PACK_FIXTURE_DIR = Path(__file__).resolve().parent / "fixtures" / "benchmark" / "public_species_packs"


def test_profile_benchmark_cli_smoke(tmp_path: Path) -> None:
    output_dir = tmp_path / "benchmark"
    command = [
        sys.executable,
        str(WORKFLOW_SCRIPT),
        "profile_benchmark",
        "--predicted-dir",
        str(FIXTURE_DIR / "predicted"),
        "--truth-dir",
        str(FIXTURE_DIR / "truth"),
        "--output",
        str(output_dir),
    ]

    completed = subprocess.run(command, check=False, capture_output=True, text=True, cwd=REPO_ROOT)

    assert completed.returncode == 0, completed.stderr
    assert "Matched samples: 2" in completed.stdout
    summary = json.loads((output_dir / "wgmlst_benchmark_summary.json").read_text(encoding="utf-8"))
    assert summary["exact_allele_id_match_count"] == 1


def test_profile_benchmark_pack_cli_smoke(tmp_path: Path) -> None:
    output_dir = tmp_path / "benchmark_packs"
    command = [
        sys.executable,
        str(WORKFLOW_SCRIPT),
        "profile_benchmark",
        "--benchmark-pack-dir",
        str(PACK_FIXTURE_DIR),
        "--output",
        str(output_dir),
    ]

    completed = subprocess.run(command, check=False, capture_output=True, text=True, cwd=REPO_ROOT)

    assert completed.returncode == 0, completed.stderr
    assert "Species packs benchmarked: 3" in completed.stdout
    summary = json.loads((output_dir / "wgmlst_benchmark_pack_summary.json").read_text(encoding="utf-8"))
    assert summary["species_count"] == 3
