#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from workflow.scripts.build_species_validation_corpus import create_corpus_package
from workflow.scripts.freeze_public_benchmark_pack import freeze_pubmlst_benchmark_pack
from workflow.scripts.import_pubmlst_scheme import DEFAULT_BASE_URL, validate_base_url


def load_config(path: str | Path) -> dict[str, object]:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def build_workspace(config: dict[str, object], output_dir: str | Path) -> dict[str, object]:
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    base_url = validate_base_url(str(config.get("base_url", DEFAULT_BASE_URL)))
    species_configs = list(config.get("species", []))
    packs_dir = output_path / "packs"
    frozen_results = []

    for species_config in species_configs:
        species_slug = str(species_config["slug"])
        pack_output_dir = packs_dir / species_slug
        result = freeze_pubmlst_benchmark_pack(
            scheme_database=str(species_config["scheme_database"]),
            isolate_database=str(species_config["isolate_database"]),
            scheme_id=int(species_config["scheme_id"]),
            isolate_ids=[int(value) for value in species_config["isolate_ids"]],
            output_dir=pack_output_dir,
            base_url=base_url,
            species=species_config.get("species"),
            source_name=str(species_config.get("source_name", "Frozen PubMLST benchmark pack")),
            selection_criteria=species_config.get("selection_criteria"),
            download_contigs=bool(species_config.get("download_contigs", False)),
            zip_name=str(species_config.get("zip_name", "frozen_benchmark_pack.zip")),
        )
        frozen_results.append(
            {
                "slug": species_slug,
                "species": species_config.get("species", species_slug),
                **result,
            }
        )

    corpus = create_corpus_package(packs_dir, output_path / "corpus", "public_species_validation_corpus.zip")
    workspace_manifest = {
        "kind": "public_benchmark_workspace",
        "base_url": base_url,
        "species_count": len(species_configs),
        "frozen_packs": frozen_results,
        "corpus": corpus,
    }
    manifest_path = output_path / "public_benchmark_workspace.json"
    manifest_path.write_text(json.dumps(workspace_manifest, indent=2), encoding="utf-8")
    return workspace_manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a frozen public benchmark workspace from a JSON config of PubMLST species packs."
    )
    parser.add_argument("--config", required=True, help="JSON config describing one or more species benchmark packs.")
    parser.add_argument("--output-dir", required=True, help="Directory for the assembled workspace.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    manifest = build_workspace(load_config(args.config), args.output_dir)
    print(f"Species processed: {manifest['species_count']}")
    print(f"Workspace manifest written to {Path(args.output_dir) / 'public_benchmark_workspace.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
