#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import sys
from datetime import date
from pathlib import Path
from urllib.request import urlopen

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from workflow.scripts.import_pubmlst_scheme import (
    DEFAULT_BASE_URL,
    REQUEST_TIMEOUT_SECONDS,
    validate_base_url,
)


def fetch_json(url: str) -> dict | list:
    with urlopen(url, timeout=REQUEST_TIMEOUT_SECONDS) as response:
        return json.load(response)


def fetch_text(url: str) -> str:
    with urlopen(url, timeout=REQUEST_TIMEOUT_SECONDS) as response:
        return response.read().decode("utf-8")


def sanitize_sample_name(raw_name: str) -> str:
    cleaned = "".join(char if char.isalnum() or char in {"_", "-"} else "_" for char in raw_name.strip())
    cleaned = cleaned.strip("._-")
    return cleaned or "sample"


def fetch_scheme_info(database: str, scheme_id: int, base_url: str) -> dict:
    return fetch_json(f"{base_url}/db/{database}/schemes/{scheme_id}")


def fetch_isolate_record(database: str, isolate_id: int, base_url: str) -> dict:
    return fetch_json(f"{base_url}/db/{database}/isolates/{isolate_id}")


def fetch_isolate_scheme_alleles(database: str, isolate_id: int, scheme_id: int, base_url: str) -> dict:
    return fetch_json(f"{base_url}/db/{database}/isolates/{isolate_id}/schemes/{scheme_id}/allele_ids")


def fetch_isolate_ids(database: str, base_url: str, limit: int) -> list[int]:
    payload = fetch_json(f"{base_url}/db/{database}/isolates?page_size={limit}")
    isolate_uris = payload.get("isolates", [])
    isolate_ids = []
    for uri in isolate_uris[:limit]:
        isolate_ids.append(int(str(uri).rstrip("/").split("/")[-1]))
    return isolate_ids


def build_truth_profile_rows(payload: dict) -> list[tuple[str, str]]:
    rows = []
    for locus, allele in sorted(payload.items()):
        if locus in {"isolate", "scheme", "sender", "curator"}:
            continue
        rows.append((locus, str(allele)))
    return rows


def write_truth_profile(output_dir: Path, sample_name: str, rows: list[tuple[str, str]]) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / f"{sample_name}_wgmlst.tsv"
    with path.open("w", encoding="utf-8") as handle:
        for locus, allele in rows:
            handle.write(f"{locus}\t{allele}\n")
    return path


def write_placeholder_predicted_dir(output_dir: Path) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    placeholder = output_dir / "README.md"
    placeholder.write_text(
        "Run Milestone on matching assemblies or reads and place predicted *_wgmlst.tsv files here.\n",
        encoding="utf-8",
    )
    return placeholder


def write_isolate_metadata(output_dir: Path, rows: list[dict[str, object]]) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / "isolate_metadata.json"
    path.write_text(json.dumps(rows, indent=2), encoding="utf-8")
    return path


def download_contigs_fasta(database: str, isolate_id: int, output_dir: Path, sample_name: str, base_url: str) -> Path:
    fasta_text = fetch_text(f"{base_url}/db/{database}/isolates/{isolate_id}/contigs_fasta")
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / f"{sample_name}.fasta"
    path.write_text(fasta_text, encoding="utf-8")
    return path


def import_benchmark_pack(
    scheme_database: str,
    isolate_database: str,
    scheme_id: int,
    isolate_ids: list[int],
    output_dir: Path,
    base_url: str = DEFAULT_BASE_URL,
    species: str | None = None,
    source_name: str = "PubMLST/BIGSdb isolate benchmark pack",
    download_contigs: bool = False,
    selection_criteria: str | None = None,
) -> list[Path]:
    scheme = fetch_scheme_info(scheme_database, scheme_id, base_url)
    truth_dir = output_dir / "truth"
    predicted_dir = output_dir / "predicted"
    assemblies_dir = output_dir / "assemblies"
    downloaded: list[Path] = []
    isolate_metadata_rows: list[dict[str, object]] = []

    for isolate_id in isolate_ids:
        isolate_record = fetch_isolate_record(isolate_database, isolate_id, base_url)
        allele_payload = fetch_isolate_scheme_alleles(isolate_database, isolate_id, scheme_id, base_url)
        sample_name = sanitize_sample_name(str(isolate_record.get("id") or f"isolate_{isolate_id}"))
        truth_rows = build_truth_profile_rows(allele_payload)
        downloaded.append(write_truth_profile(truth_dir, sample_name, truth_rows))
        isolate_metadata_rows.append(
            {
                "isolate_id": isolate_id,
                "sample_name": sample_name,
                "record_url": f"{base_url}/db/{isolate_database}/isolates/{isolate_id}",
            }
        )
        if download_contigs:
            downloaded.append(download_contigs_fasta(isolate_database, isolate_id, assemblies_dir, sample_name, base_url))

    downloaded.append(write_placeholder_predicted_dir(predicted_dir))
    downloaded.append(write_isolate_metadata(output_dir, isolate_metadata_rows))

    manifest = {
        "species": species or isolate_database,
        "schema_name": scheme.get("description") or f"{scheme_database}_scheme_{scheme_id}",
        "source_name": source_name,
        "source_url": f"{base_url}/db/{scheme_database}/schemes/{scheme_id}",
        "retrieval_date": str(date.today()),
        "source_snapshot": f"{scheme_database}:{scheme_id}",
        "selection_criteria": selection_criteria or "explicit_isolate_ids",
        "scheme_database": scheme_database,
        "isolate_database": isolate_database,
        "scheme_id": scheme_id,
        "isolate_ids": isolate_ids,
        "public_profile_count": len(isolate_ids),
        "truth_profile_count": len(isolate_ids),
    }
    manifest_path = output_dir / "benchmark_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    downloaded.append(manifest_path)
    return downloaded


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create a benchmark-pack layout from PubMLST/BIGSdb isolates and scheme allele IDs."
    )
    parser.add_argument("--database", help="Legacy shortcut: use the same database for scheme and isolate requests.")
    parser.add_argument("--scheme-database", help="PubMLST/BIGSdb sequence-definition database name.")
    parser.add_argument("--isolate-database", help="PubMLST/BIGSdb isolate database name.")
    parser.add_argument("--scheme-id", type=int, required=True, help="Scheme identifier.")
    parser.add_argument(
        "--isolate-id",
        dest="isolate_ids",
        type=int,
        action="append",
        required=False,
        help="Isolate identifier to include. Repeat for multiple isolates.",
    )
    parser.add_argument(
        "--discover-isolates",
        type=int,
        default=0,
        help="Automatically fetch the first N isolate IDs from the isolate database when explicit isolate IDs are not provided.",
    )
    parser.add_argument("--output-dir", type=Path, required=True, help="Directory for the benchmark pack.")
    parser.add_argument("--base-url", default=DEFAULT_BASE_URL, help=f"REST API base URL. Default: {DEFAULT_BASE_URL}")
    parser.add_argument("--species", help="Species label written into the benchmark manifest.")
    parser.add_argument("--source-name", default="PubMLST/BIGSdb isolate benchmark pack", help="Human-readable source label.")
    parser.add_argument("--selection-criteria", default=None, help="Free-text isolate selection criteria for the benchmark manifest.")
    parser.add_argument("--download-contigs", action="store_true", help="Also download contig FASTA for each isolate.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.base_url = validate_base_url(args.base_url)
    scheme_database = args.scheme_database or args.database
    isolate_database = args.isolate_database or args.database
    if not scheme_database or not isolate_database:
        raise SystemExit("--scheme-database and --isolate-database are required unless --database is used.")
    isolate_ids = list(args.isolate_ids or [])
    if not isolate_ids and args.discover_isolates > 0:
        isolate_ids = fetch_isolate_ids(isolate_database, args.base_url, args.discover_isolates)
    if not isolate_ids:
        raise SystemExit("Provide at least one --isolate-id or use --discover-isolates.")
    downloaded = import_benchmark_pack(
        scheme_database=scheme_database,
        isolate_database=isolate_database,
        scheme_id=args.scheme_id,
        isolate_ids=isolate_ids,
        output_dir=args.output_dir,
        base_url=args.base_url,
        species=args.species,
        source_name=args.source_name,
        download_contigs=args.download_contigs,
        selection_criteria=args.selection_criteria,
    )
    print(f"Benchmark pack written to {args.output_dir}")
    print(f"Files created: {len(downloaded)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
