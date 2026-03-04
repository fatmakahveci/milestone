#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import time
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import quote, urlparse
from urllib.request import Request, urlopen

DEFAULT_BASE_URL = "https://rest.pubmlst.org"
SUPPORTED_SCHEME_TYPES = {"wgmlst", "cgmlst"}
REQUEST_TIMEOUT_SECONDS = 30
REQUEST_RETRY_ATTEMPTS = 3
REQUEST_RETRY_BACKOFF_SECONDS = 1.0
REPO_ROOT = Path(__file__).resolve().parents[2]
DATABASE_NAME_PATTERN = re.compile(r"^[A-Za-z0-9_-]+$")

if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from workflow.scripts.schema_manifest import (
    build_schema_manifest,
    validate_schema_manifest,
    write_schema_manifest,
)


def fetch_json(url: str) -> dict | list:
    return json.loads(fetch_url(url).decode("utf-8"))


def fetch_text(url: str) -> str:
    return fetch_url(url).decode("utf-8")


def fetch_url(url: str) -> bytes:
    request = Request(url, headers={"User-Agent": "Milestone/1.0"})
    last_error: Exception | None = None
    for attempt in range(1, REQUEST_RETRY_ATTEMPTS + 1):
        try:
            with urlopen(request, timeout=REQUEST_TIMEOUT_SECONDS) as response:
                return response.read()
        except HTTPError as exc:
            last_error = exc
            if exc.code in {401, 403, 404}:
                break
        except URLError as exc:
            last_error = exc
        if attempt < REQUEST_RETRY_ATTEMPTS:
            time.sleep(REQUEST_RETRY_BACKOFF_SECONDS * attempt)
    if isinstance(last_error, HTTPError):
        raise RuntimeError(_format_remote_error(url, last_error)) from last_error
    if isinstance(last_error, URLError):
        raise RuntimeError(_format_remote_error(url, last_error)) from last_error
    raise RuntimeError(f"Remote request failed for {url}.")


def validate_base_url(base_url: str) -> str:
    parsed = urlparse(base_url)
    if parsed.scheme not in {"http", "https"}:
        raise ValueError(f"Unsupported base URL scheme: {parsed.scheme or 'none'}")
    if not parsed.netloc:
        raise ValueError("Base URL must include a network location.")
    return base_url.rstrip("/")


def validate_database_name(database: str) -> str:
    normalized = database.strip()
    if not normalized:
        raise ValueError("Database name may not be empty.")
    if not DATABASE_NAME_PATTERN.match(normalized):
        raise ValueError(
            "Database name may only contain letters, numbers, underscores, and hyphens."
        )
    return normalized


def _format_remote_error(url: str, exc: Exception) -> str:
    if isinstance(exc, HTTPError):
        if exc.code == 401:
            return f"Remote service rejected the request with HTTP 401 for {url}. Authentication may be required."
        if exc.code == 404:
            return f"Remote resource not found for {url}. The database or scheme identifier may be invalid."
        return f"Remote service returned HTTP {exc.code} for {url}."
    if isinstance(exc, URLError):
        return f"Unable to reach remote service for {url}: {exc.reason}"
    return f"Remote request failed for {url}: {exc}"


def list_schemes(database: str, base_url: str = DEFAULT_BASE_URL) -> list[dict]:
    return fetch_json(f"{base_url}/db/{database}/schemes")


def infer_scheme_type(description: str | None) -> str:
    if not description:
        return "other"
    lowered = description.lower()
    if "wgmlst" in lowered or "whole genome" in lowered:
        return "wgmlst"
    if "cgmlst" in lowered or "core genome" in lowered:
        return "cgmlst"
    return "other"


def filter_schemes_by_type(schemes: list[dict], scheme_type: str | None) -> list[dict]:
    if not scheme_type:
        return schemes
    normalized = scheme_type.lower()
    if normalized not in SUPPORTED_SCHEME_TYPES:
        raise ValueError(f"Unsupported scheme type filter: {scheme_type}")
    return [
        item
        for item in schemes
        if infer_scheme_type(item.get("description")) == normalized
    ]


def get_scheme(database: str, scheme_id: int, base_url: str = DEFAULT_BASE_URL) -> dict:
    return fetch_json(f"{base_url}/db/{database}/schemes/{scheme_id}")


def get_scheme_loci(database: str, scheme_id: int, base_url: str = DEFAULT_BASE_URL) -> list[str]:
    payload = fetch_json(f"{base_url}/db/{database}/schemes/{scheme_id}/loci")
    return payload["loci"]


def locus_name_from_uri(locus_uri: str) -> str:
    return locus_uri.rstrip("/").split("/")[-1]


def download_locus_alleles(database: str, locus: str, output_dir: Path, base_url: str = DEFAULT_BASE_URL) -> Path:
    locus_name = locus_name_from_uri(locus)
    safe_locus_name = "".join(character if character.isalnum() or character in {"_", "-", "."} else "_" for character in locus_name)
    safe_locus_name = safe_locus_name.strip("._-") or "locus"
    fasta_text = fetch_text(f"{base_url}/db/{database}/loci/{quote(locus_name, safe='')}/alleles_fasta")
    output_path = output_dir / f"{safe_locus_name}.fasta"
    output_path.write_text(fasta_text, encoding="utf-8")
    return output_path


def download_scheme(
    database: str,
    scheme_id: int,
    output_dir: Path,
    base_url: str = DEFAULT_BASE_URL,
    include_profiles: bool = False,
    expected_scheme_type: str | None = None,
    schema_name: str | None = None,
    species: str | None = None,
    schema_version: str | None = None,
) -> list[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    scheme = get_scheme(database, scheme_id, base_url=base_url)
    loci = get_scheme_loci(database, scheme_id, base_url=base_url)
    inferred_scheme_type = infer_scheme_type(scheme.get("description"))

    if expected_scheme_type and inferred_scheme_type != expected_scheme_type:
        raise ValueError(
            f"Scheme {scheme_id} in database '{database}' is inferred as '{inferred_scheme_type}', "
            f"not '{expected_scheme_type}'."
        )

    downloaded = [
        download_locus_alleles(database, locus_uri, output_dir, base_url=base_url)
        for locus_uri in loci
    ]

    metadata = {
        "database": database,
        "scheme_id": scheme_id,
        "description": scheme.get("description"),
        "scheme_type": inferred_scheme_type,
        "locus_count": scheme.get("locus_count"),
        "source": f"{base_url}/db/{database}/schemes/{scheme_id}",
    }
    (output_dir / "scheme_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    manifest = build_schema_manifest(
        output_dir,
        schema_name=schema_name or f"{database}_scheme_{scheme_id}",
        schema_type=inferred_scheme_type,
        species=species or "unspecified",
        source="PubMLST/BIGSdb",
        source_url=f"{base_url}/db/{database}/schemes/{scheme_id}",
        schema_version=schema_version or str(scheme_id),
        notes=f"Imported from database '{database}' via Milestone importer.",
    )
    manifest_path = write_schema_manifest(manifest, output_dir / "scheme_manifest.json")
    validate_schema_manifest(manifest_path, output_dir)

    if include_profiles and scheme.get("profiles_csv"):
        profiles_text = fetch_text(scheme["profiles_csv"])
        profiles_path = output_dir / "profiles.tsv"
        profiles_path.write_text(profiles_text, encoding="utf-8")
        downloaded.append(profiles_path)

    return downloaded


def build_schema_creation_command(
    schema_dir: Path,
    reference: str,
    output_dir: Path,
    threads: int,
) -> list[str]:
    workflow_script = Path(__file__).resolve().parent.parent / "milestone.py"
    return [
        sys.executable,
        str(workflow_script),
        "schema_creation",
        "--reference",
        reference,
        "--schema_name",
        str(schema_dir),
        "--output",
        str(output_dir),
        "--threads",
        str(threads),
    ]


def run_schema_creation(schema_dir: Path, reference: str, output_dir: Path, threads: int) -> int:
    command = build_schema_creation_command(schema_dir, reference, output_dir, threads)
    completed = subprocess.run(command, check=False)
    return completed.returncode


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Download public PubMLST/BIGSdb allele schemes into Milestone schema layout."
    )
    parser.add_argument("--database", required=True, help="PubMLST sequence definition database, e.g. pubmlst_neisseria_seqdef")
    parser.add_argument("--scheme-id", type=int, help="Scheme identifier to download.")
    parser.add_argument("--output-dir", type=Path, help="Directory to write locus FASTA files.")
    parser.add_argument("--base-url", default=DEFAULT_BASE_URL, help=f"REST API base URL. Default: {DEFAULT_BASE_URL}")
    parser.add_argument("--include-profiles", action="store_true", help="Also download scheme profiles TSV when available.")
    parser.add_argument("--list-schemes", action="store_true", help="List available schemes for the database and exit.")
    parser.add_argument(
        "--scheme-type",
        choices=sorted(SUPPORTED_SCHEME_TYPES),
        help="Filter listed schemes by high-level type inferred from the description.",
    )
    parser.add_argument("--reference", help="Reference prefix to use if schema_creation is triggered after download.")
    parser.add_argument("--pipeline-output-dir", type=Path, help="Output directory for chained schema_creation run.")
    parser.add_argument("--threads", type=int, default=1, help="Thread count for chained schema_creation run.")
    parser.add_argument("--schema-name", help="Human-readable schema name stored in scheme_manifest.json.")
    parser.add_argument("--species", help="Species label stored in scheme_manifest.json.")
    parser.add_argument("--schema-version", help="Schema version string stored in scheme_manifest.json.")
    parser.add_argument(
        "--run-schema-creation",
        action="store_true",
        help="Run Milestone schema_creation immediately after downloading the public scheme.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_arguments()
    args.base_url = validate_base_url(args.base_url)
    args.database = validate_database_name(args.database)

    if args.list_schemes:
        schemes = filter_schemes_by_type(
            list_schemes(args.database, base_url=args.base_url),
            args.scheme_type,
        )
        for item in schemes:
            scheme_uri = item["scheme"]
            scheme_id = scheme_uri.rstrip("/").split("/")[-1]
            print(f"{scheme_id}\t{infer_scheme_type(item.get('description'))}\t{item['description']}")
        return 0

    if args.scheme_id is None or args.output_dir is None:
        raise SystemExit("--scheme-id and --output-dir are required unless --list-schemes is used.")

    downloaded = download_scheme(
        args.database,
        args.scheme_id,
        args.output_dir,
        base_url=args.base_url,
        include_profiles=args.include_profiles,
        expected_scheme_type=args.scheme_type,
        schema_name=args.schema_name,
        species=args.species,
        schema_version=args.schema_version,
    )
    print(f"Downloaded {len(downloaded)} files to {args.output_dir}")

    if args.run_schema_creation:
        if args.reference is None or args.pipeline_output_dir is None:
            raise SystemExit("--reference and --pipeline-output-dir are required with --run-schema-creation.")
        returncode = run_schema_creation(args.output_dir, args.reference, args.pipeline_output_dir, args.threads)
        if returncode != 0:
            raise SystemExit(returncode)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
