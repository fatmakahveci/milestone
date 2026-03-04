#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import re
import time
from urllib.request import Request
from urllib.error import HTTPError, URLError
from urllib.request import urlopen


DEFAULT_BASE_URL = "https://enterobase.warwick.ac.uk"
REQUEST_TIMEOUT_SECONDS = 30
REQUEST_RETRY_ATTEMPTS = 3
REQUEST_RETRY_BACKOFF_SECONDS = 1.5
DATABASE_NAME_PATTERN = re.compile(r"^[A-Za-z0-9_-]+$")
DEFAULT_TOKEN_ENV = "ENTEROBASE_TOKEN"
DEFAULT_USERNAME_ENV = "ENTEROBASE_USERNAME"
DEFAULT_PASSWORD_ENV = "ENTEROBASE_PASSWORD"


def build_headers(
    token_env: str = DEFAULT_TOKEN_ENV,
    username_env: str = DEFAULT_USERNAME_ENV,
    password_env: str = DEFAULT_PASSWORD_ENV,
) -> dict[str, str]:
    headers = {"User-Agent": "Milestone/1.0"}
    token = os.environ.get(token_env)
    username = os.environ.get(username_env)
    password = os.environ.get(password_env)
    if token:
        headers["Authorization"] = f"Bearer {token}"
    elif username and password:
        import base64

        credential = base64.b64encode(f"{username}:{password}".encode("utf-8")).decode("ascii")
        headers["Authorization"] = f"Basic {credential}"
    return headers


def fetch_url(url: str, headers: dict[str, str] | None = None) -> bytes:
    request = Request(url, headers=headers or {})
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
        raise RuntimeError(format_remote_error(url, last_error)) from last_error
    if isinstance(last_error, URLError):
        raise RuntimeError(format_remote_error(url, last_error)) from last_error
    raise RuntimeError(f"EnteroBase request failed for {url}.")


def fetch_json(url: str, headers: dict[str, str] | None = None) -> dict | list:
    return json.loads(fetch_url(url, headers=headers).decode("utf-8"))


def fetch_bytes(url: str, headers: dict[str, str] | None = None) -> bytes:
    return fetch_url(url, headers=headers)


def validate_database_name(database: str) -> str:
    normalized = database.strip()
    if not normalized:
        raise ValueError("Database name may not be empty.")
    if not DATABASE_NAME_PATTERN.match(normalized):
        raise ValueError(
            "Database name may only contain letters, numbers, underscores, and hyphens."
        )
    return normalized


def format_remote_error(url: str, exc: Exception) -> str:
    if isinstance(exc, HTTPError):
        if exc.code == 401:
            return (
                f"EnteroBase returned HTTP 401 for {url}. The API endpoint may require authentication "
                "or the current public endpoint may no longer allow anonymous access."
            )
        if exc.code == 404:
            return f"EnteroBase resource not found for {url}. The database or scheme name may be invalid."
        return f"EnteroBase returned HTTP {exc.code} for {url}."
    if isinstance(exc, URLError):
        return f"Unable to reach EnteroBase for {url}: {exc.reason}"
    return f"EnteroBase request failed for {url}: {exc}"


def list_schemes(
    database: str,
    base_url: str = DEFAULT_BASE_URL,
    headers: dict[str, str] | None = None,
) -> list[dict]:
    payload = fetch_json(f"{base_url}/api/v2.0/{database}/schemes?limit=1000", headers=headers)
    return payload.get("Schemes", [])


def select_scheme(
    database: str,
    scheme_name: str,
    base_url: str = DEFAULT_BASE_URL,
    headers: dict[str, str] | None = None,
) -> dict:
    schemes = list_schemes(database, base_url, headers=headers)
    for scheme in schemes:
        if scheme.get("scheme_name") == scheme_name:
            return scheme
    raise ValueError(f"Scheme '{scheme_name}' not found in EnteroBase database '{database}'.")


def download_scheme_profiles(
    database: str,
    scheme_name: str,
    output_dir: Path,
    base_url: str = DEFAULT_BASE_URL,
    headers: dict[str, str] | None = None,
) -> list[Path]:
    scheme = select_scheme(database, scheme_name, base_url, headers=headers)
    output_dir.mkdir(parents=True, exist_ok=True)
    profile_link = scheme.get("download_sts_link")
    if not profile_link:
        raise ValueError(f"Scheme '{scheme_name}' does not expose a download_sts_link.")

    archive_path = output_dir / f"{scheme_name}_profiles.tar.gz"
    archive_path.write_bytes(fetch_bytes(profile_link, headers=headers))

    metadata = {
        "database": database,
        "scheme_name": scheme_name,
        "label": scheme.get("label"),
        "download_sts_link": profile_link,
        "source": f"{base_url}/api/v2.0/{database}/schemes",
        "auth_mode": (
            "bearer"
            if headers and headers.get("Authorization", "").startswith("Bearer ")
            else "basic"
            if headers and headers.get("Authorization", "").startswith("Basic ")
            else "anonymous"
        ),
    }
    metadata_path = output_dir / "enterobase_scheme_metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    return [archive_path, metadata_path]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Download an EnteroBase scheme ST-profile dump and metadata."
    )
    parser.add_argument("--database", required=True, help="EnteroBase database, e.g. senterica.")
    parser.add_argument("--scheme-name", help="Scheme name, e.g. wgMLST.")
    parser.add_argument("--output-dir", type=Path, help="Directory for the downloaded scheme dump.")
    parser.add_argument("--base-url", default=DEFAULT_BASE_URL, help=f"Base URL. Default: {DEFAULT_BASE_URL}")
    parser.add_argument("--list-schemes", action="store_true", help="List available schemes and exit.")
    parser.add_argument("--token-env", default=DEFAULT_TOKEN_ENV, help=f"Env var name for bearer token. Default: {DEFAULT_TOKEN_ENV}")
    parser.add_argument("--username-env", default=DEFAULT_USERNAME_ENV, help=f"Env var name for basic-auth username. Default: {DEFAULT_USERNAME_ENV}")
    parser.add_argument("--password-env", default=DEFAULT_PASSWORD_ENV, help=f"Env var name for basic-auth password. Default: {DEFAULT_PASSWORD_ENV}")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.database = validate_database_name(args.database)
    headers = build_headers(args.token_env, args.username_env, args.password_env)
    if args.list_schemes:
        for scheme in list_schemes(args.database, args.base_url, headers=headers):
            print(f"{scheme.get('scheme_name')}\t{scheme.get('label')}\t{scheme.get('download_sts_link')}")
        return 0
    if not args.scheme_name or args.output_dir is None:
        raise SystemExit("--scheme-name and --output-dir are required unless --list-schemes is used.")
    downloaded = download_scheme_profiles(
        args.database,
        args.scheme_name,
        args.output_dir,
        args.base_url,
        headers=headers,
    )
    print(f"EnteroBase scheme dump written to {args.output_dir}")
    print(f"Files created: {len(downloaded)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
