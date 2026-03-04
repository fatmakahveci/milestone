from __future__ import annotations

import re
import os
from pathlib import Path
from uuid import uuid4


UPLOAD_ROOT = Path(__file__).resolve().parent.parent / "webapp_uploads"
MAX_UPLOAD_SIZE_MB = int(os.environ.get("MILESTONE_MAX_UPLOAD_MB", "250"))


def ensure_upload_root() -> Path:
    UPLOAD_ROOT.mkdir(parents=True, exist_ok=True)
    return UPLOAD_ROOT


def slugify(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "-", value.strip())
    cleaned = cleaned.strip(".-")
    if cleaned in {"", ".", ".."}:
        return "upload"
    return cleaned


def create_upload_dir(prefix: str) -> Path:
    ensure_upload_root()
    target = UPLOAD_ROOT / f"{slugify(prefix)}-{uuid4().hex[:8]}"
    target.mkdir(parents=True, exist_ok=True)
    return target


def save_uploaded_file(upload_dir: Path, uploaded_file) -> Path:
    payload = uploaded_file.getbuffer()
    size_bytes = len(payload)
    if size_bytes > MAX_UPLOAD_SIZE_MB * 1024 * 1024:
        raise ValueError(
            f"Upload '{uploaded_file.name}' exceeds the {MAX_UPLOAD_SIZE_MB} MB limit."
        )
    target = upload_dir / slugify(uploaded_file.name)
    target.write_bytes(payload)
    return target


def save_many(upload_dir: Path, uploaded_files: list) -> list[Path]:
    return [save_uploaded_file(upload_dir, uploaded_file) for uploaded_file in uploaded_files]
