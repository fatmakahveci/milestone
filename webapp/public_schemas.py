from __future__ import annotations

import importlib.util
import sys
from functools import lru_cache
from pathlib import Path

MODULE_PATH = Path(__file__).resolve().parent.parent / "workflow" / "scripts" / "import_pubmlst_scheme.py"
SPEC = importlib.util.spec_from_file_location("import_pubmlst_scheme", MODULE_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"Unable to load import_pubmlst_scheme module from {MODULE_PATH}")
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules.setdefault(SPEC.name, MODULE)
SPEC.loader.exec_module(MODULE)

filter_schemes_by_type = MODULE.filter_schemes_by_type
_list_schemes = MODULE.list_schemes


@lru_cache(maxsize=32)
def list_schemes(database: str, base_url: str = MODULE.DEFAULT_BASE_URL) -> list[dict]:
    return _list_schemes(database, base_url=base_url)


def clear_schema_cache() -> None:
    list_schemes.cache_clear()
