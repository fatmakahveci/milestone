from __future__ import annotations

import sys
from pathlib import Path

import pytest

pytest.importorskip("pandas")
pytest.importorskip("streamlit")
pytest.importorskip("streamlit.testing.v1")

from streamlit.testing.v1 import AppTest

APP_PATH = Path(__file__).resolve().parent.parent / "webapp" / "app.py"


def _build_app_test() -> AppTest:
    webapp_dir = APP_PATH.parent
    webapp_dir_str = str(webapp_dir)
    if webapp_dir_str not in sys.path:
        sys.path.insert(0, webapp_dir_str)
    return AppTest.from_file(str(APP_PATH))


def test_app_renders_without_exception() -> None:
    app = _build_app_test().run()

    assert not app.exception
    assert len(app.text_input) >= 9
    assert len(app.button) >= 4


def test_schema_and_allele_forms_are_present() -> None:
    app = _build_app_test().run()

    assert not app.exception
    labels = {button.label for button in app.button}
    assert "Build demo command" in labels
    assert "Preview demo run" in labels
    assert "Preview comparison" in labels
    assert "Preview benchmark" in labels
    assert "Preview matrix" in labels
    assert "Open novel allele exports" in labels
    assert "Preview QC" in labels
    assert "Start PubMLST benchmark import" in labels
    assert "Start EnteroBase import" in labels


def test_public_schema_tab_renders() -> None:
    app = _build_app_test().run()
    assert not app.exception
    assert len(app.text_input) >= 10
