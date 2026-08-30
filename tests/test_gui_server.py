"""Flask GUI routes used by docker-compose (`/status`, `/`, `/visual`)."""
from __future__ import annotations

import importlib.util
import sys
import types
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def _load_app():
    try:
        from flask import Flask  # noqa: F401
    except ImportError:
        pytest.skip("flask is not installed")
    if importlib.util.find_spec("flask_cors") is None:
        stub = types.ModuleType("flask_cors")
        stub.CORS = lambda app, *a, **k: app
        sys.modules["flask_cors"] = stub
    from M2_data_extractor import gui_server

    return gui_server.app


@pytest.fixture(scope="module")
def client():
    app = _load_app()
    app.config["TESTING"] = True
    return app.test_client()


def test_status(client):
    response = client.get("/status")
    assert response.status_code == 200
    payload = response.get_json()
    assert payload["ok"] is True


def test_visual_serves_atom_picker(client):
    response = client.get("/visual")
    assert response.status_code == 200
    body = response.get_data(as_text=True)
    assert "3Dmol" in body
    assert "DescriPyTor" in body


def test_visual_trailing_slash(client):
    response = client.get("/visual/")
    assert response.status_code == 200
    assert "3Dmol" in response.get_data(as_text=True)


def test_packaged_picker_is_preferred():
    from M2_data_extractor.gui_server import atom_picker_html

    path = atom_picker_html()
    assert path.endswith("atom_picker.html")
    assert Path(path).is_file()


def test_root_serves_feature_gui(client):
    response = client.get("/")
    assert response.status_code == 200
    body = response.get_data(as_text=True)
    assert "html" in body.lower()
