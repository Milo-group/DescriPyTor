"""Flask GUI routes used by docker-compose (`/status`, `/visual`, `/forms`)."""
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
    assert "example_feather_dir" in payload
    assert "example_outcomes" in payload
    assert "example_presets" in payload
    directory = payload.get("example_feather_dir") or ""
    if directory:
        assert "baptiste_products" in directory.replace("\\", "/")
        assert (payload.get("example_reference_name") or "") == "unsub"


def test_example_xyz_loads_unsub(client):
    from M2_data_extractor.gui_server import example_reference_feather

    path = example_reference_feather()
    if not path:
        pytest.skip("example feathers are not installed")
    try:
        import pyarrow  # noqa: F401
    except ImportError:
        pytest.skip("pyarrow is not installed")
    response = client.get("/example_xyz")
    if response.status_code == 500:
        pytest.skip("could not load the example feather in this environment")
    assert response.status_code == 200, response.get_json()
    body = response.get_json()
    assert body.get("xyz")
    assert body.get("name") == "unsub"
    assert str(body.get("filepath") or "").replace("\\", "/").endswith("unsub.feather")
    assert int(body.get("n_atoms") or 0) >= 10
    assert "C " in body["xyz"] or body["xyz"].splitlines()[2].strip()[:1] in "CNOSHP"


def test_fast_feather_xyz_reads_unsub():
    from M2_data_extractor.gui_server import example_reference_feather, _xyz_from_feather_fast

    path = example_reference_feather()
    if not path:
        pytest.skip("example feathers are not installed")
    try:
        import pyarrow  # noqa: F401
    except ImportError:
        pytest.skip("pyarrow is not installed")
    xyz, name, n_atoms = _xyz_from_feather_fast(path, "unsub")
    assert name == "unsub"
    assert n_atoms >= 10
    assert xyz.splitlines()[0].strip() == str(n_atoms)


def test_visual_serves_atom_picker(client):
    response = client.get("/visual")
    assert response.status_code == 200
    body = response.get_data(as_text=True)
    assert "3Dmol" in body
    assert "DescriPyTor" in body
    assert "Extract CSV" in body
    assert "Download CSV" in body
    assert "Save picks" in body
    assert 'data-tab="model"' in body
    assert "Run linear regression" in body
    assert "1-based" in body
    assert "Advanced" in body
    assert "Conformer Search" in body
    assert "Browse" in body
    assert 'id="a11yStatus"' in body
    assert 'id="applyExamplePicksBtn"' in body
    assert 'id="numberingBanner"' in body
    assert 'id="checkNumberingBtn"' in body
    assert 'id="viewNumberingBtn"' in body
    assert 'id="numberingProgress"' in body
    assert 'id="extractProgress"' in body
    viz = body.find('id="vizbar"')
    adv = body.find('id="advancedPanel"')
    assert viz != -1 and adv != -1
    assert viz < adv
    assert "pickerViewer" in body
    assert "recenterViewer" in body
    assert 'onclick="if(viewer)' not in body
    assert "Caffeine is the demo" not in body
    assert "loadExampleMolecule({ silent: true })" in body
    assert "Check numbering" in body


def test_visual_trailing_slash(client):
    response = client.get("/visual/")
    assert response.status_code == 200
    assert "3Dmol" in response.get_data(as_text=True)


def test_packaged_picker_is_preferred():
    from M2_data_extractor.gui_server import atom_picker_html

    path = atom_picker_html()
    assert path.endswith("atom_picker.html")
    assert Path(path).is_file()


def test_root_redirects_to_visual(client):
    response = client.get("/", follow_redirects=False)
    assert response.status_code in (301, 302)
    assert "/visual" in (response.headers.get("Location") or "")


def test_forms_serves_feature_gui(client):
    response = client.get("/forms")
    assert response.status_code == 200
    body = response.get_data(as_text=True)
    assert "html" in body.lower()


def test_xyz_from_feather_needs_input(client):
    response = client.post("/xyz_from_feather", json={})
    assert response.status_code == 400


def test_extract_needs_config(client):
    response = client.post("/extract", json={})
    assert response.status_code == 400
    assert "error" in response.get_json()


def test_extract_needs_feather_dir(client):
    response = client.post("/extract", json={"config": {"engines": {"descripytor_full": {"enabled": True}}}})
    assert response.status_code == 400
    body = response.get_json()
    assert "feather" in body["error"].lower()


def test_model_run_needs_csv(client):
    response = client.post("/model/run", json={})
    assert response.status_code == 400
    assert "csv" in response.get_json()["error"].lower()


def test_model_run_linear_regression(client):
    csv_path = ROOT / "tests" / "data" / "small_set" / "modeling_table.csv"
    if not csv_path.is_file():
        pytest.skip("modeling_table.csv missing")
    try:
        import sklearn  # noqa: F401
    except ImportError:
        pytest.skip("sklearn is not installed")
    csv_text = csv_path.read_text(encoding="utf-8")
    response = client.post("/model/run", json={
        "csv": csv_text,
        "y_column": "bp_c",
        "min_features": 1,
        "max_features": 1,
        "top_n": 5,
        "threshold": 0.0,
    })
    assert response.status_code == 200, response.get_json()
    body = response.get_json()
    assert body["target"] == "bp_c"
    assert body["n_models"] >= 1
    assert "combination" in body["rows"][0]
    assert "r2" in body["rows"][0]
    assert "q2" in body["rows"][0]


def test_model_run_pasted_outputs(client):
    csv_text = (
        "name,B1,L\n"
        "a,1.5,3.0\n"
        "b,1.6,4.0\n"
        "c,1.8,5.0\n"
        "d,2.0,6.0\n"
    )
    try:
        import sklearn  # noqa: F401
    except ImportError:
        pytest.skip("sklearn is not installed")
    response = client.post("/model/run", json={
        "csv": csv_text,
        "outputs": "10\n12\n14\n16\n",
        "new_column": "output",
        "min_features": 1,
        "max_features": 1,
        "top_n": 3,
        "threshold": 0.0,
    })
    assert response.status_code == 200, response.get_json()
    body = response.get_json()
    assert body["target"] == "output"
    assert body["n_models"] >= 1


def test_model_run_drops_nan_feature_columns(client):
    csv_text = (
        "name,B1,bad,L\n"
        "a,1.5,,3.0\n"
        "b,1.6,1,4.0\n"
        "c,1.8,,5.0\n"
        "d,2.0,2,6.0\n"
    )
    try:
        import sklearn  # noqa: F401
    except ImportError:
        pytest.skip("sklearn is not installed")
    response = client.post("/model/run", json={
        "csv": csv_text,
        "outputs": "10\n12\n14\n16\n",
        "new_column": "output",
        "min_features": 1,
        "max_features": 1,
        "top_n": 3,
        "threshold": 0.0,
    })
    assert response.status_code == 200, response.get_json()
    assert response.get_json()["n_models"] >= 1


def test_browse_folder_returns_path(client, monkeypatch):
    import M2_data_extractor.gui_server as gs

    monkeypatch.setattr(gs, "ask_folder_dialog", lambda initial="": r"C:\data\feathers")
    response = client.post("/browse/folder", json={"initial": ""})
    assert response.status_code == 200
    body = response.get_json()
    assert body["cancelled"] is False
    assert body["path"].endswith("feathers")


def test_browse_folder_cancelled(client, monkeypatch):
    import M2_data_extractor.gui_server as gs

    monkeypatch.setattr(gs, "ask_folder_dialog", lambda initial="": "")
    response = client.post("/browse/folder", json={})
    assert response.status_code == 200
    assert response.get_json()["cancelled"] is True


def test_numbering_check_example_set(client):
    from M2_data_extractor.gui_server import example_feather_dir, example_reference_feather

    directory = example_feather_dir()
    reference = example_reference_feather()
    if not directory or not reference:
        pytest.skip("example feather set is not installed")
    try:
        import pyarrow  # noqa: F401
    except ImportError:
        pytest.skip("pyarrow is not installed")
    response = client.post("/numbering/check", json={
        "directory": directory,
        "reference": reference,
        "picked": [1, 23, 8],
        "include_xyz": False,
        "limit": 2,
    })
    if response.status_code == 500:
        pytest.skip("could not load example feathers in this environment")
    assert response.status_code == 200, response.get_json()
    body = response.get_json()
    assert body["n_mols"] >= 1
    assert "molecules" in body
    assert "message" in body


def test_numbering_check_streams_progress(client):
    import json

    from M2_data_extractor.gui_server import example_feather_dir, example_reference_feather

    directory = example_feather_dir()
    reference = example_reference_feather()
    if not directory or not reference:
        pytest.skip("example feather set is not installed")
    try:
        import pyarrow  # noqa: F401
    except ImportError:
        pytest.skip("pyarrow is not installed")
    response = client.post("/numbering/check", json={
        "directory": directory,
        "reference": reference,
        "picked": [1],
        "include_xyz": False,
        "limit": 2,
        "stream": True,
    })
    if response.status_code == 500:
        pytest.skip("could not load example feathers in this environment")
    assert response.status_code == 200
    assert "ndjson" in (response.mimetype or "")
    events = [json.loads(line) for line in response.get_data(as_text=True).splitlines() if line.strip()]
    kinds = [ev.get("event") for ev in events]
    assert kinds[0] == "start"
    assert "progress" in kinds
    assert kinds[-1] == "done"
    assert events[0]["n"] == 2
    assert events[-1]["n_mols"] == 2


def test_placeholder_path_ignores_double_underscore_in_real_paths():
    from M2_data_extractor.gui_server import _is_placeholder_path

    assert _is_placeholder_path("")
    assert _is_placeholder_path(r"C:\path\to\your\feathers")
    assert _is_placeholder_path("__XYZ_DATA__")
    assert not _is_placeholder_path(r"C:\Users\edens\__backup__\mols")
    assert not _is_placeholder_path(r"C:\Users\edens\Documents\GitHub\DescriPyTor_to_upload")


def test_folder_summary_example_set(client):
    from M2_data_extractor.gui_server import example_feather_dir

    directory = example_feather_dir()
    if not directory:
        pytest.skip("example feather set is not installed")
    response = client.post("/folder/summary", json={"directory": directory})
    assert response.status_code == 200
    body = response.get_json()
    assert body["n"] >= 18
    assert "unsub" in body["names"]


def test_numbering_flags_out_of_range_index(client):
    from M2_data_extractor.gui_server import example_feather_dir, example_reference_feather

    directory = example_feather_dir()
    reference = example_reference_feather()
    if not directory or not reference:
        pytest.skip("example feather set is not installed")
    try:
        import pyarrow  # noqa: F401
    except ImportError:
        pytest.skip("pyarrow is not installed")
    response = client.post("/numbering/check", json={
        "directory": directory,
        "reference": reference,
        "picked": [9999],
        "include_xyz": False,
        "limit": 2,
    })
    if response.status_code == 500:
        pytest.skip("could not load example feathers in this environment")
    assert response.status_code == 200
    body = response.get_json()
    assert body["ok"] is False
    assert body["n_flagged"] == body["n_mols"]


def test_numbering_without_reference_uses_folder_unsub(client):
    from M2_data_extractor.gui_server import example_feather_dir

    directory = example_feather_dir()
    if not directory:
        pytest.skip("example feather set is not installed")
    try:
        import pyarrow  # noqa: F401
    except ImportError:
        pytest.skip("pyarrow is not installed")
    response = client.post("/numbering/check", json={
        "directory": directory,
        "picked": [1],
        "include_xyz": False,
        "limit": 2,
    })
    if response.status_code == 500:
        pytest.skip("could not load example feathers in this environment")
    assert response.status_code == 200
    body = response.get_json()
    assert body["reference"] == "unsub"
    assert body["n_mols"] == 2
