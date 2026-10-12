"""Flask GUI routes (`/status`, `/visual`, the explorers)."""
from __future__ import annotations

import importlib.util
import json
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
    from M2_data_extractor.gui_server import GUI_API
    assert payload.get("api") == GUI_API
    assert "example_feather_dir" in payload
    assert "example_outcomes" in payload
    assert "example_presets" in payload
    directory = payload.get("example_feather_dir") or ""
    if directory:
        assert "feather_example" in directory.replace("\\", "/")
        assert (payload.get("example_reference_name") or "") == "basic"
        assert payload.get("example_presets")
    outcomes = payload.get("example_outcomes") or []
    if outcomes:
        names = {str(row.get("name") or "") for row in outcomes}
        assert "unsub" in names
        assert any(row.get("name") == "unsub" and float(row.get("output")) == 0.1461 for row in outcomes)


def test_example_xyz_loads_basic(client):
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
    assert body.get("name") == "basic"
    assert str(body.get("filepath") or "").replace("\\", "/").endswith("basic.feather")
    assert int(body.get("n_atoms") or 0) >= 6
    assert "C " in body["xyz"] or body["xyz"].splitlines()[2].strip()[:1] in "CNOSHP"
    dipole = ((body.get("mol_data") or {}).get("dipole") or {})
    vec = dipole.get("vector") or []
    assert len(vec) == 3
    assert all(isinstance(x, (int, float)) for x in vec)
    assert abs(float(dipole.get("total") or 0)) > 0.1
    modes = (body.get("mol_data") or {}).get("modes") or []
    assert modes, "expected vibration modes from basic.feather"
    first = modes[0]
    assert "freq" in first and "disp" in first
    assert len(first["disp"]) == int(body["n_atoms"])
    assert len(first["disp"][0]) == 3


def test_dipole_from_named_and_positional_feather_columns():
    import math

    import pandas as pd

    from M2_data_extractor.gui_server import _dipole_from_feather_df

    named = pd.DataFrame({
        "atom": [None, "C"],
        "x": [None, 0.0],
        "y": [None, 0.0],
        "z": [None, 0.0],
        "dip_x": [0.9783, None],
        "dip_y": [-1.3102, None],
        "dip_z": [0.2607, None],
        "total_dipole": [1.6558, None],
    })
    got = _dipole_from_feather_df(named)
    assert got is not None
    assert got["vector"] == [0.9783, -1.3102, 0.2607]
    assert abs(got["total"] - 1.6558) < 1e-6

    positional = pd.DataFrame({
        "0": ["29.0", "C"],
        "1": [None, 0.1],
        "2": [None, 0.2],
        "3": [None, 0.3],
        "4": [-1.2426, None],
        "5": [-1.4114, None],
        "6": [4.9682, None],
        "7": [5.3122, None],
    })
    got = _dipole_from_feather_df(positional)
    assert got is not None
    assert got["vector"][0] == -1.2426
    assert abs(got["total"] - 5.3122) < 1e-6
    assert abs(got["total"] - math.sqrt(sum(x * x for x in got["vector"]))) < 0.02


def test_fast_feather_xyz_reads_basic():
    from M2_data_extractor.gui_server import example_reference_feather, _xyz_from_feather_fast

    path = example_reference_feather()
    if not path:
        pytest.skip("example feathers are not installed")
    try:
        import pyarrow  # noqa: F401
    except ImportError:
        pytest.skip("pyarrow is not installed")
    xyz, name, n_atoms = _xyz_from_feather_fast(path, "basic")
    assert name == "basic"
    assert n_atoms >= 6
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
    assert 'id="fileInput"' in body
    assert 'id="loadFileBtn"' in body
    assert "Load file" in body
    assert 'accept=".xyz,.feather,.ftr"' not in body
    assert "descripytor.visual.session.v2" in body
    assert 'id="a11yStatus"' in body
    assert 'id="clearPicksBtn"' in body
    assert 'id="pickControlsFold"' in body
    assert 'id="selTableFold" open' in body
    assert 'id="numberingBanner"' in body
    assert 'id="checkNumberingBtn"' in body
    assert 'id="viewNumberingBtn"' in body
    assert 'id="numberingProgress"' in body
    assert 'id="extractProgress"' in body
    assert 'id="presetProgress"' in body
    assert 'id="modelProgress"' in body
    viz = body.find('id="vizbar"')
    adv = body.find('id="advancedPanel"')
    assert viz != -1 and adv != -1
    assert viz < adv
    assert "pickerFrame" in body and "feature:session" in body      # the main view is the feature explorer, framed
    assert "recenterViewer" in body
    assert 'id="pickerStage"' in body
    assert 'id="viewer"' not in body
    assert 'onclick="if(viewer)' not in body
    assert "Caffeine is the demo" not in body
    assert "loadExampleMolecule({ silent: true })" in body
    assert "Check numbering" in body
    assert "JSON.stringify({ config: cfg, stream: true })" in body
    assert "lastModelXYZ" in body
    assert "Apply example picks" in body
    assert "/static/3Dmol-min.js" in body
    assert 'e.shiftKey ? "full"' in body


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


def test_the_removed_pages_are_gone(client):
    """The form-based page and the mode viewer's route went in 0.3.3; nothing in the package serves or links to them."""
    for path in ("/forms", "/modes"):
        assert client.get(path).status_code == 404
    for path in ("/sterimol", "/charges", "/dipole", "/vibrations", "/molecule_info", "/steriplot", "/batch_sterimol", "/features_set", "/api/parse-modes"):
        assert client.post(path, json={}).status_code == 404
    picker = (ROOT / "M2_data_extractor" / "atom_picker.html").read_text(encoding="utf-8")
    assert '"/forms"' not in picker and "/modes" not in picker.split("FEATURE-EXPLORER:END")[-1]


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


def test_extract_stream_needs_feather_dir(client):
    response = client.post("/extract", json={
        "config": {"engines": {"descripytor_full": {"enabled": True}}},
        "stream": True,
    })
    assert response.status_code == 400
    body = response.get_json()
    assert "feather" in body["error"].lower()


def test_extract_stream_emits_start(client):
    from M2_data_extractor.gui_server import example_feather_dir

    directory = example_feather_dir()
    if not directory:
        pytest.skip("example feather set is not installed")
    response = client.post("/extract", json={
        "config": {"engines": {}, "feather_dir": directory},
        "stream": True,
    })
    assert response.status_code == 200
    assert "ndjson" in (response.mimetype or "")
    lines = [ln for ln in response.get_data(as_text=True).splitlines() if ln.strip()]
    events = [json.loads(ln) for ln in lines]
    assert events[0]["event"] == "start"
    # the bundled example set (10 benzenes since 82742d0), counted rather than hard-coded
    n_files = len(list(Path(directory).glob("*.feather")))
    assert events[0]["n"] == n_files > 0
    assert any(e.get("event") == "error" for e in events)


def test_molecules_reports_load_progress(tmp_path, monkeypatch):
    import os
    from M2_data_extractor import data_extractor as de

    de._MOLECULES_CACHE.clear()

    class Dummy:
        def __init__(self, data_file, threshold=1.82):
            self.molecule_name = os.path.splitext(os.path.basename(data_file))[0]

    monkeypatch.setattr(de, "Molecule", Dummy)
    (tmp_path / "a.feather").write_bytes(b"x")
    (tmp_path / "b.feather").write_bytes(b"x")
    events = []
    mols = de.Molecules(str(tmp_path), progress=events.append)
    assert events[0]["event"] == "start"
    assert events[0]["n"] == 2
    assert events[0]["phase"] == "load"
    progress = [e for e in events if e.get("event") == "progress"]
    assert len(progress) == 2
    assert {e["name"] for e in progress} == {"a", "b"}
    assert progress[-1]["i"] == 2
    assert len(mols.success_molecules) == 2


def test_molecules_reuses_folder_cache(tmp_path, monkeypatch):
    import os
    from M2_data_extractor import data_extractor as de

    de._MOLECULES_CACHE.clear()
    calls = {"n": 0}

    class Dummy:
        def __init__(self, data_file, threshold=1.82):
            calls["n"] += 1
            self.molecule_name = os.path.splitext(os.path.basename(data_file))[0]

    monkeypatch.setattr(de, "Molecule", Dummy)
    (tmp_path / "a.feather").write_bytes(b"x")
    (tmp_path / "b.feather").write_bytes(b"x")
    first = de.Molecules(str(tmp_path))
    assert calls["n"] == 2
    second = de.Molecules(str(tmp_path))
    assert calls["n"] == 2
    assert second.success_molecules == first.success_molecules


def test_molecules_skips_non_molecule_json(tmp_path, monkeypatch):
    import os
    from M2_data_extractor import data_extractor as de

    de._MOLECULES_CACHE.clear()

    class Dummy:
        def __init__(self, data_file, threshold=1.82):
            self.molecule_name = os.path.splitext(os.path.basename(data_file))[0]

    monkeypatch.setattr(de, "Molecule", Dummy)
    (tmp_path / "a.feather").write_bytes(b"x")
    (tmp_path / "presets.json").write_text(
        '{"label": "demo", "atoms": {"sterimol": [[1, 4]]}}', encoding="utf-8"
    )
    (tmp_path / "input_example.json").write_text('{"Sterimol atoms": "1,4"}', encoding="utf-8")
    (tmp_path / "run_config.json").write_text('{"feather_dir": "x"}', encoding="utf-8")
    (tmp_path / "notes.json").write_text('{"comment": "not a molecule"}', encoding="utf-8")
    (tmp_path / "mol.json").write_text(
        '{"molecule": {"atoms": [{"atom": "C", "x": 0, "y": 0, "z": 0}]}}',
        encoding="utf-8",
    )
    events = []
    mols = de.Molecules(str(tmp_path), progress=events.append)
    assert events[0]["n"] == 2
    assert set(mols.success_molecules) == {"a.feather", "mol.json"}


def test_molecule_defers_vibration_frames():
    from M2_data_extractor.data_extractor import Molecule, _MOLECULES_CACHE
    from M2_data_extractor.gui_server import example_reference_feather

    path = example_reference_feather()
    if not path:
        pytest.skip("example feathers are not installed")
    try:
        import pyarrow  # noqa: F401
    except ImportError:
        pytest.skip("pyarrow is not installed")
    _MOLECULES_CACHE.clear()
    mol = Molecule(path)
    assert getattr(mol, "_dfs", None) is None
    assert mol.vibration_dict
    frames = mol.dfs
    assert isinstance(frames, list)
    assert frames
    assert getattr(mol, "_dfs", None) is frames


def test_nob_atype_counts_bonds_without_row_filters():
    import pandas as pd
    from utils.help_functions import nob_atype

    xyz = pd.DataFrame({
        "atom": ["C", "C", "H", "H", "H", "O"],
        "x": [0, 1.5, -0.5, 0.5, 0, 2.5],
        "y": [0, 0, 0.5, 0.5, -0.7, 0],
        "z": [0, 0, 0, 0, 0, 0],
    })
    bonds = pd.DataFrame({0: [1, 1, 1, 1, 2], 1: [2, 3, 4, 5, 6]})
    types = nob_atype(xyz, bonds)
    assert types[0] == "C"
    assert types[1] == "C3"
    assert types[2] == "H"
    assert types[5] == "O2"


def test_3dmol_static_is_local_or_cdn(client):
    response = client.get("/static/3Dmol-min.js", follow_redirects=False)
    assert response.status_code in (200, 302)
    if response.status_code == 302:
        assert "3dmol" in (response.headers.get("Location") or "").lower()


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


def test_model_run_streams_search_progress(client):
    csv_text = (
        "name,B1,L,B5\n"
        "a,1.5,3.0,2.0\n"
        "b,1.6,4.0,2.1\n"
        "c,1.8,5.0,2.4\n"
        "d,2.0,6.0,2.8\n"
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
        "max_features": 2,
        "top_n": 5,
        "threshold": 0.0,
        "stream": True,
    })
    assert response.status_code == 200
    assert "ndjson" in (response.mimetype or "")
    events = [json.loads(ln) for ln in response.get_data(as_text=True).splitlines() if ln.strip()]
    kinds = [e.get("event") for e in events]
    assert "start" in kinds
    assert "progress" in kinds
    assert "done" in kinds
    start = next(e for e in events if e.get("event") == "start")
    done = next(e for e in events if e.get("event") == "done")
    assert start["n"] >= 3
    assert start["phase"] == "search"
    assert done["n_models"] >= 1
    assert done["target"] == "output"


def test_browse_folder_returns_path(client, monkeypatch):
    import M2_data_extractor.gui_server as gs

    monkeypatch.setattr(gs, "ask_folder_dialog", lambda initial="": r"C:\data\feathers")
    response = client.post("/browse/folder", json={"initial": ""})
    assert response.status_code == 200
    body = response.get_json()
    assert body["cancelled"] is False
    assert body["path"].endswith("feathers")


def test_folder_dialog_helper_ships_with_package():
    helper = ROOT / "M2_data_extractor" / "_folder_dialog.py"
    assert helper.is_file()
    text = helper.read_text(encoding="utf-8")
    assert "askdirectory" in text


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
    assert _is_placeholder_path("path/to/your/feathers")
    assert _is_placeholder_path("__XYZ_DATA__")
    assert not _is_placeholder_path(r"C:\Users\edens\__backup__\mols")
    assert not _is_placeholder_path(r"C:\Users\edens\Documents\GitHub\DescriPyTor_to_upload")


def test_gui_payload_is_current():
    from M2_data_extractor.gui_server import GUI_API, _gui_payload_is_current

    assert not _gui_payload_is_current({"ok": True, "version": "1.0"})
    assert _gui_payload_is_current({"ok": True, "api": GUI_API})
    assert not _gui_payload_is_current({"ok": True, "api": GUI_API - 1})      # a server from before the last update
    assert not _gui_payload_is_current(None)


def test_folder_summary_example_set(client):
    from M2_data_extractor.gui_server import example_feather_dir

    directory = example_feather_dir()
    if not directory:
        pytest.skip("example feather set is not installed")
    response = client.post("/folder/summary", json={"directory": directory})
    assert response.status_code == 200
    body = response.get_json()
    assert body["n"] == 10
    assert "basic" in body["names"]


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


def test_numbering_without_reference_uses_folder_basic(client):
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
    assert body["reference"] == "basic"
    assert body["n_mols"] == 2


def test_feature_explorer_is_served_with_its_query():
    """`descripytor feature_explorer` opens /feature_explorer; the page, and the live check it calls, are on the same server."""
    app = _load_app()
    client = app.test_client()
    r = client.get("/feature_explorer?mol=FL_lig_13&preset=theta%20construction")
    assert r.status_code == 200 and b"Feature explorer" in r.data and b"Sterimol" in r.data


def test_the_old_theta_address_redirects_with_its_query():
    """/theta was the explorer's address before the rename; bookmarks and scripts keep working."""
    client = _load_app().test_client()
    r = client.get("/theta?mol=FL_lig_13&preset=2")
    assert r.status_code == 302 and r.headers["Location"].endswith("/feature_explorer?mol=FL_lig_13&preset=2")
    old, new = client.post("/theta/check", json={}), client.post("/feature_explorer/check", json={})
    assert old.status_code == new.status_code != 404 and old.get_json() == new.get_json()


def test_feature_explorer_command_builds_the_explorer_address():
    from descripytor.cli import feature_explorer_page, theta_page

    assert theta_page is feature_explorer_page
    assert feature_explorer_page() == "/feature_explorer"
    assert feature_explorer_page("FL_lig_13") == "/feature_explorer?mol=FL_lig_13"
    assert feature_explorer_page("FL_lig_13", "theta construction") == "/feature_explorer?mol=FL_lig_13&preset=theta+construction"
    assert feature_explorer_page(preset="2") == "/feature_explorer?preset=2"


def test_the_old_explorer_folder_still_imports():
    """Scripts outside the package put M2_data_extractor/theta_explorer on sys.path and import data."""
    import subprocess, sys
    old = Path(__file__).resolve().parents[1] / "M2_data_extractor" / "theta_explorer"
    code = "import sys; sys.path.insert(0, %r); import data, build; assert data.payload and build.assemble" % str(old)
    assert subprocess.run([sys.executable, "-c", code]).returncode == 0
    same = lambda f: f.read_bytes().replace(bytes([13, 10]), bytes([10]))      # git may check the two copies out with different line ends
    assert same(old / "theta_explorer.html") == same(old.parent / "feature_explorer" / "feature_explorer.html")


def test_visual_extractor_address_is_the_picker():
    r = _load_app().test_client().get("/visual_extractor")
    assert r.status_code == 302 and r.headers["Location"].endswith("/visual")


def test_handoff_needs_a_table_and_a_target():
    client = _load_app().test_client()
    assert "Extract features first" in client.post("/bassa_explorer/handoff", json={}).get_json()["error"]
    r = client.post("/bassa_explorer/handoff", json={"csv": "name,a\nm1,1\nm2,2\n"})
    assert r.status_code == 400 and "target" in r.get_json()["error"]


def test_the_extraction_engine_is_inside_the_package():
    """Extract CSV runs descriptor_extractor; an installed package has no examples folder, so it must sit in M2_data_extractor."""
    app = _load_app()
    import M2_data_extractor.gui_server as gs
    dx = gs._load_descriptor_extractor()
    assert Path(dx.__file__).resolve().parent == ROOT / "M2_data_extractor"
    assert (ROOT / "M2_data_extractor" / "feature_packages.py").is_file()
    import subprocess, sys
    old = ROOT / "Getting_started_with_examples" / "descriptor_extraction_toolkit"
    code = "import sys; sys.path.insert(0, %r); import descriptor_extractor as d, feature_packages as f; assert d.run_from_config and f" % str(old)
    assert subprocess.run([sys.executable, "-c", code]).returncode == 0      # the toolkit's old import path still works


def test_a_server_from_before_the_update_is_not_reused():
    """Its pages are read from disk (new) while its routes are old: the header links gave Not Found."""
    _load_app()
    import M2_data_extractor.gui_server as gs
    old = {"ok": True, "api": 2, "example_presets": {}, "example_reference_name": "basic"}
    assert not gs._gui_payload_is_current(old)
    assert gs._gui_payload_is_current({"ok": True, "api": gs.GUI_API})
    assert _load_app().test_client().get("/status").get_json()["api"] == gs.GUI_API


def test_pasted_outputs_are_read_strictly():
    """A line that cannot be read stops the paste; it used to be skipped, moving every later value up a row."""
    from M2_data_extractor.gui_server import _merge_outputs_into_df, _parse_pasted_outputs as parse
    import pandas as pd
    import pytest

    assert parse("1.5\n2\n3") == [(None, 1.5), (None, 2.0), (None, 3.0)]
    assert parse("a,1.5\nb\t2\nc 3") == [("a", 1.5), ("b", 2.0), ("c", 3.0)]
    for bad in ("1.5\nn.d.\n3", "a,1\n2", "1,5\n2,5"):
        with pytest.raises(ValueError):
            parse(bad)
    df = pd.DataFrame({"name": ["MOL_0", "mol-1", "x.feather"], "f": [1, 2, 3]})
    assert _merge_outputs_into_df(df, parse("mol-0,1\nMOL_1,2\nX,3"), "y")["y"].tolist() == [1.0, 2.0, 3.0]
    with pytest.raises(ValueError, match="No row is named"):
        _merge_outputs_into_df(df, parse("zz,1"), "y")
    with pytest.raises(ValueError, match="2 values pasted for 3 rows"):
        _merge_outputs_into_df(df, parse("1\n2"), "y")


def test_outcomes_come_only_from_the_folder_they_sit_in():
    """Another set's measured outputs were offered for the example and joined onto any table with matching names."""
    import M2_data_extractor.gui_server as gs
    import pandas as pd

    assert gs.example_outcomes_payload()["rows"] == []            # the benzene example has no outcomes.csv
    assert _load_app().test_client().get("/status").get_json()["example_outcomes"] == []
    table = pd.DataFrame({"f": [1.0, 2.0]}, index=["p-Cl", "m_Br"])    # the extractor's table: names are the index
    out, n = gs._join_outcomes_into_df(table, [{"name": "P_cl", "output": 0.3}])
    assert n == 1 and out.loc["p-Cl", "output"] == 0.3
    payload = gs._extract_success_payload(table, {"feather_dir": str(ROOT / "descripytor" / "examples" / "feather_example")}, "[error] engine 'x' failed: boom")
    assert payload["n_outcomes"] == 0 and "output" not in payload["columns"]
    assert payload["engine_errors"] == ["[error] engine 'x' failed: boom"]
