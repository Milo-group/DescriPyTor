"""Atom renumbering with Conformaty: the adapter (M2_data_extractor/alignment.py), `descripytor renumber` and the visual
extractor's /align routes. The mapping tests run on basic.feather and a copy with its atoms shuffled, which must come
back to basic's numbering exactly (identical geometry, so Conformaty takes seconds)."""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT / "descripytor" / "examples" / "feather_example"


def _molecules(folder):
    from M2_data_extractor.data_extractor import Molecules

    return Molecules(str(folder))


def _shuffle(folder, seed=5):
    """folder/basic.feather and folder/shuffled.feather (basic's atoms in a known random order); returns the order."""
    import shutil
    from M2_data_extractor.data_extractor import Molecule

    shutil.copy(EXAMPLE / "basic.feather", folder)
    m = Molecule(str(EXAMPLE / "basic.feather"))
    perm = np.random.default_rng(seed).permutation(len(m.xyz_df)) + 1
    m.renumber_atoms({i + 1: int(p) for i, p in enumerate(perm)}, save_feather=True, verbose=False,
                     out_file=str(folder / "shuffled.feather"))
    return perm


def test_kabsch_puts_a_moved_copy_back():
    from M2_data_extractor.alignment import kabsch

    rng = np.random.default_rng(0)
    xyz = rng.normal(size=(12, 3))
    q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    q *= np.sign(np.linalg.det(q))                       # a rotation, not a reflection
    moved = xyz @ q.T + [3.0, -1.0, 2.0]
    assert np.abs(kabsch(moved, xyz, list(range(12))) - xyz).max() < 1e-9


def test_written_folder_is_renumbered_and_keeps_names(tmp_path):
    """write_renumbered applies {old: new}: atom `old` of the original is atom `new` of the copy, names unchanged."""
    from M2_data_extractor.alignment import write_renumbered

    src = tmp_path / "in"
    src.mkdir()
    import shutil
    shutil.copy(EXAMPLE / "m_Br.feather", src)
    (src / "outcomes.csv").write_text("name,output\nm_Br,1.0\n")
    mols = _molecules(src)
    n = len(mols.molecules[0].xyz_df)
    old_to_new = {i: i for i in range(1, n + 1)}
    old_to_new[1], old_to_new[2] = 2, 1
    written, skipped = write_renumbered(mols, {"m_Br": {"old_to_new": old_to_new}}, str(tmp_path / "out"), "m_Br")
    assert written == ["m_Br"] and skipped == []
    back = _molecules(tmp_path / "out").molecules[0]
    a, b = mols.molecules[0].xyz_df, back.xyz_df
    assert back.molecule_name == "m_Br" and (tmp_path / "out" / "outcomes.csv").is_file()
    assert np.allclose(b.iloc[1][["x", "y", "z"]].astype(float), a.iloc[0][["x", "y", "z"]].astype(float))
    saved = json.loads((tmp_path / "out" / "renumbering.json").read_text())
    assert saved["reference"] == "m_Br" and saved["molecules"]["m_Br"]["old_to_new"]["1"] == 2


def test_a_shuffled_molecule_is_numbered_back(tmp_path):
    pytest.importorskip("conformaty")
    from M2_data_extractor.alignment import FITS, map_folder

    perm = _shuffle(tmp_path)
    events = []
    res = map_folder(_molecules(tmp_path), "basic", progress=events.append)
    r = res["shuffled"]
    # basic's atom i+1 became the shuffled copy's atom perm[i]; numbered back, that atom is i+1 again
    assert all(r["old_to_new"][int(p)] == i + 1 for i, p in enumerate(perm))
    assert r["differs"] == [] and r["rmsd"] < 1e-3 and r["fit"] == "rmsd" == FITS[0] and r["rmsd_kabsch"] < 1e-3
    assert [e["event"] for e in events] == ["start", "molecule"] and events[1]["moved"] == len(perm) - (perm == np.arange(1, len(perm) + 1)).sum()
    with pytest.raises(ValueError):
        map_folder(_molecules(tmp_path), "basic", fit="nonsense")


def test_the_renumber_command_writes_the_folder(tmp_path):
    pytest.importorskip("conformaty")
    perm = _shuffle(tmp_path)
    run = subprocess.run([sys.executable, "-m", "descripytor.cli", "renumber", str(tmp_path), "-r", "basic", "--fit", "kabsch"],
                         capture_output=True, text=True, cwd=ROOT)
    assert run.returncode == 0, run.stderr[-2000:]
    assert "1/1  shuffled" in run.stdout and "Wrote 2 renumbered molecules" in run.stdout
    out = _molecules(tmp_path / "renumbered")
    by = {m.molecule_name: m.xyz_df for m in out.molecules}
    assert np.allclose(by["shuffled"][["x", "y", "z"]].astype(float), by["basic"][["x", "y", "z"]].astype(float))
    assert len(perm) == len(by["basic"])


@pytest.fixture()
def client():
    from M2_data_extractor.gui_server import app

    app.config["TESTING"] = True
    return app.test_client()


def test_apply_refuses_an_incomplete_numbering(client, tmp_path):
    import shutil
    shutil.copy(EXAMPLE / "basic.feather", tmp_path)
    r = client.post("/align/apply", json={"directory": str(tmp_path), "molecules": {"basic": {"1": 2, "2": 1}}})
    assert r.status_code == 400 and "each of its" in r.get_json()["error"]
    r = client.post("/align/apply", json={"directory": str(tmp_path), "molecules": {"nobody": {"1": 1}}})
    assert r.status_code == 400 and not (tmp_path / "renumbered").exists()


def test_the_map_route_streams_one_molecule_at_a_time(client, tmp_path):
    pytest.importorskip("conformaty")
    _shuffle(tmp_path)
    r = client.post("/align/map", json={"directory": str(tmp_path), "reference": "basic", "fit": "kabsch"})
    events = [json.loads(line) for line in r.get_data(as_text=True).splitlines() if line.strip()]
    assert [e["event"] for e in events] == ["start", "molecule", "done"]
    assert events[0]["names"] == ["shuffled"] and events[0]["ref_xyz"].startswith(str(events[0]["n_ref"]))
    result = events[1]["result"]
    fits = client.post("/align/fits", json={"directory": str(tmp_path), "reference": "basic", "name": "shuffled",
                                            "old_to_new": result["old_to_new"]}).get_json()["fits"]
    assert [f["fit"] for f in fits] == ["rmsd", "msd", "mixed", "mad", "sig", "kabsch"]
    assert fits[-1]["rmsd"] < 1e-3                       # the same molecule: Kabsch puts it back exactly
    wrote = client.post("/align/apply", json={"directory": str(tmp_path), "reference": "basic",
                                              "molecules": {"shuffled": result["old_to_new"]}}).get_json()
    assert wrote["written"] == ["shuffled"] and Path(wrote["out"], "renumbering.json").is_file()


def test_stop_ends_a_running_molecule(client, tmp_path):
    """Stop (or a page that was closed) used to wait for the running molecule, minutes, holding the server's one-run lock."""
    pytest.importorskip("conformaty")
    import shutil
    import threading
    import time

    for name in ("basic", "o_Cl"):                       # o_Cl takes tens of seconds
        shutil.copy(EXAMPLE / f"{name}.feather", tmp_path)
    events = []
    run = threading.Thread(target=lambda: events.append(client.post("/align/map", json={
        "directory": str(tmp_path), "reference": "basic"}).get_data(as_text=True)))
    run.start()
    time.sleep(4)
    t0 = time.perf_counter()
    assert client.post("/align/map", json={"directory": str(tmp_path)}).status_code == 409
    stopped = client.post("/align/stop").get_json()
    run.join(30)
    assert stopped["free"] and time.perf_counter() - t0 < 15 and not run.is_alive()
    assert '"event": "stopped"' in events[0]
