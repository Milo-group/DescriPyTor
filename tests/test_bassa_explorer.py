"""The BASSA model explorer's server: the page, the bundled example, the package's feature explorer, and a kept case
that lands in the cases folder (never inside the package). No sampler run: bassa_reg is optional."""
import base64
import json
import os
import socket
import subprocess
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
TOOL = ROOT / "M3_modeler" / "bassa_explorer"


@pytest.fixture()
def server(tmp_path):
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0)); port = s.getsockname()[1]
    env = dict(os.environ, BASSA_EXPLORER_CASES=str(tmp_path / "cases"))
    p = subprocess.Popen([sys.executable, str(TOOL / "server.py"), str(port)], env=env, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    url = f"http://127.0.0.1:{port}/"
    try:
        for _ in range(60):
            try: urllib.request.urlopen(url + "manifest.js", timeout=2); break
            except Exception: time.sleep(0.5)
        yield url, tmp_path / "cases"
    finally:
        p.kill()


def _get(url):
    return urllib.request.urlopen(url, timeout=10).read().decode("utf-8")


def test_page_example_and_explorer_are_served(server):
    url, cases = server
    assert "BASSA model explorer" in _get(url)
    assert _get(url + "manifest.js").strip() == 'window.CASE_FILES = ["example"];'
    assert _get(url + "example.js").startswith('(window.CASES = window.CASES || {})["example"]')
    page = _get(url + "feature_explorer.html?embed=figure")
    assert "Feature explorer" in page and "embed" in page          # the package's built explorer, not a copy
    assert page == (ROOT / "M2_data_extractor" / "feature_explorer" / "feature_explorer.html").read_text(encoding="utf-8")
    assert not cases.exists()                                        # looking writes nothing


def test_a_kept_case_goes_to_the_cases_folder(server):
    url, cases = server
    t = (TOOL / "example.js").read_text(encoding="utf-8"); case = json.loads(t[t.index("] = ") + 4:-2])
    q = urllib.request.Request(url + "keep", data=json.dumps(dict(name="Mine 1", case=case)).encode(), headers={"Content-Type": "application/json"})
    out = json.loads(urllib.request.urlopen(q, timeout=30).read())
    assert out["name"] == "mine_1" and out["cases"] == ["example", "mine_1"]
    assert (cases / "mine_1.js").exists() and not (TOOL / "mine_1.js").exists()
    assert _get(url + "manifest.js").strip() == 'window.CASE_FILES = ["example", "mine_1"];'
    assert '"mine_1"' in _get(url + "mine_1.js")[:60]


def test_a_new_table_says_what_it_needs_without_the_sampler(server):
    url, _ = server
    try:
        import bassa_reg  # noqa: F401
        pytest.skip("bassa_reg is installed here; the run itself is test_page.py's job")
    except ImportError:
        pass
    q = urllib.request.Request(url + "run", data=json.dumps(dict(csv=(TOOL / "example_small.csv").read_text(), response="ddG")).encode(), headers={"Content-Type": "application/json"})
    with pytest.raises(urllib.error.HTTPError) as e:
        urllib.request.urlopen(q, timeout=30)
    assert "bassa-reg" in json.loads(e.value.read())["error"]


def test_the_three_pages_have_one_address_each(server):
    """/visual_extractor, /feature_explorer and /bassa_explorer mean the same page on both servers."""
    url, _ = server

    class NoRedirect(urllib.request.HTTPRedirectHandler):
        def redirect_request(self, *a, **k): return None
    with pytest.raises(urllib.error.HTTPError) as e:
        urllib.request.build_opener(NoRedirect).open(url + "bassa_explorer", timeout=10)
    assert e.value.code == 302 and e.value.headers["Location"] == "/"
    page = _get(url)
    assert 'href="/visual_extractor"' in page and 'href="/feature_explorer"' in page
    picker = (ROOT / "M2_data_extractor" / "atom_picker.html").read_text(encoding="utf-8")
    assert 'href="/feature_explorer"' in picker and 'href="/bassa_explorer"' in picker
    explorer = (ROOT / "M2_data_extractor" / "feature_explorer" / "template.html").read_text(encoding="utf-8")
    assert 'href="/visual_extractor"' in explorer and 'href="/bassa_explorer"' in explorer


def test_a_table_left_by_another_page_is_taken_once(server):
    """The picker's one-click hand-off: POST /handoff leaves a table, the page opened with #handoff=id takes it."""
    url, _ = server
    body = dict(csv=(TOOL / "example_small.csv").read_text(), response="ddG", structures="", name="extracted")
    q = urllib.request.Request(url + "handoff", data=json.dumps(body).encode(), headers={"Content-Type": "application/json"})
    token = json.loads(urllib.request.urlopen(q, timeout=30).read())["id"]
    got = json.loads(_get(url + "handoff?id=" + token))
    assert got["response"] == "ddG" and got["csv"] == body["csv"]
    with pytest.raises(urllib.error.HTTPError) as e:
        _get(url + "handoff?id=" + token)
    assert e.value.code == 404
    bad = urllib.request.Request(url + "handoff", data=json.dumps(dict(csv="name,a,y\nm1,1,2\n", response="y")).encode(), headers={"Content-Type": "application/json"})
    with pytest.raises(urllib.error.HTTPError) as e:              # a table a run would refuse is refused at the door
        urllib.request.urlopen(bad, timeout=30)
    assert "at least 6 compounds" in json.loads(e.value.read())["error"]


def test_the_extractor_s_own_column_names_are_drawn():
    """A table straight from `descripytor extractor` / the visual extractor: its dipole, Sterimol and vibration columns map to drawings."""
    code = (
        "import build_case as B\n"
        "want = {'dipole_x_1-2-19': 'dipole', 'total_dipole_1-2-19': 'dipole_total', 'loc_B1_1-4': 'sterimol', 'B1_B5_azimuth_4-1': 'sterimol',\n"
        "        'B1_B5_acute_1-4': 'sterimol', 'Frequency_Stretch_1_2': 'mark', 'Frequency_Bending_5-6': 'mark', 'Cross_mag_Bending_5-6': 'mark',\n"
        "        'cross_9_0': 'mark', 'cross_angle_9_0': 'mark', 'para_9_0': 'mark', 'para_angle_9_0': 'mark', 'B1_1-4': 'sterimol',\n"
        "        '0_bond_length_1-2': 'bond_length', 'nbo_atom_3': 'charge', 'nbo_diff_3-5': 'charge_diff'}\n"
        "for c, k in want.items():\n"
        "    w = B.drawing(c); assert w and w['kind'] == k, (c, w)\n"
        "    assert B.label(c, w) != c and B.describe(c, w)\n"
        "d = B.drawing('dipole_y_20-19-27'); assert (d['origin'], d['y'], d['plane'], d['comp']) == (20, 19, 27, 'u')\n"
        "assert B.drawing('cross_9_0')['atoms'] == [9] and B.drawing('aniso') is None and B.label('aniso', None) != 'aniso'\n"
    )
    r = subprocess.run([sys.executable, "-c", code], cwd=str(TOOL), capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-600:]


@pytest.fixture()
def other_server():
    """Something else on a port: no /status of ours, as a BASSA explorer from before the check, or any other program."""
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0)); port = s.getsockname()[1]
    p = subprocess.Popen([sys.executable, "-m", "http.server", str(port), "--bind", "127.0.0.1"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    try:
        for _ in range(60):
            try: urllib.request.urlopen(f"http://127.0.0.1:{port}/", timeout=2); break
            except Exception: time.sleep(0.25)
        yield port
    finally:
        p.kill()


def _running():
    sys.path.insert(0, str(ROOT))
    from M3_modeler.bassa_explorer import running
    return running


def test_a_running_server_says_which_version_it_is(server):
    url, _ = server
    running = _running()
    st = json.loads(_get(url + "status"))
    assert st["app"] == "bassa_explorer" and st["api"] == running.API
    port = int(url.rsplit(":", 1)[1].strip("/"))
    assert running.state(port) == "current"
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0)); free = s.getsockname()[1]
    assert running.state(free) is None


def test_an_older_server_on_the_port_is_not_reused(other_server):
    """A server left from before an update would serve the new page on old routes; every way in refuses it and says how to stop it."""
    running = _running()
    assert running.state(other_server) == "old"
    assert str(other_server) in running.stop_hint(other_server, "the BASSA model explorer")
    # the command
    r = subprocess.run([sys.executable, "-m", "descripytor.cli", "bassa_explorer", "--no-browser", "--port", str(other_server)],
                       cwd=str(ROOT), capture_output=True, text=True, timeout=120)
    assert r.returncode == 1 and "older or different process" in r.stdout
    # the visual extractor's link and its one-click hand-off
    env = dict(os.environ, BASSA_EXPLORER_PORT=str(other_server))
    code = ("import sys; sys.path.insert(0, %r)\n"
            "import M2_data_extractor.gui_server as gs\n"
            "c = gs.app.test_client()\n"
            "r = c.get('/bassa_explorer'); assert r.status_code == 503 and b'older or different process' in r.data, r.data\n"
            "csv = 'name,a,b,c,d,y' + chr(10) + chr(10).join('m%%d,%%d,%%d,%%d,%%d,%%d' %% (i, i, i * i %% 7, i * 3 %% 5, i * 5 %% 11, i %% 4) for i in range(8))\n"
            "r = c.post('/bassa_explorer/handoff', json={'csv': csv, 'y_column': 'y'}); assert r.status_code == 503 and 'older' in r.get_json()['error'], r.data\n") % str(ROOT)
    r = subprocess.run([sys.executable, "-c", code], env=env, capture_output=True, text=True, timeout=120)
    assert r.returncode == 0, r.stderr[-800:]


def test_names_and_tables_from_the_review():
    code = (
        "import io, numpy as np, pandas as pd, build_case as B, server as S\n"
        "assert B.drawing('q_cm5[Cl12]')['i'] == 12 and B.drawing('q_cm5[12]')['i'] == 12 and B.drawing('nbo_atom_0') is None\n"
        "assert B.drawing('ster_B1[N5->C6]')['a'] == 5 and B.drawing('bond_length_1-2')['atoms'] == [1, 2]\n"
        "assert 'Ni' not in B.label('B1', B.drawing('B1'))\n"
        "t = pd.DataFrame({'name': ['m%d' % i for i in range(8)], ' a': range(8), 'b ': [i * i % 7 for i in range(8)], 'c': [i * 3 % 5 for i in range(8)],\n"
        "                  'd': [i * 5 % 11 for i in range(8)], 'e': [1.0, 2, 3, np.inf, 5, 6, 7, 8], ' y ': [i % 4 + 0.5 * i for i in range(8)]})\n"
        "X, y, bad, skipped, ext = S.table(dict(csv=t.to_csv(index=False), response='y'))\n"
        "assert list(X.columns) == ['a', 'b', 'c', 'd'] and bad == ['e'] and np.isfinite(X.to_numpy()).all() and ext is None\n"
        "p = pd.DataFrame({'acid': ['A1', 'A1', 'A2', 'A2', 'A3', 'A3', 'A4', 'A4'], 'amine': ['N1', 'NA'] * 4, 'y': range(8), 'a': [1, 3, 2, 5, 4, 7, 6, 9], 'b': [i % 3 for i in range(8)], 'c': [i * 7 % 5 for i in range(8)], 'd': [i * i for i in range(8)]})\n"
        "X, y, bad, skipped, ext = S.table(dict(csv=p.to_csv(index=False), response='y'))\n"
        "assert list(X.index[:2]) == ['A1-N1', 'A1-NA'] and list(X.columns) == ['a', 'b', 'c', 'd']      # two name columns: the pair, and NA is a name\n"
        "s = t.assign(Set=['train'] * 6 + ['test', 'external'])\n"
        "X, y, bad, skipped, ext = S.table(dict(csv=s.to_csv(index=False), response='y'))\n"
        "assert list(X.index) == ['m%d' % i for i in range(6)] and ext['names'] == ['m6', 'm7'] and ext['y'] == [5.0, 6.5] and list(ext['X']) == list(X.columns) and skipped == []      # the split column: test rows held out\n"
    )
    r = subprocess.run([sys.executable, "-c", code], cwd=str(TOOL), capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-800:]


def _featurize(url, case, files, headers=None):
    body = dict(case=case, files=[dict(name=Path(f).name if isinstance(f, (str, Path)) else f[0], data=base64.b64encode(Path(f).read_bytes() if isinstance(f, (str, Path)) else f[1]).decode()) for f in files])
    q = urllib.request.Request(url + "featurize", data=json.dumps(body).encode(), headers=dict({"Content-Type": "application/json"}, **(headers or {})))
    try: return 200, json.loads(urllib.request.urlopen(q, timeout=300).read())
    except urllib.error.HTTPError as e: return e.code, (json.loads(e.read() or b"{}") if e.code == 400 else {})


def test_a_structure_file_becomes_a_candidate_by_the_case_s_recipe(server):
    """POST /featurize: a measured compound's own files (a cs3_cp complex and its stored GFN2 single point), dropped in as
    new, give exactly the case's X, so the page prices it at its fitted value; what does not fit the recipe is refused, with the reason."""
    url, cases = server
    fx = ROOT / "tests" / "fixtures" / "bassa_featurize"; d = json.loads((fx / "case_cs3_cp_045.json").read_text(encoding="utf-8"))
    cases.mkdir(parents=True, exist_ok=True); (cases / "fx.js").write_text('(window.CASES = window.CASES || {})["fx"] = %s;\n' % json.dumps(d), encoding="utf-8")
    xyz, sp = fx / "045_lig_CuCl.xyz", fx / "045_lig_CuCl.txt"
    code, out = _featurize(url, "fx", [xyz, sp])
    assert code == 200, out
    assert out["x"] == {c: d["X"][c][0] for c in d["names"]}
    pred = lambda x: d["ens"]["b0"] + sum(d["ens"]["beta"].get(c, 0) * x[c] for c in d["names"])
    assert pred(out["x"]) == pred({c: d["X"][c][0] for c in d["names"]})
    # refused, with the reason
    assert "two files" in _featurize(url, "fx", [xyz])[1]["error"]
    assert "reads one .xyz" in _featurize(url, "fx", [("x.bin", b"not a structure")])[1]["error"]
    rows = xyz.read_text().splitlines(); moved = "\n".join(rows[:2] + rows[3:4] + rows[2:3] + rows[4:]) + "\n"        # metal no longer atom 1
    assert "atom 1 is" in _featurize(url, "fx", [("m.xyz", moved.encode()), sp])[1]["error"]
    d2 = dict(d); d2.pop("extraction"); (cases / "plain.js").write_text('(window.CASES = window.CASES || {})["plain"] = %s;\n' % json.dumps(d2), encoding="utf-8")
    assert "no recipe" in _featurize(url, "plain", [xyz, sp])[1]["error"]
    assert "no such case" in _featurize(url, "../fx", [xyz, sp])[1]["error"]
    assert _featurize(url, "fx", [xyz, sp], {"Origin": "http://evil.example"})[0] == 403
    q = urllib.request.Request(url + "featurize", data=b"{}", headers={"Content-Type": "text/plain"})
    with pytest.raises(urllib.error.HTTPError) as e:
        urllib.request.urlopen(q, timeout=30)
    assert e.value.code == 415
