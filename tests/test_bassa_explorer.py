"""The BASSA model explorer's server: the page, the bundled example, the package's feature explorer, and a kept case
that lands in the cases folder (never inside the package). No sampler run: bassa_reg is optional."""
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
