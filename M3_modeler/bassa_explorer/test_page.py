"""Click-through test of the explorer page, in headless Chrome: every prepared case, every view, selection, merge,
leaving a compound out, design, and that the molecule frame shows the selected compound. Then the server: a small table
through /run with a structures folder, and the single-file export.

    python test_page.py        exits 0 and prints "ok" lines, or says what failed

The /run part needs bassa_reg in this Python and is skipped without it; BASSA_EXPLORER_CS2_FEATHERS names a folder of
structures for it (the CS2 feathers), else the run is made without structures.

The page tests itself (index.html, selfTest, behind ?selftest) and writes the outcome in its title; this script starts
the server on a spare port, opens the page and reads the title.
"""
import json, re, subprocess, sys, time, urllib.request
from pathlib import Path

H = Path(__file__).resolve().parent; PORT = 8799; URL = f"http://127.0.0.1:{PORT}/"
from settings import S

CHROME = S["chrome"]


def title(url):
    r = subprocess.run([CHROME, "--headless=new", "--disable-gpu", "--virtual-time-budget=60000", "--dump-dom", url], capture_output=True, text=True, encoding="utf-8", errors="replace", timeout=300)
    m = re.search(r"<title>(.*?)</title>", r.stdout, re.S); return m.group(1) if m else "no title: " + r.stderr[-300:]


def post(path, body):
    q = urllib.request.Request(URL + path, data=json.dumps(body).encode(), headers={"Content-Type": "application/json"}); return json.loads(urllib.request.urlopen(q, timeout=600).read())


server = subprocess.Popen([sys.executable, str(H / "server.py"), str(PORT)], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL); bad = []
try:
    for _ in range(60):
        try: urllib.request.urlopen(URL + "manifest.js", timeout=2); break
        except Exception: time.sleep(1)
    t = title(URL + "?selftest=1"); print("page:", t); bad += [] if t == "SELFTEST ok" else [t]
    try: import bassa_reg
    except ImportError: print("run: skipped (no bassa_reg in this Python)"); sys.exit("FAILED: " + "; ".join(bad) if bad else 0)
    d = post("run", dict(csv=(H / "example_small.csv").read_text(), response="ddG", iters=4000, name="t", structures=S["cs2_feathers"]))
    ok = d["n"] == 17 and d["draw"][0]["kind"] == "sterimol" and (not S["cs2_feathers"] or (len(d["mols"]) == 17 and d["mols"]["L1"]["charges"])); print("run%s:" % (" with a feather folder" if S["cs2_feathers"] else ""), "ok" if ok else "FAILED", d.get("no_structure")); bad += [] if ok else ["run"]
    html = post("export", dict(cases={"t": d}))["html"]; (H / "_selftest_standalone.html").write_text(html, encoding="utf-8")
    t = title((H / "_selftest_standalone.html").as_uri() + "?selftest=1"); print("single file:", t); bad += [] if t == "SELFTEST ok" else ["standalone: " + t]
finally:
    server.kill(); (H / "_selftest_standalone.html").unlink(missing_ok=True)
sys.exit("FAILED: " + "; ".join(bad) if bad else 0)
