"""The four pictures of the README: the three pages and the hand-off between them.

Start the pages first, then run this with a Python that has playwright:

    descripytor visual_extractor --no-browser --port 7441      (with BASSA_EXPLORER_PORT=8772 and a Python that has bassa_reg
                                                                in DESCRIPYTOR_BASSA_PYTHON, for the hand-off picture)
    python docs/images/readme/make_images.py [http://127.0.0.1:7441]

It goes through the bundled set that has measured outputs (descripytor/examples/baptiste_products): its picks,
Extract CSV, its outcomes, Model in BASSA explorer.
"""
import json
import sys
from pathlib import Path

from playwright.sync_api import sync_playwright

HERE = Path(__file__).resolve().parent
EX = HERE.parents[2] / "descripytor" / "examples" / "baptiste_products"
PRESETS = json.loads((EX / "presets.json").read_text(encoding="utf-8"))
OUTCOMES = chr(10).join((EX / "outcomes.csv").read_text(encoding="utf-8").strip().splitlines()[1:])      # name,value lines
GUI = (sys.argv[1] if len(sys.argv) > 1 else "http://127.0.0.1:7441").rstrip("/")

with sync_playwright() as p:
    b = p.chromium.launch(channel="chrome")            # the installed Chrome: nothing to download
    ctx = b.new_context(viewport=dict(width=1500, height=900), device_scale_factor=1.5)
    pg = ctx.new_page()

    # 1. the feature explorer on a built-in structure
    pg.goto(GUI + "/feature_explorer?mol=FL_lig_13&preset=theta%20construction", wait_until="networkidle")
    pg.wait_for_selector("#fig svg .hit")
    pg.wait_for_timeout(800)
    pg.screenshot(path=str(HERE / "feature_explorer.png"))

    # 2. the visual extractor with the example picks drawn
    pg.goto(GUI + "/visual", wait_until="networkidle")
    pg.wait_for_function("typeof figReady !== 'undefined' && figReady && nAtoms > 0", timeout=60000)
    pg.evaluate("localStorage.setItem('descripytorTourSeen', '1')")
    # the bundled set with measured outputs: 19 substrates (examples/baptiste_products), its picks and its outcomes
    pg.evaluate("""async ([dir, atoms]) => {
        const r = await fetch('/xyz_from_feather', { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify({ filepath: dir + '/unsub.feather' }) });
        const d = await r.json(); if (!d.xyz) throw new Error(JSON.stringify(d).slice(0, 200));
        showMolecule(d.xyz, 'unsub.feather', d.mol_data || d.molData || d.data || null, { silent: true });
        document.getElementById('featherDirInput').value = dir;
        clearAllPicks(); Object.entries(atoms).forEach(([k, groups]) => { state[k] = groups.map((g) => g.slice()); }); refresh();
        ['vizSterimol'].forEach((id) => { const c = document.getElementById(id); if (c) c.checked = true; });
        buildSterPairOptions(); highlight(); }""", [str(EX).replace(chr(92), "/"), PRESETS["atoms"]])
    pg.wait_for_timeout(1500)
    pg.screenshot(path=str(HERE / "visual_extractor.png"))

    # 3. extract, give the outputs, hand over
    pg.evaluate("extractFeatures()")
    pg.wait_for_function("lastExtractCsv.length > 0", timeout=300000)
    pg.evaluate("(t) => { showOutTab('model'); document.getElementById('lrOutputs').value = t; }", OUTCOMES)
    pg.wait_for_timeout(500)
    print("outputs given:", len(pg.evaluate("document.getElementById('lrOutputs').value").splitlines()), "lines")
    pg.locator("#bassaBtn").scroll_into_view_if_needed()
    pg.screenshot(path=str(HERE / "handoff.png"))
    with ctx.expect_page() as new:
        pg.click("#bassaBtn")
    bp = new.value
    bp.wait_for_selector('#newstatus', state='attached', timeout=180000)      # the tab opens blank, then goes to the BASSA page
    bp.wait_for_function("typeof C !== 'undefined' && C && C.uploaded || /Could not/.test(document.getElementById('newstatus').textContent)", timeout=900000)
    print("picker says:", pg.evaluate("document.getElementById('lrStatus').textContent"))
    print("BASSA page:", bp.evaluate("C && C.uploaded ? C.title + ' | ' + C.n + ' compounds, ' + C.names.length + ' descriptors, ' + Object.keys(C.mols).length + ' structures' : document.getElementById('newstatus').textContent"))

    # 4. the BASSA model explorer on that run
    bp.evaluate("() => { try { document.getElementById('guide').close(); } catch (e) {} }")
    bp.wait_for_timeout(2500)
    bp.screenshot(path=str(HERE / "bassa_explorer.png"))
    b.close()
print("wrote", sorted(f.name for f in HERE.glob("*.png")))
