"""Regenerate the screenshots of docs/THETA_EXPLORER.md.

    pip install playwright pillow      # Chrome is used from the machine (channel="chrome")
    python docs/images/theta_explorer/make_images.py [name ...]

Serves M2_data_extractor/theta_explorer on a free local port, drives the page, adds numbered callouts
at the real element positions (so a layout change moves them with it) and writes 2x PNGs next to this
file. With names, only those images are rebuilt. The numbering of the callouts is the numbering of
the lists in the manual: edit both together.
"""
import socket
import subprocess
import sys
import time
from pathlib import Path
from urllib.parse import quote

from PIL import Image, ImageDraw, ImageFont
from playwright.sync_api import sync_playwright

HERE = Path(__file__).resolve().parent
EXPLORER = HERE.parents[2] / "M2_data_extractor" / "theta_explorer"
W, H, DPR = 1480, 1000, 2


def _free_port():
    with socket.socket() as k:
        k.bind(("127.0.0.1", 0))
        return k.getsockname()[1]


PORT = _free_port()
URL = f"http://127.0.0.1:{PORT}/theta_explorer.html"

BADGE_JS = """
(args) => {
  const [items, ink] = args;
  const mk = (n, x, y) => {
    const d = document.createElement('div');
    d.className = 'mbadge';
    d.textContent = n;
    d.style.cssText = `position:absolute;left:${x - 11}px;top:${y - 11}px;width:22px;height:22px;border-radius:50%;` +
      `background:${ink};color:#fff;font:700 13px/22px system-ui,sans-serif;text-align:center;z-index:99999;` +
      `box-shadow:0 0 0 2px #fff;pointer-events:none`;
    document.body.appendChild(d);
  };
  for (const [n, sel, anchor, dx, dy, re] of items) {
    let el = null;
    if (re) el = [...document.querySelectorAll(sel)].find((e) => new RegExp(re).test(e.textContent.replace(/\\s+/g, '')));
    else el = document.querySelector(sel);
    if (!el) { console.warn('badge target missing', n, sel, re); continue; }
    const r = el.getBoundingClientRect();
    const x = anchor.includes('r') ? r.right : anchor.includes('l') ? r.left : r.left + r.width / 2;
    const y = anchor.includes('b') ? r.bottom : anchor.includes('t') ? r.top : r.top + r.height / 2;
    mk(n, x + scrollX + (dx || 0), y + scrollY + (dy || 0));
  }
}
"""


class Shots:
    def __init__(self, pw):
        self.b = pw.chromium.launch(channel="chrome", headless=True)

    def page(self, query="", width=W, height=H, wait=True):
        pg = self.b.new_page(viewport={"width": width, "height": height}, device_scale_factor=DPR)
        pg.on("console", lambda m: print("   console:", m.text[:110]) if m.type == "warning" and "badge" in m.text else None)
        pg.goto(URL + query, wait_until="domcontentloaded")
        if wait:
            pg.wait_for_selector("#fig circle.hit", state="attached", timeout=20000)
        pg.wait_for_timeout(1200)
        return pg

    @staticmethod
    def chips(pg, **state):
        """Switch Show chips and construction toggles on or off: chips(pg, sterimol=True, dipole=False)."""
        for key, want in state.items():
            el = pg.locator(f'button[data-chip="{key}"], button[data-flag="{key}"]').first
            if (el.get_attribute("aria-pressed") == "true") != want:
                el.click()
        pg.wait_for_timeout(400)

    @staticmethod
    def settle(pg):
        """Top of the page, toolbar un-stuck: sticky bars would otherwise lie across a clipped section."""
        pg.add_style_tag(content=".bar{position:static!important}")
        pg.evaluate("window.scrollTo(0, 0)")
        pg.wait_for_timeout(300)

    def badges(self, pg, items, ink="#c0392b"):
        self.settle(pg)
        pg.evaluate(BADGE_JS, [items, ink])

    @staticmethod
    def region(pg, *selectors, pad=10):
        r = pg.evaluate("""(sels) => { const b = sels.map((s) => document.querySelector(s).getBoundingClientRect());
            return [Math.min(...b.map(x => x.left)) + scrollX, Math.min(...b.map(x => x.top)) + scrollY,
                    Math.max(...b.map(x => x.right)) + scrollX, Math.max(...b.map(x => x.bottom)) + scrollY]; }""", list(selectors))
        return dict(x=max(0, r[0] - pad), y=max(0, r[1] - pad), width=r[2] - r[0] + 2 * pad, height=r[3] - r[1] + 2 * pad)

    @staticmethod
    def done(path, max_w=1900):
        im = Image.open(path)
        if im.width > max_w:                                    # 2x captures are crisper than a manual page needs
            im = im.resize((max_w, round(im.height * max_w / im.width)), Image.LANCZOS)
        im.save(path, optimize=True)
        print(f"  {path.name}  {im.size[0]}x{im.size[1]}  {path.stat().st_size // 1024} KB")

    def save(self, pg, name, locator=None, clip=None):
        self.settle(pg)
        path = HERE / f"{name}.png"
        if locator is not None:
            pg.locator(locator).first.screenshot(path=path)
        else:
            pg.screenshot(path=path, clip=clip, full_page=True)
        self.done(path)

    @staticmethod
    def atom(pg, label_re, nth=0):
        """Viewport centre of the nth atom whose title matches label_re."""
        return pg.evaluate("""([re, nth]) => { const c = [...document.querySelectorAll('#fig circle.hit')].filter((e) => new RegExp(re).test(e.textContent))[nth];
            if (!c) return null; const r = c.getBoundingClientRect(); return [r.left + r.width / 2, r.top + r.height / 2]; }""", [label_re, nth])

    @staticmethod
    def open_section(pg, title):
        """Close every side-panel section but `title` and put it at the top of the panel."""
        pg.evaluate("""(t) => { const A = document.querySelector('aside#layers');
            A.querySelectorAll('details').forEach((d) => { d.open = d.querySelector('summary').textContent.trim().startsWith(t); });
            A.scrollTop = 0; window.scrollTo(0, 0); }""", title)
        pg.wait_for_timeout(500)

    @staticmethod
    def pick(pg, label_text, option_text):
        """Choose an option in the select that belongs to the label containing label_text."""
        pg.evaluate("""([l, o]) => { const lab = [...document.querySelectorAll('label')].find((x) => x.textContent.trim().startsWith(l));
            const sel = lab.querySelector('select') || document.getElementById(lab.htmlFor);
            const opt = [...sel.options].find((x) => x.textContent.startsWith(o)); sel.value = opt.value;
            sel.dispatchEvent(new Event('change', { bubbles: true })); }""", [label_text, option_text])
        pg.wait_for_timeout(500)

    def theta_view(self, mol="FL_lig_13", values=False):
        """The explorer's own teaching preset, with values on the labels."""
        pg = self.page(f"?mol={mol}&preset=" + quote("theta construction"))
        if values:
            self.pick(pg, "feature labels", "symbol and value")
        return pg


FONT = next((f for f in ("C:/Windows/Fonts/segoeui.ttf", "C:/Windows/Fonts/arial.ttf") if Path(f).exists()), None)


def grid(paths, captions, out, cols=2, cap_h=70, gap=24):
    ims = [Image.open(p).convert("RGB") for p in paths]
    w, h = max(i.width for i in ims), max(i.height for i in ims)
    rows = (len(ims) + cols - 1) // cols
    sheet = Image.new("RGB", (cols * w + (cols + 1) * gap, rows * (h + cap_h) + (rows + 1) * gap), "white")
    d = ImageDraw.Draw(sheet)
    font = ImageFont.truetype(FONT, 40) if FONT else ImageFont.load_default()
    for k, (im, cap) in enumerate(zip(ims, captions)):
        x, y = gap + (k % cols) * (w + gap), gap + (k // cols) * (h + cap_h + gap)
        d.text((x + 8, y + 8), cap, fill="#222", font=font)
        sheet.paste(im, (x + (w - im.width) // 2, y + cap_h))
    sheet.save(out, optimize=True)
    Shots.done(out)


# ---------------------------------------------------------------------------------------- images
def overview(s):
    pg = s.theta_view()
    s.badges(pg, [(1, "#mol", "tl", -2, -4, None), (2, "label[for=a]", "tl", -4, -4, None), (3, "#orient", "tl", -4, -4, None),
                  (4, "[data-chip=sterimol]", "tl", -4, -4, None), (5, "[data-flag=axis]", "tl", -4, -4, None),
                  (6, "#fig", "tl", 30, 40, None), (7, "#fig svg text", "tl", -10, -8, r"^end-on"),
                  (8, "#layers", "tl", 4, 4, None), (9, "#vals", "tl", 0, 0, None), (10, ".exports", "tl", -4, -4, None)])
    s.save(pg, "01_overview", clip=s.region(pg, "main", pad=16))
    pg.close()


def toolbar(s):
    pg = s.theta_view()
    s.badges(pg, [(1, "#mol", "tl", -2, -4, None), (2, ".filebtn", "tl", -2, -4, None), (3, "label[for=a]", "tl", -4, -4, None),
                  (4, "#swap", "tl", -4, -4, None), (5, "label[for=b]", "tl", -4, -4, None), (6, "#orient", "tl", -4, -4, None),
                  (7, "[data-chip=sterimol]", "tl", -4, -4, None), (8, "[data-flag=axis]", "tl", -4, -4, None),
                  (9, "[data-flag=profile]", "tl", -4, -4, None), (10, "details.paste summary", "tl", -4, -4, None)])
    s.save(pg, "02_toolbar", clip=s.region(pg, ".bar", pad=18))
    pg.close()


def anatomy(s):
    pg = s.theta_view()
    s.badges(pg, [(1, "#fig svg text", "l", -12, 0, r"^L$"), (2, "#fig svg text", "l", -12, 0, r"^B1$"),
                  (3, "#fig svg text", "l", -12, 0, r"^B5$"), (4, "#fig svg text", "l", -12, 0, r"^locB5"),
                  (5, "#fig svg text", "l", -12, 0, r"^θ"), (6, "#fig svg text", "tl", -10, -8, r"^end-on")])
    s.save(pg, "03_anatomy", locator="#figwrap")
    pg.close()


def popup(s):
    pg = s.theta_view()
    xy = s.atom(pg, r"^C\d+$", 2)
    pg.mouse.click(*xy)
    pg.wait_for_selector("#pop", state="visible")
    pg.wait_for_timeout(400)
    s.badges(pg, [(1, "#pop button[data-role=a]", "l", -12, 0, None), (2, "#pop button[data-role=b]", "l", -12, 0, None),
                  (3, "#pop button[data-act=axis]", "l", -12, 0, None), (4, "#pop button[data-act=q]", "l", -12, 0, None),
                  (5, "#pop button[data-act=fade]", "l", -12, 0, None)])
    s.save(pg, "04_atom_popup", locator="#figwrap")
    pg.close()


def show_chips(s):
    tmp = []
    cases = [("Dipole", "?mol=m-CN", dict(sterimol=False, dipole=True), None),
             ("Charges", "?mol=m-CN", dict(sterimol=False, dipole=False, charges=True), None),
             ("Cone", "?mol=Cu-PHOX", dict(sterimol=False, cone=True), "Cu"),
             ("%Vbur", "?mol=Cu-PHOX", dict(sterimol=False, vbur=True), "Cu")]
    for name, q, state, metal in cases:
        pg = s.page(q)
        if metal:
            lab = pg.evaluate("[...document.querySelectorAll('#fig circle.hit')].map(c => c.textContent).find(t => /^Cu\\d+$/.test(t))")
            pg.fill("#a", lab[2:])
            pg.press("#a", "Tab")
            pg.wait_for_timeout(300)
        s.chips(pg, **state)
        path = HERE / f"_tmp_{name[:3]}.png"
        pg.locator("#figwrap").screenshot(path=path)
        tmp.append(path)
        pg.close()
    grid(tmp, ["Dipole", "Charges", "Cone", "%Vbur"], HERE / "05_show_chips.png")
    for p in tmp:
        p.unlink()


def b1_scan(s):
    pg = s.theta_view()
    s.chips(pg, profile=True)
    s.save(pg, "06_b1_scan", locator="#figwrap")
    pg.close()


def side_panel(s):
    pg = s.page("?mol=m-CN")
    pg.evaluate("""() => { document.querySelectorAll('aside#layers details').forEach((d) => d.open = true);
        const a = document.querySelector('#layers'); a.style.maxHeight = 'none'; a.style.height = 'auto'; a.style.overflow = 'visible'; }""")
    pg.wait_for_timeout(500)
    tall = HERE / "_tmp_panel.png"
    pg.locator("#layers").screenshot(path=tall)
    pg.close()
    im = Image.open(tall).convert("RGB")
    half = im.height // 2
    cut = half
    px = im.load()
    for y in range(half, min(im.height, half + 400)):           # cut in a blank row between two sections
        if all(px[x, y] == px[0, y] for x in range(0, im.width, 7)):
            cut = y
            break
    a, b = im.crop((0, 0, im.width, cut)), im.crop((0, cut, im.width, im.height))
    out = Image.new("RGB", (im.width * 2 + 40, max(a.height, b.height)), "white")
    out.paste(a, (0, 0))
    out.paste(b, (im.width + 40, 0))
    path = HERE / "07_side_panel.png"
    out.save(path, optimize=True)
    s.done(path)
    tall.unlink()


def stack(s):
    pg = s.theta_view()
    s.open_section(pg, "Structures")
    pg.select_option("#stAddSel", label="FL_lig_7")
    pg.click("#stAdd")
    pg.wait_for_timeout(900)
    pg.select_option("#alMode", "atoms")
    pg.fill("#alAtoms", "6, 7, 4")
    pg.press("#alAtoms", "Tab")
    pg.wait_for_timeout(900)
    s.open_section(pg, "Structures")
    s.badges(pg, [(1, "#stAddSel", "tl", -4, -4, None), (2, "#stAdd", "tl", -4, -4, None), (3, "#stList", "tl", -4, -4, None),
                  (4, "#alMode", "tl", -4, -4, None), (5, "#alAtoms", "tl", -4, -4, None), (6, "#spread", "tl", -4, -4, None)])
    s.save(pg, "08_stack", locator="section.stage")
    pg.close()


def sheet(s):
    q = "?sheet=" + quote("m-CN::Case Study 1;FL_lig_13::Case Study 2;045_lig::Case Study 3") + "&cols=3"
    pg = s.page(q)
    s.open_section(pg, "Sheet")
    s.badges(pg, [(1, "#shAdd", "tl", -4, -4, None), (2, "#shList", "tl", -4, -4, None), (3, "#shCols", "tl", -4, -4, None),
                  (4, "#shSvg", "tl", -4, -4, None)])
    s.save(pg, "09_sheet", locator="section.stage")
    pg.close()


def figure_mode(s):
    pg = s.page("?mol=m-CN&preset=" + quote("Case Study 1") + "&figure=1", width=900, height=700, wait=False)
    pg.wait_for_selector("svg", state="attached")
    pg.screenshot(path=HERE / "10_figure_mode.png", full_page=True)
    s.done(HERE / "10_figure_mode.png")
    pg.close()


def export_row(s):
    pg = s.page("?mol=m-CN")
    s.badges(pg, [(1, "#exW", "tl", -4, -4, None), (2, "#exDpi", "tl", -4, -4, None), (3, "#dlsvg", "tl", -4, -4, None),
                  (4, "#dlpng", "tl", -4, -4, None), (5, "#copylink", "tl", -4, -4, None), (6, "#savesession", "tl", -4, -4, None),
                  (7, "#copypy", "tl", -4, -4, None)])
    s.save(pg, "11_export_row", clip=s.region(pg, ".exports", pad=24))
    pg.close()


def paste(s):
    pg = s.page("?mol=m-CN")
    pg.evaluate("document.querySelector('details.paste').open = true")
    pg.fill("#paste", "3\nwater\nO   0.000   0.000   0.117\nH   0.000   0.757  -0.470\nH   0.000  -0.757  -0.470")
    pg.wait_for_timeout(300)
    s.badges(pg, [(1, "details.paste summary", "tl", -4, -4, None), (2, "#paste", "tl", -4, -4, None), (3, "#load", "tl", -4, -4, None),
                  (4, ".filebtn", "tl", -4, -4, None)])
    s.save(pg, "12_paste", clip=s.region(pg, ".bar", pad=18))
    pg.close()


def substructure(s):
    pg = s.page("?mol=045_lig")
    s.open_section(pg, "Substructure")
    s.save(pg, "13_substructure", locator="#layers")
    pg.close()


def live_check(s):
    """Served by `descripytor visual`, the page re-checks its five numbers against the Python on every pick."""
    import urllib.request
    port = _free_port()
    root = HERE.parents[2]
    srv = subprocess.Popen([sys.executable, "-m", "descripytor.cli", "visual", "--no-browser", "--port", str(port)],
                           cwd=root, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    try:
        for _ in range(80):
            try:
                urllib.request.urlopen(f"http://127.0.0.1:{port}/status", timeout=0.5)
                break
            except Exception:
                time.sleep(0.5)
        else:
            raise RuntimeError("descripytor visual did not start")
        for mol, name in (("FL_lig_13", "14_live_check"), ("p-OTf", "15_live_check_otf")):
            pg = s.b.new_page(viewport={"width": W, "height": H}, device_scale_factor=DPR)
            pg.goto(f"http://127.0.0.1:{port}/theta?mol={mol}&preset=" + quote("theta construction"), wait_until="domcontentloaded")
            pg.wait_for_selector("#fig circle.hit", state="attached", timeout=20000)
            pg.wait_for_function("/same five values|different/.test(document.querySelector('#status').innerText)", timeout=30000)
            pg.wait_for_timeout(600)
            s.save(pg, name, clip=s.region(pg, "#vals", ".foot", pad=14))
            pg.close()
    finally:
        srv.terminate()


STEPS = {f.__name__: f for f in (overview, toolbar, anatomy, popup, show_chips, b1_scan, side_panel, stack, sheet,
                                 figure_mode, export_row, paste, substructure, live_check)}


def main(names):
    server = subprocess.Popen([sys.executable, "-m", "http.server", str(PORT), "--bind", "127.0.0.1", "--directory", str(EXPLORER)],
                              stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    time.sleep(1.5)
    try:
        with sync_playwright() as pw:
            s = Shots(pw)
            for name in names or STEPS:
                print(name)
                try:
                    STEPS[name](s)
                except Exception as e:                       # keep going: one broken image should not hide the rest
                    print(f"  FAILED {name}: {type(e).__name__}: {str(e)[:220]}")
            s.b.close()
    finally:
        server.terminate()


if __name__ == "__main__":
    main(sys.argv[1:])
