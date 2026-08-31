"""Record a longer picker demo: cmd -> example feathers -> click atoms -> CSV -> model."""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = Path(__file__).resolve().parent
FEATHER = ROOT / "descripytor" / "examples" / "feather_example" / "basic.feather"
GUI = os.environ.get("DESCRIPYTOR_GUI", "http://127.0.0.1:7433/visual")
W, H = 1280, 720


def ffmpeg_exe() -> str:
    import imageio_ffmpeg

    return imageio_ffmpeg.get_ffmpeg_exe()


def font(size: int):
    for name in (
        r"C:\Windows\Fonts\cascadiamono.ttf",
        r"C:\Windows\Fonts\CascadiaMono.ttf",
        r"C:\Windows\Fonts\consola.ttf",
        r"C:\Windows\Fonts\lucon.ttf",
    ):
        if Path(name).is_file():
            return ImageFont.truetype(name, size)
    return ImageFont.load_default()


def draw_terminal(lines, cursor=None):
    img = Image.new("RGB", (W, H), (12, 12, 14))
    d = ImageDraw.Draw(img)
    d.rectangle((0, 0, W, 36), fill=(40, 40, 44))
    d.ellipse((16, 12, 28, 24), fill=(255, 95, 86))
    d.ellipse((36, 12, 48, 24), fill=(255, 189, 46))
    d.ellipse((56, 12, 68, 24), fill=(39, 201, 63))
    d.text((88, 10), "Command Prompt  —  DescriPyTor", fill=(200, 200, 204), font=font(16))
    y = 56
    f = font(22)
    for line in lines:
        color = (80, 200, 120) if line.startswith("C:\\") or line.startswith("PS ") else (220, 220, 220)
        if line.strip().startswith("DescriPyTor") or line.strip().startswith("Open:"):
            color = (90, 200, 255)
        d.text((28, y), line, fill=color, font=f)
        y += 30
        if y > H - 40:
            break
    if cursor is not None:
        d.rectangle((28 + cursor[0], 56 + cursor[1] * 30, 40 + cursor[0], 78 + cursor[1] * 30), fill=(200, 200, 200))
    return img


def write_terminal_clip(path: Path) -> None:
    prompt = r"C:\Users\edens\Documents\GitHub\DescriPyTor_to_upload> "
    cmd = "descripytor visual"
    frames_dir = path.parent / "_term_frames"
    frames_dir.mkdir(parents=True, exist_ok=True)
    n = 0

    def save(im, hold=2):
        nonlocal n
        for _ in range(hold):
            im.save(frames_dir / f"{n:04d}.png")
            n += 1

    save(draw_terminal([prompt], cursor=(len(prompt) * 11, 0)), 8)
    typed = ""
    for ch in cmd:
        typed += ch
        save(draw_terminal([prompt + typed], cursor=((len(prompt) + len(typed)) * 11, 0)), 2)
    save(draw_terminal([prompt + cmd], cursor=None), 6)
    banner = [
        prompt + cmd,
        "",
        "  DescriPyTor GUI",
        "  Open:  http://127.0.0.1:7432/visual",
        "  Forms: http://127.0.0.1:7432/forms",
        "  Press Ctrl+C to stop",
        "",
        " * Serving Flask app 'M2_data_extractor.gui_server'",
        " * Running on http://127.0.0.1:7432",
    ]
    shown = []
    for line in banner:
        shown.append(line)
        save(draw_terminal(shown), 3)
    save(draw_terminal(shown), 50)

    pattern = str(frames_dir / "%04d.png").replace("\\", "/")
    subprocess.check_call(
        [
            ffmpeg_exe(), "-y", "-framerate", "8", "-i", pattern,
            "-r", "25", "-c:v", "libx264", "-pix_fmt", "yuv420p", str(path),
        ],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    shutil.rmtree(frames_dir, ignore_errors=True)


def atom_click_xy(page, serial_1based: int):
    try:
        return page.evaluate(
            """(serial) => {
              const v = window.pickerViewer;
              if (!v || typeof v.getModel !== 'function') return null;
              const model = v.getModel();
              const atoms = (model && typeof model.selectedAtoms === 'function')
                ? (model.selectedAtoms({}) || [])
                : ((model && model.atoms) || []);
              const a = (atoms || []).find(x => x && (x.serial + 1) === serial);
              if (!a || typeof v.modelToScreen !== 'function') return null;
              const s = v.modelToScreen({x: a.x, y: a.y, z: a.z});
              const canvas = document.querySelector('#pickerStage canvas');
              if (!canvas || !s) return null;
              const r = canvas.getBoundingClientRect();
              return { x: r.left + s.x, y: r.top + s.y };
            }""",
            serial_1based,
        )
    except Exception:
        return None


def record_browser(path: Path) -> None:
    from playwright.sync_api import sync_playwright

    if not FEATHER.is_file():
        raise SystemExit(f"missing {FEATHER}")

    with sync_playwright() as p:
        browser = p.chromium.launch(
            headless=True,
            slow_mo=280,
            args=["--use-gl=angle", "--use-angle=swiftshader", "--ignore-gpu-blocklist"],
        )
        context = browser.new_context(
            viewport={"width": W, "height": H},
            record_video_dir=str(path.parent / "_pw"),
            record_video_size={"width": W, "height": H},
        )
        page = context.new_page()
        page.set_default_timeout(240_000)
        video = page.video
        try:
            page.goto(GUI, wait_until="domcontentloaded")
            page.wait_for_timeout(4500)

            page.locator("#useExampleFeathersBtn").click()
            page.wait_for_timeout(3200)

            page.set_input_files("#fileInput", str(FEATHER))
            page.wait_for_function(
                "() => (document.getElementById('molInfoBadge')||{}).innerText && document.getElementById('molInfoBadge').innerText.toLowerCase().includes('basic')"
            )
            page.wait_for_timeout(4000)

            page.locator("#fieldSelect").select_option("sterimol")
            page.wait_for_timeout(1500)
            for serial in (1, 4):
                xy = atom_click_xy(page, serial)
                if xy and xy.get("x"):
                    page.mouse.move(xy["x"], xy["y"], steps=12)
                    page.wait_for_timeout(350)
                    page.mouse.click(xy["x"], xy["y"])
                    page.wait_for_timeout(400)
                else:
                    page.fill("#typedAtoms", str(serial))
                    page.locator("#addTypedAtomsBtn").click()
                page.wait_for_timeout(2200)

            page.locator("#extractBtn").evaluate("el => el.scrollIntoView({block:'center'})")
            page.wait_for_timeout(2500)
            page.locator("#extractBtn").click()
            page.wait_for_function(
                "() => (document.getElementById('extractStatus')||{}).innerText && document.getElementById('extractStatus').innerText.includes('molecules')",
                timeout=240_000,
            )
            page.wait_for_timeout(5500)

            page.locator('button[data-tab="model"]').click()
            page.wait_for_timeout(2800)
            n = page.evaluate(
                """() => {
                  const t = document.getElementById('extractTable');
                  return t ? Math.max(0, t.querySelectorAll('tbody tr').length) : 0;
                }"""
            )
            n = n or 10
            outputs = "\n".join(f"{0.15 + i * 0.07:.2f}" for i in range(n))
            page.fill("#lrOutputs", outputs)
            page.fill("#lrThr", "0")
            page.wait_for_timeout(2200)
            page.locator("#lrRunBtn").click()
            page.wait_for_function(
                """() => {
                  const el = document.getElementById('lrStatus');
                  return el && /models for target/i.test(el.textContent || '');
                }""",
                timeout=60_000,
            )
            page.wait_for_timeout(800)
            page.evaluate(
                """() => {
                  const el = document.getElementById('lrTableWrap');
                  if (el && el.style.display !== 'none') el.scrollIntoView({block:'center'});
                }"""
            )
            page.wait_for_timeout(7000)
        finally:
            context.close()
            browser.close()
            src = Path(video.path()) if video else None
            if src and src.is_file():
                shutil.move(str(src), str(path))
            shutil.rmtree(path.parent / "_pw", ignore_errors=True)


def concat(clips, dest: Path) -> None:
    lst = dest.parent / "_concat.txt"
    lines = []
    for c in clips:
        lines.append(f"file '{c.resolve().as_posix()}'")
    lst.write_text("\n".join(lines), encoding="utf-8")
    subprocess.check_call(
        [
            ffmpeg_exe(), "-y", "-f", "concat", "-safe", "0", "-i", str(lst),
            "-c:v", "libx264", "-pix_fmt", "yuv420p", "-movflags", "+faststart",
            str(dest),
        ],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    lst.unlink(missing_ok=True)


def main() -> None:
    term = OUT_DIR / "_term.mp4"
    gui = OUT_DIR / "_gui.webm"
    dest = OUT_DIR / "picker-full-walkthrough.mp4"
    print("1/3 terminal clip")
    write_terminal_clip(term)
    print("2/3 browser recording (extract can take a few minutes)")
    record_browser(gui)
    print("3/3 concatenate")
    gui_mp4 = OUT_DIR / "_gui.mp4"
    subprocess.check_call(
        [
            ffmpeg_exe(), "-y", "-i", str(gui),
            "-r", "25", "-c:v", "libx264", "-pix_fmt", "yuv420p", str(gui_mp4),
        ],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    concat([term, gui_mp4], dest)
    for p in (term, gui, gui_mp4):
        p.unlink(missing_ok=True)
    print("wrote", dest, dest.stat().st_size)


if __name__ == "__main__":
    os.chdir(ROOT)
    main()
