#!/usr/bin/env python3
"""py3Dmol / 3Dmol.js viewer for structures and ORCA normal modes.

Reuse:
  parse_hess          utils/orca_compare/.../vibration_tokenizer.py
  animation frames    utils/orca_compare/.../vibration_py3dmol.py
  Gallarati windows   ts_generator/parse_orca_imag.py
  GUI 3Dmol styling   M2_data_extractor/feature_extraction_gui.html

Usage:
  python utils/mode_viewer.py path/to/TSRE_R.hess --open
  python utils/mode_viewer.py path/to/run_dir --open
  python utils/mode_viewer.py --serve          # drop .hess / .out / .xyz in the browser

A legal TS is exactly one imaginary mode (optionally inside the TSRC/TSRE window).
"""
from __future__ import annotations

import argparse
import json
import math
import os
import re
import sys
import tempfile
import webbrowser
from dataclasses import dataclass, field
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Optional
from urllib.parse import urlparse

import numpy as np

HERE = Path(__file__).resolve().parent
DATASET = HERE / "orca_compare" / "validation_outputs" / "dataset"
HTML_PATH = HERE / "mode_viewer.html"
THREEDMOL = DATASET / "vendor" / "3Dmol-min.js"
if str(DATASET) not in sys.path:
    sys.path.insert(0, str(DATASET))

from vibration_tokenizer import HessData, parse_hess  # noqa: E402

WINDOWS = {
    "TSRC": (-110.0, -25.0),
    "TSRE": (-230.0, -35.0),
    "any": None,
}
IMAG_CUT = -1.0
BOHR = 0.529177210903


@dataclass
class ViewMol:
    name: str
    source: str
    symbols: list
    coordinates: np.ndarray
    frequencies: np.ndarray = field(default_factory=lambda: np.zeros(0))
    modes: Optional[np.ndarray] = None  # (3N, n_modes)
    intensities: Optional[np.ndarray] = None
    trajectory_xyz: str = ""
    note: str = ""


def parse_xyz_text(text: str, name: str = "xyz") -> ViewMol:
    lines = [ln.rstrip("\n") for ln in text.splitlines()]
    frames = []
    i = 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        try:
            n = int(lines[i].split()[0])
        except (ValueError, IndexError) as exc:
            raise ValueError(f"{name}: not XYZ (line {i + 1})") from exc
        comment = lines[i + 1] if i + 1 < len(lines) else ""
        body = lines[i + 2:i + 2 + n]
        if len(body) < n:
            raise ValueError(f"{name}: expected {n} atoms, found {len(body)}")
        atoms = []
        for row in body:
            p = row.split()
            atoms.append((p[0], float(p[1]), float(p[2]), float(p[3])))
        frames.append((comment, atoms))
        i += 2 + n
    if not frames:
        raise ValueError(f"{name}: empty XYZ")
    symbols = [a[0] for a in frames[0][1]]
    coords = np.array([a[1:] for a in frames[0][1]], dtype=float)
    traj = ""
    if len(frames) > 1:
        blocks = []
        for comment, atoms in frames:
            block = [str(len(atoms)), comment]
            block.extend(f"{el:<2s} {x: .8f} {y: .8f} {z: .8f}" for el, x, y, z in atoms)
            blocks.append("\n".join(block))
        traj = "\n".join(blocks) + "\n"
    return ViewMol(
        name=name, source=name, symbols=symbols, coordinates=coords,
        trajectory_xyz=traj,
        note=f"{len(frames)} XYZ frame(s)",
    )


def _hess_to_view(data: HessData) -> ViewMol:
    return ViewMol(
        name=data.molecule_id,
        source=data.source,
        symbols=list(data.symbols),
        coordinates=np.asarray(data.coordinates, dtype=float),
        frequencies=np.asarray(data.frequencies, dtype=float),
        modes=np.asarray(data.modes, dtype=float),
        intensities=np.asarray(data.intensities, dtype=float),
        note="ORCA .hess (analytic/numeric Hessian + normal modes)",
    )


def _last_header(lines: list[str], title: str) -> int:
    hit = -1
    t = title.lower()
    for i, line in enumerate(lines):
        if t in line.lower():
            hit = i
    return hit


def _orca_last_xyz(lines: list[str]) -> tuple[list[str], np.ndarray]:
    idx = _last_header(lines, "CARTESIAN COORDINATES (ANGSTROEM)")
    if idx < 0:
        idx = _last_header(lines, "CARTESIAN COORDINATES (A.U.)")
        au = True
    else:
        au = "A.U." in lines[idx].upper()
    if idx < 0:
        raise ValueError("No CARTESIAN COORDINATES block in ORCA output")
    symbols, coords = [], []
    for line in lines[idx + 1:]:
        s = line.strip()
        if not s or s.startswith("-"):
            if symbols:
                break
            continue
        if s.lower().startswith("cartesian") or s.lower().startswith("internal"):
            break
        parts = s.split()
        if len(parts) < 4:
            if symbols:
                break
            continue
        try:
            xyz = [float(parts[1]), float(parts[2]), float(parts[3])]
        except ValueError:
            if symbols:
                break
            continue
        symbols.append(parts[0])
        coords.append(xyz)
    if not symbols:
        raise ValueError("Empty coordinate block in ORCA output")
    arr = np.asarray(coords, dtype=float)
    if au:
        arr = arr * BOHR
    return symbols, arr


def _orca_frequencies(lines: list[str]) -> tuple[np.ndarray, list[bool]]:
    idx = _last_header(lines, "VIBRATIONAL FREQUENCIES")
    if idx < 0:
        return np.zeros(0), []
    freqs, flags = [], []
    pat = re.compile(r"^\s*(\d+):\s+(-?\d+\.\d+)\s+cm")
    for line in lines[idx + 1:]:
        s = line.strip()
        if s.startswith("NORMAL") or s.startswith("IR SPECTRUM") or s.startswith("RAMAN"):
            break
        m = pat.match(line)
        if not m:
            continue
        v = float(m.group(2))
        imag = "imaginary" in line.lower()
        if imag and v > 0:
            v = -v
        freqs.append(v)
        flags.append(imag or v < IMAG_CUT)
    return np.asarray(freqs, dtype=float), flags


def _orca_normal_modes(lines: list[str], n_cart: int, n_modes: int) -> Optional[np.ndarray]:
    idx = _last_header(lines, "NORMAL MODES")
    if idx < 0 or n_modes <= 0:
        return None
    matrix = np.zeros((n_cart, n_modes), dtype=float)
    columns: list[int] = []
    started = False
    for line in lines[idx + 1:]:
        s = line.strip()
        if not s:
            continue
        if s.startswith("IR SPECTRUM") or s.startswith("RAMAN") or s.startswith("VIBRATIONAL"):
            break
        if s.startswith("-") or s.lower().startswith("these values"):
            continue
        parts = s.split()
        if not parts:
            continue
        header = all("." not in p and p.lstrip("+-").isdigit() for p in parts)
        if header:
            columns = [int(p) for p in parts]
            started = True
            continue
        if not started or not columns:
            continue
        if not parts[0].lstrip("+-").isdigit():
            break
        row = int(parts[0])
        vals = [float(p) for p in parts[1:]]
        if row >= n_cart or len(vals) != len(columns):
            continue
        for col, val in zip(columns, vals):
            if 0 <= col < n_modes:
                matrix[row, col] = val
    if not started:
        return None
    return matrix


def _orca_ir(lines: list[str], n_modes: int) -> np.ndarray:
    idx = _last_header(lines, "IR SPECTRUM")
    out = np.full(n_modes, math.nan)
    if idx < 0 or n_modes <= 0:
        return out
    pat = re.compile(r"^\s*(\d+):\s+")
    for line in lines[idx + 1:]:
        s = line.strip()
        if not s:
            if np.isfinite(out).any():
                break
            continue
        if s.startswith("RAMAN") or s.startswith("NORMAL") or s.startswith("THE"):
            break
        m = pat.match(line)
        if not m:
            continue
        parts = s.replace(":", " ").split()
        try:
            k = int(parts[0])
            # ORCA: Mode freq eps Int ...  Int is km/mol, typically 4th number
            nums = [float(p) for p in parts[1:] if re.match(r"^-?\d", p)]
            if k < n_modes and len(nums) >= 3:
                out[k] = nums[2]
            elif k < n_modes and len(nums) >= 1:
                out[k] = nums[-1]
        except (ValueError, IndexError):
            continue
    return out


def parse_orca_out_text(text: str, name: str = "orca") -> ViewMol:
    lines = text.splitlines()
    symbols, coords = _orca_last_xyz(lines)
    freqs, _flags = _orca_frequencies(lines)
    n_cart = 3 * len(symbols)
    n_modes = int(len(freqs)) if len(freqs) else n_cart
    modes = _orca_normal_modes(lines, n_cart, n_modes)
    intensities = _orca_ir(lines, n_modes) if n_modes else None
    note = "ORCA output"
    if freqs.size:
        n_imag = int((freqs < IMAG_CUT).sum())
        note += f" · {n_modes} frequencies · {n_imag} imaginary"
        if modes is None:
            note += " · no NORMAL MODES block (cannot animate)"
    else:
        note += " · no VIBRATIONAL FREQUENCIES (OptTS did not print a Hessian)"
    return ViewMol(
        name=name, source=name, symbols=symbols, coordinates=coords,
        frequencies=freqs, modes=modes, intensities=intensities, note=note,
    )


def parse_text(filename: str, text: str) -> ViewMol:
    lower = filename.lower()
    stem = Path(filename).stem
    if lower.endswith(".hess") or "$orca_hessian_file" in text[:400].lower() or "$vibrational_frequencies" in text[:8000].lower():
        tmp = Path(tempfile.mkdtemp()) / (stem + ".hess")
        tmp.write_text(text, encoding="utf-8")
        try:
            return _hess_to_view(parse_hess(tmp))
        finally:
            try:
                tmp.unlink()
                tmp.parent.rmdir()
            except OSError:
                pass
    if lower.endswith(".xyz"):
        return parse_xyz_text(text, name=stem)
    if "O   R   C   A" in text or "CARTESIAN COORDINATES" in text or lower.endswith((".out", ".log", ".optts.out")):
        return parse_orca_out_text(text, name=stem)
    if text.lstrip()[:1].isdigit():
        return parse_xyz_text(text, name=stem)
    raise ValueError(f"Unrecognized file: {filename}")


def parse_path(path: Path) -> ViewMol:
    path = path.resolve()
    text = path.read_text(encoding="utf-8", errors="replace")
    mol = parse_text(path.name, text)
    mol.source = str(path)
    if mol.name in ("xyz", "orca"):
        mol.name = path.stem
    return mol


def collect_paths(inputs: list[Path]) -> list[Path]:
    out: list[Path] = []
    suffixes = {".hess", ".xyz", ".out", ".log"}
    for item in inputs:
        if item.is_dir():
            found = []
            for p in sorted(item.rglob("*")):
                if not p.is_file():
                    continue
                if p.suffix.lower() in suffixes or p.name.lower().endswith(".optts.out"):
                    found.append(p)
            # Prefer hess over out over xyz for the same stem
            by_stem: dict[str, Path] = {}
            rank = {".hess": 0, ".out": 1, ".log": 2, ".xyz": 3}
            for p in found:
                suf = ".out" if p.name.lower().endswith(".out") else p.suffix.lower()
                prev = by_stem.get(p.stem)
                if prev is None or rank.get(suf, 9) < rank.get(prev.suffix.lower(), 9):
                    by_stem[p.stem] = p
            out.extend(by_stem[k] for k in sorted(by_stem))
        elif item.is_file():
            out.append(item)
        else:
            raise FileNotFoundError(item)
    return out


def _normalized_mode(modes: np.ndarray, k: int) -> np.ndarray:
    disp = modes[:, k].reshape(-1, 3).copy()
    largest = float(np.linalg.norm(disp, axis=1).max())
    if largest > 1e-15:
        disp = disp / largest
    return disp


def default_mode_index(frequencies: np.ndarray) -> int:
    if frequencies.size == 0:
        return 0
    imag = np.where(frequencies < IMAG_CUT)[0]
    if len(imag):
        return int(imag[0])
    vib = np.where(np.abs(frequencies) > 10.0)[0]
    return int(vib[0]) if len(vib) else 0


def window_label(freq: float) -> str:
    if freq >= IMAG_CUT:
        return ""
    for name, win in WINDOWS.items():
        if win is None:
            continue
        lo, hi = win
        if lo <= freq <= hi:
            return name
    return "out-of-window"


def molecule_payload(mol: ViewMol, max_modes: int = 400) -> dict:
    n = len(mol.symbols)
    freqs = mol.frequencies if mol.frequencies is not None else np.zeros(0)
    n_modes = int(freqs.size)
    modes_out = []
    if mol.modes is not None and n_modes:
        take = min(n_modes, max_modes)
        for k in range(take):
            freq = float(freqs[k])
            imag = freq < IMAG_CUT
            disp = _normalized_mode(mol.modes, k)
            amps = np.linalg.norm(disp, axis=1)
            order = np.argsort(-amps)[:8]
            top = [
                {"i": int(i) + 1, "elem": mol.symbols[i], "amp": round(float(amps[i]), 3)}
                for i in order if amps[i] > 0.05
            ]
            ir = None
            if mol.intensities is not None and k < len(mol.intensities) and np.isfinite(mol.intensities[k]):
                ir = round(float(mol.intensities[k]), 3)
            modes_out.append({
                "n": k,
                "f": round(freq, 2),
                "imag": imag,
                "window": window_label(freq),
                "ir": ir,
                "v": np.round(disp, 5).tolist(),
                "top": top,
            })
    n_imag = sum(1 for m in modes_out if m["imag"]) if modes_out else int((freqs < IMAG_CUT).sum()) if n_modes else 0
    if modes_out:
        if n_imag == 1:
            w = next(m["window"] for m in modes_out if m["imag"])
            verdict = f"one imaginary mode ({w or 'no Gallarati window'})"
            ok = True
        elif n_imag == 0:
            verdict = "no imaginary mode — minimum, not a TS"
            ok = False
        else:
            verdict = f"{n_imag} imaginary modes — not a first-order saddle"
            ok = False
    elif n_modes:
        verdict = f"{n_modes} frequencies listed, no displacement vectors"
        ok = n_imag == 1
    else:
        verdict = "geometry only — drop a .hess or NUMFREQ .out to see vibrations"
        ok = False
    return {
        "name": mol.name,
        "source": mol.source,
        "note": mol.note,
        "symbols": mol.symbols,
        "coords": np.round(mol.coordinates, 6).tolist(),
        "n_atoms": n,
        "n_imag": n_imag,
        "ok": ok,
        "verdict": verdict,
        "default_mode": default_mode_index(freqs) if n_modes else None,
        "modes": modes_out,
        "traj": mol.trajectory_xyz,
    }


def xyz_of(mol: ViewMol) -> str:
    lines = [str(len(mol.symbols)), mol.name]
    for s, (x, y, z) in zip(mol.symbols, mol.coordinates):
        lines.append(f" {s:<2s} {x:12.8f} {y:12.8f} {z:12.8f}")
    return "\n".join(lines) + "\n"


def py3dmol_view(mol: ViewMol, mode: Optional[int] = None, amplitude: float = 0.45,
                 frames: int = 24, interval: int = 70, arrows: bool = True,
                 width: int = 720, height: int = 520):
    """Notebook helper: animated py3Dmol view of one mode."""
    import py3Dmol
    from vibration_py3dmol import build_view

    if mol.modes is None or mol.frequencies.size == 0:
        view = py3Dmol.view(width=width, height=height)
        view.addModel(xyz_of(mol), "xyz")
        view.setStyle({}, {"stick": {"radius": 0.12}, "sphere": {"scale": 0.25}})
        view.zoomTo()
        return view
    k = default_mode_index(mol.frequencies) if mode is None else mode
    data = HessData(
        molecule_id=mol.name, source=mol.source, symbols=mol.symbols,
        masses=np.ones(len(mol.symbols)), coordinates=mol.coordinates,
        frequencies=mol.frequencies, modes=mol.modes,
        intensities=mol.intensities if mol.intensities is not None else np.full(mol.frequencies.shape, np.nan),
    )
    view = build_view(data, k, amplitude, frames, interval, arrows)
    view._width = width
    view._height = height
    return view


def render_html(payloads: list[dict]) -> str:
    template = HTML_PATH.read_text(encoding="utf-8")
    blob = json.dumps({"molecules": payloads}, separators=(",", ":"))
    marker = '<script type="application/json" id="boot">__DATA__</script>'
    if marker not in template:
        raise RuntimeError("mode_viewer.html missing boot placeholder")
    return template.replace(marker, f'<script type="application/json" id="boot">{blob}</script>', 1)


def write_html(payloads: list[dict], dest: Path) -> Path:
    dest.write_text(render_html(payloads), encoding="utf-8")
    return dest


class _Handler(BaseHTTPRequestHandler):
    def log_message(self, fmt, *args):
        sys.stderr.write("mode_viewer: " + (fmt % args) + "\n")

    def _send(self, code: int, body: bytes, ctype: str):
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):
        path = urlparse(self.path).path
        if path in ("/", "/modes", "/mode_viewer.html"):
            html = HTML_PATH.read_text(encoding="utf-8")
            marker = '<script type="application/json" id="boot">__DATA__</script>'
            html = html.replace(marker, '<script type="application/json" id="boot">null</script>', 1)
            self._send(200, html.encode("utf-8"), "text/html; charset=utf-8")
            return
        if path.endswith("3Dmol-min.js") and THREEDMOL.is_file():
            self._send(200, THREEDMOL.read_bytes(), "application/javascript")
            return
        self._send(404, b"not found", "text/plain")

    def do_POST(self):
        path = urlparse(self.path).path
        if path not in ("/api/parse-modes", "/parse"):
            self._send(404, b"not found", "text/plain")
            return
        n = int(self.headers.get("Content-Length") or 0)
        raw = self.rfile.read(n)
        try:
            msg = json.loads(raw.decode("utf-8"))
            mol = parse_text(msg.get("name") or "upload", msg.get("text") or "")
            payload = molecule_payload(mol)
            body = json.dumps(payload).encode("utf-8")
            self._send(200, body, "application/json")
        except Exception as exc:
            body = json.dumps({"error": str(exc)}).encode("utf-8")
            self._send(400, body, "application/json")


def serve(port: int = 7433, host: str = "127.0.0.1"):
    httpd = ThreadingHTTPServer((host, port), _Handler)
    url = f"http://{host}:{port}/"
    print(f"Mode viewer  {url}", flush=True)
    print("Drop a .hess, ORCA .out / .optts.out, or .xyz", flush=True)
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print("\nstopped")


def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("inputs", nargs="*", type=Path, help=".hess / .out / .xyz / directory")
    ap.add_argument("--open", action="store_true", help="write HTML and open the browser")
    ap.add_argument("--serve", action="store_true", help="local drop-file server")
    ap.add_argument("--port", type=int, default=int(os.environ.get("MODE_VIEWER_PORT", "7433")))
    ap.add_argument("--output", type=Path, default=None)
    ap.add_argument("--window", choices=sorted(WINDOWS), default="any")
    args = ap.parse_args(argv)

    if args.serve and not args.inputs:
        serve(port=args.port)
        return 0

    if not args.inputs:
        ap.print_help()
        print("\nNo files given — starting drop server. Ctrl+C to stop.")
        serve(port=args.port)
        return 0

    paths = collect_paths(args.inputs)
    if not paths:
        raise SystemExit("No .hess / .out / .xyz files found")
    payloads = []
    for p in paths:
        mol = parse_path(p)
        payload = molecule_payload(mol)
        if args.window != "any" and payload["modes"]:
            lo, hi = WINDOWS[args.window]
            for m in payload["modes"]:
                if m["imag"]:
                    m["in_window"] = lo <= m["f"] <= hi
        payloads.append(payload)
        print(f"{payload['name']}\t{payload['verdict']}\tn_imag={payload['n_imag']}\t{p}")

    dest = args.output or (Path(tempfile.gettempdir()) / "mode_viewer.html")
    write_html(payloads, dest)
    print("wrote", dest.resolve())
    if args.open or args.serve:
        webbrowser.open(dest.resolve().as_uri())
    if args.serve:
        serve(port=args.port)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
