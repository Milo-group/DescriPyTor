"""
gui_server.py — local backend for the browser GUIs.

Preferred:
    descripytor visual

Or:
    python M2_data_extractor/gui_server.py

Then open http://localhost:7432  (redirects to /visual)
"""

import sys
import os
import json
import traceback
import math

# ── Flask ─────────────────────────────────────────────────────
try:
    from flask import Flask, request, jsonify, send_file, redirect, send_from_directory
    from flask_cors import CORS
except ImportError:
    print("Missing dependencies. Run:  pip install flask flask-cors")
    sys.exit(1)

app = Flask(__name__, static_folder=None)
# Every route reads or writes files on this machine, so only pages served from this machine may call it. Other ports here
# are fine (the feature explorer can be served by a plain http.server and call this one); any other site is refused.
_LOCAL_HOSTS = {"127.0.0.1", "localhost", "::1"}
CORS(app, origins=[r"^https?://(127\.0\.0\.1|localhost|\[::1\])(:\d+)?$"])


@app.before_request
def _local_pages_only():
    """CORS only stops a foreign page reading the answer; this stops the request. A Host that is not this machine's
    name is DNS rebinding; an Origin from elsewhere is another site's script (Origin "null" too: sandboxed or file://);
    a cross-site request without Origin is an <img> or <form> (opening a page by a link is still allowed)."""
    from urllib.parse import urlsplit

    origin = request.headers.get("Origin")
    refused = (
        urlsplit("//" + request.host).hostname not in _LOCAL_HOSTS
        or (origin is not None and urlsplit(origin).hostname not in _LOCAL_HOSTS)
        or (origin is None and request.headers.get("Sec-Fetch-Site") == "cross-site"
            and request.headers.get("Sec-Fetch-Mode") != "navigate")
    )
    if refused:
        return jsonify({"error": "This server only answers pages opened from this computer."}), 403

PORT = int(os.environ.get("GUI_PORT", "7432"))

# ── path helper ───────────────────────────────────────────────
def ensure_path(root: str):
    """Add DescriPyTor root to sys.path so imports work."""
    if root and root not in sys.path:
        sys.path.insert(0, root)
    # also try the directory containing this file as fallback
    here = os.path.dirname(os.path.abspath(__file__))
    parent = os.path.dirname(here)
    for p in (here, parent):
        if p not in sys.path:
            sys.path.insert(0, p)


# ── routes ────────────────────────────────────────────────────

def benzene_example_dir():
    """GUI default example: bundled substituted benzenes (basic.feather)."""
    try:
        from descripytor.examples import feather_example_dir as _feather_dir

        path = _feather_dir()
        if path.is_dir():
            return str(path)
    except Exception:
        pass
    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.dirname(here)
    for candidate in (
        os.path.join(root, "descripytor", "examples", "feather_example"),
        os.path.join(root, "Getting_started_with_examples", "feather_example"),
    ):
        if os.path.isdir(candidate) and os.path.isfile(os.path.join(candidate, "basic.feather")):
            return candidate
    return ""


def example_feather_dir():
    """GUI default example: the bundled substituted benzenes."""
    return benzene_example_dir()


def example_reference_feather():
    """3D reference molecule: basic.feather, else unsub.feather in the example folder."""
    directory = example_feather_dir()
    if not directory:
        return ""
    for name in ("basic.feather", "unsub.feather"):
        path = os.path.join(directory, name)
        if os.path.isfile(path):
            return path
    files = _list_feather_files(directory)
    return files[0] if files else ""


def example_basic_feather():
    """Back-compat alias for the example 3D reference molecule."""
    return example_reference_feather()


def _read_example_json(filename, default=None):
    directory = example_feather_dir()
    if not directory:
        return default
    path = os.path.join(directory, filename)
    if not os.path.isfile(path):
        return default
    try:
        with open(path, encoding="utf-8") as handle:
            return json.load(handle)
    except Exception:
        return default


def example_presets():
    data = _read_example_json("presets.json") or {}
    atoms = data.get("atoms") if isinstance(data, dict) else {}
    return atoms if isinstance(atoms, dict) else {}


def _norm_mol_name(name):
    text = str(name or "").strip()
    lower = text.lower()
    if lower.endswith(".feather"):
        text = text[:-8]
    elif lower.endswith(".ftr"):
        text = text[:-4]
    return text.strip().lower().replace("_", "-")


def _read_outcomes_csv(path):
    import csv

    rows = []
    if not path or not os.path.isfile(path):
        return {"rows": [], "y_column": "output", "name_column": "name", "path": path or ""}
    with open(path, newline="", encoding="utf-8-sig") as handle:      # Excel's "CSV UTF-8" starts with a BOM
        reader = csv.DictReader(handle)
        for row in reader:
            name = (row.get("name") or "").strip()
            raw = row.get("output")
            if not name:
                continue
            try:
                value = float(raw)
            except (TypeError, ValueError):
                continue
            rows.append({"name": name, "output": value})
    return {"rows": rows, "y_column": "output", "name_column": "name", "path": path}


def example_outcomes_payload():
    """Measured outputs of the example set the page starts on: the outcomes.csv in that folder, if it has one.
    Never another set's: the Baptiste outcomes used to be offered for the benzene set, and joined onto any table
    whose compound names happened to match (p-Cl, p-NO2, ...)."""
    empty = {"rows": [], "y_column": "output", "name_column": "name"}
    example_dir = example_feather_dir()
    if not example_dir:
        return empty
    payload = _read_outcomes_csv(os.path.join(example_dir, "outcomes.csv"))
    return payload if payload.get("rows") else empty


def _join_outcomes_into_df(df, rows, y_column="output"):
    """Add y_column from name/output rows. Returns (df, n_matched)."""
    if df is None or getattr(df, "empty", True) or not rows:
        return df, 0
    name_col = next(
        (
            c
            for c in df.columns
            if str(c).lower() in (
                "name", "names", "molecule", "molecule_name", "mol", "compound",
            )
        ),
        None,
    )
    by_norm = {}
    for row in rows:
        key = _norm_mol_name(row.get("name"))
        if key:
            by_norm[key] = row.get("output")
    # the extractor's table carries the compound names as its index, not as a column
    names = df[name_col] if name_col is not None else df.index.to_series()
    mapped = names.map(lambda value: by_norm.get(_norm_mol_name(value)))
    n_matched = int(mapped.notna().sum())
    if n_matched <= 0:
        return df, 0
    out = df.copy()
    out[y_column] = mapped.to_numpy()
    return out, n_matched


def _is_placeholder_path(path):
    if not path or not str(path).strip():
        return True
    text = str(path).replace("/", "\\").lower()
    if "path\\to\\your" in text or "path\\to\\logs" in text:
        return True
    stripped = str(path).strip()
    return stripped.startswith("__") and stripped.endswith("__")


def _load_descriptor_extractor():
    here = os.path.dirname(os.path.abspath(__file__))       # the engine is in this folder, so an installed package has it
    if here not in sys.path:
        sys.path.insert(0, here)
    import descriptor_extractor as dx  # noqa: WPS433

    return dx


@app.route("/status")
def status():
    outcomes = example_outcomes_payload()
    return jsonify({
        "ok": True,
        "version": "1.0",
        "api": GUI_API,
        "example_feather_dir": example_feather_dir(),
        "example_basic_feather": example_basic_feather(),
        "example_reference_feather": example_reference_feather(),
        "example_reference_name": os.path.splitext(
            os.path.basename(example_reference_feather() or "")
        )[0],
        "example_presets": example_presets(),
        "example_outcomes": outcomes.get("rows") or [],
        "example_y_column": outcomes.get("y_column") or "output",
        "cpu_count": os.cpu_count() or 1,
        "bassa_available": _module_available("bassa_reg"),
    })


@app.route("/")
def root():
    """First-run URL: the 3D picker."""
    return redirect("/visual")


@app.route("/visual_extractor")
def visual_extractor():
    """The picker under its command's name: the three pages link to each other by /visual_extractor, /feature_explorer, /bassa_explorer."""
    return redirect("/visual")


def _listening(port):
    import socket
    with socket.socket() as s:
        s.settimeout(0.3)
        return s.connect_ex(("127.0.0.1", port)) == 0


def _bassa_running():
    root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    if root_dir not in sys.path:
        sys.path.insert(0, root_dir)
    from M3_modeler.bassa_explorer import running
    return running


def _bassa_url():
    """The BASSA model explorer's address. It has its own server (it may run under another Python): started here if it
    is not up. None when it does not come up, or when the port is held by a server from before an update."""
    import subprocess
    import time
    global BASSA_DOWN
    port = int(os.environ.get("BASSA_EXPLORER_PORT", "8765"))
    BASSA_DOWN = BASSA_NOT_STARTED
    if _bassa_running().state(port) == "old":
        BASSA_DOWN = _bassa_running().stop_hint(port, "the BASSA model explorer")
        return None
    if not _listening(port):
        root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        env = dict(os.environ, PYTHONPATH=root_dir + os.pathsep + os.environ.get("PYTHONPATH", ""))
        subprocess.Popen([sys.executable, "-m", "descripytor.cli", "bassa_explorer", "--no-browser", "--port", str(port)], env=env)
        for _ in range(120):
            if _listening(port):
                break
            time.sleep(0.25)
        else:
            return None
    return f"http://127.0.0.1:{port}/"


BASSA_NOT_STARTED = "The BASSA model explorer did not start; run `descripytor bassa_explorer` to see why."
BASSA_DOWN = BASSA_NOT_STARTED            # what the last failed _bassa_url() has to say


@app.route("/bassa_explorer")
def bassa_explorer():
    """Go to the BASSA model explorer, starting its server if needed."""
    url = _bassa_url()
    return redirect(url) if url else (BASSA_DOWN, 503)


@app.route("/bassa_explorer/handoff", methods=["POST"])
def bassa_handoff():
    """One click from the picker: the extracted table with its target, and the folder its structures are in, are left
    with the BASSA explorer; the answer is the address that opens it and starts the run."""
    import json
    import urllib.error
    import urllib.request

    data = request.json or {}
    try:
        df, y_name = _table_with_target(data)
    except ValueError as exc:
        return jsonify({"error": str(exc)}), 400
    url = _bassa_url()
    if not url:
        return jsonify({"error": BASSA_DOWN}), 503
    body = dict(csv=df.to_csv(index=False), response=y_name, structures=str(data.get("structures") or "").strip(),
                name=str(data.get("name") or "extracted"))
    req = urllib.request.Request(url + "handoff", data=json.dumps(body).encode("utf-8"), headers={"Content-Type": "application/json"})
    try:
        token = json.loads(urllib.request.urlopen(req, timeout=30).read())["id"]
    except urllib.error.HTTPError as exc:                    # its own words: too few compounds, a constant response, ...
        try:
            why = json.loads(exc.read())["error"]
        except Exception:
            why = str(exc)
        return jsonify({"error": why}), 400
    except Exception as exc:
        return jsonify({"error": f"The BASSA model explorer did not take the table: {exc}"}), 502
    return jsonify({"url": f"{url}#handoff={token}", "rows": int(df.shape[0]), "columns": int(df.shape[1]) - 2, "response": y_name})


def atom_picker_html():
    """The picker page, from this folder."""
    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.dirname(here)
    candidates = (
        os.path.join(here, "atom_picker.html"),
    )
    for path in candidates:
        if os.path.isfile(path):
            return path
    raise FileNotFoundError(
        "atom_picker.html not found. Reinstall descripytor or run from a clone."
    )


_3DMOL_URL = "https://3dmol.org/build/3Dmol-min.js"


def _ensure_3dmol_js(folder):
    """Cache 3Dmol next to the GUI after the first successful fetch."""
    os.makedirs(folder, exist_ok=True)
    dest = os.path.join(folder, "3Dmol-min.js")
    if os.path.isfile(dest) and os.path.getsize(dest) > 1000:
        return dest
    try:
        import urllib.request

        req = urllib.request.Request(_3DMOL_URL, headers={"User-Agent": "descripytor"})
        with urllib.request.urlopen(req, timeout=20) as response:
            data = response.read()
        if not data or len(data) < 1000:
            return None
        with open(dest, "wb") as handle:
            handle.write(data)
        return dest
    except Exception:
        return None


@app.route("/visual", strict_slashes=False)
def visual_gui():
    """Serve the combined atom-picker, extraction, and modeling workflow."""
    path = atom_picker_html()
    response = send_file(path, mimetype="text/html")
    response.headers["Cache-Control"] = "no-cache"       # revalidate: after an update the browser kept showing the old page
    return response


@app.route("/static/<path:filename>")
def visual_static(filename):
    """Serve bundled viewer assets next to the picker HTML."""
    import threading

    folder = os.path.join(os.path.dirname(os.path.abspath(__file__)), "static")
    local = os.path.join(folder, filename)
    if os.path.isfile(local) and os.path.getsize(local) > 1000:
        return send_from_directory(folder, filename)
    if filename == "3Dmol-min.js":
        threading.Thread(target=_ensure_3dmol_js, args=(folder,), daemon=True).start()
        return redirect(_3DMOL_URL)
    return jsonify({"error": "not found"}), 404


def _xyz_text_from_coords(df, name):
    label = name or "molecule"
    lines = [str(len(df)), str(label)]
    for _, row in df.iterrows():
        lines.append(
            f"{str(row['atom']):2s} {float(row['x']):14.8f} "
            f"{float(row['y']):14.8f} {float(row['z']):14.8f}"
        )
    return "\n".join(lines), str(label), int(len(df))


def _dipole_from_feather_df(data):
    """Gaussian dipole from named columns or the legacy 4–7 positional block."""
    import math

    import pandas as pd

    if data is None or getattr(data, "empty", True):
        return None
    cols = [str(c).strip() for c in data.columns]
    frame = data.copy()
    frame.columns = cols
    dx = dy = dz = total = None
    named_xyz = all(c in frame.columns for c in ("atom", "x", "y", "z"))
    if all(c in frame.columns for c in ("dip_x", "dip_y", "dip_z")):
        xyz_names = ("dip_x", "dip_y", "dip_z")
    elif all(c in frame.columns for c in ("dipole_x", "dipole_y", "dipole_z")):
        xyz_names = ("dipole_x", "dipole_y", "dipole_z")
    else:
        xyz_names = None
    if xyz_names:
        block = frame[list(xyz_names)].apply(pd.to_numeric, errors="coerce").dropna()
        if not block.empty:
            row = block.iloc[0]
            dx, dy, dz = float(row[xyz_names[0]]), float(row[xyz_names[1]]), float(row[xyz_names[2]])
        for total_name in ("total_dipole", "total"):
            if total_name in frame.columns:
                totals = pd.to_numeric(frame[total_name], errors="coerce").dropna()
                if len(totals):
                    total = float(totals.iloc[0])
                    break
    elif (not named_xyz) and frame.shape[1] >= 8:
        block = frame.iloc[:, 4:8].apply(pd.to_numeric, errors="coerce").dropna(how="any")
        if not block.empty:
            row = block.iloc[0]
            dx, dy, dz = float(row.iloc[0]), float(row.iloc[1]), float(row.iloc[2])
            total = float(row.iloc[3])
    if dx is None:
        return None
    if total is None:
        total = math.sqrt(dx * dx + dy * dy + dz * dz)
    return {"vector": [dx, dy, dz], "total": total}


def _viewer_payload_from_feather(path, display_name="molecule"):
    """XYZ plus dipole overlay from one feather read."""
    import pandas as pd

    data = pd.read_feather(path)
    data.columns = [str(c).strip() for c in data.columns]
    needed = ["atom", "x", "y", "z"]
    if all(c in data.columns for c in needed):
        df = data[needed].copy()
    else:
        xyz = data.iloc[:, 0:4].copy()
        xyz.columns = needed
        df = xyz
    df["x"] = pd.to_numeric(df["x"], errors="coerce")
    df["y"] = pd.to_numeric(df["y"], errors="coerce")
    df["z"] = pd.to_numeric(df["z"], errors="coerce")
    df = df.dropna(subset=["atom", "x", "y", "z"]).reset_index(drop=True)
    if df.empty:
        raise ValueError("No atom coordinates in this feather file")
    label = display_name or os.path.splitext(os.path.basename(path))[0]
    xyz, name, n_atoms = _xyz_text_from_coords(df, label)
    overlay = {"dipole": _dipole_from_feather_df(data), "modes": _modes_from_feather_df(data, n_atoms)}
    return xyz, name, n_atoms, overlay


def _modes_from_feather_df(data, n_atoms, max_modes=250):
    """Normal-mode displacements for the vibration overlay (rows = modes)."""
    import numpy as np
    import pandas as pd

    if data is None or getattr(data, "empty", True) or not n_atoms:
        return []
    try:
        n_atoms = int(n_atoms)
        cols = [str(c).strip() for c in data.columns]
        frame = data.copy()
        frame.columns = cols
        named_xyz = all(c in frame.columns for c in ("atom", "x", "y", "z"))
        freq = ir = vectors = None
        if "Frequency" in frame.columns:
            freq = pd.to_numeric(frame["Frequency"], errors="coerce")
            ir = pd.to_numeric(frame["IR"], errors="coerce") if "IR" in frame.columns else None
            start = (
                list(frame.columns).index("IR") + 1
                if "IR" in frame.columns
                else list(frame.columns).index("Frequency") + 1
            )
            if start < frame.shape[1]:
                vectors = frame.iloc[:, start:].apply(pd.to_numeric, errors="coerce")
        elif (not named_xyz) and frame.shape[1] > 15:
            freq = pd.to_numeric(frame.iloc[:, 13], errors="coerce")
            ir = pd.to_numeric(frame.iloc[:, 14], errors="coerce")
            vectors = frame.iloc[:, 15:].apply(pd.to_numeric, errors="coerce")
        if freq is None or vectors is None or vectors.shape[1] < 3:
            return []
        mask = freq.notna()
        freq_v = freq[mask].to_numpy(dtype=float)
        if freq_v.size == 0:
            return []
        ir_v = ir[mask].to_numpy(dtype=float) if ir is not None else None
        arr = vectors.loc[mask].to_numpy(dtype=float)
        n_vib_atoms = arr.shape[1] // 3
        if n_vib_atoms < 1:
            return []
        use = min(n_vib_atoms, n_atoms)
        modes = []
        for i, fval in enumerate(freq_v):
            if i >= max_modes:
                break
            row = arr[i]
            disp = []
            for atom_i in range(n_atoms):
                if atom_i < use:
                    sl = row[atom_i * 3:(atom_i + 1) * 3]
                    if sl.shape[0] == 3 and np.isfinite(sl).all():
                        disp.append([round(float(sl[0]), 4), round(float(sl[1]), 4), round(float(sl[2]), 4)])
                    else:
                        disp.append([0.0, 0.0, 0.0])
                else:
                    disp.append([0.0, 0.0, 0.0])
            ir_val = None
            if ir_v is not None and i < len(ir_v) and np.isfinite(ir_v[i]):
                ir_val = round(float(ir_v[i]), 4)
            modes.append({"freq": round(float(fval), 4), "ir": ir_val, "disp": disp})
        return modes
    except Exception:
        return []


def _mol_data_from_molecule(mol):
    import math

    overlay = {"dipole": None, "modes": []}
    try:
        ddf = getattr(mol, "gauss_dipole_df", None)
        if ddf is not None and len(ddf) > 0:
            row = ddf.iloc[0]
            dx = float(row.get("dip_x", row.iloc[0]))
            dy = float(row.get("dip_y", row.iloc[1]))
            dz = float(row.get("dip_z", row.iloc[2]))
            total = row.get("total_dipole", row.get("total", None))
            if total is None and len(row) > 3:
                total = row.iloc[3]
            overlay["dipole"] = {
                "vector": [dx, dy, dz],
                "total": float(total) if total is not None else math.sqrt(dx * dx + dy * dy + dz * dz),
            }
    except Exception:
        pass
    try:
        vmd = getattr(mol, "vibration_mode_dict", None) or {}
        info = getattr(mol, "info_df", None)
        ir_for = {}
        if info is not None and "Frequency" in getattr(info, "columns", []):
            ir_col = "IR" if "IR" in info.columns else None
            for _, row in info.iterrows():
                try:
                    ir_for[round(float(row["Frequency"]), 2)] = (
                        float(row[ir_col]) if ir_col else None
                    )
                except Exception:
                    pass
        modes = []
        for freq, disp in vmd.items():
            try:
                fval = float(freq)
            except (TypeError, ValueError):
                continue
            try:
                rows = [[round(float(x), 4) for x in vec] for vec in disp]
            except Exception:
                continue
            if not rows or len(rows[0]) != 3:
                continue
            ir_val = ir_for.get(round(fval, 2))
            modes.append({
                "freq": round(fval, 4),
                "ir": None if ir_val is None else round(float(ir_val), 4),
                "disp": rows,
            })
            if len(modes) >= 250:
                break
        modes.sort(key=lambda m: m["freq"])
        overlay["modes"] = modes
    except Exception:
        pass
    return overlay


def _xyz_text_from_molecule(mol, name=None):
    df = mol.xyz_df[["atom", "x", "y", "z"]]
    label = name or getattr(mol, "molecule_name", "molecule")
    return _xyz_text_from_coords(df, label)


def _molecule_view_from_path(path, display_name="molecule"):
    """XYZ, dipole, and vibration overlays for the 3D viewer."""
    try:
        return _viewer_payload_from_feather(path, display_name)
    except Exception:
        pass
    here = os.path.dirname(os.path.abspath(__file__))
    ensure_path(here)
    ensure_path(os.path.dirname(here))
    try:
        from M2_data_extractor.data_extractor import Molecule
    except ImportError:
        from data_extractor import Molecule
    mol = Molecule(path)
    xyz, name, n_atoms = _xyz_text_from_molecule(mol, display_name)
    return xyz, name, n_atoms, _mol_data_from_molecule(mol)


def _molecule_xyz_from_path(path, display_name="molecule"):
    """XYZ for the 3D viewer. Prefer a fast feather read; fall back to Molecule()."""
    xyz, name, n_atoms, _overlay = _molecule_view_from_path(path, display_name)
    return xyz, name, n_atoms


@app.route("/xyz_from_feather", methods=["POST"])
def xyz_from_feather():
    """Turn a uploaded or on-disk .feather into an XYZ block for the 3D viewer."""
    import tempfile

    path = None
    cleanup = False
    display_name = "molecule"
    try:
        uploaded = request.files.get("file")
        if uploaded and uploaded.filename:
            display_name = os.path.splitext(os.path.basename(uploaded.filename))[0]
            suffix = os.path.splitext(uploaded.filename)[1] or ".feather"
            handle = tempfile.NamedTemporaryFile(delete=False, suffix=suffix)
            uploaded.save(handle.name)
            handle.close()
            path = handle.name
            cleanup = True
        else:
            data = request.json or {}
            path = data.get("filepath")
            if path:
                display_name = data.get("name") or os.path.splitext(os.path.basename(path))[0]
                ensure_path("")
        if not path:
            return jsonify({"error": "Upload a .feather file or pass filepath"}), 400
        if not os.path.isfile(path):
            return jsonify({"error": "File not found: " + str(path)}), 404
        xyz, name, n_atoms, overlay = _molecule_view_from_path(path, display_name)
        payload = {
            "xyz": xyz, "name": name, "n_atoms": n_atoms, "mol_data": overlay,
        }
        if not cleanup:
            payload["filepath"] = path
        return jsonify(payload)
    except Exception as exc:
        traceback.print_exc()
        return jsonify({"error": str(exc)}), 500
    finally:
        if cleanup and path and os.path.isfile(path):
            try:
                os.remove(path)
            except OSError:
                pass


@app.route("/example_xyz", methods=["GET", "POST"])
def example_xyz():
    """XYZ for the bundled example reference molecule (basic.feather)."""
    path = example_reference_feather()
    if not path or not os.path.isfile(path):
        return jsonify({
            "error": (
                "No example molecule was found. Expected basic.feather in the "
                "benzene set (descripytor/examples/feather_example)."
            ),
            "filepath": path or "",
        }), 404
    display = os.path.splitext(os.path.basename(path))[0]
    try:
        xyz, name, n_atoms, overlay = _molecule_view_from_path(path, display)
        return jsonify({
            "xyz": xyz, "name": name, "n_atoms": n_atoms, "filepath": path,
            "mol_data": overlay,
        })
    except Exception as exc:
        traceback.print_exc()
        return jsonify({"error": str(exc), "filepath": path}), 500


def ask_folder_dialog(initial=""):
    """Native folder picker on the machine running the GUI (not the browser sandbox).

    Tk must own the process main thread. Flask ``threaded=True`` handles
    ``/browse/folder`` on a worker, so the dialog runs in a short-lived
    subprocess instead of in-process Tk.
    """
    import subprocess
    import sys

    helper = os.path.join(os.path.dirname(os.path.abspath(__file__)), "_folder_dialog.py")
    args = [sys.executable, helper]
    if initial and os.path.isdir(initial):
        args.append(initial)
    try:
        completed = subprocess.run(args, capture_output=True, timeout=600, check=False)
    except subprocess.TimeoutExpired:
        return ""
    if completed.returncode not in (0, None):
        err = (completed.stderr or b"").decode("utf-8", errors="replace").strip()
        raise RuntimeError(err or f"the dialog helper exited with code {completed.returncode}")
    return (completed.stdout or b"").decode("utf-8", errors="replace").strip()


@app.route("/browse/folder", methods=["POST"])
def browse_folder():
    data = request.json or {}
    initial = data.get("initial") or example_feather_dir() or ""
    try:
        path = ask_folder_dialog(initial)
    except Exception as exc:
        traceback.print_exc()
        return jsonify({
            "error": "Could not open a folder dialog. Paste the folder path instead. (" + str(exc) + ")",
        }), 500
    if not path:
        return jsonify({"cancelled": True, "path": ""})
    return jsonify({"path": str(path), "cancelled": False})


def _list_feather_files(directory):
    """The molecule files of a folder, one per name: .feather, structured .json, Gaussian .log, ORCA .out, .xyz."""
    try:
        ensure_path("")
        from data_extractor import _molecules_folder_files

        path, files = _molecules_folder_files(directory)
    except Exception:
        return []
    return [os.path.join(path, name) for name, _stamp, _size in files]


def _flatten_picked_indices(picked):
    found = []

    def walk(node):
        if node is None or node is False:
            return
        if isinstance(node, dict):
            for value in node.values():
                walk(value)
            return
        if isinstance(node, (list, tuple)):
            if node and all(isinstance(x, (int, float)) for x in node):
                for x in node:
                    found.append(int(x))
                return
            for value in node:
                walk(value)
            return
        if isinstance(node, (int, float)) and not isinstance(node, bool):
            found.append(int(node))

    walk(picked)
    seen = set()
    unique = []
    for idx in found:
        if idx > 0 and idx not in seen:
            seen.add(idx)
            unique.append(idx)
    return unique


def _elements_from_xyz_text(xyz):
    lines = [ln for ln in str(xyz or "").splitlines() if ln.strip()]
    if len(lines) < 3:
        return []
    try:
        n_atoms = int(lines[0].split()[0])
    except (TypeError, ValueError):
        n_atoms = max(0, len(lines) - 2)
    elements = []
    for line in lines[2:2 + n_atoms]:
        parts = line.split()
        if parts:
            elements.append(parts[0])
    return elements


@app.route("/folder/summary", methods=["POST"])
def folder_summary():
    """The molecules of a folder (any kind the package reads); with inspect, what their files carry, so the page offers
    only the descriptors every molecule can give."""
    data = request.json or {}
    directory = (data.get("directory") or "").strip()
    if _is_placeholder_path(directory) or not os.path.isdir(directory):
        return jsonify({"error": "Set a folder of molecule files first."}), 400
    files = _list_feather_files(directory)
    names = [os.path.splitext(os.path.basename(p))[0] for p in files]
    payload = {"n": len(files), "names": names, "files": files}   # files: the 3D view loads one by its path
    if data.get("inspect"):
        try:
            from data_extractor import folder_inputs

            payload["inputs"] = folder_inputs(directory)
        except Exception as exc:
            payload["inputs_error"] = str(exc)
    return jsonify(payload)


def _iter_numbering_check(files, reference, picked, include_xyz, limit_xyz, reference_xyz=None, reference_name=""):
    """Yield start/progress/done events while comparing feathers to the reference: the xyz the picks were made on
    when the page sends it, else a reference file."""
    if reference_xyz:
        ref_xyz, ref_name = reference_xyz, reference_name or "the picked molecule"
        ref_n = len(_elements_from_xyz_text(ref_xyz))
    else:
        ref_xyz, ref_name, ref_n = _molecule_xyz_from_path(
            reference, os.path.splitext(os.path.basename(reference))[0]
        )
    ref_elements = _elements_from_xyz_text(ref_xyz)
    n_files = len(files)
    yield {"event": "start", "n": n_files, "reference": ref_name, "i": 0}

    molecules = []
    flagged = []
    for i, path in enumerate(files, 1):
        name = os.path.splitext(os.path.basename(path))[0]
        try:
            xyz, _, n_atoms = _molecule_xyz_from_path(path, name)
            elements = _elements_from_xyz_text(xyz)
        except Exception as exc:
            item = {
                "name": name,
                "path": path,
                "n_atoms": 0,
                "ok": False,
                "error": str(exc),
                "picked": [],
            }
            molecules.append(item)
            flagged.append(item)
            yield {"event": "progress", "i": i, "n": n_files, "name": name, "ok": False}
            continue
        picked_rows = []
        ok = True
        for idx in picked:
            if idx > n_atoms:
                picked_rows.append({
                    "index": idx, "ok": False, "ref": None, "got": None,
                    "detail": "index is past the end of this molecule",
                })
                ok = False
                continue
            got = elements[idx - 1]
            ref_el = ref_elements[idx - 1] if idx <= len(ref_elements) else None
            same = ref_el is None or got == ref_el
            picked_rows.append({
                "index": idx, "ok": same, "ref": ref_el, "got": got,
                "detail": None if same else f"reference is {ref_el}, this file is {got}",
            })
            if not same:
                ok = False
        item = {
            "name": name,
            "path": path,
            "n_atoms": n_atoms,
            "ok": ok,
            "picked": picked_rows,
        }
        if include_xyz:
            item["xyz"] = xyz
            item["elements"] = elements
        molecules.append(item)
        if not ok:
            flagged.append(item)
        yield {"event": "progress", "i": i, "n": n_files, "name": name, "ok": ok}

    overlay = []
    if include_xyz:
        overlay_names = {ref_name}
        overlay.append({
            "name": ref_name, "xyz": ref_xyz, "n_atoms": ref_n, "reference": True,
        })
        for item in flagged:
            if item.get("xyz") and item["name"] not in overlay_names:
                overlay.append({
                    "name": item["name"], "xyz": item["xyz"],
                    "n_atoms": item["n_atoms"], "reference": False,
                })
                overlay_names.add(item["name"])
            if len(overlay) >= limit_xyz:
                break
        if len(overlay) < min(4, limit_xyz):
            for item in molecules:
                if item["name"] in overlay_names or not item.get("xyz"):
                    continue
                overlay.append({
                    "name": item["name"], "xyz": item["xyz"],
                    "n_atoms": item["n_atoms"], "reference": False,
                })
                overlay_names.add(item["name"])
                if len(overlay) >= min(4, limit_xyz):
                    break

    n_ok = sum(1 for m in molecules if m.get("ok"))
    done = {
        "event": "done",
        "ok": n_ok == len(molecules),
        "n_mols": len(molecules),
        "n_ok": n_ok,
        "n_flagged": len(molecules) - n_ok,
        "reference": ref_name,
        "reference_n_atoms": ref_n,
        "picked": picked,
        "molecules": [
            {k: v for k, v in m.items() if k != "xyz"}
            for m in molecules
        ],
        "overlay": overlay,
        "message": (
            "Picked atoms match the reference element in every file."
            if n_ok == len(molecules)
            else (
                f"{len(molecules) - n_ok} of {len(molecules)} molecules differ "
                f"from {ref_name} at a picked atom index. Align & renumber (below) numbers them like it; "
                "Conformer Viewer overlays them."
            )
        ),
    }
    yield done


@app.route("/numbering/check", methods=["POST"])
def numbering_check():
    """Compare atom elements at shared indices vs a reference molecule."""
    from flask import Response, stream_with_context

    data = request.json or {}
    directory = (data.get("directory") or "").strip()
    if _is_placeholder_path(directory) or not os.path.isdir(directory):
        return jsonify({"error": "Set a folder of .feather files first."}), 400
    files = _list_feather_files(directory)
    if not files:
        return jsonify({"error": "No .feather files in that folder."}), 400
    limit = data.get("limit")
    if limit:
        try:
            files = files[: max(1, int(limit))]
        except (TypeError, ValueError):
            pass

    reference = (data.get("reference") or "").strip()
    if not reference or not os.path.isfile(reference):
        preferred = os.path.join(directory, "basic.feather")
        if not os.path.isfile(preferred):
            preferred = os.path.join(directory, "unsub.feather")
        reference = preferred if os.path.isfile(preferred) else files[0]
    picked = _flatten_picked_indices(data.get("picked") or data.get("atoms") or [])
    ref_args = dict(reference_xyz=(data.get("reference_xyz") or "").strip() or None,
                    reference_name=str(data.get("reference_name") or ""))
    include_xyz = bool(data.get("include_xyz"))
    try:
        limit_xyz = max(1, min(int(data.get("limit_xyz") or 8), 18))
    except (TypeError, ValueError):
        limit_xyz = 8
    stream = bool(data.get("stream"))

    if not stream:
        try:
            payload = None
            for ev in _iter_numbering_check(files, reference, picked, include_xyz, limit_xyz, **ref_args):
                if ev.get("event") == "done":
                    payload = {k: v for k, v in ev.items() if k != "event"}
            return jsonify(payload or {"error": "Numbering check produced no result."})
        except Exception as exc:
            traceback.print_exc()
            return jsonify({"error": f"Could not check numbering: {exc}"}), 500

    def generate():
        try:
            for ev in _iter_numbering_check(files, reference, picked, include_xyz, limit_xyz, **ref_args):
                yield json.dumps(ev) + "\n"
        except Exception as exc:
            traceback.print_exc()
            yield json.dumps({"event": "error", "error": str(exc)}) + "\n"

    return Response(
        stream_with_context(generate()),
        mimetype="application/x-ndjson",
        headers={"Cache-Control": "no-cache", "X-Accel-Buffering": "no"},
    )


@app.route("/theta")
def theta_redirect():
    """The explorer's address before it was renamed; the query goes along."""
    q = request.query_string.decode()
    return redirect("/feature_explorer" + ("?" + q if q else ""))


@app.route("/feature_explorer")
def feature_explorer():
    """Serve the live feature explorer (sources in feature_explorer/)."""
    ensure_path("")
    from feature_explorer.build import assemble
    return assemble(), 200, {"Content-Type": "text/html; charset=utf-8"}


@app.route("/feature_explorer/check", methods=["POST"])
def theta_check():
    """The Python Sterimol on the explorer's coordinates and axis, so the page
    can confirm its in-browser port against the package on every pick."""
    data = request.json or {}
    try:
        ensure_path("")
        import io
        import pandas as pd
        from data_extractor import extract_connectivity, get_sterimol_df

        lines = (data.get("xyz") or "").strip().splitlines()
        if lines and lines[0].strip().isdigit():
            lines = lines[2:]
        xyz = pd.read_csv(io.StringIO("\n".join(lines)), sep=r"\s+", header=None,
                          names=["atom", "x", "y", "z"], usecols=range(4))
        bonds = extract_connectivity(xyz)
        r = get_sterimol_df(xyz, bonds, [int(data["a"]), int(data["b"])], None, radii="CPK").iloc[0]
        return jsonify({k: float(r[k]) for k in ("B1", "B5", "L", "loc_B5", "B1_B5_angle")})
    except Exception as e:
        traceback.print_exc()
        return jsonify({"error": str(e)}), 400


@app.route("/figure/load", methods=["POST"])
def figure_load():
    """Geometry, dipole and partial charges for the explorer's Descriptors view,
    read by the package from an uploaded .feather or Gaussian .log (or a local path)."""
    import tempfile
    try:
        ensure_path("")
        from feature_explorer.data import payload
        up = request.files.get("file")
        if up is not None:
            name, ext = os.path.splitext(os.path.basename(up.filename or "upload"))
            if ext.lower() not in (".feather", ".log", ".out"):
                return jsonify({"error": "expected a .feather or Gaussian .log file"}), 400
            with tempfile.TemporaryDirectory() as tmp:
                path = os.path.join(tmp, name + ext.lower().replace(".out", ".log"))
                up.save(path)
                return jsonify(payload(path, name=name))
        data = request.json or {}
        return jsonify(payload(data["filepath"]))
    except Exception as e:
        traceback.print_exc()
        return jsonify({"error": str(e)}), 400


def _module_available(name: str) -> bool:
    import importlib.util
    return importlib.util.find_spec(name) is not None


def _extract_empty_message(log_text):
    lines = (log_text or "").splitlines()
    engine_errors = [
        line.split("[error]", 1)[-1].strip()
        for line in lines
        if "[error]" in line
    ]
    load_errors = [
        line.strip()
        for line in lines
        if line.strip().startswith("Error:") or "could not be processed" in line
    ]
    if engine_errors:
        msg = "Feature extraction failed: " + "; ".join(engine_errors[:4])
        if load_errors:
            msg += " | " + " | ".join(load_errors[:3])
        return msg
    if load_errors:
        return (
            "Extraction produced no columns. "
            + " | ".join(load_errors[:4])
        )
    return (
        "Extraction produced no columns. Pick atoms, check the feather "
        "folder, and look at Advanced engines."
    )


def _failed_descriptors(log_text, names):
    """A descriptor that fails on a molecule is printed by the package and its column left out of the table, which the
    page then reported as a plain success. One line per message, the molecule names folded into a count."""
    import re

    names = sorted({str(n) for n in names if n}, key=len, reverse=True)
    name_re = re.compile(r"(?<![\w-])(" + "|".join(map(re.escape, names)) + r")(?![\w-])") if names else None
    counts = {}
    for line in str(log_text or "").splitlines():
        line = line.strip()
        if not line.startswith("Error"):
            continue
        key = " ".join(name_re.sub("", line).split()) if name_re else line
        counts[key] = counts.get(key, 0) + 1
    return [f"{k} ({c} of {len(names)} molecules)" if name_re and c > 1 else k for k, c in counts.items()]


def _extract_success_payload(df, cfg, log_text):
    feather_dir = (cfg or {}).get("feather_dir") or ""
    names = [os.path.splitext(os.path.basename(f))[0] for f in _list_feather_files(feather_dir)] if feather_dir else []
    warnings = _failed_descriptors(log_text, names)
    local_outcomes = os.path.join(feather_dir, "outcomes.csv") if feather_dir else ""
    rows = []
    if local_outcomes and os.path.isfile(local_outcomes):
        try:
            rows = _read_outcomes_csv(local_outcomes).get("rows") or []
        except Exception as exc:          # a cp1252 file used to lose the finished table to a 500
            warnings.append(f"outcomes.csv was not read ({exc}); save it as CSV UTF-8 with the columns name,output")
    df, n_outcomes = _join_outcomes_into_df(df, rows)      # only this folder's own outcomes.csv
    if os.path.isfile(local_outcomes or "") and not n_outcomes and not any("outcomes.csv" in w for w in warnings):
        warnings.append("outcomes.csv in this folder matched no molecule: it needs the columns name,output, "
                        "one row per .feather name")
    engine_errors = [ln.strip() for ln in str(log_text or "").splitlines() if ln.strip().startswith("[error]")]
    return {
        "n_mols": int(df.shape[0]),
        "n_features": int(df.shape[1]),
        "n_files": len(_list_feather_files(cfg["feather_dir"])),
        "n_outcomes": int(n_outcomes),
        "y_column": "output" if n_outcomes else "",
        "engine_errors": engine_errors,       # an engine that failed gave no columns: the page says so beside the count
        "warnings": warnings,                 # descriptors that failed on some molecules, an unread outcomes.csv
        "columns": [str(c) for c in df.columns],
        "index": [str(i) for i in df.index],
        "csv": df.to_csv(),
        "filename": "features.csv",
        "log": log_text,
    }


def _run_extract(cfg, progress=None):
    import contextlib
    import io

    log = io.StringIO()
    try:
        dx = _load_descriptor_extractor()
        with contextlib.redirect_stdout(log):
            df = dx.run_from_config(cfg, progress=progress)
        return df, log.getvalue(), None
    except Exception as exc:
        return None, log.getvalue(), exc


@app.route("/extract", methods=["POST"])
def extract():
    """Run the picker's run_config and return a features CSV for download."""
    data = request.json or {}
    cfg = data.get("config") if isinstance(data.get("config"), dict) else data
    if not cfg or not isinstance(cfg, dict) or cfg.keys() <= {"config"}:
        return jsonify({"error": "Send the picker run_config as JSON ({config: {...}})"}), 400

    cfg = dict(cfg)
    feather_dir = (data.get("feather_dir") or cfg.get("feather_dir") or "").strip()
    if _is_placeholder_path(feather_dir):
        return jsonify({
            "error": "Set the folder of .feather files (or click Use example set), then Extract CSV."
        }), 400
    if not os.path.isdir(feather_dir):
        return jsonify({"error": f"feather folder not found: {feather_dir}"}), 400

    cfg["feather_dir"] = os.path.abspath(feather_dir)
    cfg["root_dir"] = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))      # never a folder the request names: it goes first on sys.path
    cfg["output_csv"] = None
    engines = cfg.get("engines")
    if isinstance(engines, dict) and isinstance(engines.get("descripytor_full"), dict):
        engines["descripytor_full"] = dict(engines["descripytor_full"])
        engines["descripytor_full"]["save_as"] = False
    if data.get("limit") is not None:
        try:
            cfg["file_limit"] = max(1, int(data["limit"]))
        except (TypeError, ValueError):
            pass

    ensure_path(cfg["root_dir"])
    stream = bool(data.get("stream"))
    n_files = len(_list_feather_files(cfg["feather_dir"]))
    if cfg.get("file_limit"):
        n_files = min(n_files, int(cfg["file_limit"]))

    if not stream:
        df, log_text, exc = _run_extract(cfg)
        if exc is not None:
            traceback.print_exception(type(exc), exc, exc.__traceback__)
            return jsonify({
                "error": f"Feature extraction failed: {exc}",
                "log": log_text,
            }), 500
        if df is None or getattr(df, "empty", True):
            return jsonify({
                "error": _extract_empty_message(log_text),
                "log": log_text,
            }), 500
        return jsonify(_extract_success_payload(df, cfg, log_text))

    from flask import Response, stream_with_context
    from queue import Queue
    import threading

    def generate():
        q = Queue()

        def progress(ev):
            if ev:
                q.put(ev)

        def worker():
            try:
                df, log_text, exc = _run_extract(cfg, progress=progress)
                if exc is not None:
                    traceback.print_exception(type(exc), exc, exc.__traceback__)
                    q.put({
                        "event": "error",
                        "error": f"Feature extraction failed: {exc}",
                        "log": log_text,
                    })
                elif df is None or getattr(df, "empty", True):
                    q.put({
                        "event": "error",
                        "error": _extract_empty_message(log_text),
                        "log": log_text,
                    })
                else:
                    payload = _extract_success_payload(df, cfg, log_text)
                    payload["event"] = "done"
                    q.put(payload)
            except Exception as exc:
                traceback.print_exc()
                q.put({
                    "event": "error",
                    "error": f"Feature extraction failed: {exc}",
                })
            finally:
                q.put(None)

        q.put({"event": "start", "n": n_files, "phase": "load", "i": 0})
        threading.Thread(target=worker, daemon=True).start()
        while True:
            ev = q.get()
            if ev is None:
                break
            yield json.dumps(ev) + "\n"

    return Response(
        stream_with_context(generate()),
        mimetype="application/x-ndjson",
        headers={"Cache-Control": "no-cache", "X-Accel-Buffering": "no"},
    )


import threading

_ALIGN_BUSY = threading.Lock()          # one renumbering at a time
_ALIGN_STOP = {"event": None}           # the running one's stop, for /align/stop (a refreshed page lost its request)


def _folder_molecules(directory):
    """The folder's molecules, loaded as the extraction loads them (same module, same cache)."""
    ensure_path("")
    from data_extractor import Molecules

    mols = Molecules(directory)
    if not mols.molecules:
        raise ValueError("No molecules could be read in that folder.")
    return mols


@app.route("/align/map", methods=["POST"])
def align_map():
    """Number every molecule of a folder like a reference (Conformaty), streamed one molecule at a time as ndjson:
    start (names, reference and its xyz), molecule (one per molecule, as it finishes), done or error. Closing the
    request stops it after the molecules already running."""
    from flask import Response, stream_with_context
    from queue import Queue

    data = request.json or {}
    directory = (data.get("directory") or "").strip()
    if _is_placeholder_path(directory) or not os.path.isdir(directory):
        return jsonify({"error": "Set a folder of .feather files first."}), 400
    try:
        from M2_data_extractor.alignment import _conformaty, map_folder

        _conformaty()
        mols = _folder_molecules(directory)
    except (ImportError, ValueError) as exc:
        return jsonify({"error": str(exc)}), 400
    sizes = {m.molecule_name: len(m.xyz_df) for m in mols.molecules}
    reference = (data.get("reference") or "").strip() or min(sizes, key=lambda k: (sizes[k], k))
    if reference not in sizes:
        return jsonify({"error": f"No molecule named {reference!r} in that folder."}), 400
    try:
        anchors = {str(name): [(int(r), int(a)) for r, a in pairs] for name, pairs in (data.get("anchors") or {}).items()}
        workers = int(data["workers"]) if data.get("workers") else None
    except (TypeError, ValueError):
        return jsonify({"error": "Anchors are {molecule: [[reference atom, atom], ...]}, workers a number."}), 400
    stages = [x for x in (data.get("stages") or []) if x in ("complete_rings", "extend_mcs", "coords_proximity_any")]
    fit = data.get("fit") or "rmsd"
    try:
        sigma = float(data.get("sigma") or 30)
        fit_atoms = [int(a) for a in data.get("fit_atoms") or []] or None
    except (TypeError, ValueError):
        return jsonify({"error": "sigma is a number and fit_atoms a list of atom numbers."}), 400
    if not _ALIGN_BUSY.acquire(blocking=False):
        return jsonify({"error": "A renumbering is still running.", "busy": True}), 409
    ref_mol = next(m for m in mols.molecules if m.molecule_name == reference)
    ref_xyz = _xyz_text_from_coords(ref_mol.xyz_df, reference)[0]
    stop = threading.Event()
    _ALIGN_STOP["event"] = stop
    q = Queue()

    def progress(ev):
        if ev.get("event") == "start":
            ev = dict(ev, ref_xyz=ref_xyz, n_ref=sizes[reference])
        q.put(ev)

    def worker():
        import time

        t0 = time.perf_counter()
        try:
            map_folder(mols, reference, anchors=anchors, stages=stages or None, match_charge=bool(data.get("match_charge")),
                       workers=workers, progress=progress, cancelled=stop.is_set, fit=fit, sigma=sigma, fit_atoms=fit_atoms,
                       isolate=True, only=[str(n) for n in data["only"]] if data.get("only") else None)
            q.put({"event": "stopped" if stop.is_set() else "done", "seconds": time.perf_counter() - t0})
        except Exception as exc:
            traceback.print_exc()
            q.put({"event": "error", "error": f"{type(exc).__name__}: {exc}"})
        finally:
            _ALIGN_BUSY.release()
            q.put(None)

    def generate():
        threading.Thread(target=worker, daemon=True).start()
        try:
            while True:
                ev = q.get()
                if ev is None:
                    break
                yield json.dumps(ev) + "\n"
        finally:                       # the page closed the request (Stop, or left): no new molecules are started
            stop.set()

    return Response(stream_with_context(generate()), mimetype="application/x-ndjson",
                    headers={"Cache-Control": "no-cache", "X-Accel-Buffering": "no"})


@app.route("/align/stop", methods=["POST"])
def align_stop():
    """Stop the running renumbering (its molecule's process is ended); answers once the server is free again."""
    import time

    event = _ALIGN_STOP["event"]
    if event is not None:
        event.set()
    for _ in range(40):
        if not _ALIGN_BUSY.locked():
            return jsonify({"stopped": event is not None, "free": True})
        time.sleep(0.25)
    return jsonify({"stopped": event is not None, "free": False}), 409


@app.route("/align/warm", methods=["POST"])
def align_warm():
    """Start the renumbering worker (its process and its torch/Conformaty imports, ~5 s) before the first run."""
    if _ALIGN_BUSY.locked():          # a run with more workers would lose its pool
        return jsonify({"warm": False})
    try:
        from M2_data_extractor.alignment import _conformaty, _pool

        _conformaty()
        _pool(1).submit(int)
    except ImportError as exc:
        return jsonify({"warm": False, "error": str(exc)})
    return jsonify({"warm": True})


@app.route("/align/example", methods=["POST"])
def align_example():
    """Conformaty's test case 03 (a methyl-oxane numbered onto cyclohexane) copied to a fresh temporary folder, so the
    renumbered copies are written there and never into the package."""
    import shutil
    import tempfile

    if _ALIGN_BUSY.locked():
        return jsonify({"error": "A renumbering is still running.", "busy": True}), 409
    try:
        from descripytor.examples import conformaty_rings_dir

        source = conformaty_rings_dir()
    except (ImportError, FileNotFoundError) as exc:
        return jsonify({"error": f"The example is not installed with this copy: {exc}"}), 404
    dest = os.path.join(tempfile.gettempdir(), "descripytor_conformaty_rings")
    shutil.rmtree(dest, ignore_errors=True)
    shutil.copytree(source, dest)
    return jsonify({"directory": dest, "reference": "stage_order_ref",
                    "about": "Conformaty's test case 03: a methyl-oxane numbered onto cyclohexane. With the default steps the "
                             "ring O takes the reference's carbon 6; with only 'extend the common part' the methyl carbon "
                             "takes 6 and the O is numbered last. Reorder the steps under More options and run it again to see the numbering change."})


@app.route("/align/fits", methods=["POST"])
def align_fits():
    """One molecule, in the numbering the page holds (as mapped, or with its swaps), fitted onto the reference every
    way: Conformaty's five losses and Kabsch, side by side."""
    from M2_data_extractor.alignment import compare_fits

    data = request.json or {}
    directory = (data.get("directory") or "").strip()
    if _is_placeholder_path(directory) or not os.path.isdir(directory):
        return jsonify({"error": "Set a folder of .feather files first."}), 400
    try:
        by = {m.molecule_name: m for m in _folder_molecules(directory).molecules}
        ref, mol = by[str(data.get("reference"))], by[str(data.get("name"))]
        old_to_new = {int(k): int(v) for k, v in (data.get("old_to_new") or {}).items()}
        if sorted(old_to_new.values()) != list(range(1, len(mol.xyz_df) + 1)):
            raise ValueError("the numbering does not give each atom one new number")
        order = [o - 1 for o, _n in sorted(old_to_new.items(), key=lambda kv: kv[1])]
        syms = [str(mol.xyz_df["atom"].iloc[i]) for i in order]
        coords = mol.xyz_df[["x", "y", "z"]].to_numpy(float)[order]
        refa = (ref.xyz_df["atom"].astype(str).tolist(), ref.xyz_df[["x", "y", "z"]].to_numpy(float))
        return jsonify({"fits": compare_fits(refa, syms, coords, float(data.get("sigma") or 30))})
    except (KeyError, ValueError, TypeError) as exc:
        return jsonify({"error": f"Cannot compare the fits: {exc}"}), 400


@app.route("/align/apply", methods=["POST"])
def align_apply():
    """Write the folder renumbered: <folder>/renumbered/<name>.feather from the page's {name: {old: new}} dicts (as
    mapped, or edited on the page), with renumbering.json. Every dict must be a complete renumbering of its molecule;
    the originals are never changed."""
    from M2_data_extractor.alignment import write_renumbered

    data = request.json or {}
    directory = (data.get("directory") or "").strip()
    if _is_placeholder_path(directory) or not os.path.isdir(directory):
        return jsonify({"error": "Set a folder of .feather files first."}), 400
    try:
        mols = _folder_molecules(directory)
        sizes = {m.molecule_name: len(m.xyz_df) for m in mols.molecules}
        results = {}
        for name, mapping in (data.get("molecules") or {}).items():
            if name not in sizes:
                raise ValueError(f"No molecule named {name!r} in that folder.")
            old_to_new = {int(k): int(v) for k, v in mapping.items()}
            full = set(range(1, sizes[name] + 1))
            if set(old_to_new) != full or set(old_to_new.values()) != full:
                raise ValueError(f"{name}: the numbering must give each of its {sizes[name]} atoms one new number.")
            results[name] = {"old_to_new": old_to_new}
        if not results:
            raise ValueError("Nothing to write: number the folder first.")
    except (ValueError, TypeError, AttributeError) as exc:
        return jsonify({"error": str(exc)}), 400
    out = os.path.join(os.path.abspath(directory), "renumbered")
    written, skipped = write_renumbered(mols, results, out, str(data.get("reference") or ""))
    return jsonify({"out": out, "written": written, "skipped": skipped})


def _parse_pasted_outputs(text):
    """Pasted targets as [(name or None, value)]: one number per line (or one line of numbers), or name and value
    per line separated by a comma, semicolon, tab or spaces. Raises ValueError naming the line it cannot read:
    a skipped line used to move every later value up a row."""
    import re

    lines = [ln.strip() for ln in (text or "").splitlines() if ln.strip()]
    if not lines:
        return []

    def number(token):
        try:
            return float(token)
        except (TypeError, ValueError):
            return None

    if len(lines) == 1:
        tokens = [t for t in re.split(r"[\s,;]+", lines[0]) if t]
        if len(tokens) > 1 and all(number(t) is not None for t in tokens):
            return [(None, float(t)) for t in tokens]
    pairs = []
    for k, line in enumerate(lines, 1):
        value = number(line)
        if value is not None:
            pairs.append((None, value))
            continue
        m = re.match(r"^(.*?)\s*[\t,;]\s*(\S+)$", line) or re.match(r"^(.*\S)\s+(\S+)$", line)
        value = number(m.group(2)) if m else None
        if value is None or not m.group(1).strip():
            raise ValueError(f"Line {k} of the pasted values is neither a number nor name,value: {line!r}")
        pairs.append((m.group(1).strip(), value))
    named = [p for p in pairs if p[0] is not None]
    if named and len(named) != len(pairs):
        raise ValueError("Pasted values mix lines with a name and lines without: give every line a name, or none.")
    if named and all(number(n) is not None for n, _ in named):
        raise ValueError("Every pasted name is a number (1,5 was read as name 1, value 5). Use a decimal point, "
                         "one value per line.")
    return pairs


def _merge_outputs_into_df(df, pairs, new_col_name):
    """Add or overwrite new_col_name from (name, value) or (None, value) pairs.
    Unnamed values must be one per row. Named values are matched as outcomes.csv is (case, .feather, - and _ ignored);
    a pasted name that matches no row is an error, a row without a value is left empty."""
    out = df.copy()
    if not pairs:
        return out
    if all(p[0] is not None for p in pairs):
        name_col = next(
            (
                c
                for c in out.columns
                if str(c).lower() in (
                    "name", "names", "molecule", "molecule_name", "mol", "compound",
                )
            ),
            out.columns[0],
        )
        value_map = {_norm_mol_name(n): v for n, v in pairs}
        keys = out[name_col].map(_norm_mol_name)
        unknown = [n for n, _ in pairs if _norm_mol_name(n) not in set(keys)]
        if unknown:
            raise ValueError("No row is named " + ", ".join(map(str, unknown[:8])) + ("..." if len(unknown) > 8 else "")
                             + f" ({len(unknown)} of {len(pairs)} pasted names). Names are matched to the first column.")
        out[new_col_name] = keys.map(value_map).to_numpy()
    else:
        values = [v for _, v in pairs]
        if len(values) != len(out):
            raise ValueError(f"{len(values)} values pasted for {len(out)} rows: give one per row, in the table's order, "
                             "or name,value lines.")
        out[new_col_name] = values
    return out


def _linear_combo_search(df, y_name, min_features, max_features, top_n, threshold, progress=None):
    """OLS feature-combination search with leave-one-out Q². Returns a results DataFrame."""
    from itertools import combinations

    import numpy as np
    import pandas as pd
    from sklearn.linear_model import LinearRegression
    from sklearn.metrics import r2_score
    from sklearn.model_selection import LeaveOneOut
    from sklearn.preprocessing import StandardScaler

    if y_name not in df.columns:
        raise ValueError(f"No target column '{y_name}'. Choose a CSV column or paste values.")

    y = pd.to_numeric(df[y_name], errors="coerce")
    drop_names = {
        str(c)
        for c in df.columns
        if str(c).lower() in (
            "name", "names", "smiles", "id", "molecule", "molecule_name", "mol", "compound",
        )
        or str(c).startswith("Unnamed")
    }
    drop_names.add(y_name)
    X = df.drop(columns=[c for c in df.columns if c in drop_names or c == y_name], errors="ignore")
    X = X.select_dtypes(include=["number", "bool"])
    mask = y.notna()
    X = X.loc[mask]
    y = y.loc[mask]
    X = X.replace([np.inf, -np.inf], np.nan)
    X = X.dropna(axis=1, how="all")
    X = X.dropna(axis=1, thresh=max(3, int(len(X) * 0.7)))
    row_ok = X.notna().all(axis=1)
    X = X.loc[row_ok]
    y = y.loc[row_ok]
    X = X.fillna(X.median(numeric_only=True))
    X = X.loc[:, X.nunique(dropna=True) > 1]
    if X.shape[0] < 3:
        raise ValueError("Need at least 3 rows with a numeric target.")
    if X.shape[1] < 1:
        raise ValueError("No numeric feature columns left after dropping the target.")

    n_feat = int(X.shape[1])
    min_features = max(1, min(int(min_features), n_feat, 3))
    max_features = max(min_features, min(int(max_features), n_feat, 3))
    n_combos = sum(
        1
        for k in range(min_features, max_features + 1)
        for _ in combinations(range(n_feat), k)
    )
    if n_combos > 15000:
        raise ValueError(
            f"{n_combos:,} combinations is too many for the picker. "
            "Lower max features (1–2) or drop extra columns."
        )

    scaler = StandardScaler()
    Xs = scaler.fit_transform(X.to_numpy(dtype=float))
    ys = y.to_numpy(dtype=float)
    n = len(ys)
    cols = list(X.columns)
    loo = LeaveOneOut()
    rows = []
    if progress:
        progress({
            "event": "start",
            "n": n_combos,
            "i": 0,
            "phase": "search",
            "n_features": n_feat,
            "min_features": min_features,
            "max_features": max_features,
        })

    def _q2(Xm, ym):
        pred = np.zeros(len(ym))
        for train, test in loo.split(Xm):
            model = LinearRegression()
            model.fit(Xm[train], ym[train])
            pred[test] = model.predict(Xm[test])
        ss_res = float(np.sum((ym - pred) ** 2))
        ss_tot = float(np.sum((ym - ym.mean()) ** 2))
        if ss_tot <= 0:
            return float("nan")
        return 1.0 - ss_res / ss_tot

    done = 0
    step = max(1, n_combos // 40)
    for k in range(min_features, max_features + 1):
        for combo in combinations(range(n_feat), k):
            Xm = Xs[:, combo]
            model = LinearRegression()
            model.fit(Xm, ys)
            pred = model.predict(Xm)
            r2 = float(r2_score(ys, pred))
            if r2 >= threshold:
                p = k
                adj = 1.0 - (1.0 - r2) * (n - 1) / max(n - p - 1, 1)
                rows.append({
                    "combination": ", ".join(cols[i] for i in combo),
                    "n_features": k,
                    "r2": round(r2, 4),
                    "adj_r2": round(float(adj), 4),
                    "q2": round(float(_q2(Xm, ys)), 4),
                })
            done += 1
            if progress and (done == 1 or done == n_combos or done % step == 0):
                progress({
                    "event": "progress",
                    "i": done,
                    "n": n_combos,
                    "phase": "search",
                    "kept": len(rows),
                    "k": k,
                })

    if not rows:
        raise ValueError(
            f"No models reached R² ≥ {threshold}. Lower the threshold or paste a numeric target."
        )
    out = pd.DataFrame(rows).sort_values(["q2", "r2"], ascending=False).head(int(top_n))
    return out.reset_index(drop=True)


def _model_result_payload(results, y_name, n_rows):
    return {
        "n_models": int(len(results)),
        "target": y_name,
        "n_rows": int(n_rows),
        "columns": [str(c) for c in results.columns],
        "rows": results.to_dict(orient="records"),
        "csv": results.to_csv(index=False),
    }


def _table_with_target(data):
    """The picker's extracted CSV with its target column: one already in the table, or pasted values merged in.
    Returns (table, target name); raises ValueError with the message the page shows."""
    import io

    import pandas as pd

    csv_text = data.get("csv") or ""
    if not str(csv_text).strip():
        raise ValueError("Extract features first (no CSV on the page).")
    try:
        df = pd.read_csv(io.StringIO(csv_text))
    except Exception as exc:
        raise ValueError(f"Could not parse CSV: {exc}")
    if df.empty:
        raise ValueError("CSV is empty.")
    pasted = data.get("outputs") or data.get("pasted") or ""
    y_name = (data.get("y_column") or data.get("target") or "").strip()
    new_col = (data.get("new_column") or "output").strip() or "output"
    if str(pasted).strip():
        pairs = _parse_pasted_outputs(pasted)
        if not pairs:
            raise ValueError("Could not parse pasted values. Use one number per line or name,value.")
        y_name = new_col
        df = _merge_outputs_into_df(df, pairs, y_name)
        for c in list(df.columns):            # an older output column with another name would otherwise be searched as a descriptor
            if c != y_name and str(c).strip().lower() == "output":
                df = df.drop(columns=[c])
    elif not y_name:
        raise ValueError("Choose a target column from the CSV, or paste output values.")
    return df, y_name


@app.route("/model/run", methods=["POST"])
def model_run():
    """Fit linear regression on the extracted CSV (chosen or pasted target)."""
    import io

    import pandas as pd

    data = request.json or {}
    try:
        df, y_name = _table_with_target(data)
    except ValueError as exc:
        return jsonify({"error": str(exc)}), 400

    try:
        min_f = int(data.get("min_features", 1) or 1)
        max_f = int(data.get("max_features", 2) or 2)
        top_n = int(data.get("top_n", 8) or 8)
        threshold = float(data.get("threshold", 0.2) or 0.0)
    except (TypeError, ValueError) as exc:
        return jsonify({"error": str(exc)}), 400

    stream = bool(data.get("stream"))
    if not stream:
        try:
            results = _linear_combo_search(df, y_name, min_f, max_f, top_n, threshold)
        except ValueError as exc:
            return jsonify({"error": str(exc)}), 400
        except Exception as exc:
            traceback.print_exc()
            return jsonify({"error": f"Linear regression failed: {exc}"}), 500
        return jsonify(_model_result_payload(results, y_name, df.shape[0]))

    from flask import Response, stream_with_context
    from queue import Queue
    import threading

    def generate():
        q = Queue()

        def progress(ev):
            if ev:
                q.put(ev)

        def worker():
            try:
                results = _linear_combo_search(
                    df, y_name, min_f, max_f, top_n, threshold, progress=progress
                )
                payload = _model_result_payload(results, y_name, df.shape[0])
                payload["event"] = "done"
                q.put(payload)
            except ValueError as exc:
                q.put({"event": "error", "error": str(exc)})
            except Exception as exc:
                traceback.print_exc()
                q.put({"event": "error", "error": f"Linear regression failed: {exc}"})
            finally:
                q.put(None)

        threading.Thread(target=worker, daemon=True).start()
        while True:
            ev = q.get()
            if ev is None:
                break
            yield json.dumps(ev) + "\n"

    return Response(
        stream_with_context(generate()),
        mimetype="application/x-ndjson",
        headers={"Cache-Control": "no-cache", "X-Accel-Buffering": "no"},
    )


# Raised whenever the pages need a route an older server does not have. A server still running from before then answers
# /status with a lower number, and a new start tells the user to stop it instead of reusing it (its pages are read from
# disk, so they would be new while its routes are old: the links to /feature_explorer and /bassa_explorer gave Not Found).
GUI_API = 4      # 4: pages from other sites refused, the removed modeling panel's routes gone


def _tcp_port_open(port):
    import socket

    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    sock.settimeout(0.4)
    try:
        return sock.connect_ex(("127.0.0.1", int(port))) == 0
    finally:
        sock.close()


def _status_from_running_gui(port):
    import urllib.request

    try:
        with urllib.request.urlopen(f"http://127.0.0.1:{int(port)}/status", timeout=1.2) as resp:
            return json.loads(resp.read().decode("utf-8") or "{}")
    except Exception:
        return None


def _gui_payload_is_current(payload):
    if not isinstance(payload, dict) or not payload.get("ok"):
        return False
    return payload.get("api") == GUI_API


def _pids_listening_on_port(port):
    """PIDs with a LISTEN socket on this TCP port (Windows netstat / Unix lsof)."""
    import subprocess

    port = int(port)
    pids = []
    suffix = ":" + str(port)
    try:
        if os.name == "nt":
            completed = subprocess.run(
                ["netstat", "-ano", "-p", "tcp"],
                capture_output=True, text=True, timeout=4, check=False,
            )
            for line in (completed.stdout or "").splitlines():
                parts = line.split()
                if len(parts) < 5 or parts[0].upper() != "TCP":
                    continue
                if "LISTEN" not in parts[3].upper():
                    continue
                if not parts[1].endswith(suffix):
                    continue
                pid = parts[-1].strip()
                if pid.isdigit() and pid not in pids:
                    pids.append(pid)
        else:
            completed = subprocess.run(
                ["lsof", f"-iTCP:{port}", "-sTCP:LISTEN", "-n", "-P", "-t"],
                capture_output=True, text=True, timeout=4, check=False,
            )
            for pid in (completed.stdout or "").splitlines():
                pid = pid.strip()
                if pid.isdigit() and pid not in pids:
                    pids.append(pid)
    except Exception:
        return []
    return pids


def _reuse_or_refuse_existing_gui(host, port, open_browser, page="/visual"):
    """If 7432 is already taken, reuse a current GUI or tell the user to stop the old one."""
    import webbrowser

    url = f"http://127.0.0.1:{port}{page}"
    payload = _status_from_running_gui(port)
    if _gui_payload_is_current(payload):
        print(f"\n  DescriPyTor GUI is already running.")
        print(f"  Open:  {url}")
        print(f"  This process will not start a second server.\n")
        if open_browser:
            try:
                webbrowser.open(url)
            except Exception:
                pass
        return True
    if payload is not None or _tcp_port_open(port):
        pids = _pids_listening_on_port(port)
        print(f"\n  Port {port} is already in use by an older or different process.")
        print("  It was started before the package was updated, so it lacks the newer pages. Stop it, then run the command again.")
        if pids:
            print("  PID(s): " + ", ".join(pids))
            if os.name == "nt":
                kills = "  &  ".join("taskkill /PID %s /F" % pid for pid in pids)
                print("  Stop them:  " + kills)
            else:
                print("  Stop them:  kill " + " ".join(pids))
        else:
            print("  Stop the old `descripytor visual` terminal with Ctrl+C, then retry.")
        print()
        return False
    return None


def serve(host=None, port=None, open_browser=True, page="/visual"):
    """Start the Flask GUI. Used by `descripytor visual_extractor`, `descripytor feature_explorer` and `__main__`.

    ``page`` is the address the browser opens: "/visual" (the atom picker) or "/feature_explorer?mol=..." (the explorer).
    """
    import threading
    import time
    import urllib.request
    import webbrowser

    host = host or os.environ.get("GUI_HOST", "127.0.0.1")
    port = int(port or os.environ.get("GUI_PORT", str(PORT)))
    if host in ("0.0.0.0", "::", ""):
        print("  Listening on every interface, but pages are answered only at 127.0.0.1 / localhost. To open it from "
              "another computer, start with --host <this computer's address>.")
    else:
        _LOCAL_HOSTS.add(host.strip("[]"))
    os.environ["GUI_PORT"] = str(port)                        # a BASSA explorer started from a page links back to this port
    existing = _reuse_or_refuse_existing_gui(host, port, open_browser, page)
    if existing is True:
        return
    if existing is False:
        sys.exit(1)
    url = f"http://127.0.0.1:{port}{page}"
    print(f"\n  DescriPyTor GUI")
    print(f"  Open:  {url}")
    print(f"  Feature explorer: http://127.0.0.1:{port}/feature_explorer")
    print(f"  Press Ctrl+C to stop\n")
    if open_browser:
        def _open():
            status = f"http://127.0.0.1:{port}/status"
            for _ in range(50):
                try:
                    urllib.request.urlopen(status, timeout=0.4)
                    webbrowser.open(url)
                    return
                except Exception:
                    time.sleep(0.2)

        threading.Thread(target=_open, daemon=True).start()
    try:
        app.run(host=host, port=port, debug=False, use_reloader=False, threaded=True)
    except OSError as exc:
        print(f"\n  Could not bind port {port}: {exc}")
        print("  Stop the process already using that port, then retry.\n")
        sys.exit(1)


# ── run ───────────────────────────────────────────────────────
if __name__ == "__main__":
    serve()
