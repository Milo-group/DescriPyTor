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

# ── Flask ─────────────────────────────────────────────────────
try:
    from flask import Flask, request, jsonify, send_file, redirect, send_from_directory
    from flask_cors import CORS
except ImportError:
    print("Missing dependencies. Run:  pip install flask flask-cors")
    sys.exit(1)

app = Flask(__name__, static_folder=None)
CORS(app)

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


def load_molecule(filepath: str, root: str = ""):
    ensure_path(root)
    from data_extractor import Molecule          # noqa: E402
    return Molecule(filepath)


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


def baptiste_example_dir():
    """Bundled Baptiste product set, if present."""
    try:
        from descripytor.examples import baptiste_example_dir as _baptiste_dir

        path = _baptiste_dir()
        if path.is_dir():
            return str(path)
    except Exception:
        pass
    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.dirname(here)
    for candidate in (
        os.path.join(root, "descripytor", "examples", "baptiste_products"),
        os.path.join(root, "Getting_started_with_examples", "baptiste_products"),
    ):
        if os.path.isdir(candidate) and os.path.isfile(os.path.join(candidate, "unsub.feather")):
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


def example_outcomes_payload():
    directory = example_feather_dir()
    if not directory:
        return {"rows": [], "y_column": "output", "name_column": "name"}
    path = os.path.join(directory, "outcomes.csv")
    if not os.path.isfile(path):
        return {"rows": [], "y_column": "output", "name_column": "name"}
    import csv

    rows = []
    with open(path, newline="", encoding="utf-8") as handle:
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


def _is_placeholder_path(path):
    if not path or not str(path).strip():
        return True
    text = str(path).replace("/", "\\").lower()
    if "path\\to\\your" in text or "path\\to\\logs" in text:
        return True
    stripped = str(path).strip()
    return stripped.startswith("__") and stripped.endswith("__")


def _load_descriptor_extractor():
    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.dirname(here)
    toolkit = os.path.join(
        root, "Getting_started_with_examples", "descriptor_extraction_toolkit"
    )
    if os.path.isdir(toolkit) and toolkit not in sys.path:
        sys.path.insert(0, toolkit)
    import descriptor_extractor as dx  # noqa: WPS433

    return dx


@app.route("/status")
def status():
    outcomes = example_outcomes_payload()
    return jsonify({
        "ok": True,
        "version": "1.0",
        "api": 2,
        "example_feather_dir": example_feather_dir(),
        "example_basic_feather": example_basic_feather(),
        "example_reference_feather": example_reference_feather(),
        "example_reference_name": os.path.splitext(
            os.path.basename(example_reference_feather() or "")
        )[0],
        "example_presets": example_presets(),
        "example_outcomes": outcomes.get("rows") or [],
        "example_y_column": outcomes.get("y_column") or "output",
    })


@app.route("/")
def root():
    """First-run URL: the 3D picker."""
    return redirect("/visual")


@app.route("/forms")
def gui():
    """Form-based feature GUI (same engine, no 3D picker)."""
    return send_file(
        os.path.join(os.path.dirname(__file__), "feature_extraction_gui.html"),
        mimetype="text/html",
    )


def atom_picker_html():
    """Picker page: packaged copy first, then the clone path used by make_picker."""
    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.dirname(here)
    candidates = (
        os.path.join(here, "atom_picker.html"),
        os.path.join(
            root,
            "Getting_started_with_examples",
            "descriptor_extraction_toolkit",
            "atom_picker.html",
        ),
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
    return send_file(path, mimetype="text/html")


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


def _xyz_from_feather_fast(path, display_name="molecule"):
    """Read only atom coordinates from a .feather — no full Molecule() load."""
    import pandas as pd

    data = pd.read_feather(path)
    data.columns = [str(c).strip() for c in data.columns]
    needed = ["atom", "x", "y", "z"]
    if all(c in data.columns for c in needed):
        df = data[needed].copy()
    else:
        # Older benzene example feathers store xyz in the first four columns.
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
    return _xyz_text_from_coords(df, label)


def _xyz_text_from_molecule(mol, name=None):
    df = mol.xyz_df[["atom", "x", "y", "z"]]
    label = name or getattr(mol, "molecule_name", "molecule")
    return _xyz_text_from_coords(df, label)


def _molecule_xyz_from_path(path, display_name="molecule"):
    """XYZ for the 3D viewer. Prefer a fast feather read; fall back to Molecule()."""
    try:
        return _xyz_from_feather_fast(path, display_name)
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
    return _xyz_text_from_molecule(mol, display_name)


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
                ensure_path(data.get("root", ""))
        if not path:
            return jsonify({"error": "Upload a .feather file or pass filepath"}), 400
        if not os.path.isfile(path):
            return jsonify({"error": "File not found: " + str(path)}), 404
        xyz, name, n_atoms = _molecule_xyz_from_path(path, display_name)
        return jsonify({"xyz": xyz, "name": name, "n_atoms": n_atoms})
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
        xyz, name, n_atoms = _molecule_xyz_from_path(path, display)
        return jsonify({
            "xyz": xyz, "name": name, "n_atoms": n_atoms, "filepath": path,
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
        completed = subprocess.run(
            args,
            capture_output=True,
            timeout=600,
            check=False,
        )
    except Exception:
        return _ask_folder_dialog_inline(initial)
    if completed.returncode not in (0, None):
        err = (completed.stderr or b"").decode("utf-8", errors="replace").strip()
        if err:
            raise RuntimeError(err)
        return _ask_folder_dialog_inline(initial)
    return (completed.stdout or b"").decode("utf-8", errors="replace").strip()


def _ask_folder_dialog_inline(initial=""):
    """Last-resort in-process Tk dialog (works if this is the main thread)."""
    import tkinter as tk
    from tkinter import filedialog

    root = tk.Tk()
    root.withdraw()
    try:
        root.wm_attributes("-topmost", True)
    except tk.TclError:
        pass
    kwargs = {"title": "Choose a folder of .feather files"}
    if initial:
        kwargs["initialdir"] = initial
    path = filedialog.askdirectory(**kwargs)
    try:
        root.destroy()
    except tk.TclError:
        pass
    return path or ""


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
    try:
        names = os.listdir(directory)
    except OSError:
        return []
    out = []
    for name in sorted(names):
        if name.lower().endswith((".feather", ".ftr")):
            out.append(os.path.join(directory, name))
    return out


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
    """Fast count of .feather files — used for extract progress text."""
    data = request.json or {}
    directory = (data.get("directory") or "").strip()
    if _is_placeholder_path(directory) or not os.path.isdir(directory):
        return jsonify({"error": "Set a folder of .feather files first."}), 400
    files = _list_feather_files(directory)
    names = [os.path.splitext(os.path.basename(p))[0] for p in files]
    return jsonify({"n": len(files), "names": names})


def _iter_numbering_check(files, reference, picked, include_xyz, limit_xyz):
    """Yield start/progress/done events while comparing feathers to the reference."""
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
                f"from {ref_name} at a picked atom index. Open Conformer Viewer "
                "to overlay them."
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
    include_xyz = bool(data.get("include_xyz"))
    try:
        limit_xyz = max(1, min(int(data.get("limit_xyz") or 8), 18))
    except (TypeError, ValueError):
        limit_xyz = 8
    stream = bool(data.get("stream"))

    if not stream:
        try:
            payload = None
            for ev in _iter_numbering_check(files, reference, picked, include_xyz, limit_xyz):
                if ev.get("event") == "done":
                    payload = {k: v for k, v in ev.items() if k != "event"}
            return jsonify(payload or {"error": "Numbering check produced no result."})
        except Exception as exc:
            traceback.print_exc()
            return jsonify({"error": f"Could not check numbering: {exc}"}), 500

    def generate():
        try:
            for ev in _iter_numbering_check(files, reference, picked, include_xyz, limit_xyz):
                yield json.dumps(ev) + "\n"
        except Exception as exc:
            traceback.print_exc()
            yield json.dumps({"event": "error", "error": str(exc)}) + "\n"

    return Response(
        stream_with_context(generate()),
        mimetype="application/x-ndjson",
        headers={"Cache-Control": "no-cache", "X-Accel-Buffering": "no"},
    )


@app.route("/sterimol", methods=["POST"])
def sterimol():
    data = request.json or {}
    try:
        mol = load_molecule(data["filepath"], data.get("root", ""))

        base_atoms = data.get("base_atoms")
        if not base_atoms or len(base_atoms) != 3:
            return jsonify({"error": "base_atoms must be [origin, direction, from_dir]"}), 400

        drop = data.get("drop_atoms") or None
        df   = mol.get_sterimol(
            base_atoms    = base_atoms,
            radii         = data.get("radii", "CPK"),
            sub_structure = data.get("sub_structure", True),
            drop_atoms    = drop if drop else None,
            mode          = data.get("mode", "all"),
        )
        # flatten the DataFrame to a plain dict  {B1: 2.34, B5: ...}
        result = {k: float(v) for k, v in df.iloc[:, 0].items()}
        return jsonify({"result": result})

    except Exception as e:
        traceback.print_exc()
        return jsonify({"error": str(e)}), 500


@app.route("/charges", methods=["POST"])
def charges():
    data = request.json or {}
    try:
        mol = load_molecule(data["filepath"], data.get("root", ""))

        indices     = data.get("atom_indices") or None
        charge_type = data.get("charge_type", "all")

        df = mol.get_charge_df(atoms_indices=indices, type=charge_type)
        result = json.loads(df.to_json())
        return jsonify({"result": result})

    except Exception as e:
        traceback.print_exc()
        return jsonify({"error": str(e)}), 500


@app.route("/dipole", methods=["POST"])
def dipole():
    data = request.json or {}
    try:
        mol = load_molecule(data["filepath"], data.get("root", ""))
        df  = mol.gauss_dipole_df
        if df is None or len(df) == 0:
            return jsonify({"error": "No dipole data in this file"}), 404
        result = {k: float(df[k].iloc[0]) for k in df.columns}
        return jsonify({"result": result})

    except Exception as e:
        traceback.print_exc()
        return jsonify({"error": str(e)}), 500


@app.route("/vibrations", methods=["POST"])
def vibrations():
    data = request.json or {}
    try:
        mol = load_molecule(data["filepath"], data.get("root", ""))
        df  = mol.info_df
        if df is None or len(df) == 0:
            return jsonify({"error": "No vibrational data in this file"}), 404
        result = json.loads(df[["Frequency", "IR"]].to_json())
        return jsonify({"result": result})

    except Exception as e:
        traceback.print_exc()
        return jsonify({"error": str(e)}), 500


@app.route("/molecule_info", methods=["POST"])
def molecule_info():
    """Return XYZ + all available metadata for a single molecule."""
    data = request.json or {}
    try:
        mol    = load_molecule(data["filepath"], data.get("root", ""))
        atoms  = mol.xyz_df
        n      = len(atoms)

        xyz_lines = [str(n), mol.molecule_name]
        for _, row in atoms.iterrows():
            xyz_lines.append(
                f"{row['atom']}  {float(row['x']):.6f}"
                f"  {float(row['y']):.6f}  {float(row['z']):.6f}"
            )

        energy = None
        try:
            ev = mol.energy_value
            if ev is not None and len(ev) > 0:
                energy = float(ev.iloc[0, 0])
        except Exception:
            pass

        dipole_vals = None
        try:
            d = mol.gauss_dipole_df
            if d is not None and len(d) > 0:
                dipole_vals = {k: float(d[k].iloc[0]) for k in d.columns}
        except Exception:
            pass

        charge_types = []
        try:
            cd = mol.charge_dict
            if cd:
                for ct in ("nbo", "hirshfeld", "cm5"):
                    if cd.get(ct) is not None and len(cd[ct]) > 0:
                        charge_types.append(ct)
        except Exception:
            pass

        return jsonify({
            "name":         mol.molecule_name,
            "n_atoms":      n,
            "xyz":          "\n".join(xyz_lines),
            "energy":       energy,
            "has_dipole":   dipole_vals is not None,
            "dipole":       dipole_vals,
            "charge_types": charge_types,
        })

    except Exception as e:
        traceback.print_exc()
        return jsonify({"error": str(e)}), 500


def _fig_to_b64(fig, dpi: int = 130) -> str:
    """Render a matplotlib Figure to a base64-encoded PNG string."""
    import io, base64
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight")
    buf.seek(0)
    b64 = base64.b64encode(buf.read()).decode("utf-8")
    return b64


def _sterimol_plots(mol, base_atoms, radii, sub_structure, drop_atoms, mode,
                    n_points, dpi, endon_title, side_title):
    """
    Compute Sterimol, generate both steriplots, return
    (sterimol_dict, endon_b64, side_b64).
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from extractor_utils.sterimol_utils import (
        STERIMOL_CIRCLE_POINTS,
        get_sterimol_plot_data,
    )
    from utils.visualize import plot_b1_visualization, plot_L_B5_plane

    st_df = mol.get_sterimol(
        base_atoms    = base_atoms,
        radii         = radii,
        sub_structure = sub_structure,
        drop_atoms    = drop_atoms or None,
        mode          = mode,
    )
    st_result = {k: float(v) for k, v in st_df.iloc[:, 0].items()}

    # Same pipeline get_sterimol_df runs, so the plots show the geometry the
    # numbers above came from.
    rotated_df, rotated_plane = get_sterimol_plot_data(
        mol.coordinates_df, mol.bonds_df, base_atoms,
        radii=radii, sub_structure=sub_structure,
        drop_atoms=drop_atoms or None, mode=mode,
    )

    fig_endon = plot_b1_visualization(
        rotated_plane, rotated_df,
        sterimol_df=st_df,
        # the plane was sampled at this rate; plot_b1_visualization slices it
        # back into circles, so it must not use the caller's n_points here
        n_points=STERIMOL_CIRCLE_POINTS, title=endon_title,
    )
    endon_b64 = _fig_to_b64(fig_endon, dpi)
    plt.close(fig_endon)

    fig_side = plot_L_B5_plane(
        rotated_df, st_df, n_points=n_points, title=side_title,
    )
    side_b64 = _fig_to_b64(fig_side, dpi)
    plt.close(fig_side)

    return st_result, endon_b64, side_b64


@app.route("/steriplot", methods=["POST"])
def steriplot():
    """Compute Sterimol + generate both steriplots for a single molecule."""
    data = request.json or {}
    try:
        base_atoms = data.get("base_atoms")
        if not base_atoms or len(base_atoms) != 3:
            return jsonify({"error": "base_atoms must be [origin, direction, from_dir]"}), 400

        mol = load_molecule(data["filepath"], data.get("root", ""))

        st_result, endon_b64, side_b64 = _sterimol_plots(
            mol         = mol,
            base_atoms  = base_atoms,
            radii       = data.get("radii", "CPK"),
            sub_structure = data.get("sub_structure", True),
            drop_atoms  = data.get("drop_atoms") or None,
            mode        = data.get("mode", "all"),
            n_points    = int(data.get("n_points", 100)),
            dpi         = int(data.get("dpi", 130)),
            endon_title = data.get("endon_title", "XZ plane — End-on view"),
            side_title  = data.get("side_title",  "YZ plane — Side view"),
        )
        return jsonify({
            "name":      os.path.splitext(os.path.basename(data["filepath"]))[0],
            "result":    st_result,
            "endon_img": endon_b64,
            "side_img":  side_b64,
        })

    except Exception as e:
        traceback.print_exc()
        return jsonify({"error": str(e)}), 500


@app.route("/batch_sterimol", methods=["POST"])
def batch_sterimol():
    """
    Compute Sterimol + steriplots for a list of molecules.

    Body:
      filepaths   – list of absolute feather file paths
      base_atoms  – [origin, direction, from_dir] (same for all)
      radii, sub_structure, drop_atoms, mode, n_points, dpi (optional)
    """
    data = request.json or {}
    filepaths  = data.get("filepaths", [])
    base_atoms = data.get("base_atoms")

    if not base_atoms or len(base_atoms) != 3:
        return jsonify({"error": "base_atoms must be [origin, direction, from_dir]"}), 400
    if not filepaths:
        return jsonify({"error": "filepaths list is empty"}), 400

    results = []
    for fp in filepaths:
        entry = {"filepath": fp,
                 "name": os.path.splitext(os.path.basename(fp))[0]}
        try:
            mol = load_molecule(fp, data.get("root", ""))
            st_result, endon_b64, side_b64 = _sterimol_plots(
                mol         = mol,
                base_atoms  = base_atoms,
                radii       = data.get("radii", "CPK"),
                sub_structure = data.get("sub_structure", True),
                drop_atoms  = data.get("drop_atoms") or None,
                mode        = data.get("mode", "all"),
                n_points    = int(data.get("n_points", 80)),
                dpi         = int(data.get("dpi", 110)),
                endon_title = data.get("endon_title", "XZ"),
                side_title  = data.get("side_title",  "YZ"),
            )
            entry.update({"result": st_result,
                          "endon_img": endon_b64, "side_img": side_b64})
        except Exception as e:
            traceback.print_exc()
            entry["error"] = str(e)
        results.append(entry)

    return jsonify({"results": results})


@app.route("/features_set", methods=["POST"])
def features_set():
    """
    Run Molecules.get_molecules_features_set() over the loaded dataset.

    Body (JSON):
      dir_path        – absolute path to the directory of .feather files
      root            – DescriPyTor root to add to sys.path (optional)
      entry_widgets   – dict of string values: {ring, stretching, stretch,
                        upper_stretch, bending, bend, dipole,
                        charges, charge_diff, sterimol, drop_atoms,
                        bond_angle, bond_length}
      parameters      – {Radii, Isotropic}  (optional, defaults used if absent)
      selected_names  – list of molecule names to keep (optional, all if absent)
      save_as         – bool
      csv_file_name   – str
      corr_thresh     – float, default 0.8
    """
    import numpy as np

    data = request.json or {}
    dir_path = data.get("dir_path", "")
    if not dir_path or not os.path.isdir(dir_path):
        return jsonify({"error": f"dir_path not found: {dir_path!r}"}), 400

    ensure_path(data.get("root", ""))

    try:
        from data_extractor import Molecules
    except ImportError as e:
        return jsonify({"error": f"Cannot import Molecules: {e}"}), 500

    try:
        mols = Molecules(dir_path)
    except Exception as e:
        traceback.print_exc()
        return jsonify({"error": f"Molecules load failed: {e}"}), 500

    # optional filtering by selected molecule names
    selected_names = data.get("selected_names") or None
    if selected_names:
        selected_set = set(selected_names)
        mols.molecules = [m for m in mols.molecules if m.molecule_name in selected_set]
        if not mols.molecules:
            return jsonify({"error": "No molecules matched selected_names"}), 400

    entry_widgets = data.get("entry_widgets", {})
    parameters    = data.get("parameters", {"Radii": "CPK", "Isotropic": True})
    save_as       = bool(data.get("save_as", False))
    csv_file_name = data.get("csv_file_name", "features_output")
    corr_thresh   = float(data.get("corr_thresh", 0.8))

    try:
        res_df = mols.get_molecules_features_set(
            entry_widgets = entry_widgets,
            parameters    = parameters,
            save_as       = save_as,
            csv_file_name = csv_file_name,
        )
    except Exception as e:
        traceback.print_exc()
        return jsonify({"error": f"Feature extraction failed: {e}"}), 500

    if res_df is None or res_df.empty:
        return jsonify({"error": "Feature extraction returned an empty DataFrame"}), 500

    # ── diagnostics ──────────────────────────────────────────
    n_mols     = len(res_df)
    n_features = len(res_df.columns)

    # NaN summary
    nan_info = {}
    for col in res_df.columns:
        pct = float(res_df[col].isna().mean() * 100)
        if pct > 0:
            nan_info[col] = round(pct, 1)

    # correlation pairs above threshold
    corr_pairs = []
    try:
        numeric_df = res_df.select_dtypes(include=[np.number])
        if len(numeric_df.columns) > 1:
            corr_mat = numeric_df.corr().abs()
            arr      = corr_mat.to_numpy(copy=True)
            arr[np.tril_indices_from(arr)] = np.nan
            corr_mat2 = _pd_from_numpy(arr, corr_mat.index, corr_mat.columns)
            pairs_idx = corr_mat2.stack()[corr_mat2.stack() >= corr_thresh].index.tolist()
            corr_pairs = [
                {"a": a, "b": b, "r": round(float(corr_mat.loc[a, b]), 4)}
                for a, b in pairs_idx
            ]
    except Exception:
        pass

    # serialise DataFrame: {col: {mol_name: value}}
    records = json.loads(res_df.to_json())

    return jsonify({
        "n_mols":     n_mols,
        "n_features": n_features,
        "columns":    list(res_df.columns),
        "index":      list(res_df.index),
        "data":       records,
        "nan_info":   nan_info,
        "corr_pairs": corr_pairs,
        "saved":      save_as,
    })


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


def _extract_success_payload(df, cfg, log_text):
    return {
        "n_mols": int(df.shape[0]),
        "n_features": int(df.shape[1]),
        "n_files": len(_list_feather_files(cfg["feather_dir"])),
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
    cfg["root_dir"] = cfg.get("root_dir") or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
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


def _parse_pasted_outputs(text):
    """Parse pasted targets: one number per line, comma-separated numbers, or name,value."""
    text = (text or "").strip()
    if not text:
        return []
    lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
    pairs = []
    if len(lines) == 1 and "," in lines[0]:
        parts = [p.strip() for p in lines[0].split(",") if p.strip()]
        try:
            return [(None, float(p)) for p in parts]
        except ValueError:
            pass
    for line in lines:
        if "," in line:
            name, _, rest = line.partition(",")
            try:
                pairs.append((name.strip(), float(rest.strip())))
            except ValueError:
                continue
        else:
            try:
                pairs.append((None, float(line)))
            except ValueError:
                continue
    return pairs


def _merge_outputs_into_df(df, pairs, new_col_name):
    """Add or overwrite new_col_name from (name, value) or (None, value) pairs."""
    import pandas as pd

    out = df.copy()
    if not pairs:
        return out
    by_name = all(p[0] is not None for p in pairs)
    if by_name:
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
        value_map = {str(n): v for n, v in pairs}
        out[new_col_name] = out[name_col].astype(str).map(value_map)
    else:
        values = [v for _, v in pairs]
        series = pd.Series(values[: len(out)], index=out.index[: len(values)])
        out[new_col_name] = series.reindex(out.index).values
    return out


def _linear_combo_search(df, y_name, min_features, max_features, top_n, threshold):
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

    for k in range(min_features, max_features + 1):
        for combo in combinations(range(n_feat), k):
            Xm = Xs[:, combo]
            model = LinearRegression()
            model.fit(Xm, ys)
            pred = model.predict(Xm)
            r2 = float(r2_score(ys, pred))
            if r2 < threshold:
                continue
            p = k
            adj = 1.0 - (1.0 - r2) * (n - 1) / max(n - p - 1, 1)
            rows.append({
                "combination": ", ".join(cols[i] for i in combo),
                "n_features": k,
                "r2": round(r2, 4),
                "adj_r2": round(float(adj), 4),
                "q2": round(float(_q2(Xm, ys)), 4),
            })

    if not rows:
        raise ValueError(
            f"No models reached R² ≥ {threshold}. Lower the threshold or paste a numeric target."
        )
    out = pd.DataFrame(rows).sort_values(["q2", "r2"], ascending=False).head(int(top_n))
    return out.reset_index(drop=True)


@app.route("/model/run", methods=["POST"])
def model_run():
    """Fit linear regression on the extracted CSV (chosen or pasted target)."""
    import io

    import pandas as pd

    data = request.json or {}
    csv_text = data.get("csv") or ""
    if not str(csv_text).strip():
        return jsonify({"error": "Extract features first (no CSV on the page)."}), 400

    try:
        df = pd.read_csv(io.StringIO(csv_text))
    except Exception as exc:
        return jsonify({"error": f"Could not parse CSV: {exc}"}), 400
    if df.empty:
        return jsonify({"error": "CSV is empty."}), 400

    pasted = data.get("outputs") or data.get("pasted") or ""
    y_name = (data.get("y_column") or data.get("target") or "").strip()
    new_col = (data.get("new_column") or "output").strip() or "output"
    if str(pasted).strip():
        pairs = _parse_pasted_outputs(pasted)
        if not pairs:
            return jsonify({
                "error": "Could not parse pasted values. Use one number per line or name,value.",
            }), 400
        y_name = new_col
        df = _merge_outputs_into_df(df, pairs, y_name)
    elif not y_name:
        return jsonify({
            "error": "Choose a target column from the CSV, or paste output values.",
        }), 400

    try:
        min_f = int(data.get("min_features", 1) or 1)
        max_f = int(data.get("max_features", 2) or 2)
        top_n = int(data.get("top_n", 8) or 8)
        threshold = float(data.get("threshold", 0.2) or 0.0)
        results = _linear_combo_search(df, y_name, min_f, max_f, top_n, threshold)
    except ValueError as exc:
        return jsonify({"error": str(exc)}), 400
    except Exception as exc:
        traceback.print_exc()
        return jsonify({"error": f"Linear regression failed: {exc}"}), 500

    return jsonify({
        "n_models": int(len(results)),
        "target": y_name,
        "n_rows": int(df.shape[0]),
        "columns": [str(c) for c in results.columns],
        "rows": results.to_dict(orient="records"),
        "csv": results.to_csv(index=False),
    })


def _pd_from_numpy(arr, index, columns):
    """Helper: rebuild DataFrame from numpy array (avoids pandas import at top)."""
    import pandas as pd
    return pd.DataFrame(arr, index=index, columns=columns)


GUI_API = 2


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
    if payload.get("api") == GUI_API:
        return True
    return "example_presets" in payload and "example_reference_name" in payload


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


def _pid_listening_on_port(port):
    pids = _pids_listening_on_port(port)
    return pids[0] if pids else None


def _reuse_or_refuse_existing_gui(host, port, open_browser):
    """If 7432 is already taken, reuse a current GUI or tell the user to stop the old one."""
    import webbrowser

    url = f"http://127.0.0.1:{port}/visual"
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
        print("  Stop every listener, then run descripytor visual again.")
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


def serve(host=None, port=None, open_browser=True):
    """Start the Flask GUI. Used by `descripytor visual` and `__main__`."""
    import threading
    import time
    import urllib.request
    import webbrowser

    host = host or os.environ.get("GUI_HOST", "127.0.0.1")
    port = int(port or os.environ.get("GUI_PORT", str(PORT)))
    existing = _reuse_or_refuse_existing_gui(host, port, open_browser)
    if existing is True:
        return
    if existing is False:
        sys.exit(1)
    url = f"http://127.0.0.1:{port}/visual"
    print(f"\n  DescriPyTor GUI")
    print(f"  Open:  {url}")
    print(f"  Forms: http://127.0.0.1:{port}/forms")
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
