"""
gui_server.py — local backend for feature_extraction_gui.html

Run once:
    python gui_server.py

Then open feature_extraction_gui.html in your browser.
The GUI detects the server and enables live computation buttons.

Dependencies:
    pip install flask flask-cors
    (all other deps come from DescriPytor itself)
"""

import sys
import os
import json
import traceback
import math

# ── Flask ─────────────────────────────────────────────────────
try:
    from flask import Flask, request, jsonify, send_file
    from flask_cors import CORS
except ImportError:
    print("Missing dependencies. Run:  pip install flask flask-cors")
    sys.exit(1)

app = Flask(__name__)
CORS(app)

PORT = int(os.environ.get("GUI_PORT", "7432"))

# ── path helper ───────────────────────────────────────────────
def ensure_path(root: str):
    """Add MolFeatures root to sys.path so DescriPytor imports work."""
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

@app.route("/status")
def status():
    return jsonify({
        "ok": True,
        "version": "1.1",
        "cpu_count": os.cpu_count() or 1,
        "bassa_available": _module_available("bassa_reg"),
    })


@app.route("/")
def gui():
    """Serve the GUI from the same origin as the API."""
    return send_file(os.path.join(os.path.dirname(__file__), "feature_extraction_gui.html"))


@app.route("/modes")
def modes_viewer():
    """Serve the py3Dmol TS / vibration checker."""
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    return send_file(os.path.join(root, "utils", "mode_viewer.html"))


@app.route("/api/parse-modes", methods=["POST"])
def parse_modes_api():
    """Parse an uploaded .hess / ORCA .out / .xyz into a 3Dmol payload."""
    data = request.json or {}
    try:
        root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        if root not in sys.path:
            sys.path.insert(0, root)
        from utils.mode_viewer import molecule_payload, parse_text
        mol = parse_text(data.get("name") or "upload", data.get("text") or "")
        return jsonify(molecule_payload(mol))
    except Exception as e:
        traceback.print_exc()
        return jsonify({"error": str(e)}), 400


@app.route("/visual")
def visual_gui():
    """Serve the combined atom-picker, extraction, and modeling workflow."""
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    toolkit = os.path.join(
        root, "Getting_started_with_examples", "descriptor_extraction_toolkit",
    )
    loaded = os.path.join(toolkit, "atom_picker_loaded.html")
    template = os.path.join(toolkit, "atom_picker.html")
    return send_file(loaded if os.path.isfile(loaded) else template)


@app.route("/theta")
def theta_explorer():
    """Serve the live Sterimol theta explorer (sources in theta_explorer/)."""
    ensure_path("")
    from theta_explorer.build import assemble
    return assemble(), 200, {"Content-Type": "text/html; charset=utf-8"}


@app.route("/theta/check", methods=["POST"])
def theta_check():
    """The Python Sterimol on the explorer's coordinates and axis, so the page
    can confirm its in-browser port against the package on every pick."""
    data = request.json or {}
    try:
        ensure_path(data.get("root", ""))
        import io
        import pandas as pd
        from data_extractor import extract_connectivity, get_sterimol_df

        lines = (data.get("xyz") or "").strip().splitlines()
        if lines and lines[0].strip().isdigit():
            lines = lines[2:]
        xyz = pd.read_csv(io.StringIO("\n".join(lines)), sep=r"\s+", header=None,
                          names=["atom", "x", "y", "z"], usecols=range(4))
        bonds = extract_connectivity(xyz, threshold_distance=1.82)
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
        from theta_explorer.data import payload
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


def _json_value(value):
    """Convert pandas/NumPy values into strict JSON-compatible values."""
    import numpy as np
    import pandas as pd

    if value is None or value is pd.NA:
        return None
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        return float(value) if math.isfinite(float(value)) else None
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, (list, tuple, set, np.ndarray)):
        return [_json_value(item) for item in value]
    if isinstance(value, dict):
        return {str(key): _json_value(item) for key, item in value.items()}
    return str(value) if not isinstance(value, (str, int)) else value


def _frame_records(frame, limit: int = 100):
    if frame is None:
        return []
    return [
        {str(key): _json_value(value) for key, value in row.items()}
        for row in frame.head(limit).to_dict(orient="records")
    ]


@app.route("/model/merge-outputs", methods=["POST"])
def model_merge_outputs():
    """Merge pasted target values into a descriptor CSV and save a new CSV."""
    import re
    import pandas as pd

    data = request.json or {}
    source = str(data.get("features_csv", "")).strip()
    target = str(data.get("target", "output")).strip() or "output"
    raw = str(data.get("outputs", "")).strip()
    mode = str(data.get("mode", "row_order")).strip()
    name_column = str(data.get("name_column", "")).strip()

    if not source or not os.path.isfile(source):
        return jsonify({"error": f"CSV file not found: {source!r}"}), 400
    if not raw:
        return jsonify({"error": "Paste at least one output value."}), 400

    frame = pd.read_csv(source)
    try:
        if mode == "name_value":
            if not name_column or name_column not in frame.columns:
                return jsonify({"error": f"Name column {name_column!r} was not found."}), 400
            mapping = {}
            for line_number, line in enumerate(raw.splitlines(), start=1):
                if not line.strip():
                    continue
                parts = [item.strip() for item in re.split(r"[\t,;]", line, maxsplit=1)]
                if len(parts) != 2 or not parts[0]:
                    return jsonify({"error": f"Line {line_number} must be: name,value"}), 400
                mapping[parts[0]] = float(parts[1])
            names = frame[name_column].astype(str)
            missing = names[~names.isin(mapping)].tolist()
            if missing:
                return jsonify({
                    "error": f"No output supplied for {len(missing)} rows.",
                    "missing_names": missing[:20],
                }), 400
            frame[target] = names.map(mapping).astype(float)
        else:
            tokens = [item for item in re.split(r"[\s,;]+", raw) if item]
            values = [float(item) for item in tokens]
            if len(values) != len(frame):
                return jsonify({
                    "error": f"Received {len(values)} outputs for {len(frame)} CSV rows."
                }), 400
            frame[target] = values
    except ValueError as exc:
        return jsonify({"error": f"Outputs must be numeric: {exc}"}), 400

    requested_output = str(data.get("output_csv", "")).strip()
    if requested_output:
        output = os.path.abspath(requested_output)
    else:
        stem, extension = os.path.splitext(os.path.abspath(source))
        output = f"{stem}_with_{target}{extension or '.csv'}"
    os.makedirs(os.path.dirname(output), exist_ok=True)
    frame.to_csv(output, index=False)
    return jsonify({
        "ok": True,
        "output_csv": output,
        "rows": len(frame),
        "target": target,
        "preview": _frame_records(frame[[target]], 5),
    })


@app.route("/model/dataset", methods=["POST"])
def model_dataset():
    """Inspect a modeling CSV before starting an expensive search."""
    import pandas as pd

    data = request.json or {}
    path = data.get("features_csv", "")
    target = data.get("target", "")
    if not path or not os.path.isfile(path):
        return jsonify({"error": f"CSV file not found: {path!r}"}), 400

    try:
        frame = pd.read_csv(path)
        numeric = frame.select_dtypes(include="number").columns.tolist()
        candidate_features = [column for column in numeric if column != target]
        target_found = bool(target and target in frame.columns)
        target_numeric = bool(target_found and target in numeric)
        target_missing = int(frame[target].isna().sum()) if target_found else None
        duplicate_columns = frame.columns[frame.columns.duplicated()].tolist()
        ready = target_found and target_numeric and bool(candidate_features) and not duplicate_columns
        return jsonify({
            "rows": len(frame),
            "columns": frame.columns.tolist(),
            "numeric_columns": numeric,
            "feature_count": len(candidate_features),
            "target": target,
            "target_found": target_found,
            "target_numeric": target_numeric,
            "target_missing": target_missing,
            "duplicate_columns": duplicate_columns,
            "ready": ready,
            "preview": _frame_records(frame, 5),
        })
    except Exception as exc:
        return jsonify({"error": str(exc)}), 400


@app.route("/model/run", methods=["POST"])
def model_run():
    """Run M3 model search or BASSA from the browser GUI."""
    data = request.json or {}
    engine = str(data.get("engine", "m3")).lower()
    features_csv = data.get("features_csv", "")
    target_csv = data.get("target_csv") or ""
    target = str(data.get("target", "output")).strip()

    if not features_csv or not os.path.isfile(features_csv):
        return jsonify({"error": f"CSV file not found: {features_csv!r}"}), 400
    if not target:
        return jsonify({"error": "A target column is required."}), 400

    try:
        ensure_path(data.get("root", ""))
        if engine == "bassa":
            return _run_bassa(data, features_csv, target)
        return _run_m3(data, features_csv, target_csv, target)
    except Exception as exc:
        traceback.print_exc()
        return jsonify({"error": str(exc)}), 500


def _run_m3(data, features_csv: str, target_csv: str, target: str):
    try:
        from M3_modeler.modeling import ClassificationModel, LinearRegressionModel
    except ImportError:
        from MolFeatures.M3_modeler.modeling import ClassificationModel, LinearRegressionModel

    task = str(data.get("task", "regression")).lower()
    process_method = "two csvs" if target_csv else "one csv"
    paths = {
        "features_csv_filepath": features_csv,
        "target_csv_filepath": target_csv or None,
    }
    common = dict(
        csv_filepaths=paths,
        process_method=process_method,
        y_value=target,
        names_column=data.get("names_column") or None,
        leave_out=data.get("leave_out") or None,
        min_features_num=int(data.get("min_features", 1)),
        max_features_num=int(data.get("max_features", 3)),
        n_splits=int(data.get("n_splits", 5)),
        db_path=data.get("db_path") or "results",
    )
    if task == "classification":
        model = ClassificationModel(
            **common,
            ordinal=bool(data.get("ordinal", False)),
        )
    else:
        model = LinearRegressionModel(
            **common,
            model_type=data.get("regression_type", "linear"),
            alpha=float(data.get("alpha", 1.0)),
            scale=bool(data.get("scale", True)),
            seed=int(data.get("seed", 42)),
        )

    results = model.search_models(
        top_n=int(data.get("top_n", 20)),
        n_jobs=int(data.get("n_jobs", 1)),
        threshold=float(data.get("threshold", 0.7)),
        bool_parallel=bool(data.get("parallel", False)),
        required_features=data.get("required_features") or None,
    )
    return jsonify({
        "engine": "m3",
        "task": task,
        "rows": len(results),
        "columns": list(results.columns),
        "results": _frame_records(results),
        "run_directory": str(model.paths.root),
    })


def _run_bassa(data, features_csv: str, target: str):
    if not _module_available("bassa_reg"):
        return jsonify({
            "error": "BASSA is not installed in this Python environment. Run: pip install bassa-reg"
        }), 400

    import pandas as pd
    from bassa_reg import Bassa
    from bassa_reg.spike_and_slab.spike_and_slab import (
        SpikeAndSlabConfigurations,
        SpikeAndSlabRegression,
    )
    from bassa_reg.spike_and_slab.spike_and_slab_util_models import SpikeAndSlabPriors

    frame = pd.read_csv(features_csv)
    if target not in frame.columns:
        return jsonify({"error": f"Target column {target!r} was not found."}), 400

    excluded = {target, data.get("names_column") or ""}
    requested = data.get("feature_columns") or []
    if requested:
        feature_columns = [column for column in requested if column in frame.columns]
    else:
        feature_columns = [
            column for column in frame.select_dtypes(include="number").columns
            if column not in excluded
        ]
    if not feature_columns:
        return jsonify({"error": "No numeric feature columns were found."}), 400

    clean = frame[feature_columns + [target]].dropna()
    project_path = os.path.abspath(data.get("project_path") or "bassa_runs")
    os.makedirs(project_path, exist_ok=True)
    config = SpikeAndSlabConfigurations(
        sampler_iterations=int(data.get("sampler_iterations", 5000))
    )
    regression = SpikeAndSlabRegression(
        x=clean[feature_columns],
        y=clean[target],
        priors=SpikeAndSlabPriors(),
        config=config,
        project_path=project_path,
        experiment_name=data.get("experiment_name") or "molfeatures_gui",
    )
    regression.run()
    Bassa(model=regression).run()
    return jsonify({
        "engine": "bassa",
        "rows": len(clean),
        "features": feature_columns,
        "project_path": project_path,
        "message": "BASSA completed. Plots and summary files were written to the project path.",
    })


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
        get_extended_df_for_sterimol,
        preform_coordination_transformation,
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

    extended_df = get_extended_df_for_sterimol(
        mol.coordinates_df, mol.bonds_df, radii=radii
    )
    rotated_df, rotated_plane = preform_coordination_transformation(
        extended_df, indices=base_atoms
    )

    fig_endon = plot_b1_visualization(
        rotated_plane, rotated_df,
        sterimol_df=st_df, n_points=n_points, title=endon_title,
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
      root            – MolFeatures root to add to sys.path (optional)
      entry_widgets   – dict of string values: {ring, stretching, stretch,
                        upper_stretch, bending, bend, sub_atoms, npa, dipole,
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


def _pd_from_numpy(arr, index, columns):
    """Helper: rebuild DataFrame from numpy array (avoids pandas import at top)."""
    import pandas as pd
    return pd.DataFrame(arr, index=index, columns=columns)


# ── run ───────────────────────────────────────────────────────
if __name__ == "__main__":
    print(f"\n  DescriPytor GUI server")
    print(f"  Listening on http://localhost:{PORT}")
    print(f"  Open feature_extraction_gui.html in your browser")
    print(f"  TS / vibration viewer: http://localhost:{PORT}/modes")
    print(f"  Press Ctrl+C to stop\n")
    app.run(host=os.environ.get("GUI_HOST", "127.0.0.1"), port=PORT, debug=False)
