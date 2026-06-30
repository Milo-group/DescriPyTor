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

# ── Flask ─────────────────────────────────────────────────────
try:
    from flask import Flask, request, jsonify
    from flask_cors import CORS
except ImportError:
    print("Missing dependencies. Run:  pip install flask flask-cors")
    sys.exit(1)

app = Flask(__name__)
CORS(app)

PORT = 7432

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
    return jsonify({"ok": True, "version": "1.0"})


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
    print(f"  Press Ctrl+C to stop\n")
    app.run(host="127.0.0.1", port=PORT, debug=False)
