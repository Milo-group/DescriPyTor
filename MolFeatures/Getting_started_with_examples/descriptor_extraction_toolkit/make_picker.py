"""
make_picker.py
==============

Embed a *real* molecule from your set into `atom_picker.html` and open it in
the browser, so you click atoms on the actual structure you'll extract from.

The picker then exports a complete `run_config.json` (dataset paths baked in +
your atom selections) that `descriptor_extractor.py` runs as-is.

Examples
--------
    # From a feather set (DescriPytor), pick which molecule to display:
    python make_picker.py --feather-dir path/to/feathers --index 0

    # From one .xyz file or a folder of xyz:
    python make_picker.py --xyz path/to/mol.xyz
    python make_picker.py --xyz-dir path/to/xyz --name conf_2

From a notebook where you already have `mols` loaded:
    from make_picker import launch_from_molecule
    launch_from_molecule(mols.molecules[0])      # writes + opens picker
"""

from __future__ import annotations

import argparse
import os
import sys
import webbrowser
from glob import glob
from pathlib import Path

HERE = Path(__file__).resolve().parent
TEMPLATE = HERE / "atom_picker.html"


# --------------------------------------------------------------------------- #
# xyz helpers
# --------------------------------------------------------------------------- #
def molecule_to_xyz(mol) -> str:
    """
    Turn a DescriPytor Molecule (with .xyz_df: atom,x,y,z) into an XYZ block.
    1-indexed atom order is preserved, so picker indices == DescriPytor indices.
    """
    df = mol.xyz_df[["atom", "x", "y", "z"]]
    name = getattr(mol, "molecule_name", "molecule")
    lines = [str(len(df)), str(name)]
    for _, r in df.iterrows():
        lines.append(f"{r['atom']:2s} {float(r['x']):14.8f} {float(r['y']):14.8f} {float(r['z']):14.8f}")
    return "\n".join(lines)


def export_molecules_xyz(mols, out_dir: str) -> str:
    """Export every Molecule.xyz_df in a Molecules set to .xyz files."""
    import re

    out = Path(out_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    count = 0
    for i, mol in enumerate(mols.molecules):
        name = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(getattr(mol, "molecule_name", f"mol_{i}")))
        (out / f"{name}.xyz").write_text(molecule_to_xyz(mol), encoding="utf-8")
        count += 1
    print(f"[picker] exported {count} xyz files -> {out}")
    return str(out)


def read_xyz_text(path: str) -> str:
    with open(path, "r", encoding="utf-8") as f:
        return f.read().strip()


def molecule_payload(mol) -> dict:
    """
    Pull dipole vector + vibrational modes off a DescriPytor Molecule so the
    picker can draw the dipole arrow and animate normal modes.

    Returns {"dipole": [x,y,z] | None,
             "modes": [{"freq": float, "ir": float, "disp": [[dx,dy,dz], ...]}, ...]}
    Anything missing degrades gracefully to None / [].
    """
    import numpy as np

    payload = {"dipole": None, "modes": []}

    # -- dipole vector (gauss_dipole_df: dip_x, dip_y, dip_z) ----------------
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
            payload["dipole"] = {
                "vector": [dx, dy, dz],
                "total": float(total) if total is not None else float(np.linalg.norm([dx, dy, dz])),
            }
    except Exception as e:
        print(f"[picker] no dipole data: {e}")

    # -- vibrational normal modes -------------------------------------------
    try:
        vmd = getattr(mol, "vibration_mode_dict", None) or {}
        info = getattr(mol, "info_df", None)
        ir_for = {}
        if info is not None and "Frequency" in getattr(info, "columns", []):
            ir_col = "IR" if "IR" in info.columns else None
            for _, r in info.iterrows():
                try:
                    ir_for[round(float(r["Frequency"]), 2)] = float(r[ir_col]) if ir_col else None
                except Exception:
                    pass
        for freq, disp in vmd.items():
            arr = np.asarray(disp, dtype=float)
            if arr.ndim != 2 or arr.shape[1] != 3:
                continue
            try:
                f = float(freq)
            except Exception:
                continue
            payload["modes"].append({
                "freq": f,
                "ir": ir_for.get(round(f, 2)),
                "disp": [[float(x) for x in v] for v in arr],
            })
        payload["modes"].sort(key=lambda m: m["freq"])
    except Exception as e:
        print(f"[picker] no vibration data: {e}")

    return payload


def _js_escape(s: str) -> str:
    """Make a string safe inside a JS backtick template literal."""
    return s.replace("\\", "\\\\").replace("`", "\\`").replace("$", "\\$")


def write_picker(xyz_text, out_html=None, open_browser=True, paths=None, mol_data=None):
    """
    Inject xyz_text (and, optionally, dataset paths + molecule data) into the
    template and write a ready-to-open HTML file. `mol_data` (from
    molecule_payload) enables the dipole arrow and vibration animation.
    """
    import json as _json
    template = TEMPLATE.read_text(encoding="utf-8")
    html = template.replace("__XYZ_DATA__", _js_escape(xyz_text))
    html = html.replace("__MOL_DATA__", _js_escape(_json.dumps(mol_data or {"dipole": None, "modes": []})))

    paths = paths or {}
    for placeholder, key in (
        ("__ROOT_DIR__", "root_dir"),
        ("__FEATHER_DIR__", "feather_dir"),
        ("__XYZ_DIR__", "xyz_dir"),
        ("__OUTPUT_CSV__", "output_csv"),
    ):
        v = paths.get(key)
        if v:  # only replace when we actually have a value; else keep default
            html = html.replace(placeholder, _js_escape(str(v)))

    out_html = out_html or str(HERE / "atom_picker_loaded.html")
    Path(out_html).write_text(html, encoding="utf-8")
    print(f"[picker] wrote {out_html}")
    if open_browser:
        webbrowser.open(Path(out_html).resolve().as_uri())
    return out_html


# --------------------------------------------------------------------------- #
# notebook-friendly entry points
# --------------------------------------------------------------------------- #
def launch_from_molecule(mol, out_html=None, open_browser=True, paths=None):
    return write_picker(molecule_to_xyz(mol), out_html, open_browser, paths,
                        mol_data=molecule_payload(mol))


def launch_from_xyz(path, out_html=None, open_browser=True, paths=None):
    return write_picker(read_xyz_text(path), out_html, open_browser, paths)


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def _resolve_root(root_dir):
    root = root_dir or os.environ.get(
        "DESCRIPYTOR_ROOT",
        r"C:\Users\edens\Desktop\DescriPytor\DescriPyTor-main\MolFeatures",
    )
    for sub in ("", "M2_data_extractor", "MolAlign"):
        p = os.path.join(root, sub) if sub else root
        if os.path.isdir(p) and p not in sys.path:
            sys.path.append(p)
    return root


def _default_csv(folder):
    return os.path.join(os.path.abspath(folder), "merged_features.csv")


def _main(argv=None):
    ap = argparse.ArgumentParser(description="Embed a molecule into the 3D atom picker and open it.")
    ap.add_argument("--xyz", help="Single .xyz file to load.")
    ap.add_argument("--xyz-dir", help="Folder of .xyz files.")
    ap.add_argument("--name", help="With --xyz-dir: file stem to pick (e.g. conf_2). Default: first.")
    ap.add_argument("--feather-dir", help="DescriPytor feather directory.")
    ap.add_argument("--index", type=int, default=0, help="With --feather-dir: molecule index (default 0).")
    ap.add_argument("--export-xyz-dir", help="With --feather-dir: also export all feather molecules to this xyz directory.")
    ap.add_argument("--root-dir", help="MolFeatures root (for feather loading).")
    ap.add_argument("--out", help="Output HTML path.")
    ap.add_argument("--no-open", action="store_true", help="Write the HTML but don't open a browser.")
    args = ap.parse_args(argv)

    open_browser = not args.no_open
    root = _resolve_root(args.root_dir)

    if args.xyz:
        folder = os.path.dirname(os.path.abspath(args.xyz))
        paths = {"root_dir": root, "xyz_dir": folder, "output_csv": _default_csv(folder)}
        launch_from_xyz(args.xyz, args.out, open_browser, paths)
        return

    if args.xyz_dir:
        files = sorted(glob(os.path.join(args.xyz_dir, "*.xyz")))
        if not files:
            ap.error(f"no .xyz files in {args.xyz_dir}")
        chosen = files[0]
        if args.name:
            hit = [f for f in files if Path(f).stem == args.name]
            if not hit:
                ap.error(f"no file named {args.name}.xyz in {args.xyz_dir}")
            chosen = hit[0]
        print(f"[picker] using {chosen}")
        paths = {"root_dir": root, "xyz_dir": os.path.abspath(args.xyz_dir),
                 "output_csv": _default_csv(args.xyz_dir)}
        launch_from_xyz(chosen, args.out, open_browser, paths)
        return

    if args.feather_dir:
        from data_extractor import Molecules
        mols = Molecules(args.feather_dir)
        xyz_dir = None
        if args.export_xyz_dir:
            xyz_dir = export_molecules_xyz(mols, args.export_xyz_dir)
        mol = mols.molecules[args.index]
        print(f"[picker] using molecule[{args.index}] = {getattr(mol, 'molecule_name', '?')}")
        paths = {"root_dir": root, "feather_dir": os.path.abspath(args.feather_dir),
                 "output_csv": _default_csv(args.feather_dir)}
        if xyz_dir:
            paths["xyz_dir"] = xyz_dir
        launch_from_molecule(mol, args.out, open_browser, paths)
        return

    ap.error("provide one of --xyz / --xyz-dir / --feather-dir")


if __name__ == "__main__":
    _main()
