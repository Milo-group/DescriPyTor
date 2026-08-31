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


def atom_picker_template() -> Path:
    """Canonical picker HTML: packaged M2 copy, else this toolkit folder."""
    m2 = HERE.parents[2] / "M2_data_extractor" / "atom_picker.html"
    local = HERE / "atom_picker.html"
    if m2.is_file():
        return m2
    if local.is_file():
        return local
    raise FileNotFoundError(
        "atom_picker.html not found in M2_data_extractor or this toolkit folder."
    )


TEMPLATE = atom_picker_template()


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


# --------------------------------------------------------------------------- #
# 2D bond-graph (igraph) - an atom-numbered diagram to read off indices from,
# meant to sit next to the 3D picker rather than replace it.
# --------------------------------------------------------------------------- #
COVALENT_RADII = {  # Angstrom, Cordero et al. (rounded) - good enough for a bond-perception cutoff
    "H": 0.31, "HE": 0.28, "LI": 1.28, "BE": 0.96, "B": 0.84, "C": 0.76, "N": 0.71,
    "O": 0.66, "F": 0.57, "NE": 0.58, "NA": 1.66, "MG": 1.41, "AL": 1.21, "SI": 1.11,
    "P": 1.07, "S": 1.05, "CL": 1.02, "AR": 1.06, "K": 2.03, "CA": 1.76, "SC": 1.70,
    "TI": 1.60, "V": 1.53, "CR": 1.39, "MN": 1.39, "FE": 1.32, "CO": 1.26, "NI": 1.24,
    "CU": 1.32, "ZN": 1.22, "GA": 1.22, "GE": 1.20, "AS": 1.19, "SE": 1.20, "BR": 1.20,
    "KR": 1.16, "RB": 2.20, "SR": 1.95, "PD": 1.39, "AG": 1.45, "CD": 1.44, "SN": 1.39,
    "SB": 1.39, "I": 1.39, "XE": 1.40, "CS": 2.44, "BA": 2.15, "PT": 1.36, "AU": 1.36,
    "HG": 1.32, "PB": 1.46,
}
DEFAULT_COVALENT_RADIUS = 0.75


def _parse_xyz_atoms(xyz_text: str):
    """Return (elements, coords) from an XYZ-block string."""
    lines = [ln for ln in xyz_text.splitlines() if ln.strip()]
    if not lines:
        return [], []
    try:
        n = int(lines[0].split()[0])
    except (ValueError, IndexError):
        n = None
    body = lines[2:] if n is not None else lines
    elements, coords = [], []
    for ln in body:
        parts = ln.split()
        if len(parts) < 4:
            continue
        try:
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
        except ValueError:
            continue
        elements.append(parts[0])
        coords.append((x, y, z))
        if n is not None and len(elements) >= n:
            break
    return elements, coords


def _infer_bonds(elements, coords, tolerance: float = 1.3):
    """Pairwise distance + covalent-radius cutoff bond perception (0-indexed pairs)."""
    bonds = []
    radii = [COVALENT_RADII.get(el.strip().upper(), DEFAULT_COVALENT_RADIUS) for el in elements]
    n = len(coords)
    for i in range(n):
        xi, yi, zi = coords[i]
        for j in range(i + 1, n):
            xj, yj, zj = coords[j]
            d2 = (xi - xj) ** 2 + (yi - yj) ** 2 + (zi - zj) ** 2
            cutoff = (radii[i] + radii[j]) * tolerance
            if d2 <= cutoff * cutoff:
                bonds.append((i, j))
    return bonds


def graph2d_png_base64(xyz_text: str) -> str:
    """
    Build an atom-numbered 2D bond-graph from an XYZ block using igraph for the
    graph/layout and matplotlib for rendering, so it can sit next to the 3D
    picker as a quick "which number is which atom" reference.

    Returns a bare base64 PNG string (no `data:image/...` prefix), or "" if
    igraph/matplotlib aren't available or the xyz can't be parsed.
    """
    try:
        import base64
        import io
        import networkx as nx
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:  # noqa
        print(f"[picker] 2D graph view unavailable ({e}); needs networkx + matplotlib.")
        return ""

    elements, coords = _parse_xyz_atoms(xyz_text)
    if not elements:
        return ""

    bonds = _infer_bonds(elements, coords)
    g = nx.Graph()
    g.add_nodes_from(range(len(elements)))
    g.add_edges_from(bonds)
    pos = nx.kamada_kawai_layout(g) if bonds else nx.circular_layout(g)
    xs = [pos[i][0] for i in range(len(elements))]
    ys = [pos[i][1] for i in range(len(elements))]

    fig, ax = plt.subplots(figsize=(5.2, 5.2), dpi=120)
    for (i, j) in bonds:
        ax.plot([xs[i], xs[j]], [ys[i], ys[j]], color="#8a94a6", linewidth=1.6, zorder=1)
    ax.scatter(xs, ys, s=380, color="#eef2f8", edgecolors="#2a3644", linewidths=1.3, zorder=2)
    for idx, (x, y, el) in enumerate(zip(xs, ys, elements), start=1):
        ax.annotate(f"{idx}\n{el}", (x, y), ha="center", va="center", fontsize=8, fontweight="bold", zorder=3)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.tight_layout()

    buf = io.BytesIO()
    fig.savefig(buf, format="png", transparent=True)
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode("ascii")


# --------------------------------------------------------------------------- #
# RDKit structure + fingerprint helpers (2D depiction, scope grids, Morgan bits).
# Shared by the webapp's Tab 1 "structure & fingerprint viewer" and "Build scope
# page" button. Bond perception reuses the same rdDetermineBonds approach as
# feature_packages.extract_smiles, but works straight off an xyz-block string
# (no ase/file round-trip needed) and returns an actual Mol, not just a SMILES.
# --------------------------------------------------------------------------- #
def rdkit_mol_from_xyz(xyz_text: str, charges=(0, 1, -1, 2, -2), remove_hs: bool = True):
    """
    Perceive bonds/valences for an xyz-block string and return a sanitized
    RDKit Mol with 2D depiction coordinates, or None if perception fails or
    RDKit isn't installed.
    """
    try:
        from rdkit import Chem, RDLogger
        from rdkit.Chem import AllChem, rdDetermineBonds
        RDLogger.DisableLog("rdApp.*")
    except Exception as e:  # noqa
        print(f"[picker] rdkit unavailable ({e})")
        return None

    try:
        base = Chem.MolFromXYZBlock(xyz_text)
    except Exception as e:  # noqa
        print(f"[picker] MolFromXYZBlock failed: {e}")
        return None
    if base is None:
        return None

    mol = None
    last_err = None
    for charge in charges:
        try:
            import copy
            candidate = copy.deepcopy(base)
            rdDetermineBonds.DetermineBonds(candidate, charge=charge)
            Chem.SanitizeMol(candidate)
            mol = candidate
            break
        except Exception as e:  # noqa
            last_err = e
    if mol is None:
        print(f"[picker] bond perception failed for all charge guesses: {last_err}")
        return None

    if remove_hs:
        try:
            mol = Chem.RemoveHs(mol)
        except Exception:
            pass
    try:
        AllChem.Compute2DCoords(mol)
    except Exception as e:  # noqa
        print(f"[picker] 2D coordinate generation failed: {e}")
    return mol


def mol_structure_png_base64(mol, size=(320, 280), highlight_atoms=None) -> str:
    """2D structure depiction of an RDKit Mol as a bare base64 PNG string."""
    try:
        import base64
        import io
        from rdkit.Chem import Draw
    except Exception as e:  # noqa
        print(f"[picker] structure image unavailable ({e})")
        return ""
    try:
        img = Draw.MolToImage(mol, size=size, highlightAtoms=highlight_atoms or [])
        buf = io.BytesIO()
        img.save(buf, format="PNG")
        return base64.b64encode(buf.getvalue()).decode("ascii")
    except Exception as e:  # noqa
        print(f"[picker] structure image render failed: {e}")
        return ""


def morgan_fingerprint_bits(mol, n_bits: int = 256, radius: int = 2):
    """
    Compute an ECFP-style Morgan fingerprint and return
    (sorted list of "on" bit indices, bitInfo dict) for visualization.
    bitInfo maps bit_id -> [(atom_idx, radius), ...] as required by DrawMorganBit.
    """
    from rdkit.Chem import rdMolDescriptors

    bit_info = {}
    fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(
        mol, radius, nBits=n_bits, bitInfo=bit_info)
    on_bits = sorted(fp.GetOnBits())
    return on_bits, bit_info


def fingerprint_bit_png_base64(mol, bit_id: int, bit_info: dict, size=(300, 300)) -> str:
    """Render the substructure responsible for one Morgan fingerprint bit."""
    try:
        import base64
        import io
        from rdkit.Chem import Draw
    except Exception as e:  # noqa
        print(f"[picker] fingerprint bit image unavailable ({e})")
        return ""
    try:
        img = Draw.DrawMorganBit(mol, bit_id, bit_info, useSVG=False)
        if hasattr(img, "save"):  # PIL Image
            buf = io.BytesIO()
            img.save(buf, format="PNG")
            data = buf.getvalue()
        else:  # some rdkit versions return raw PNG bytes directly
            data = img
        return base64.b64encode(data).decode("ascii")
    except Exception as e:  # noqa
        print(f"[picker] fingerprint bit render failed for bit {bit_id}: {e}")
        return ""


def scope_grid_png_bytes(entries, mols_per_row: int = 4, sub_img_size=(260, 230)) -> bytes:
    """
    Build a "scope page"-style grid image: one 2D structure + legend per cell.

    entries: list of (mol, legend_str). Mols that are None are skipped.
    Returns PNG bytes, or b"" if nothing could be rendered.
    """
    from rdkit.Chem import Draw

    pairs = [(m, lg) for m, lg in entries if m is not None]
    if not pairs:
        return b""
    mols = [m for m, _ in pairs]
    legends = [lg for _, lg in pairs]
    img = Draw.MolsToGridImage(
        mols, molsPerRow=mols_per_row, subImgSize=sub_img_size,
        legends=legends, returnPNG=False,
    )
    import io
    buf = io.BytesIO()
    img.save(buf, format="PNG")
    return buf.getvalue()


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
    graph2d_b64 = graph2d_png_base64(xyz_text) if xyz_text else ""
    html = html.replace("__GRAPH2D_IMG__", graph2d_b64)
    html = html.replace("__GRAPH2D_DISPLAY__", "" if graph2d_b64 else "none")

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
        str(Path(__file__).resolve().parents[2]),
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
    ap.add_argument("--root-dir", help="DescriPyTor root (for feather loading).")
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
