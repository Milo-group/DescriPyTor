"""Extract geometry descriptors from the small alcohol / Ni XYZ set.

Run from the repository root::

    python tests/data/small_set/extract.py

Writes ``modeling_table.csv`` next to this file (name, Sterimol, bonds, %Vbur,
optional topology, boiling point). That table is what the modeling tests load.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from M2_data_extractor.metal_complex import (
    PKG_CPK,
    MetalComplex,
    adjacency,
    buried_volume,
    nob_types,
    sterimol,
)
from M2_data_extractor.xyz_io import XYZEnsemble

HERE = Path(__file__).resolve().parent
XYZ_DIR = HERE / "xyz"
TARGETS = HERE / "targets.csv"
TABLE = HERE / "modeling_table.csv"


def _angle(a, b, c) -> float:
    v1, v2 = a - b, c - b
    cos = np.clip(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)), -1, 1)
    return float(np.degrees(np.arccos(cos)))


def extract_alcohol(path: Path) -> dict:
    ens = XYZEnsemble(path)
    frame = ens.frames[0]
    symbols, coords = frame.symbols, frame.coords
    if symbols[0] != "C" or symbols[1] != "O" or symbols[2] != "H":
        raise ValueError(f"{path.name}: expected atom order C, O, H(O), ...")
    adj = adjacency(symbols, coords)
    radii = [PKG_CPK.get(t, 1.92) for t in nob_types(symbols, adj)]
    s = sterimol(symbols, coords, 2, 1, radii)
    if not s:
        raise RuntimeError(f"Sterimol failed for {path.name}")
    return {
        "B1": s["B1"],
        "B5": s["B5"],
        "L": s["L"],
        "sterimol_angle": s["angle"],
        "CO": float(np.linalg.norm(coords[0] - coords[1])),
        "HOC": _angle(coords[2], coords[1], coords[0]),
        "vbur_3.5": buried_volume(symbols, coords, coords[1], radius=3.5),
        "n_atoms": len(symbols),
        "n_heavy": int(sum(sym != "H" for sym in symbols)),
    }


def extract_table(with_topology: bool = True) -> pd.DataFrame:
    targets = pd.read_csv(TARGETS)
    rows = []
    for rec in targets.itertuples(index=False):
        path = XYZ_DIR / f"{rec.name}.xyz"
        row = {"name": rec.name, "smiles": rec.smiles, "bp_c": rec.bp_c}
        row.update(extract_alcohol(path))
        rows.append(row)
    table = pd.DataFrame(rows).set_index("name")
    if with_topology:
        try:
            from M2_data_extractor.ligand_topology import LigandTopology
        except ImportError:
            return table
        topo_rows = {}
        for name, smiles in table["smiles"].items():
            try:
                topo_rows[name] = LigandTopology.from_smiles(smiles, name=name).size_normalized_features()
            except Exception:
                continue
        if topo_rows:
            topo = pd.DataFrame.from_dict(topo_rows, orient="index")
            keep = [c for c in ("tf_kappa1", "tf_kappa2", "tf_zagreb1", "tf_chi1") if c in topo.columns]
            table = table.join(topo[keep])
    return table


def extract_tiny_ni() -> dict:
    path = XYZ_DIR / "tiny_nih2n2.xyz"
    mc = MetalComplex.from_xyz(path)
    geom = mc.geometric_features()
    geom["name"] = "tiny_nih2n2"
    geom["n_atoms"] = len(mc.symbols)
    return geom


def main():
    table = extract_table(with_topology=True)
    table.to_csv(TABLE)
    print(f"wrote {TABLE}  shape={table.shape}")
    print(table.round(3).to_string())
    ni = extract_tiny_ni()
    print(
        "tiny_nih2n2 bite={:.2f} MD_mean={:.3f} vbur_3.5={:.2f}".format(
            ni["bite"], ni["MD_mean"], ni["vbur_3.5"]
        )
    )


if __name__ == "__main__":
    main()
