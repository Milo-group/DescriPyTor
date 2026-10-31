"""Write the small-set XYZ files from Z-matrices (no RDKit).

Atom order is the same in every alcohol XYZ:

    1  C_alpha   (carbon bound to OH)
    2  O
    3  H of OH
    4+ remainder (first extra heavy atom is C_beta when present)

Sterimol of the alkyl group is then origin=O (2), direction=C (1).
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
XYZ_DIR = HERE / "xyz"

# Experimental boiling points (°C), 1 atm.
ALCOHOLS = [
    ("methanol", "CO", 64.7),
    ("ethanol", "CCO", 78.4),
    ("n_propanol", "CCCO", 97.2),
    ("isopropanol", "CC(C)O", 82.3),
    ("n_butanol", "CCCCO", 117.7),
    ("isobutanol", "CC(C)CO", 108.0),
    ("sec_butanol", "CCC(C)O", 99.5),
    ("tert_butanol", "CC(C)(C)O", 82.4),
]


def _unit(v):
    n = np.linalg.norm(v)
    if n < 1e-12:
        raise ValueError("zero vector")
    return v / n


def _place(coords, bond_i, dist, ang_i, ang_deg, dih_i=None, dih_deg=None):
    """Append one atom using a Z-matrix row (0-based references)."""
    coords = [np.asarray(c, dtype=float) for c in coords]
    n = len(coords)
    if n == 0:
        return [np.zeros(3)]
    if n == 1:
        return coords + [np.array([float(dist), 0.0, 0.0])]

    a = math.radians(ang_deg)
    A = coords[bond_i]
    B = coords[ang_i]
    ax = _unit(B - A)

    if n == 2 or dih_i is None:
        helper = np.array([0.0, 0.0, 1.0])
        if abs(np.dot(ax, helper)) > 0.9:
            helper = np.array([0.0, 1.0, 0.0])
        perp = _unit(np.cross(ax, helper))
        pos = A + dist * (math.cos(a) * ax + math.sin(a) * perp)
        return coords + [pos]

    C = coords[dih_i]
    d = math.radians(dih_deg)
    nrm = np.cross(ax, C - A)
    if np.linalg.norm(nrm) < 1e-8:
        nrm = np.cross(ax, np.array([0.0, 0.0, 1.0]))
        if np.linalg.norm(nrm) < 1e-8:
            nrm = np.cross(ax, np.array([0.0, 1.0, 0.0]))
    nrm = _unit(nrm)
    perp = _unit(np.cross(nrm, ax))
    ad = dist * (
        math.cos(a) * ax
        + math.sin(a) * (math.cos(d) * perp + math.sin(d) * nrm)
    )
    return coords + [A + ad]


def _h_on(c_idx, heavy_neighbor, existing, n_h, dist=1.09, angle=109.5):
    """Add tetrahedral hydrogens on carbon ``c_idx``."""
    coords = list(existing)
    # Seed dihedrals so successive H atoms stagger.
    dihedrals = [180.0, 60.0, -60.0, 0.0]
    added = 0
    ref = 0 if c_idx != 0 else 1
    for dih in dihedrals:
        if added >= n_h:
            break
        if ref == c_idx:
            ref = heavy_neighbor
        coords = _place(coords, c_idx, dist, heavy_neighbor, angle, ref, dih)
        added += 1
    return coords


def methanol():
    # C, O, H(O), H, H, H
    c = _place([], 0, 0, 0, 0)
    c = _place(c, 0, 1.43, 0, 0)
    c = _place(c, 1, 0.96, 0, 108.5)
    c = _place(c, 0, 1.09, 1, 109.5, 2, 180)
    c = _place(c, 0, 1.09, 1, 109.5, 2, 60)
    c = _place(c, 0, 1.09, 1, 109.5, 2, -60)
    return ["C", "O", "H", "H", "H", "H"], c


def ethanol():
    c = methanol()[1][:3]  # C, O, H(O)
    c = _place(c, 0, 1.52, 1, 109.5, 2, 180)  # C_beta
    c = _place(c, 0, 1.09, 1, 109.5, 3, 60)   # H_alpha
    c = _place(c, 0, 1.09, 1, 109.5, 3, -60)
    c = _place(c, 3, 1.09, 0, 109.5, 1, 180)  # H_beta
    c = _place(c, 3, 1.09, 0, 109.5, 1, 60)
    c = _place(c, 3, 1.09, 0, 109.5, 1, -60)
    return ["C", "O", "H", "C", "H", "H", "H", "H", "H"], c


def n_propanol():
    symbols, c = ethanol()
    # replace nothing: add C_gamma on C_beta (index 3)
    c = c[:4]  # C O H C_beta
    # rebuild H set after adding C_gamma
    c = _place(c, 3, 1.53, 0, 109.5, 1, 180)  # C_gamma
    c = _place(c, 0, 1.09, 1, 109.5, 3, 60)
    c = _place(c, 0, 1.09, 1, 109.5, 3, -60)
    c = _place(c, 3, 1.09, 0, 109.5, 4, 60)
    c = _place(c, 3, 1.09, 0, 109.5, 4, -60)
    c = _place(c, 4, 1.09, 3, 109.5, 0, 180)
    c = _place(c, 4, 1.09, 3, 109.5, 0, 60)
    c = _place(c, 4, 1.09, 3, 109.5, 0, -60)
    return ["C", "O", "H", "C", "C", "H", "H", "H", "H", "H", "H", "H"], c


def isopropanol():
    c = methanol()[1][:3]
    c = _place(c, 0, 1.52, 1, 109.5, 2, 180)   # C_a
    c = _place(c, 0, 1.52, 1, 109.5, 2, 60)    # C_b
    c = _place(c, 0, 1.09, 1, 109.5, 2, -60)   # H_alpha
    for carbon, other in ((3, 4), (4, 3)):
        c = _place(c, carbon, 1.09, 0, 109.5, other, 180)
        c = _place(c, carbon, 1.09, 0, 109.5, other, 60)
        c = _place(c, carbon, 1.09, 0, 109.5, other, -60)
    return ["C", "O", "H", "C", "C", "H", "H", "H", "H", "H", "H", "H"], c


def n_butanol():
    symbols, c = n_propanol()
    c = c[:5]  # C O H C_beta C_gamma
    c = _place(c, 4, 1.53, 3, 109.5, 0, 180)  # C_delta
    c = _place(c, 0, 1.09, 1, 109.5, 3, 60)
    c = _place(c, 0, 1.09, 1, 109.5, 3, -60)
    c = _place(c, 3, 1.09, 0, 109.5, 4, 60)
    c = _place(c, 3, 1.09, 0, 109.5, 4, -60)
    c = _place(c, 4, 1.09, 3, 109.5, 5, 60)
    c = _place(c, 4, 1.09, 3, 109.5, 5, -60)
    c = _place(c, 5, 1.09, 4, 109.5, 3, 180)
    c = _place(c, 5, 1.09, 4, 109.5, 3, 60)
    c = _place(c, 5, 1.09, 4, 109.5, 3, -60)
    sym = ["C", "O", "H", "C", "C", "C"] + ["H"] * 9
    return sym, c


def isobutanol():
    c = methanol()[1][:3]
    c = _place(c, 0, 1.52, 1, 109.5, 2, 180)  # C_beta
    c = _place(c, 3, 1.53, 0, 109.5, 1, 180)  # C_a
    c = _place(c, 3, 1.53, 0, 109.5, 1, 60)   # C_b
    c = _place(c, 0, 1.09, 1, 109.5, 3, 60)
    c = _place(c, 0, 1.09, 1, 109.5, 3, -60)
    c = _place(c, 3, 1.09, 0, 109.5, 4, -60)
    for carbon, other in ((4, 5), (5, 4)):
        c = _place(c, carbon, 1.09, 3, 109.5, other, 180)
        c = _place(c, carbon, 1.09, 3, 109.5, other, 60)
        c = _place(c, carbon, 1.09, 3, 109.5, other, -60)
    # C4H10O = 15 atoms
    return ["C", "O", "H", "C", "C", "C"] + ["H"] * 9, c


def sec_butanol():
    c = methanol()[1][:3]
    c = _place(c, 0, 1.52, 1, 109.5, 2, 180)  # C_methyl
    c = _place(c, 0, 1.52, 1, 109.5, 2, 60)   # C_ethyl
    c = _place(c, 4, 1.53, 0, 109.5, 1, 180)  # C_ethyl terminal
    c = _place(c, 0, 1.09, 1, 109.5, 2, -60)
    for carbon, other in ((3, 4),):
        c = _place(c, carbon, 1.09, 0, 109.5, other, 180)
        c = _place(c, carbon, 1.09, 0, 109.5, other, 60)
        c = _place(c, carbon, 1.09, 0, 109.5, other, -60)
    c = _place(c, 4, 1.09, 0, 109.5, 5, 60)
    c = _place(c, 4, 1.09, 0, 109.5, 5, -60)
    c = _place(c, 5, 1.09, 4, 109.5, 0, 180)
    c = _place(c, 5, 1.09, 4, 109.5, 0, 60)
    c = _place(c, 5, 1.09, 4, 109.5, 0, -60)
    # C4H10O = 15 atoms
    return ["C", "O", "H", "C", "C", "C"] + ["H"] * 9, c


def tert_butanol():
    c = methanol()[1][:3]
    c = _place(c, 0, 1.52, 1, 109.5, 2, 180)
    c = _place(c, 0, 1.52, 1, 109.5, 2, 60)
    c = _place(c, 0, 1.52, 1, 109.5, 2, -60)
    for carbon, other in ((3, 4), (4, 5), (5, 3)):
        c = _place(c, carbon, 1.09, 0, 109.5, other, 180)
        c = _place(c, carbon, 1.09, 0, 109.5, other, 60)
        c = _place(c, carbon, 1.09, 0, 109.5, other, -60)
    sym = ["C", "O", "H", "C", "C", "C"] + ["H"] * 9
    return sym, c


def tiny_nih2n2():
    """Minimal Ni(H)2L geometry: metal, two hydrides, two N donors."""
    symbols = ["Ni", "H", "H", "N", "N"]
    coords = [
        np.array([0.0, 0.0, 0.0]),
        np.array([0.0, 0.0, 1.45]),
        np.array([0.0, 0.0, -1.45]),
        np.array([1.90, 0.0, 0.0]),
        np.array([-1.90, 0.0, 0.0]),
    ]
    return symbols, coords


BUILDERS = {
    "methanol": methanol,
    "ethanol": ethanol,
    "n_propanol": n_propanol,
    "isopropanol": isopropanol,
    "n_butanol": n_butanol,
    "isobutanol": isobutanol,
    "sec_butanol": sec_butanol,
    "tert_butanol": tert_butanol,
    "tiny_nih2n2": tiny_nih2n2,
}


def write_xyz(path: Path, symbols, coords, comment: str):
    coords = np.asarray(coords, dtype=float)
    if coords.shape[0] != len(symbols):
        raise ValueError(
            f"{comment}: {len(symbols)} symbols vs {coords.shape[0]} coordinates"
        )
    lines = [str(len(symbols)), comment]
    for sym, xyz in zip(symbols, coords):
        lines.append(f"{sym:2s} {xyz[0]:12.8f} {xyz[1]:12.8f} {xyz[2]:12.8f}")
    path.write_text("\n".join(lines) + "\n", encoding="ascii")


def main():
    XYZ_DIR.mkdir(parents=True, exist_ok=True)
    for name, builder in BUILDERS.items():
        symbols, coords = builder()
        write_xyz(XYZ_DIR / f"{name}.xyz", symbols, coords, name)
        print(f"wrote {name}.xyz ({len(symbols)} atoms)")
    targets = HERE / "targets.csv"
    rows = ["name,smiles,bp_c"]
    for name, smiles, bp in ALCOHOLS:
        rows.append(f"{name},{smiles},{bp}")
    targets.write_text("\n".join(rows) + "\n", encoding="ascii")
    print(f"wrote {targets}")


if __name__ == "__main__":
    main()
