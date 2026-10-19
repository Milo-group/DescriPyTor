"""Usage: python make_guide.py  -- writes the "Align & renumber" picture guide (SVG panels + guide.json) next to itself."""
import json
import logging
import os
import sys

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdDepictor, rdDetermineBonds
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Geometry import Point2D

sys.path.insert(0, "C:/Users/edens/Documents/GitHub/DescriPyTor_to_upload")
import conformaty as cf  # noqa: E402
from M2_data_extractor.alignment import _atoms, fit_onto, map_atoms  # noqa: E402
from M2_data_extractor.data_extractor import Molecule  # noqa: E402

logging.disable(logging.WARNING)
OUT = os.path.dirname(os.path.abspath(__file__))
CASES = "C:/Users/edens/Documents/GitHub/Conformaty/test_cases/"
FEATHER = "C:/Users/edens/Documents/GitHub/DescriPyTor_to_upload/descripytor/examples/feather_example/"
GREEN, AMBER, BLUE, GREY = "#2F7A4F", "#F59E0B", "#3B82F6", "#C8CCD2"
W, H = 520, 270


def rgb(hex_, a=1.0):
    return tuple(int(hex_[i:i + 2], 16) / 255 for i in (1, 3, 5)) + (a,)


def mol_from(syms, xyz):
    """RDKit mol from symbols and coordinates; each atom keeps its 0-based input index as 'i'."""
    block = f"{len(syms)}\n\n" + "\n".join(f"{s} {x:.6f} {y:.6f} {z:.6f}" for s, (x, y, z) in zip(syms, xyz))
    m = Chem.MolFromXYZBlock(block)
    try:
        rdDetermineBonds.DetermineBonds(m, charge=0)
    except Exception:  # bond orders can fail (odd valence): connectivity is enough for a picture
        m = Chem.MolFromXYZBlock(block)
        rdDetermineBonds.DetermineConnectivity(m)
    Chem.RemoveStereochemistry(m)          # no wedges in the pictures
    for a in m.GetAtoms():
        a.SetIntProp("i", a.GetIdx())
    return m


def frame(df):
    x = df.iloc[2:]
    return x["Symbol"].astype(str).tolist(), x[["X", "Y", "Z"]].to_numpy(float)


def heavy(m):
    h = Chem.RemoveHs(m, sanitize=False)
    return h, {a.GetIntProp("i"): a.GetIdx() for a in h.GetAtoms()}


def draw_pair(name, ref, qry, ref_notes, qry_notes, pairs, ref_hl, qry_hl, legends, stack=False):
    """Reference and query (heavy atoms) side by side (stack: one above the other); the query laid out like the reference over `pairs` {q: r}
    (0-based full indices). notes: {full index: label}; hl: {full index: colour}."""
    (r, rmap), (q, qmap) = heavy(ref), heavy(qry)
    rdDepictor.Compute2DCoords(r)
    rc = r.GetConformer()
    coord = {qmap[k]: Point2D(rc.GetAtomPosition(rmap[v]).x, rc.GetAtomPosition(rmap[v]).y)
             for k, v in pairs.items() if k in qmap and v in rmap}
    rdDepictor.Compute2DCoords(q, coordMap=coord)
    w, h, top = (W, 140, 22) if stack else (W // 2, H - 22, 22)
    d = rdMolDraw2D.MolDraw2DSVG(W, 2 * (h + top) if stack else h + top, w, h, noFreetype=True)
    o = d.drawOptions()
    o.setBackgroundColour((1, 1, 1, 1))
    o.annotationColour = (0.07, 0.07, 0.07, 1)
    o.annotationFontScale, o.fixedBondLength, o.padding = 0.9, 32, 0.08
    labels = []
    for k, (m, mp, notes, hl, legend) in enumerate(((r, rmap, ref_notes, ref_hl, legends[0]),
                                                    (q, qmap, qry_notes, qry_hl, legends[1]))):
        for i, v in notes.items():
            if i in mp:
                m.GetAtomWithIdx(mp[i]).SetProp("atomNote", str(v))
        c = {mp[i]: rgb(v, 0.7) for i, v in hl.items() if i in mp}
        b = {x.GetIdx(): c[x.GetBeginAtomIdx()] for x in m.GetBonds()
             if c.get(x.GetBeginAtomIdx(), 0) == c.get(x.GetEndAtomIdx(), 1)}
        x0, y0 = (0, k * (h + top)) if stack else (k * w, 0)
        d.SetOffset(x0, y0 + top)          # RDKit's own legends are letter-spaced without FreeType: drawn by hand
        d.DrawMolecule(m, highlightAtoms=list(c), highlightAtomColors=c, highlightBonds=list(b), highlightBondColors=b)
        labels.append(f'<text x="{x0 + 10}" y="{y0 + 17}" font-family="Arial, Helvetica, sans-serif" font-size="14" '
                      f'fill="#222">{legend}</text>')
    d.FinishDrawing()
    write(name, d.GetDrawingText().replace("</svg>", "\n".join(labels) + "\n</svg>"))


def proj_svg(name, panels, basis, centre):
    """Hand-drawn 2D projection of 3D layers. panels: [(caption, [(syms, xyz, bonds, colour, width, notes, fills)])]."""
    pts = [(np.asarray(xyz) - centre) @ basis.T for _c, layers in panels for _s, xyz, *_ in layers]
    allp = np.vstack(pts)
    lo, hi = allp.min(0), allp.max(0)
    pw, pad = W / len(panels), 22
    s = min((pw - 2 * pad) / (hi[0] - lo[0]), (H - 2 * pad - 24) / (hi[1] - lo[1]))
    h = int(s * (hi[1] - lo[1]) + 2 * pad + 30)              # as tall as the drawing needs
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{h}" viewBox="0 0 {W} {h}" '
           'font-family="Arial, Helvetica, sans-serif">', f'<rect width="{W}" height="{h}" fill="#FFFFFF"/>']
    k = 0
    for j, (caption, layers) in enumerate(panels):
        ox = j * pw + (pw - s * (hi[0] - lo[0])) / 2

        def P(p):
            return ox + s * (p[0] - lo[0]), pad + s * (hi[1] - p[1])
        for syms, _xyz, bonds, colour, width, notes, fills in layers:
            p = [P(v) for v in pts[k]]
            k += 1
            for a, b in bonds:
                out.append(f'<line x1="{p[a][0]:.1f}" y1="{p[a][1]:.1f}" x2="{p[b][0]:.1f}" y2="{p[b][1]:.1f}" '
                           f'stroke="{colour}" stroke-width="{width}" stroke-linecap="round"/>')
            for i in sorted(notes, key=lambda i: syms[i] != "H"):  # H first, heavy atoms on top
                x, y = p[i]
                r = 7.5 if syms[i] == "H" else 9.5
                out.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r}" fill="{fills.get(i, "#FFFFFF")}" '
                           f'stroke="{colour}" stroke-width="1"/>')
                out.append(f'<text x="{x:.1f}" y="{y + 3.5:.1f}" font-size="{9 if syms[i] == "H" else 10}" '
                           f'text-anchor="middle" fill="{"#FFF" if i in fills else "#111"}">{notes[i]}</text>')
        out.append(f'<text x="{j * pw + pw / 2:.1f}" y="{h - 10}" font-size="14" text-anchor="middle" '
                   f'fill="#222">{caption}</text>')
    write(name, "\n".join(out + ["</svg>"]))


def write(name, svg):
    with open(os.path.join(OUT, name), "w", encoding="utf-8") as f:
        f.write(svg)


def basis_of(xyz, tilt=0.0):
    """Best-plane axes (PCA) of xyz, the view tilted by `tilt` radians about the long axis."""
    c = xyz.mean(0)
    vt = np.linalg.svd(xyz - c)[2]
    up = np.cos(tilt) * vt[1] + np.sin(tilt) * vt[2]
    return np.vstack([vt[0], up]), c


def bonds_of(m):
    return [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in m.GetBonds()]


guide = []

# 1-6 data. p-F onto the unsubstituted parent, DescriPyTor's own example set.
basic, pf = _atoms(Molecule(FEATHER + "basic.feather")), _atoms(Molecule(FEATHER + "p_F.feather"))
res = map_atoms(basic, pf)
new = {o - 1: n - 1 for o, n in res["old_to_new"].items()}          # query index -> new index (0-based)
seed = [o for o, n in new.items() if res["stage"][n] == "init_mcs"]
rmol, qmol = mol_from(*basic), mol_from(*pf)
f_atom = pf[0].index("F")
draw_pair("01_common_part.svg", rmol, qmol, {i: i + 1 for i in range(len(basic[0]))},
          {o: n + 1 for o, n in new.items()}, {o: new[o] for o in seed},
          {new[o]: GREEN for o in seed}, {o: GREEN for o in seed}, ["reference (parent)", "query (p-F)"], stack=True)
guide.append(("01_common_part.svg", "The common part",
              "Conformaty first finds the largest common substructure, the biggest set of bonded atoms both molecules "
              f"share (green), and gives those atoms the reference's numbers. The F, outside it, takes number "
              f"{new[f_atom] + 1}, the number of the H it replaces."))

# 2-4: one stage each on test case 03 (cyclohexane reference, 2-methyl-tetrahydropyran query).
ref03, q03 = cf.xyz_to_df(CASES + "03_ring_numbering/stage_order_ref.xyz"), cf.xyz_to_df(CASES + "03_ring_numbering/stage_order_query.xyz")
rs, rx = frame(ref03)
qs, qx = frame(q03)
r3, q3 = mol_from(rs, rx), mol_from(qs, qx)
o_at = qs.index("O")
me = next(i for i, s in enumerate(qs) if s == "C" and not q3.GetAtomWithIdx(i).IsInRing())
runs = {}
for st in ("complete_rings", "extend_mcs", "coords_proximity_any"):
    state = cf.map_frames(ref03.copy(), q03.copy(), stage_order=[st])
    runs[st] = (state.mapping, {q: k for k, q in enumerate(cf.query_order(state))})
ring_pairs = {q: r for q, r in runs["complete_rings"][0].items() if qs[q] != "H"}   # one layout for all three
six = {st: next(q for q, k in num.items() if k == 5) for st, (_m, num) in runs.items()}  # who takes reference atom 6
assert six == {"complete_rings": o_at, "extend_mcs": me, "coords_proximity_any": o_at}, six
text = {
    "complete_rings": ("complete rings", "Once part of a ring matches, the rest of the ring is matched too, whatever "
                       "the element: the ring O (amber) takes reference atom 6, so the two rings correspond atom for atom."),
    "extend_mcs": ("extend the common part", "The common part grows outward one bonded atom at a time, element for "
                   "element: the methyl carbon, bonded to atom 1 as reference carbon 6 is, takes number 6 (amber), and "
                   f"the ring O, left over, goes to the end ({runs['extend_mcs'][1][o_at] + 1})."),
    "coords_proximity_any": ("any element by position", "After the fit, an atom lying within 0.2 Å of a reference "
                             "atom takes its number, whatever its element: the ring O sits on reference carbon 6, so it "
                             "takes 6 (amber), the same answer as complete rings for a different reason."),
}
for n, (st, (_m, num)) in enumerate(runs.items(), start=2):
    title, cap = text[st]
    fname = f"0{n}_{st}.svg"
    draw_pair(fname, r3, q3, {i: i + 1 for i in range(len(rs))}, {q: k + 1 for q, k in num.items()}, ring_pairs,
              {5: AMBER}, {six[st]: AMBER}, ["reference (cyclohexane)", f"query, {title}"])
    guide.append((fname, title.capitalize(), cap))

# 5: equivalent H atoms by position, test case 02 (the smaller molecule is the reference, as Conformaty's pair does).
ref02, q02 = cf.xyz_to_df(CASES + "02_prochirality/cyclohexane_boat.xyz"), cf.xyz_to_df(CASES + "02_prochirality/cyclohexanol.xyz")
state = cf.map_frames(ref02, q02)
rs2, rx2 = frame(state.df_ref)            # map_frames centres the reference with the query
qs2, qx2 = frame(state.df_query)          # fitted onto the reference and already in the new order
ring = rx2[[i for i, s in enumerate(rs2) if s != "H"]]


def crowding(t):                          # closest pair of projected reference atoms: larger is more readable
    b, c0 = basis_of(ring, t)
    p = (rx2 - c0) @ b.T
    return min(np.linalg.norm(p[i] - p[j]) for i in range(len(p)) for j in range(i))


basis, c = basis_of(ring, max(np.linspace(0, np.pi, 37), key=crowding))
byp = {i: BLUE for i, s in enumerate(qs2) if s == "H" and i < len(rs2) and rs2[i] == "H"}
o2 = qs2.index("O")
proj_svg("05_equivalent_atoms.svg", [
    ("reference (cyclohexane)", [(rs2, rx2, bonds_of(mol_from(rs2, rx2)), "#555", 2,
                                  {i: i + 1 for i in range(len(rs2))}, {})]),
    ("query (cyclohexanol), after the fit", [(qs2, qx2, bonds_of(mol_from(qs2, qx2)), "#555", 2,
                                              {i: (f"O{i + 1}" if i == o2 else i + 1) for i in range(len(qs2))}, byp)]),
], basis, c)
guide.append(("05_equivalent_atoms.svg", "Equivalent atoms by position",
              "The two H on a CH2 are the same by bonding, so after the fit each one takes the number of the "
              "reference H nearest to it (blue); swapping them would look the same on paper and be wrong in 3D. The "
              f"OH oxygen takes the number of the H it replaces ({o2 + 1}); the reference is a boat and the query a "
              "chair, so only carbon 4 and its H sit apart."))

# 6: the fit, rmsd loss (the default in DescriPyTor) against Conformaty's mixed loss, p-F onto the parent.
order = sorted(new, key=new.get)
syms = [pf[0][o] for o in order]
moving = np.asarray(pf[1], float)[order]
bx = np.asarray(basic[1], float)
same = [i for i in range(len(basic[0])) if syms[i] == basic[0][i]]
fits = {f: fit_onto(basic, syms, moving, f) for f in ("rmsd", "mixed")}
dev = {f: float(np.sqrt(((x[same] - bx[same]) ** 2).sum(1).mean())) for f, x in fits.items()}
basis, c = basis_of(bx)
rb, qb = bonds_of(mol_from(*basic)), bonds_of(mol_from(syms, moving))
proj_svg("06_the_fit.svg", [
    (f"rmsd loss: {dev[f]:.2f} Å" if f == "rmsd" else f"mixed loss: {dev[f]:.1f} Å",
     [(basic[0], bx, rb, GREY, 7, {}, {}), (syms, fits[f], qb, GREEN, 2, {}, {})]) for f in ("rmsd", "mixed")
], basis, c)
guide.append(("06_the_fit.svg", "The fit",
              "To overlay the query on the reference (light grey), Conformaty moves it as a rigid body. With the rmsd "
              f"loss, the default here, p-F lands on the parent at {dev['rmsd']:.2f} Å; with Conformaty's mixed loss "
              f"the same fit stops {dev['mixed']:.1f} Å away and never finds the overlay."))
assert dev["rmsd"] < 1 < dev["mixed"], dev

with open(os.path.join(OUT, "guide.json"), "w", encoding="utf-8") as f:
    json.dump([{"file": a, "title": b, "caption": t} for a, b, t in guide], f, indent=2)
for a, *_ in guide:
    print(a, os.path.getsize(os.path.join(OUT, a)), "bytes")
