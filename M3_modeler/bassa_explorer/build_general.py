"""Build M(H)2L complexes from metal-free bidentate ligand geometries.

The existing build_complexes.py superposes a same-scaffold template, which is
exact but only works when a template of that scaffold exists. The Corminboeuf
deposit is not one scaffold: 75 of its 100 ligands are oxazoline-type, but 11
are N-only (bipyridine / imidazoline / diamine) and 6 contain phosphorus
(phosphinooxazoline). So the metal is placed geometrically instead:

    1. find the bidentate donor pair from the bond graph alone -- two lone-pair
       heteroatoms separated by 3 or 4 bonds, i.e. a 5- or 6-membered chelate
       once the metal closes the ring
    2. put M where both donor lone pairs point, at the right M-D distances
    3. add two hydrides trans to the donors, completing square-planar d8

Validated against our DFT-optimised Doyle complexes before use: strip Ni and
both H from a known-good structure, rebuild, and compare.
"""
import argparse
import glob
import os
from pathlib import Path
from collections import deque

import numpy as np
from scipy.optimize import minimize

COV = {"H": .31, "C": .76, "N": .71, "O": .66, "F": .57, "P": 1.07,
       "S": 1.05, "Cl": 1.02, "Br": 1.20, "Ni": 1.24, "Cu": 1.32}
# M-donor distances. Ni-N is read off our own optimised set; the Cu values are
# standard crystallographic ranges (Cu(I) is longer and softer than Cu(II)).
MD_BY_METAL = {
    "Ni": {"N": 2.09, "P": 2.20, "S": 2.25, "O": 2.05},
    "Cu": {"N": 2.00, "P": 2.20, "S": 2.30, "O": 2.00},
}
MD = MD_BY_METAL["Ni"]                       # default, overridden per build
# M-ancillary distances
M_ANC = {("Ni", "H"): 1.48, ("Cu", "Cl"): 2.20, ("Ni", "Cl"): 2.24}
MH = 1.48
DONORS = ("N", "P", "S")


def read_xyz(path):
    L = open(path).read().splitlines()
    n = int(L[0].split()[0])
    sym = [x.split()[0] for x in L[2:2 + n]]
    xyz = np.array([[float(v) for v in x.split()[1:4]] for x in L[2:2 + n]])
    return sym, xyz


def write_xyz(path, sym, xyz, note=""):
    with open(path, "w", newline="\n") as f:
        f.write(f"{len(sym)}\n{note}\n")
        for a, (x, y, z) in zip(sym, xyz):
            f.write(f"{a:<3s} {x:18.10f} {y:18.10f} {z:18.10f}\n")


def adjacency(sym, xyz, scale=1.3):
    n = len(sym)
    adj = [[] for _ in range(n)]
    d = np.linalg.norm(xyz[:, None] - xyz[None], axis=-1)
    for i in range(n):
        for j in range(i + 1, n):
            if d[i, j] < scale * (COV.get(sym[i], .77) + COV.get(sym[j], .77)):
                adj[i].append(j)
                adj[j].append(i)
    return adj


def path_len(adj, a, b, cap=6):
    """Shortest bond-path length between two atoms, or None beyond cap."""
    q, seen = deque([(a, 0)]), {a}
    while q:
        u, d = q.popleft()
        if u == b:
            return d
        if d >= cap:
            continue
        for v in adj[u]:
            if v not in seen:
                seen.add(v)
                q.append((v, d + 1))
    return None


def lone_pair(i, sym, xyz, adj):
    """Direction the lone pair points: away from the mean of the neighbours."""
    nb = adj[i]
    if not nb:
        return None
    v = -sum((xyz[j] - xyz[i]) / np.linalg.norm(xyz[j] - xyz[i]) for j in nb)
    nv = np.linalg.norm(v)
    return v / nv if nv > 1e-6 else None


def donor_candidates(sym, xyz, adj):
    """Rank every plausible bidentate pair; the caller rejects on geometry.

    3 bonds apart -> 5-membered chelate (bi-oxazoline, bipyridine);
    4 -> 6-membered (bis-oxazoline with a bridge).

    Ranking prefers (a) fewer heavy neighbours, because a two-coordinate imine
    N is a donor and a three-coordinate N-aryl N is not -- this is what
    separates the real donors in bis(imidazoline) from the N-aryl nitrogens --
    and (b) lone pairs that converge on a common point rather than diverging.
    """
    cand = [i for i, s in enumerate(sym)
            if s in DONORS and len([j for j in adj[i] if sym[j] != "H"]) <= 3]
    out = []
    for a in range(len(cand)):
        for b in range(a + 1, len(cand)):
            i, j = cand[a], cand[b]
            L = path_len(adj, i, j)
            if L not in (3, 4):
                continue
            li, lj = lone_pair(i, sym, xyz, adj), lone_pair(j, sym, xyz, adj)
            if li is None or lj is None:
                continue
            axis = xyz[j] - xyz[i]
            axis /= np.linalg.norm(axis)
            conv = np.dot(li, axis) - np.dot(lj, axis)
            heavy = (len([k for k in adj[i] if sym[k] != "H"])
                     + len([k for k in adj[j] if sym[k] != "H"]))
            out.append((-heavy, conv, i, j, L))
    out.sort(reverse=True)
    return [(i, j, L) for _, _, i, j, L in out]


def place_metal(sym, xyz, adj, i, j, metal="Ni"):
    md = MD_BY_METAL.get(metal, MD)
    li, lj = lone_pair(i, sym, xyz, adj), lone_pair(j, sym, xyz, adj)
    ri, rj = md.get(sym[i], 2.1), md.get(sym[j], 2.1)

    def cost(m):
        di, dj = m - xyz[i], m - xyz[j]
        ni, nj = np.linalg.norm(di), np.linalg.norm(dj)
        # Right distance first, lone-pair alignment second. At weight 2.0 the
        # angular term dominated and pushed M outward along the bisector,
        # giving M-D 2.16 against a DFT 2.07 and a bite 4 deg too narrow.
        return (3.0 * ((ni - ri) ** 2 + (nj - rj) ** 2)
                + 0.5 * ((1 - np.dot(di / ni, li)) + (1 - np.dot(dj / nj, lj))))

    start = (xyz[i] + xyz[j]) / 2 + (li + lj) / np.linalg.norm(li + lj) * 1.6
    res = minimize(cost, start, method="Nelder-Mead",
                   options=dict(xatol=1e-4, fatol=1e-8, maxiter=4000))
    return res.x


def ancillary_one(m, d1, d2, ligand, dist):
    """A single ancillary for a 3-coordinate metal (Cu(I) with a bidentate L).

    Cu(I) is d10 and prefers trigonal, not square planar, so one ligand sits on
    the external bisector of the chelate rather than trans to a donor.
    """
    u1 = (d1 - m) / np.linalg.norm(d1 - m)
    u2 = (d2 - m) / np.linalg.norm(d2 - m)
    v = -(u1 + u2)
    n = np.linalg.norm(v)
    v = v / n if n > 1e-6 else np.cross(u1, [0, 0, 1.0])
    if ligand is None:
        return [m + v * dist]
    best, bestd = v, -1.0
    perp = np.cross(np.cross(u1, u2), v)
    pn = np.linalg.norm(perp)
    perp = perp / pn if pn > 1e-6 else np.array([0., 0., 1.])
    for tilt in np.radians([0, 10, 20, 30]):
        for phi in np.linspace(0, 2 * np.pi, 12, endpoint=False):
            axis = np.cos(phi) * perp + np.sin(phi) * np.cross(v, perp)
            w = v * np.cos(tilt) + np.cross(axis, v) * np.sin(tilt)
            w /= np.linalg.norm(w)
            clear = np.linalg.norm(ligand - (m + w * dist), axis=1).min()
            score = min(clear, 3.0) - 0.15 * tilt
            if score > bestd:
                best, bestd = w, score
    return [m + best * dist]


def hydrides(m, d1, d2, ligand=None, dist=None, target=2.0):
    """Each hydride trans to a donor. Ideally in the chelate plane.

    On 14 of the 100 Corminboeuf ligands a substituent already occupies the
    in-plane trans position, putting the hydride 1.1-1.4 A from a ligand atom.
    The metal placement is right in every one of those cases -- it is the
    hydride that is in the way. A real hydride bends out of the way, so when
    the ideal position clashes we tilt it off-plane, up to 30 degrees, taking
    the direction with the most clearance. GOAT relaxes the rest.
    """
    D = MH if dist is None else dist
    u1 = (d1 - m) / np.linalg.norm(d1 - m)
    u2 = (d2 - m) / np.linalg.norm(d2 - m)
    ideal = [-u1, -u2]
    if ligand is None:
        return m + ideal[0] * D, m + ideal[1] * D
    n = np.cross(u1, u2)
    n = n / np.linalg.norm(n) if np.linalg.norm(n) > 1e-6 else np.array([0., 0., 1.])
    out = []
    for u in ideal:
        best, bestd = u, -1.0
        for tilt in np.radians([0, 10, 20, 30]):
            perp = np.cross(n, u)
            perp /= np.linalg.norm(perp)
            for phi in np.linspace(0, 2 * np.pi, 12, endpoint=False):
                axis = np.cos(phi) * n + np.sin(phi) * perp
                v = u * np.cos(tilt) + np.cross(axis, u) * np.sin(tilt)
                v /= np.linalg.norm(v)
                clear = np.linalg.norm(ligand - (m + v * D), axis=1).min()
                # prefer clearance, but do not wander from trans for nothing
                # target scales with the ancillary: 2.0 A is fine for a
                # hydride but far too low for chloride, which left a
                # 1.70 A Cl...ligand contact on the Cu(II) builds.
                score = min(clear, target) - 0.15 * tilt
                if score > bestd:
                    best, bestd = v, score
            # Scale the early exit with the target, or raising `target` does
            # nothing: at 1.8 the loop stopped at tilt=0 and never tried the
            # larger tilts a chloride needs.
            if bestd > 0.9 * target:
                break
        out.append(m + best * D)
    return out[0], out[1]


def build(path, metal="Ni", anc="H", n_anc=2):
    """Try candidate donor pairs in rank order; accept the first sane geometry.

    A pair is rejected if the metal or a hydride lands on top of the ligand, or
    if the bite angle is nowhere near a real chelate. Ranking alone put the
    metal 0.67 A from a ligand atom on the bis(imidazoline) ligands, so the
    geometric check -- not the heuristic -- is what makes this reliable.
    """
    sym, xyz = read_xyz(path)
    adj = adjacency(sym, xyz)
    cands = donor_candidates(sym, xyz, adj)
    if not cands:
        return None, "no bidentate donor pair found"
    rejected = []
    for i, j, L in cands:
        m = place_metal(sym, xyz, adj, i, j, metal)
        dist = M_ANC.get((metal, anc), 2.2)
        if n_anc == 2:
            tgt = 2.0 if anc == "H" else 3.0
            extra = list(hydrides(m, xyz[i], xyz[j], ligand=xyz,
                                  dist=dist, target=tgt))
        elif n_anc == 1:
            extra = ancillary_one(m, xyz[i], xyz[j], xyz, dist)
        else:
            extra = []
        new = np.array([m] + extra)
        dmin = np.linalg.norm(xyz[:, None] - new[None], axis=-1).min()
        bite = np.degrees(np.arccos(np.clip(np.dot(
            (xyz[i] - m) / np.linalg.norm(xyz[i] - m),
            (xyz[j] - m) / np.linalg.norm(xyz[j] - m)), -1, 1)))
        d1, d2 = np.linalg.norm(m - xyz[i]), np.linalg.norm(m - xyz[j])
        if dmin < 1.45 or not (65 < bite < 105):
            rejected.append((sym[i], sym[j], round(bite, 1), round(dmin, 2)))
            continue
        info = dict(donors=(sym[i], sym[j]), ring=L + 2, d1=d1, d2=d2,
                    bite=bite, closest=dmin, tried=len(rejected) + 1)
        order = [i, j] + [k for k in range(len(sym)) if k not in (i, j)]
        return ([metal] + [anc] * len(extra) + [sym[k] for k in order],
                np.vstack([new, xyz[order]])), info
    return None, f"all {len(cands)} candidate pairs rejected: {rejected[:3]}"


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("indir")
    ap.add_argument("outdir")
    ap.add_argument("--metal", default="Ni")
    ap.add_argument("--anc", default="H", help="ancillary element")
    ap.add_argument("--n-anc", type=int, default=2)
    ap.add_argument("--only", type=Path, help="file listing ligand stems to build")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    ok = bad = 0
    keep = None
    if a.only:
        keep = {l.strip() for l in a.only.read_text().split() if l.strip()}
    for p in sorted(glob.glob(os.path.join(a.indir, "*.xyz"))):
        if keep is not None and os.path.basename(p)[:-4] not in keep:
            continue
        res, info = build(p, a.metal, a.anc, a.n_anc)
        name = os.path.basename(p)
        if res is None:
            print(f"  FAIL {name}: {info}")
            bad += 1
            continue
        sym, xyz = res
        write_xyz(os.path.join(a.outdir, name), sym, xyz,
                  f"{name[:-4]} {a.metal}({a.anc}){a.n_anc}L  donors {info['donors']} "
                  f"{info['ring']}-ring bite {info['bite']:.1f}")
        flag = "" if (1.9 < info["d1"] < 2.4 and 1.9 < info["d2"] < 2.4
                      and info["closest"] > 1.3) else "  <-- CHECK"
        print(f"  {name[:-4]:<14s} {info['donors'][0]}/{info['donors'][1]} "
              f"{info['ring']}-ring  M-D {info['d1']:.2f}/{info['d2']:.2f}  "
              f"bite {info['bite']:5.1f}  closest {info['closest']:.2f}{flag}")
        ok += 1
    print(f"\nbuilt {ok}, failed {bad}")
