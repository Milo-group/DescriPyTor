"""Metal-complex geometric and electronic descriptors from XYZ + xTB.

This is the Case Study 3 extractor, turned into classes. Canonical atom order
(the layout ``build_general.py`` writes, and that GOAT/xTB then preserve)::

    0 = metal
    1 .. n_ancillary-1 = H / Cl / F / Br
    then the two chelate donors

No Gaussian feathers are required. Sterimol is the Verloop CPK implementation
used for the CS3 tables (not morfeus). Buried volume is a Fibonacci-grid
Cavallo %Vbur with the CS3 Bondi table (Ni 1.63 Å, Cu 1.40 Å).

Typical use::

    from M2_data_extractor import (
        MetalComplex, MetalComplexEnsemble, XtbSinglePoint, parse_props_file,
    )

    ens = MetalComplexEnsemble.from_xyz("081_lig.finalensemble.xyz")
    geom = ens.geometric_features()          # Boltzmann-averaged at 298.15 K

    sp = XtbSinglePoint.from_xtb_dir("chg_ni", "081_lig", props=props["081_lig"])
    mc = MetalComplex.from_xyz("081_lig.xyz")
    elec = mc.electronic_features(sp)
"""

from __future__ import annotations

import glob
import os
import warnings
from collections import deque

import numpy as np
import pandas as pd

from .xyz_io import XYZEnsemble, XYZFrame, boltzmann_average, parse_xyz_ensemble
from .xtb_singlepoint import XtbSinglePoint

# Connectivity radii (build_general.py). Used for donor/stereocentre walks.
CONNECT_COV = {
    "H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "P": 1.07,
    "S": 1.05, "Cl": 1.02, "Br": 1.20, "Ni": 1.24, "Cu": 1.32,
}

# Sterimol fragment adjacency (regen_b1.py). Intentionally has Ni but not Cu —
# that is how the CS3 tables were built.
STERIMOL_COV = {
    "H": 0.31, "B": 0.84, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57,
    "Si": 1.11, "P": 1.07, "S": 1.05, "Cl": 1.02, "Br": 1.20, "I": 1.39,
    "Ni": 1.24,
}

# Verloop CPK types (approach_b1.py / sterimol_standalone.CPK_RADII).
PKG_CPK = {
    "C": 1.50, "C3": 1.60, "C6/N6": 1.70, "H": 1.00, "N": 1.50,
    "N4": 1.45, "O": 1.35, "O2": 1.35, "P": 1.40, "S": 1.70, "S1": 1.00,
    "F": 1.35, "Cl": 1.80, "S4": 1.40, "Br": 1.95, "I": 2.15, "X": 1.92,
    "Ni": 1.97, "Pd": 2.10, "Cu": 1.96, "Zn": 2.01, "Fe": 2.04, "Co": 2.00,
}

# Cavallo %Vbur Bondi radii as used in extract_cb.py (not Alvarez).
VDW_BURIAL = {
    "H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47, "P": 1.80,
    "S": 1.80, "Cl": 1.75, "Br": 1.85, "I": 1.98, "Ni": 1.63, "Cu": 1.40,
}

ANCILLARY_ELEMENTS = ("H", "Cl", "F", "Br")
# Per-arm Sterimol values written as <frame>_<key>_sym/_asym. "theta" is the angle
# between the vector to the B5 atom and the B1 plane; it is defined for every arm
# (0 for a lone H on the axis). "angle" is phi, the in-plane azimuth, which has no
# defined value when the B5 atom sits on the axis; add it to rebuild the deposited
# CS3 tables. Both follow STERIMOL_THETA_RULE.
STERIMOL_KEYS = ("B1", "B5", "L", "theta")

# Coarse step of the B1 rotation scan, in degrees. The best coarse angle is
# refined at 1° within one step either side, so 1 tries every degree. The
# deposited CS3 tables were built with 18, which misses the narrowest plane on
# some arms; set it to 18 to rebuild them.
STERIMOL_SCAN_STEP = 1

# Where the rotation scan starts. "fragment" (the default) takes the fragment atom
# furthest off the axis, so the result does not depend on how the file is oriented.
# "lab" is the frame every CS3 table up to tag paper-v3 was built in (cross product
# of the axis with the lab z-axis); set it to rebuild those tables.
STERIMOL_FRAME = "fragment"

# How theta is read when several B1 directions nearly tie (a tBu group has three within
# 0.01 A, giving theta 4.7, 38.8 or 45.5 deg depending on which the scan meets first).
# "soft" (the default) averages the tilt over every direction, weighted by
# exp(-(width - B1) / THETA_SOFT_TAU), so equal directions count equally and the value moves
# continuously with the geometry. "scan" is the direction the rotation scan picks, which is
# how every table up to tag paper-v3 was built.
STERIMOL_THETA_RULE = "soft"
THETA_SOFT_TAU = 0.02          # A
THETA_SOFT_STEP = 0.05         # deg, direction grid for the soft rule

# Where the substituent fragment (C* -> R) stops. "donor_ring" (the default) blocks every atom
# of the donor's own ring, so a substituent fused to that ring cannot walk into the chelate
# backbone and the other arm (7 of 60 arms in the CS3 cyclopropanation set did). "donor" blocks
# only the donor and the metal, as tag paper-v3 did.
SUB_FRAGMENT_BOUND = "donor_ring"

# Which (C*, R) pair stands for an arm's substituent when several qualify. "alpha" (the default):
# the largest substituent on the donor's ring; among equals the one on the carbon bonded to the
# donor, then the one nearest the metal, so the choice does not depend on the atom order of the
# file. A donor outside any ring whose only candidates are bonds inside a ring (ipso -> ortho of a
# P-aryl group) has no substituent axis and the arm is skipped. "first" keeps the first of the
# largest in atom order, as every release up to 0.2.1 did.
STEREOCENTRE_RULE = "alpha"

# mu_outofplane and mu_desym are components along axes whose sign follows which donor is listed
# first in the file (the plane normal is D1 x D2). True reports their magnitudes, the only
# meaningful part; False reports the signed values of tag paper-v3.
DIPOLE_ABS = True

# build_general.py rejects a metal placement whose bite falls outside this window,
# but only before relaxation; a donor can still come off during it (009_lig in the
# CS3 oxy-alkynylation set: 87 deg placed, 58 deg relaxed). geometric_features
# re-checks the structure it is given.
BITE_WINDOW = (65.0, 105.0)
ELECTRONIC_COLUMNS = (
    "mu_bisector", "mu_outofplane", "mu_desym", "homo", "lumo", "gap",
    "q_metal", "q_donor_sym", "q_donor_asym", "q_anc_sum", "q_absmax",
    "q_spread", "wbo_MD_sym", "wbo_MD_asym", "wbo_metal_tot",
)


def _sphere_grid(n=4000):
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    theta = np.pi * (1 + 5 ** 0.5) * i
    return np.c_[np.cos(theta) * np.sin(phi),
                 np.sin(theta) * np.sin(phi),
                 np.cos(phi)]


_SHELLS = np.linspace(0.15, 1.0, 12)
_DIRS = _sphere_grid()
# Unit-sphere sample points and shell weights, reused for every %Vbur call.
_UNIT_PTS = (_SHELLS[:, None, None] * _DIRS[None, :, :]).reshape(-1, 3)
_VBUR_WEIGHTS = np.repeat(_SHELLS ** 2, len(_DIRS))


def adjacency(symbols, coords, scale=1.3, cov=None):
    """Undirected bond list from covalent radii × ``scale``."""
    cov = cov or CONNECT_COV
    n = len(symbols)
    adj = [[] for _ in range(n)]
    dist = np.linalg.norm(coords[:, None] - coords[None], axis=-1)
    for i in range(n):
        for j in range(i + 1, n):
            cutoff = scale * (cov.get(symbols[i], 0.77) + cov.get(symbols[j], 0.77))
            if dist[i, j] < cutoff:
                adj[i].append(j)
                adj[j].append(i)
    return adj


def nob_types(symbols, adj):
    """Verloop CPK atom types from coordination number (``nob_atype``)."""
    out = []
    for i, symbol in enumerate(symbols):
        nob = len(adj[i])
        if symbol in ("H", "F", "P", "Cl", "Br", "I"):
            atom_type = symbol
        elif symbol == "O":
            atom_type = "O2" if nob < 1.5 else "O"
        elif symbol == "S":
            atom_type = "S" if nob < 2.5 else ("S4" if nob < 5.5 else "S1")
        elif symbol == "N":
            atom_type = "C6/N6" if nob < 2.5 else "N"
        elif symbol == "C":
            atom_type = "C3" if nob < 2.5 else ("C6/N6" if nob < 3.5 else "C")
        elif symbol in PKG_CPK:
            atom_type = symbol
        else:
            atom_type = "X"
        out.append(atom_type)
    return out


def fragment(symbols, coords, a, b, block=()):
    """Atoms on b's side of the a–b bond, never traversing ``block``.

    ``a`` / ``b`` / ``block`` are 1-based, matching the original Sterimol call.
    Bond finding matches ``regen_b1.py`` pairwise (not a vectorized distance
    matrix) so near-cutoff contacts stay bit-identical.
    """
    n = len(symbols)
    adj = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i + 1, n):
            cutoff = (STERIMOL_COV.get(symbols[i], 0.8) + STERIMOL_COV.get(symbols[j], 0.8)) * 1.3
            if np.linalg.norm(coords[i] - coords[j]) < cutoff:
                adj[i].add(j)
                adj[j].add(i)
    seen = {a - 1, b - 1} | {k - 1 for k in block}
    out = {b - 1}
    queue = deque([b - 1])
    while queue:
        x = queue.popleft()
        for y in adj[x]:
            if y not in seen:
                seen.add(y)
                out.add(y)
                queue.append(y)
    return sorted(out)


def _theta_soft(coords, idx, origin, axis, radii):
    """(theta, phi) under the soft tie rule: the tilt of the B5 vector out of the plane normal to
    each in-plane direction, and the in-plane angle between that direction and the B5 atom, both
    averaged over directions with weights exp(-(width - B1) / tau).

    ``radii`` is aligned with ``idx``. The B5 atom is the one furthest off the axis (centre
    distance plus radius). A one-atom fragment, or one lying on the axis, has theta = phi = 0.
    """
    V = coords[idx] - origin
    perp = V - np.outer(V @ axis, axis)
    if len(idx) == 1 or np.linalg.norm(perp, axis=1).max() < 1e-6:
        return 0.0, 0.0
    e1 = perp[int(np.argmax((perp ** 2).sum(1)))]
    e1 = e1 / np.linalg.norm(e1)
    e2 = np.cross(axis, e1)
    P = np.c_[V @ e1, V @ e2]
    r = np.asarray(radii, float)
    phis = np.radians(np.arange(0, 360, THETA_SOFT_STEP))
    dirs = np.c_[np.cos(phis), np.sin(phis)]
    width = (P @ dirs.T + r[:, None]).max(0)
    j5 = int(np.argmax(np.linalg.norm(P, axis=1) + r))
    v5 = V[j5]
    tilt = np.degrees(np.arcsin(np.minimum(1.0, np.abs((dirs @ np.c_[e1, e2].T) @ v5) / np.linalg.norm(v5))))
    p5 = np.linalg.norm(P[j5])
    phi = np.degrees(np.arccos(np.clip(dirs @ P[j5] / p5, -1, 1))) if p5 > 1e-6 else np.zeros(len(dirs))
    w = np.exp(-(width - width.min()) / THETA_SOFT_TAU)
    return float((w * tilt).sum() / w.sum()), float((w * phi).sum() / w.sum())


def donor_ring(adj, donor, metal=0):
    """Atoms of the smallest ring through ``donor`` that does not pass through the metal
    (0-based), or ``{donor}`` when the donor is in no such ring."""
    best = None
    nb = [x for x in adj[donor] if x != metal]
    for i, s in enumerate(nb):
        for t in nb[i + 1:]:
            prev = {s: None}
            queue = deque([s])
            while queue:
                x = queue.popleft()
                for y in adj[x]:
                    if y in (donor, metal) or y in prev:
                        continue
                    prev[y] = x
                    queue.append(y)
            if t in prev:
                path = [t]
                while prev[path[-1]] is not None:
                    path.append(prev[path[-1]])
                ring = set(path) | {donor}
                if best is None or len(ring) < len(best):
                    best = ring
    return best or {donor}


def sterimol(symbols, coords, a, b, radii, block=()):
    """Verloop B1/B5/L and the two B1–B5 angles for the fragment on the a→b axis.

    ``angle`` is phi, measured in the plane perpendicular to the axis; ``theta``
    is the angle between the full vector to the B5 atom and the B1 plane, so
    sin(theta) = |cos(phi)| * rho with rho the in-plane share of that vector.

    Besides the numbers, the result carries the frame they were measured in --
    ``origin``, ``axis``, ``normal`` (the unit normal of the B1 plane, pointing at
    the supporting plane), ``b5_atom`` and ``atoms`` -- so a caller can draw the
    construction without running the rotation scan again.
    """
    idx = fragment(symbols, coords, a, b, block)
    if not idx:
        return None
    origin = coords[a - 1]
    axis = coords[b - 1] - origin
    axis = axis / np.linalg.norm(axis)
    # The rotation scan starts at e1, so e1 has to turn with the molecule. Taking it from the
    # lab frame made every scan-derived angle depend on how the input file happened to be
    # oriented (theta moved up to 25 deg and the azimuth up to 83 deg under a pure rotation);
    # the platform's own path builds its frame from three atoms and is invariant. Use the
    # fragment atom furthest off the axis: deterministic, and never near-degenerate.
    perp = (coords[idx] - origin) - np.outer((coords[idx] - origin) @ axis, axis)
    e1 = perp[int(np.argmax((perp ** 2).sum(1)))] if STERIMOL_FRAME == "fragment" else np.zeros(3)
    if np.linalg.norm(e1) < 1e-8:                     # "lab", or every fragment atom on the axis
        e1 = np.cross(axis, [0, 0, 1.0])
        if np.linalg.norm(e1) < 1e-6:
            e1 = np.cross(axis, [0, 1.0, 0])
    e1 = e1 / np.linalg.norm(e1)
    e2 = np.cross(axis, e1)
    P, R, along = [], [], []
    for i in idx:
        vec = coords[i] - origin
        proj = vec - np.dot(vec, axis) * axis
        P.append([np.dot(proj, e1), np.dot(proj, e2)])
        R.append(radii[i])
        along.append(np.dot(vec, axis))
    P, R, along = np.array(P), np.array(R), np.array(along)
    theta = np.linspace(0, 2 * np.pi, 100)
    # Per-atom vstack matches regen_b1.py. A raveled broadcast can flip B1
    # direction ties and jump the B1–B5 angle by ~90° with B1/B5 unchanged.
    cloud = np.vstack([
        np.column_stack((x + r * np.cos(theta), y + r * np.sin(theta)))
        for (x, y), r in zip(P, R)
    ])

    def at(deg):
        ang = np.radians(deg)
        c, s = np.cos(ang), np.sin(ang)
        t = cloud @ np.array([[c, -s], [s, c]]).T
        ev = [t[:, 0].max(), t[:, 0].min(), t[:, 1].max(), t[:, 1].min()]
        k = int(np.argmin(np.abs(ev)))
        b1c = [(ev[0], 0.0), (ev[1], 0.0), (0.0, ev[2]), (0.0, ev[3])][k]
        j = int(np.argmax((t ** 2).sum(1)))
        a1 = np.arctan2(b1c[1], b1c[0]) % (2 * np.pi)
        a5 = np.arctan2(t[j, 1], t[j, 0]) % (2 * np.pi)
        delta = abs(a5 - a1)
        if delta > np.pi:
            delta = 2 * np.pi - delta
        return abs(ev[k]), float(np.degrees(delta)), float(np.linalg.norm(t[j])), (c, s, k, j)

    step = STERIMOL_SCAN_STEP
    coarse = [(at(x)[0], x) for x in range(18, 108, step)]
    best_deg = min(coarse)[1]
    B1, ang, B5, (c, s, k, j) = min((at(x) for x in range(best_deg - step, best_deg + step + 1)),
                                    key=lambda q: q[0])
    # The scan rotates the cloud, t = R p, so the B1 axis maps back with R^T:
    # rotated x is [c, -s] and rotated y is [s, c] in the (e1, e2) frame.
    n1, n2 = (c, -s) if k in (0, 1) else (s, c)
    b5 = idx[j // len(theta)]                             # cloud rows are per atom, idx order
    v5 = coords[b5] - origin
    tilt = float(np.degrees(np.arcsin(min(1.0, abs(v5 @ (n1 * e1 + n2 * e2)) / np.linalg.norm(v5)))))
    if STERIMOL_THETA_RULE == "soft":
        tilt, ang = _theta_soft(coords, idx, origin, axis, R)
    # Orient the normal towards the supporting plane, so origin + B1 * normal lies on it.
    d = np.array([n1, n2])
    if abs((P @ d + R).max() - B1) > abs((P @ -d + R).max() - B1):
        d = -d
    return dict(B1=B1, B5=B5, angle=ang, theta=tilt, L=float(np.max(along + R)), n=len(idx),
                origin=origin, axis=axis, normal=d[0] * e1 + d[1] * e2, b5_atom=b5, atoms=list(idx))


def centre_index_guard(centre, coords):
    dist = np.linalg.norm(coords - centre, axis=1)
    j = int(np.argmin(dist))
    return j if dist[j] < 1e-6 else -1


def _occupy_numpy(pts, atoms, r2):
    buried = np.zeros(len(pts), dtype=bool)
    px, py, pz = pts[:, 0], pts[:, 1], pts[:, 2]
    for pos, rr in zip(atoms, r2):
        dx = px - pos[0]
        dy = py - pos[1]
        dz = pz - pos[2]
        buried |= (dx * dx + dy * dy + dz * dz) < rr
    return buried


try:
    from numba import njit

    @njit(cache=True)
    def _occupy_numba(pts, atoms, r2):
        n = pts.shape[0]
        buried = np.zeros(n, dtype=np.bool_)
        for i in range(n):
            px = pts[i, 0]
            py = pts[i, 1]
            pz = pts[i, 2]
            hit = False
            for j in range(atoms.shape[0]):
                dx = px - atoms[j, 0]
                dy = py - atoms[j, 1]
                dz = pz - atoms[j, 2]
                if dx * dx + dy * dy + dz * dz < r2[j]:
                    hit = True
                    break
            buried[i] = hit
        return buried

    _occupy = _occupy_numba
except Exception:  # pragma: no cover - optional accelerator
    _occupy = _occupy_numpy


def buried_volume(symbols, coords, centre, radius=3.5, scale=1.17, include_h=True):
    """Percent of a sphere around ``centre`` filled by scaled Bondi spheres.

    Same Fibonacci 4000 × 12 shell grid as the CS3 scratchpad. Occupancy uses
    squared distances (equivalent to ``norm < r`` at float64 for this grid).
    Numba is used when installed; otherwise a NumPy fallback.
    """
    pts = np.ascontiguousarray(centre + radius * _UNIT_PTS, dtype=np.float64)
    skip = centre_index_guard(centre, coords)
    keep = []
    radii = []
    for i, (symbol, pos) in enumerate(zip(symbols, coords)):
        if i == skip:
            continue
        if not include_h and symbol == "H":
            continue
        r = VDW_BURIAL.get(symbol, 1.7) * scale
        if np.linalg.norm(pos - centre) > radius + r:
            continue
        keep.append(pos)
        radii.append(r)
    if not keep:
        return 0.0
    atoms = np.asarray(keep, dtype=np.float64)
    r2 = np.asarray(radii, dtype=np.float64) ** 2
    buried = _occupy(pts, atoms, r2)
    return 100.0 * float((_VBUR_WEIGHTS * buried).sum() / _VBUR_WEIGHTS.sum())


def donor_indices(symbols):
    """``(M, D1, D2)`` from the canonical metal-then-ancillaries-then-donors order."""
    n = 0
    while 1 + n < len(symbols) and symbols[1 + n] in ANCILLARY_ELEMENTS:
        n += 1
    return 0, 1 + n, 2 + n


def backbone(adj, d1, d2, metal=0):
    """Atoms on the shortest d1→d2 path that does not go through the metal."""
    queue, seen = deque([(d1, [d1])]), {d1, metal}
    while queue:
        u, path = queue.popleft()
        if u == d2:
            return set(path)
        for v in adj[u]:
            if v not in seen:
                seen.add(v)
                queue.append((v, path + [v]))
    return {d1, d2}


def ring_through(adj, a, b, metal=0, cap=8):
    """Members of the smallest ring containing bond a–b, or ``{a, b}``."""
    queue = deque([(b, [a, b])])
    while queue:
        u, path = queue.popleft()
        if len(path) > cap:
            continue
        for v in adj[u]:
            if v == a and len(path) > 2:
                return set(path)
            if v in path or v == metal:
                continue
            queue.append((v, path + [v]))
    return {a, b}


def stereocentre(symbols, adj, donor, other, metal=0, coords=None):
    """``(C*, R)`` for one chelate arm: substituted ring carbon and substituent.

    ``coords`` lets the "alpha" rule break its last tie by distance to the metal; without them
    the first candidate in atom order is kept.
    """
    bb = backbone(adj, donor, other, metal)
    if STEREOCENTRE_RULE == "alpha":
        cand = []
        for carbon in adj[donor]:
            if carbon == metal or symbols[carbon] == "H" or carbon in bb:
                continue
            ring = ring_through(adj, donor, carbon, metal)
            for cc in sorted(ring):
                if cc in bb or symbols[cc] != "C":
                    continue
                for r in adj[cc]:
                    if r in (donor, metal) or r in ring:
                        continue
                    if len(ring) == 2 and len(ring_through(adj, cc, r, metal)) > 2:
                        continue                      # acyclic donor, bond inside a ring: not a substituent axis
                    size = -1 if symbols[r] == "H" else len([k for k in adj[r] if symbols[k] != "H"])
                    far = 0.0 if coords is None else float(np.linalg.norm(coords[r] - coords[metal]))
                    cand.append((-size, cc not in adj[donor], far, len(cand), cc, r))
        return None if not cand else min(cand)[4:]
    best = None
    for carbon in adj[donor]:
        if carbon == metal or symbols[carbon] == "H" or carbon in bb:
            continue
        ring = ring_through(adj, donor, carbon, metal)
        for cc in sorted(ring):
            if cc in bb or symbols[cc] != "C":
                continue
            for r in adj[cc]:
                if r in (donor, metal) or r in ring:
                    continue
                if symbols[r] == "H":
                    if best is None:
                        best = (-1, cc, r)
                    continue
                size = len([k for k in adj[r] if symbols[k] != "H"])
                if best is None or size > best[0]:
                    best = (size, cc, r)
    return None if best is None else (best[1], best[2])


def _donor_frame(coords, metal, d1, d2):
    """Local axes: bisector in the M–D1–D2 plane, normal, desymmetrizing axis."""
    u = coords[d1] - coords[metal]
    v = coords[d2] - coords[metal]
    u = u / np.linalg.norm(u)
    v = v / np.linalg.norm(v)
    bis = u + v
    bis = bis / np.linalg.norm(bis)
    nrm = np.cross(u, v)
    nrm = nrm / np.linalg.norm(nrm)
    des = np.cross(nrm, bis)
    return bis, nrm, des


class MetalComplex:
    """One metal-complex geometry in the CS3 canonical atom order.

    Parameters
    ----------
    symbols : sequence of str
    coords : (n, 3) array
    name : str, optional
    energy : float, optional
        Hartree. Used only when this frame is part of an ensemble average.
    """

    def __init__(self, symbols, coords, name=None, energy=None):
        self.symbols = list(symbols)
        self.coords = np.asarray(coords, dtype=float)
        if self.coords.shape != (len(self.symbols), 3):
            raise ValueError("coords must be (n_atoms, 3) matching symbols")
        self.name = name
        self.energy = energy
        self.adj = adjacency(self.symbols, self.coords)
        self.metal, self.donor_1, self.donor_2 = donor_indices(self.symbols)

    @classmethod
    def from_xyz(cls, filepath, energy_convention: str = "last") -> "MetalComplex":
        """Load the first (or only) frame of an XYZ file."""
        frames = parse_xyz_ensemble(filepath, energy_convention=energy_convention)
        return cls.from_frame(frames[0], name=_stem(filepath))

    @classmethod
    def from_frame(cls, frame: XYZFrame, name=None) -> "MetalComplex":
        energy = None if np.isnan(frame.energy) else float(frame.energy)
        return cls(frame.symbols, frame.coords, name=name, energy=energy)

    @property
    def ancillary_indices(self):
        """Co-ligands sitting between the metal and the two donors (0-based)."""
        return list(range(1, min(self.donor_1, self.donor_2)))

    def geometric_features(self) -> dict:
        """Metal-referenced Sterimol, arm Sterimol, bite, M–D, %Vbur."""
        radii = [PKG_CPK.get(t, 1.92) for t in nob_types(self.symbols, self.adj)]
        metal, d1, d2 = self.metal, self.donor_1, self.donor_2
        xyz = self.coords
        symbols = self.symbols
        out = {}

        for tag, donor, other in (("a", d1, d2), ("b", d2, d1)):
            s = sterimol(symbols, xyz, metal + 1, donor + 1, radii, block=(other + 1,))
            if s:
                for key in STERIMOL_KEYS:
                    out[f"_fromM_{key}_{tag}"] = s[key]
        for radius in (3.0, 3.5, 5.0):
            out[f"vbur_{radius:g}"] = buried_volume(symbols, xyz, xyz[metal], radius=radius)
        out["vbur_noH_3.5"] = buried_volume(
            symbols, xyz, xyz[metal], radius=3.5, include_h=False,
        )
        v1, v2 = xyz[d1] - xyz[metal], xyz[d2] - xyz[metal]
        out["bite"] = np.degrees(np.arccos(np.clip(
            np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)), -1, 1)))
        out["MD_mean"] = (np.linalg.norm(v1) + np.linalg.norm(v2)) / 2
        out["MD_asym"] = abs(np.linalg.norm(v1) - np.linalg.norm(v2))
        lo, hi = BITE_WINDOW
        if not lo <= out["bite"] <= hi:
            warnings.warn(
                f"{self.name or 'complex'}: bite {out['bite']:.1f} deg is outside {lo:g}-{hi:g} deg "
                f"(M-D {np.linalg.norm(v1):.2f} / {np.linalg.norm(v2):.2f} A); a donor has probably "
                "come off the metal, so the sym/asym descriptors compare a bound arm with a free one",
                stacklevel=2,
            )

        for tag, donor in (("a", d1), ("b", d2)):
            other = d2 if donor == d1 else d1
            sc = stereocentre(symbols, self.adj, donor, other, metal=metal, coords=xyz)
            if sc is None:
                continue
            stereo, subst = sc
            if SUB_FRAGMENT_BOUND == "donor_ring":
                block = tuple(sorted({x + 1 for x in donor_ring(self.adj, donor, metal) | {metal}} - {stereo + 1}))
            else:
                block = (donor + 1, metal + 1)
            s = sterimol(symbols, xyz, stereo + 1, subst + 1, radii, block=block)
            if s:
                for key in STERIMOL_KEYS:
                    out[f"_sub_{key}_{tag}"] = s[key]
            u, w = xyz[subst] - xyz[stereo], xyz[donor] - xyz[stereo]
            out[f"_a_R_Cstereo_N_{tag}"] = np.degrees(np.arccos(np.clip(
                np.dot(u, w) / (np.linalg.norm(u) * np.linalg.norm(w)), -1, 1)))

        final = {k: v for k, v in out.items() if not k.startswith("_")}
        bases = {k[1:-2] for k in out if k.startswith("_")}
        for base in sorted(bases):
            pa, pb = out.get(f"_{base}_a"), out.get(f"_{base}_b")
            if pa is None or pb is None:
                continue
            final[f"{base}_sym"] = (pa + pb) / 2
            final[f"{base}_asym"] = abs(pa - pb)
        return final

    def electronic_features(self, xtb: XtbSinglePoint) -> dict:
        """Map an xTB single point onto the metal–donor frame.

        Column definitions match the reverse-engineered ``elec_*.csv`` tables:
        ancillary charges are the co-ligands between metal and donors, not every
        hydrogen in the molecule.
        """
        if len(xtb.charges) != len(self.symbols):
            raise ValueError(
                f"Charge vector length {len(xtb.charges)} != n_atoms {len(self.symbols)}"
            )
        metal, d1, d2 = self.metal, self.donor_1, self.donor_2
        q = xtb.charges
        wbo = xtb.wbo or {}
        dip = xtb.dipole_debye
        bis, nrm, des = _donor_frame(self.coords, metal, d1, d2)
        anc = self.ancillary_indices
        w1 = wbo.get((min(metal, d1), max(metal, d1)), 0.0)
        w2 = wbo.get((min(metal, d2), max(metal, d2)), 0.0)
        w_metal = sum(v for (i, j), v in wbo.items() if metal in (i, j))
        homo = xtb.homo
        lumo = xtb.lumo
        gap = None if homo is None or lumo is None else lumo - homo
        return {
            "mu_bisector": float(dip @ bis),
            "mu_outofplane": float(abs(dip @ nrm) if DIPOLE_ABS else dip @ nrm),
            "mu_desym": float(abs(dip @ des) if DIPOLE_ABS else dip @ des),
            "homo": homo,
            "lumo": lumo,
            "gap": gap,
            "q_metal": float(q[metal]),
            "q_donor_sym": float((q[d1] + q[d2]) / 2),
            "q_donor_asym": float(abs(q[d1] - q[d2])),
            "q_anc_sum": float(q[anc].sum()) if anc else 0.0,
            "q_absmax": float(np.abs(q).max()),
            "q_spread": float(q.std()),
            "wbo_MD_sym": float((w1 + w2) / 2),
            "wbo_MD_asym": float(abs(w1 - w2)),
            "wbo_metal_tot": float(w_metal),
        }

    def features(self, xtb: XtbSinglePoint | None = None) -> dict:
        """Geometric block, plus electronics when ``xtb`` is given."""
        out = self.geometric_features()
        if xtb is not None:
            out.update(self.electronic_features(xtb))
        return out


class MetalComplexEnsemble:
    """GOAT/CREST multi-XYZ of a metal complex, Boltzmann-averaged at 298.15 K.

    Parameters
    ----------
    filepath : str
        ``*.ens.xyz``, ``*.finalensemble.xyz``, CREST ``*_conformers.xyz``, or
        a single-structure XYZ.
    energy_convention : {'last', 'first', 'auto'}
        Default ``last`` (GOAT-safe). Pass ``first`` for CREST.
    relative_zero_threshold : float or None
        If set to ``1.0``, CREST relative-zero dumps are marked unusable.
    temperature : float
    name : str, optional
    """

    def __init__(
        self,
        filepath,
        energy_convention: str = "last",
        relative_zero_threshold: float | None = None,
        temperature: float = 298.15,
        name: str | None = None,
    ):
        self.ensemble = XYZEnsemble(
            filepath,
            energy_convention=energy_convention,
            temperature=temperature,
            relative_zero_threshold=relative_zero_threshold,
            molecule_name=name,
        )
        self.name = self.ensemble.molecule_name
        self.complexes = [
            MetalComplex.from_frame(frame, name=f"{self.name}_conf{i + 1}")
            for i, frame in enumerate(self.ensemble.frames)
        ]

    @classmethod
    def from_xyz(cls, filepath, **kwargs) -> "MetalComplexEnsemble":
        return cls(filepath, **kwargs)

    def __len__(self) -> int:
        return len(self.complexes)

    @property
    def n_conformers(self) -> int:
        return len(self.complexes)

    def lowest(self) -> MetalComplex:
        return self.complexes[self.ensemble.lowest_index]

    def geometric_features(self) -> dict:
        """Per-frame geometry, Boltzmann-averaged. Adds ``n_conformers``."""
        usable = []
        for complex_ in self.complexes:
            try:
                row = complex_.geometric_features()
            except Exception:
                continue
            if row:
                usable.append((row, complex_.energy))
        if not usable:
            return {}
        if len(usable) == 1:
            out = dict(usable[0][0])
            out["n_conformers"] = 1
            return out
        rows, energies = [], []
        for row, energy in usable:
            if energy is None or (isinstance(energy, (float, np.floating)) and np.isnan(energy)):
                continue
            rows.append(row)
            energies.append(energy)
        averaged = boltzmann_average(rows, energies, temperature=self.ensemble.temperature)
        if not averaged:
            return {}
        averaged = dict(averaged)
        averaged["n_conformers"] = len(rows)
        return averaged

    def electronic_features(self, xtb_by_frame: dict) -> dict:
        """Boltzmann-average electronics. ``xtb_by_frame`` maps frame index → SP.

        Frames without a matching single point are skipped.
        """
        rows, energies = [], []
        for i, complex_ in enumerate(self.complexes):
            sp = xtb_by_frame.get(i)
            if sp is None:
                continue
            try:
                row = complex_.electronic_features(sp)
            except Exception:
                continue
            rows.append(row)
            energies.append(complex_.energy if complex_.energy is not None else np.nan)
        return boltzmann_average(rows, energies, temperature=self.ensemble.temperature)

    def features(self, xtb_by_frame: dict | None = None) -> dict:
        out = self.geometric_features()
        if xtb_by_frame:
            out.update(self.electronic_features(xtb_by_frame))
        return out


class MetalComplexSet:
    """Batch of :class:`MetalComplex` / :class:`MetalComplexEnsemble` objects."""

    def __init__(self, members: dict):
        self.members = dict(members)

    def __len__(self) -> int:
        return len(self.members)

    @classmethod
    def from_xyz_dir(
        cls,
        directory,
        pattern: str = "*.xyz",
        energy_convention: str = "last",
        relative_zero_threshold: float | None = None,
        temperature: float = 298.15,
    ) -> "MetalComplexSet":
        """Load every matching XYZ as an ensemble (n=1 is fine)."""
        members = {}
        for path in sorted(glob.glob(os.path.join(directory, pattern))):
            ens = MetalComplexEnsemble(
                path,
                energy_convention=energy_convention,
                relative_zero_threshold=relative_zero_threshold,
                temperature=temperature,
            )
            members[ens.name] = ens
        return cls(members)

    @classmethod
    def from_paths(cls, paths, **kwargs) -> "MetalComplexSet":
        members = {}
        for path in paths:
            ens = MetalComplexEnsemble(path, **kwargs)
            members[ens.name] = ens
        return cls(members)

    def geometric_dataframe(self) -> pd.DataFrame:
        """One row per complex; index is the ligand name."""
        rows = {}
        for name, member in self.members.items():
            if isinstance(member, MetalComplexEnsemble):
                row = member.geometric_features()
            else:
                row = member.geometric_features()
                row = dict(row)
                row["n_conformers"] = 1
            if row:
                rows[name] = row
        if not rows:
            return pd.DataFrame()
        frame = pd.DataFrame.from_dict(rows, orient="index")
        if "n_conformers" in frame.columns:
            cols = ["n_conformers"] + [c for c in frame.columns if c != "n_conformers"]
            frame = frame[cols]
        return frame

    def electronic_dataframe(self, xtb_by_name: dict) -> pd.DataFrame:
        """``xtb_by_name`` maps ligand name → :class:`XtbSinglePoint` (n=1 SP)."""
        rows = {}
        for name, member in self.members.items():
            sp = xtb_by_name.get(name)
            if sp is None:
                continue
            complex_ = member.lowest() if isinstance(member, MetalComplexEnsemble) else member
            try:
                rows[name] = complex_.electronic_features(sp)
            except Exception:
                continue
        if not rows:
            return pd.DataFrame()
        return pd.DataFrame.from_dict(rows, orient="index")


def _stem(filepath) -> str:
    return XYZEnsemble._name_from_path(filepath)
