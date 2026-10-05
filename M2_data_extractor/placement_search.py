"""Placement search: where on a series of molecules a Sterimol or dipole descriptor should be read, with the controls that
keep the answer honest.

A *placement* is how one descriptor is anchored in one molecule: an origin atom, an axis atom and the atoms the substituent
fragment may not cross (Sterimol), or an axis / ring frame (dipole). The same placement is resolved in every molecule, either
from a shared atom numbering (organic series) or from chemical roles (chelate arms of metal complexes).

    ds  = Dataset.from_xyz(files)                        # or Dataset.from_feathers(files) to carry the Gaussian dipole
    roles = audit_roles(ds, [1, 3, 8, 23])               # same element and neighbourhood in every molecule?
    P   = bond_placements(ds) + ring_placements(ds, anchor=23)      # or chelate_placements()
    F, flags = compute(ds, P)                            # molecules x features, and per-placement sanity flags
    res = rank_placements(F, y, base=F0, candidates=[c for c in F if c.endswith('_theta')])

Conventions (they differ from the older extractor paths on purpose):
  bonding  covalent radii x 1.3 (a fixed 1.82 A cutoff drops e.g. the S-CF3 bond of a triflate, 1.84 A);
  fragment everything on the axis atom's side, never crossing the origin, its declared blocks, or any ring that contains the
           origin but not the axis atom (a fused substituent would otherwise walk back through the ring into the scaffold);
  angles   theta (tilt of the B5 vector out of the B1 plane) and phi (in-plane B1-B5 angle) are ill-defined when several B1
           directions are within a few hundredths of an angstrom (tBu has three); tie='soft' weights every direction and every
           B5 atom by exp(-excess/0.02 A), tie='mean' averages the near-tied local minima; a lone-atom fragment gives 0.
Validation (rank_placements): each candidate is scored in the slot of a declared base model; the score of the best candidate
is compared with the best of the same number of random vectors, with a selection null (the whole choice repeated on permuted
responses), and with nested leave-one-out (the choice repeated inside every fold). A candidate declared before looking is
compared with single random vectors instead.
"""
from __future__ import annotations
import itertools
from collections import deque
from dataclasses import dataclass, field
from typing import Callable, Optional
import numpy as np
import pandas as pd
try:
    from .metal_complex import CONNECT_COV, PKG_CPK, nob_types, stereocentre
except ImportError:                                           # run as a plain module
    from metal_complex import CONNECT_COV, PKG_CPK, nob_types, stereocentre

TOL, TAU = 0.02, 0.02

# ------------------------------------------------------------------ molecules
@dataclass
class Mol:
    name: str
    symbols: list
    coords: np.ndarray
    dipole: Optional[np.ndarray] = None           # Debye, same frame as coords
    adj: list = field(default=None, repr=False)
    radii: list = field(default=None, repr=False)

    def __post_init__(self):
        X = self.coords; n = len(self.symbols); D = np.linalg.norm(X[:, None] - X[None], axis=-1)
        cov = np.array([CONNECT_COV.get(s, 0.77) for s in self.symbols])
        B = (D < 1.3 * (cov[:, None] + cov[None])) & ~np.eye(n, dtype=bool)
        self.adj = [list(np.where(B[i])[0]) for i in range(n)]
        self.radii = [PKG_CPK.get(t, 1.92) for t in nob_types(self.symbols, self.adj)]


class Dataset(list):
    @classmethod
    def from_xyz(cls, files, names=None):
        out = cls()
        for k, f in enumerate(files):
            L = open(f).read().split("\n"); n = int(L[0].split()[0])
            S = [l.split()[0] for l in L[2:2 + n]]; X = np.array([[float(v) for v in l.split()[1:4]] for l in L[2:2 + n]])
            out.append(Mol(names[k] if names else str(f).split("/")[-1].rsplit(".", 1)[0], S, X))
        return out

    @classmethod
    def from_feathers(cls, files, names=None):
        """Gaussian output via DescriPyTor's Molecule: coordinates and the dipole vector in the same frame."""
        import io, contextlib
        try:
            from .data_extractor import Molecule
        except ImportError:
            from data_extractor import Molecule
        out = cls()
        for k, f in enumerate(files):
            with contextlib.redirect_stdout(io.StringIO()):
                m = Molecule(str(f))
            S = list(m.xyz_df["atom"]); X = m.xyz_df[["x", "y", "z"]].to_numpy(float)
            dp = m.gauss_dipole_df; mu = dp[["dip_x", "dip_y", "dip_z"]].to_numpy(float)[0] if len(dp) else None
            out.append(Mol(names[k] if names else str(f).split("/")[-1].rsplit(".", 1)[0], S, X, mu))
        return out

# ------------------------------------------------------------------ roles
def audit_roles(ds, atoms):
    """1-based atoms -> table of element + sorted neighbour elements per molecule; `consistent` marks the majority signature."""
    rows = []
    for m in ds:
        r = {"name": m.name}
        for a in atoms:
            i = a - 1
            r[a] = "?" if i >= len(m.symbols) else m.symbols[i] + "(" + "".join(sorted(m.symbols[j] for j in m.adj[i])) + ")"
        rows.append(r)
    T = pd.DataFrame(rows).set_index("name")
    for a in atoms:
        T[f"{a}_ok"] = T[a] == T[a].mode()[0]
    T["consistent"] = T[[f"{a}_ok" for a in atoms]].all(axis=1)
    return T

def ring_of(m, anchor, avoid=()):
    """smallest ring through 1-based `anchor` (heavy atoms, not through `avoid`), 1-based, in ring order from the anchor"""
    a = anchor - 1; avoid = {x - 1 for x in avoid}; H = lambda i: [j for j in m.adj[i] if m.symbols[j] != "H" and j not in avoid]
    best = None
    for s in H(a):
        prev = {s: None}; q = deque([s])
        while q:
            x = q.popleft()
            for y in H(x):
                if y == a and x != s:
                    p = [x]
                    while prev[p[-1]] is not None: p.append(prev[p[-1]])
                    ring = [a] + p[::-1]
                    if best is None or len(ring) < len(best): best = ring
                elif y not in prev and y != a: prev[y] = x; q.append(y)
    return None if best is None else [i + 1 for i in best]

# ------------------------------------------------------------------ placements
@dataclass
class Placement:
    name: str
    kind: str                       # 'sterimol' | 'dipole_axis' | 'dipole_ring'
    resolve: Callable               # Mol -> tuple of 1-based atoms (or None when the molecule has no such placement)

def core_atoms(ds):
    n = min(len(m.symbols) for m in ds)
    return [i + 1 for i in range(n) if len({m.symbols[i] for m in ds}) == 1 and ds[0].symbols[i] != "H"]

def core_bonds(ds):
    C = core_atoms(ds)
    return [(a, b) for a, b in itertools.combinations(C, 2) if all((b - 1) in m.adj[a - 1] for m in ds)]

def bond_placements(ds, sterimol=True, dipole=True):
    """both directions of every heavy-atom bond present in all molecules (shared numbering)"""
    P = []
    for a, b in core_bonds(ds):
        if sterimol:
            P += [Placement(f"{a}->{b}", "sterimol", lambda m, a=a, b=b: (a, b)), Placement(f"{b}->{a}", "sterimol", lambda m, a=a, b=b: (b, a))]
        if dipole:
            P.append(Placement(f"mu[{a}->{b}]", "dipole_axis", lambda m, a=a, b=b: (a, b)))
    return P

def ring_placements(ds, anchor, plane_ref=None):
    """dipole ring frame on the ring through `anchor`, found per molecule: origin = ring centroid, u = centroid -> anchor,
    v in the ring plane towards the ring neighbour of the anchor nearest to `plane_ref` (1-based; default: first neighbour)"""
    def res(m):
        r = ring_of(m, anchor)
        if r is None: return None
        nb = [r[1], r[-1]]
        pa = nb[0] if plane_ref is None else min(nb, key=lambda k: np.linalg.norm(m.coords[k - 1] - m.coords[plane_ref - 1]))
        return (tuple(r), anchor, pa)
    return [Placement(f"ring@{anchor}", "dipole_ring", res)]

def chelate_placements():
    """seven Sterimol placements on each chelate arm of a metal complex, resolved by role (metal, donors, stereocentre C*,
    substituent R, next atom R'); names end in _a / _b for the two arms"""
    def roles(m, arm):
        M = next(i for i, s in enumerate(m.symbols) if s in ("Cu", "Ni", "Pd", "Co", "Fe", "Rh", "Ir", "Zn", "Ag", "Au", "Pt", "Ru"))
        dn = sorted([j for j in m.adj[M] if m.symbols[j] in ("N", "P", "O", "S")], key=lambda j: np.linalg.norm(m.coords[j] - m.coords[M]))[:2]
        if len(dn) < 2: return None
        dn = sorted(dn); D, Do = (dn[0], dn[1]) if arm == "a" else (dn[1], dn[0])
        sc = stereocentre(m.symbols, m.adj, D, Do, metal=M, coords=m.coords)
        if sc is None: return None
        C, R = sc; ring = ring_of(m, D + 1, avoid=(M + 1,)) or [D + 1]; ring0 = {x - 1 for x in ring}
        heavy = [k for k in m.adj[R] if k != C and m.symbols[k] != "H" and k not in ring0]
        Rn = max(heavy, key=lambda k: len(m.adj[k])) if heavy else None
        imine = [x for x in m.adj[D] if x in ring0 and x != C]
        return dict(M=M, D=D, Do=Do, C=C, R=R, Rn=Rn, ring=ring0, imine=imine)
    spec = {"CsR": lambda r: (r["C"], r["R"], r["ring"] | {r["M"]}), "DCs": lambda r: (r["D"], r["C"], {r["M"], r["Do"], *r["imine"]}),
            "MD": lambda r: (r["M"], r["D"], {r["Do"]}), "MCs": lambda r: (r["M"], r["C"], {r["D"], r["Do"]}),
            "MR": lambda r: (r["M"], r["R"], r["ring"] | {r["C"]}), "DR": lambda r: (r["D"], r["R"], r["ring"] | {r["M"], r["C"]}),
            "RRp": lambda r: None if r["Rn"] is None else (r["R"], r["Rn"], r["ring"] | {r["C"]})}
    P = []
    for key, f in spec.items():
        for arm in ("a", "b"):
            def res(m, f=f, arm=arm):
                r = roles(m, arm)
                if r is None: return None
                t = f(r)
                return None if t is None else (t[0] + 1, t[1] + 1, tuple(sorted(x + 1 for x in t[2] - {t[0], t[1]})))
            P.append(Placement(f"{key}_{arm}", "sterimol", res))
    return P

# ------------------------------------------------------------------ features
def fragment(m, o, a, block=()):
    """0-based atoms on a's side of o->a, never crossing o, `block`, or a ring that holds o but not a"""
    o, a = o - 1, a - 1; stop = {o} | {b - 1 for b in block}
    for s in m.adj[o]:                                           # rings through o that avoid a
        if s == a: continue
        r = ring_of(m, o + 1, avoid=(a + 1,))
        if r: stop |= {x - 1 for x in r} - {o}
        break
    out = {a}; q = deque([a])
    while q:
        x = q.popleft()
        for y in m.adj[x]:
            if y not in stop and y not in out: out.add(y); q.append(y)
    return sorted(out)

def sterimol_tie(m, o, a, idx, tie="soft"):
    X = m.coords; origin = X[o - 1]; u = X[a - 1] - origin; u /= np.linalg.norm(u)
    V = X[idx] - origin; perp = V - np.outer(V @ u, u); r = np.array([m.radii[i] for i in idx])
    e1 = perp[int(np.argmax((perp ** 2).sum(1)))]
    if np.linalg.norm(e1) < 1e-8: e1 = np.cross(u, [0, 0, 1.0])
    e1 /= np.linalg.norm(e1); e2 = np.cross(u, e1); P = np.c_[V @ e1, V @ e2]
    th = np.radians(np.arange(0, 360, 0.1)); dirs = np.c_[np.cos(th), np.sin(th)]
    ext = (P @ dirs.T + r[:, None]).max(0); B1 = float(ext.min())
    b5 = np.linalg.norm(P, axis=1) + r; k5 = int(b5.argmax()); B5 = float(b5[k5]); L = float(((V @ u) + r).max()); loc = float(V[k5] @ u)
    res = dict(B1=B1, B5=B5, L=L, loc_B5=loc, theta=0.0, phi=0.0, n_tied=1, n_frag=len(idx))
    if len(idx) == 1 or np.linalg.norm(P, axis=1).max() < 1e-6: return res
    Pn = P / np.maximum(np.linalg.norm(P, axis=1, keepdims=True), 1e-12)
    phi = np.degrees(np.arccos(np.clip(dirs @ Pn.T, -1, 1)))
    theta = np.degrees(np.arcsin(np.clip(np.abs(dirs @ P.T) / np.linalg.norm(V, axis=1), 0, 1)))
    n = len(ext); mins = np.array([i for i in range(n) if ext[i] <= ext[i - 1] and ext[i] <= ext[(i + 1) % n] and ext[i] < B1 + TOL])
    near5 = np.where(b5 > B5 - TOL)[0]; res["n_tied"] = int(len(mins))
    if tie == "mean":
        res["theta"] = float(theta[np.ix_(mins, near5)].mean()); res["phi"] = float(phi[np.ix_(mins, near5)].mean())
    else:
        W = np.outer(np.exp(-(ext - B1) / TAU), np.exp(-(B5 - b5) / TAU)); res["theta"] = float((W * theta).sum() / W.sum()); res["phi"] = float((W * phi).sum() / W.sum())
    return res

def compute(ds, placements, tie="soft"):
    """features (molecules x '<placement>_<quantity>') and per-placement flags"""
    rows, flag = {}, {}
    for p in placements:
        vals, sizes, ties, missing = {}, [], [], 0
        for m in ds:
            t = p.resolve(m)
            if t is None: missing += 1; continue
            if p.kind == "sterimol":
                o, a, *blk = t; idx = fragment(m, o, a, blk[0] if blk else ())
                s = sterimol_tie(m, o, a, idx, tie); sizes.append(s["n_frag"]); ties.append(s["n_tied"])
                for q in ("B1", "B5", "L", "loc_B5", "theta", "phi"): vals.setdefault(q, {})[m.name] = s[q]
            elif m.dipole is not None:
                if p.kind == "dipole_axis":
                    a, b = t; u = m.coords[b - 1] - m.coords[a - 1]; vals.setdefault("mu", {})[m.name] = float(m.dipole @ u / np.linalg.norm(u))
                else:
                    ring, a, pa = t; c = m.coords[[i - 1 for i in ring]].mean(0); u = m.coords[a - 1] - c; u /= np.linalg.norm(u)
                    w = np.cross(u, m.coords[pa - 1] - c); w /= np.linalg.norm(w); v = np.cross(w, u)
                    for q, e in (("mu_u", u), ("mu_v", v), ("mu_w", -w)): vals.setdefault(q, {})[m.name] = float(m.dipole @ e)
        for q, d in vals.items(): rows[f"{p.name}_{q}" if not q.startswith("mu") or p.kind == "dipole_ring" else p.name] = d
        sz = np.array(sizes) if sizes else np.array([0])
        flag[p.name] = dict(kind=p.kind, missing=missing, frag_min=int(sz.min()), frag_max=int(sz.max()),
                            frag_spread=float(sz.max() / max(np.median(sz), 1)), tied_B1_share=float(np.mean(np.array(ties) > 1)) if ties else 0.0)
    return pd.DataFrame(rows).loc[[m.name for m in ds]], pd.DataFrame(flag).T

def aggregate(F, pairs, how=("sym", "asym"), fill=0.0):
    """combine equivalent placements (e.g. the two arms: pairs = [('CsR_a','CsR_b'), ...]) into sym / asym columns.
    A placement a molecule does not have (an arm whose substituent is a lone H has no R') contributes `fill`."""
    out = {}
    for a, b in pairs:
        stem = a.rsplit("_", 1)[0]
        for q in [c[len(a) + 1:] for c in F.columns if c.startswith(a + "_")]:
            x, z = F[f"{a}_{q}"].fillna(fill), F[f"{b}_{q}"].fillna(fill)
            if "sym" in how: out[f"{stem}_{q}_sym"] = (x + z) / 2
            if "asym" in how: out[f"{stem}_{q}_asym"] = (x - z).abs()
    return pd.DataFrame(out)

# ------------------------------------------------------------------ ranking with controls
def _score(Xc, y, ridge):
    n = len(y); sd = Xc.std(0, ddof=1); sd[sd < 1e-12] = 1; Z = np.c_[np.ones(n), (Xc - Xc.mean(0)) / sd]
    R = np.eye(Z.shape[1]) * ridge; R[0, 0] = 0; P = Z @ np.linalg.inv(Z.T @ Z + R) @ Z.T
    p = y - (y - P @ y) / (1 - np.diag(P)); q = 1 - ((y - p) ** 2).sum() / ((y - y.mean()) ** 2).sum(); k = Z.shape[1] - 1
    return 1 - (1 - q) * (n - 1) / (n - k - 1)

def _fitpred(Xc, y, tr, te, ridge):
    A = Xc[tr]; mu, sd = A.mean(0), A.std(0, ddof=1); sd[sd < 1e-12] = 1; Z = np.c_[np.ones(len(tr)), (A - mu) / sd]
    R = np.eye(Z.shape[1]) * ridge; R[0, 0] = 0; w = np.linalg.solve(Z.T @ Z + R, Z.T @ y[tr])
    return np.c_[np.ones(len(te)), (Xc[te] - mu) / sd] @ w

def rank_placements(F, y, candidates, base=None, ridge=1.0, n_random=2000, n_perm=200, declared=None, seed=0):
    """Score each candidate column in the slot of `base` (DataFrame of fixed columns, may be None) by adjusted Q2_LOO
    (ridge on standardized columns; ridge=0 gives OLS). Returns (table, summary):
      random_best_of_N  does the chosen placement add more than the best of N random vectors would (N = candidates)?
                        -> the test of the placement itself;
      selection_null    does the whole equation, placement chosen, beat chance?  (base + best candidate on permuted y)
      nested            is the choice stable when it is made without the held-out molecule?
      declared          for a placement named before looking: its rank and its percentile against single random vectors."""
    y = np.asarray(y, float); n = len(y); B = np.zeros((n, 0)) if base is None else np.asarray(base, float)
    cand = [c for c in candidates if F[c].std() > 1e-9]
    s = {c: _score(np.c_[B, F[c].to_numpy(float)], y, ridge) for c in cand}
    T = pd.DataFrame({"adjQ2": s}).sort_values("adjQ2", ascending=False); T["rank"] = np.arange(1, len(T) + 1)
    rng = np.random.default_rng(seed); N = len(cand); best = float(T.adjQ2.iloc[0])
    rv_best = np.array([max(_score(np.c_[B, rng.standard_normal(n)], y, ridge) for _ in range(N)) for _ in range(max(1, n_random // max(N, 1)))])
    Xall = F[cand].to_numpy(float)
    null = np.array([max(_score(np.c_[B, Xall[:, j]], yp, ridge) for j in range(N)) for yp in (rng.permutation(y) for _ in range(n_perm))])
    nested, picks = np.empty(n), []
    for i in range(n):
        tr = np.delete(np.arange(n), i)
        j = int(np.argmax([_score(np.c_[B[tr], Xall[tr, k]], y[tr], ridge) for k in range(N)])); picks.append(cand[j])
        nested[i] = _fitpred(np.c_[B, Xall[:, j]], y, tr, np.array([i]), ridge)[0]
    base_only = _score(B, y, ridge) if B.shape[1] else np.nan
    summ = dict(n=n, n_candidates=N, best=T.index[0], best_adjQ2=best, base_adjQ2=float(base_only) if B.shape[1] else None,
                random_best_of_N=dict(median=float(np.median(rv_best)), p95=float(np.percentile(rv_best, 95)), percentile_of_best=float((rv_best < best).mean() * 100)),
                selection_null=dict(p95=float(np.percentile(null, 95)), p=float((null >= best).mean())),
                nested=dict(R2=float(1 - ((y - nested) ** 2).sum() / ((y - y.mean()) ** 2).sum()), same_as_full=f"{sum(p == T.index[0] for p in picks)}/{n}",
                            picks=pd.Series(picks).value_counts().to_dict()))
    if declared is not None:
        one = np.array([_score(np.c_[B, rng.standard_normal(n)], y, ridge) for _ in range(n_random)])
        summ["declared"] = dict(name=declared, adjQ2=float(s[declared]), rank=int(T.loc[declared, "rank"]), random_percentile=float((one < s[declared]).mean() * 100))
    return T, summ

# ------------------------------------------------------------------ self-check
if __name__ == "__main__":
    import sys, glob, os, json
    arc = sys.argv[1]                                        # the paper's archive root
    # 1. roles: the ring frame the v4 CS1 deposit declared as atoms 23-28 is not the ring in three products
    feathers = sorted(glob.glob(f"{arc}/raw/case_study_1_mazet/*.feather"))
    ds1 = Dataset.from_feathers(feathers)
    bad = [m.name for m in ds1 if set(ring_of(m, 23)) != set(range(23, 29))]
    assert sorted(bad) == ["pyridine-m-OMe", "thiophene", "unsub"], bad
    # 2. dipole ring frame reproduces the rebuilt CS1 matrix
    Fr, _ = compute(ds1, ring_placements(ds1, anchor=23, plane_ref=8))
    M = pd.read_csv(f"{arc}/feature_matrices/mazet_features.csv").set_index("name").loc[Fr.index]
    assert np.allclose(Fr["ring@23_mu_v"], M["dipole_v_ring"], atol=1e-5) and np.allclose(Fr["ring@23_mu_u"], M["dipole_u_ring"], atol=1e-5)
    # 3. chelate arms reproduce the CS3 corrected placements (fused rings bounded: 13 atoms on 024)
    files = sorted(glob.glob(f"{arc}/geometries/case_study_3_corminboeuf/single_structure_xtb/CuCl/*_CuCl_xtbopt.xyz"))
    cp = pd.read_csv(f"{arc}/feature_matrices/single_structure/cp_angles_soft.csv").set_index("name")
    ds3 = Dataset.from_xyz([f for f in files if os.path.basename(f).split("_CuCl")[0] in cp.index], None)
    for m in ds3: m.name = m.name.split("_CuCl")[0]
    F3, flags3 = compute(ds3, chelate_placements(), tie="soft")
    A3 = aggregate(F3, [(f"{k}_a", f"{k}_b") for k in ("CsR", "DCs", "MD", "MCs", "MR", "DR", "RRp")]).loc[cp.index]
    assert np.allclose(A3["CsR_B1_sym"], cp["CsR_B1_sym"], atol=2e-3) and np.allclose(A3["RRp_theta_sym"], cp["RRp_theta_sym"], atol=0.5)
    assert flags3.loc[["CsR_a", "CsR_b"], "frag_max"].max() == 29  # the quinoline arms (037-039) are still flagged by size
    # 4. ranking with controls: theta over the seven placements in the slot of B1 + HOMO
    base = pd.concat([A3["CsR_B1_sym"], cp["homo"]], axis=1)
    cand = [c for c in A3.columns if c.endswith(("_theta_sym", "_phi_sym", "_theta_asym", "_phi_asym"))]
    T, S = rank_placements(A3, cp["ddg"].to_numpy(float), cand, base=base, n_random=2000, n_perm=100)
    print(T.head(5).round(3).to_string()); print(json.dumps(S, indent=1, default=str))
    print("self-check passed")
