"""A new compound's descriptors by the recipe that made a prepared case's columns, read under the environment that has
the DescriPyTor readers (the server may run under another one).

    <reader python> featurize.py              JSON on stdin {extraction, paths}; JSON on stdout {x: {column: value}} or {error}
    <reader python> featurize.py check <case> <folder> [names...]
                                              re-extract measured compounds of a case from their own files and compare with the case's X

`extraction` is the field of a case file that says which recipe made its columns: its `kind`, the files it `accepts`, the
`companion` file it needs beside the structure (if any), `settings`, how the case's X was rounded (`round`), the scripts
it runs with the hashes it was checked with (`code`), and the source `table` the check compares with. Nothing here builds or optimises a geometry: the user supplies the structure, computed as the set's were.
A structure whose atoms cannot play the roles the columns assume is refused, with the reason.
"""
import json, re, sys
from pathlib import Path
import numpy as np

from settings import S

sys.path.insert(0, S["descripytor"])
unit = lambda v: v / np.linalg.norm(v)


class Refused(ValueError):
    """the structure does not fit the recipe: the reason goes back to the page"""


# ------------------------------------------------------------------ cs1_acids: benzoic acids, Gaussian opt+freq log
# notes/cs1_validation/benzoic_acids/extract_acids.py (table explorer_in/tables/cs1_acids.csv) and extract_frames.py; their
# log reader (aryl_bromides/glog.py) and role finder in rule; Sterimol is the package's metal_complex.sterimol.
EL = {1: "H", 5: "B", 6: "C", 7: "N", 8: "O", 9: "F", 14: "Si", 15: "P", 16: "S", 17: "Cl", 35: "Br", 53: "I"}; EV = 27.211386


def _glog(path):
    """the last Standard orientation (Input orientation for a NoSymm log), the last field-independent dipole (Debye), the last
    NPA and Hirshfeld blocks, HOMO/LUMO of the last alpha eigenvalues in eV, the last isotropic polarizability (a.u.)"""
    L = Path(path).read_text(errors="replace").splitlines()
    def last(key):
        k = [i for i, l in enumerate(L) if key in l]
        if not k: raise Refused("the log has no '%s' block (the recipe needs an opt+freq log with pop=NBO, Hirshfeld and polar)" % key.strip())
        return k[-1]
    k = last("Standard orientation:" if any("Standard orientation:" in l for l in L) else "Input orientation:") + 5; sym, X = [], []
    while not L[k].startswith(" ---"): t = L[k].split(); sym.append(EL[int(t[1])]); X.append([float(v) for v in t[3:6]]); k += 1
    n = len(sym); mu = np.array([float(v) for v in re.findall(r"[XYZ]=\s*(-?\d+\.\d+)", L[last("Dipole moment (field-independent basis, Debye)") + 1])])
    k = last("Summary of Natural Population Analysis:") + 6; nbo = [float(L[k + i].split()[2]) for i in range(n)]; assert [L[k + i].split()[0] for i in range(n)] == sym
    k = last("Hirshfeld charges, spin densities, dipoles, and CM5 charges") + 2; hir = [float(L[k + i].split()[2]) for i in range(n)]; assert [L[k + i].split()[1] for i in range(n)] == sym
    k = last("Alpha  occ. eigenvalues"); homo = float(L[k].split("--")[1].split()[-1]) * EV; lumo = float(next(l for l in L[k:] if "Alpha virt. eigenvalues" in l).split("--")[1].split()[0]) * EV
    iso = [float(l.split()[-2]) for l in L if "Isotropic polarizability for W=" in l]
    if not iso: raise Refused("the log has no isotropic polarizability (polar keyword)")
    return dict(sym=sym, X=np.array(X), mu=mu, nbo=nbo, hirshfeld=hir, homo=homo, lumo=lumo, iso=iso[-1]), L


def _modes(L, n):
    """frequencies and Cartesian displacements of the last 3n-6 modes printed"""
    F, D = [], []; i = 0
    while i < len(L):
        if L[i].startswith(" Frequencies --"):
            f = [float(v) for v in L[i].split("--")[1].split()]; j = next(k for k in range(i, i + 12) if L[k].lstrip().startswith("Atom  AN")); blk = np.array([[float(v) for v in L[j + 1 + a].split()[2:]] for a in range(n)])
            for m in range(len(f)): F.append(f[m]); D.append(blk[:, 3 * m:3 * m + 3])
            i = j + n
        i += 1
    if len(F) < 3 * n - 6: raise Refused("the log has %d vibrational modes, not 3N-6 = %d: a frequency job is needed" % (len(F), 3 * n - 6))
    return np.array(F[-(3 * n - 6):]), np.array(D[-(3 * n - 6):])


def _acid_roles(sym, X):
    """0-based C_acid, C_ipso, O=, O-H, H_acid by distance bonds (1.25 A with an H, else 1.75 A): the case's atoms 1 to 5"""
    n = len(sym); nb = lambda i: [j for j in range(n) if j != i and np.linalg.norm(X[i] - X[j]) < (1.25 if "H" in (sym[i], sym[j]) else 1.75)]
    acid = [i for i in range(n) if sym[i] == "C" and sum(sym[j] == "O" for j in nb(i)) == 2 and any(sym[j] == "O" and any(sym[k] == "H" for k in nb(j)) for j in nb(i))]
    if len(acid) != 1: raise Refused("%d carboxylic-acid carbons (a C with two O, one of them O-H); the columns need exactly one" % len(acid))
    cc = acid[0]; oh = next(j for j in nb(cc) if sym[j] == "O" and any(sym[k] == "H" for k in nb(j))); oc = next(j for j in nb(cc) if sym[j] == "O" and j != oh)
    ha = next(k for k in nb(oh) if sym[k] == "H"); ci = [j for j in nb(cc) if sym[j] == "C"]
    if not ci: raise Refused("the acid carbon has no carbon neighbour (atom 2 of the columns is the ring carbon it sits on)")
    if sum(sym[j] == "C" for j in nb(ci[0])) < 2: raise Refused("the carbon under the acid group (atom %d) has fewer than two carbon neighbours: not an aryl acid" % (ci[0] + 1))
    return cc, ci[0], oc, oh, ha, nb


def _ring(sym, nb, cc, ci):
    def walk(p):                                          # extract_frames.ring: the smallest heavy-atom ring through C_ipso, 5 or 6 atoms
        best = None
        for j in nb(p[-1]):
            if sym[j] == "H" or j == cc: continue
            if j == ci and len(p) >= 5: best = min([best, p], key=lambda r: 99 if r is None else len(r))
            elif j not in p and len(p) < 6: best = min([best, walk(p + [j])], key=lambda r: 99 if r is None else len(r))
        return best
    return walk([ci])


def acid_log(path, settings, extra):
    from M2_data_extractor.metal_complex import adjacency, nob_types, sterimol, PKG_CPK
    text = Path(path).read_text(errors="replace")
    if text.count("Normal termination") < 2: raise Refused("the log does not end two jobs normally (the set's logs are opt+freq, both terminated normally)")
    if re.search(r"\d+ imaginary frequencies", text): raise Refused("the structure has an imaginary frequency: not a minimum")
    g, L = _glog(path); sym, X, mu = g["sym"], g["X"], g["mu"]; cc, ci, oc, oh, ha, nb = _acid_roles(sym, X)
    F, D = _modes(L, len(sym)); b = unit(X[oc] - X[cc]); win = np.where((F > 1650) & (F < 1950))[0]           # C=O stretch: in 1650-1950, the mode moving O= against C most along C=O
    if not len(win): raise Refused("no mode between 1650 and 1950 cm-1 for the C=O stretch")
    nu = float(F[win[int(np.argmax([abs((D[m][oc] - D[m][cc]) @ b) for m in win]))]])
    u = unit(X[ci] - X[cc]); par = float(mu @ u)
    s = sterimol(sym, X, cc + 1, ci + 1, [PKG_CPK.get(t, 1.92) for t in nob_types(sym, adjacency(sym, X))])     # every atom, no block
    out = {"nbo_atom_2": g["nbo"][ci], "nbo_atom_1": g["nbo"][cc], "nbo_atom_3": g["nbo"][oc], "nbo_atom_5": g["nbo"][ha],
           "hirshfeld_atom_2": g["hirshfeld"][ci], "hirshfeld_atom_1": g["hirshfeld"][cc], "hirshfeld_atom_3": g["hirshfeld"][oc], "hirshfeld_atom_5": g["hirshfeld"][ha],
           "0_bond_length_1-2": float(np.linalg.norm(X[ci] - X[cc])), "0_bond_length_1-3": float(np.linalg.norm(X[oc] - X[cc])), "nu_CO": nu,
           "mu_total": float(np.linalg.norm(mu)), "mu_y[C1->C2|O3]": par, "mu_perp": float(np.linalg.norm(mu - par * u)),
           "homo": g["homo"], "lumo": g["lumo"], "gap": g["lumo"] - g["homo"], "B1_1-2": s["B1"], "B5_1-2": s["B5"], "L_1-2": s["L"], "iso": g["iso"]}
    if not settings.get("frames"): return out
    R = _ring(sym, nb, cc, ci)
    if not R: raise Refused("no five- or six-membered ring through the carbon under the acid group")
    y = b; x = unit((X[ci] - X[cc]) - ((X[ci] - X[cc]) @ y) * y)                                            # carbonyl frame
    c = X[R].mean(0); nrm = np.linalg.svd(X[R] - c)[2][2]; yr = unit((X[ci] - c) - ((X[ci] - c) @ nrm) * nrm); xr = np.cross(yr, nrm); xr = xr if (X[oc] - c) @ xr > 0 else -xr   # ring frame, signed towards O=
    return {**out, "mu_y[C1->O3|C2]": float(mu @ y), "mu_x[C1->O3|C2]": float(mu @ x), "mu_x[ring->C2]": float(mu @ xr)}


# ------------------------------------------------------------------ lucas_pdh2: (N,N)PdH2 complexes, Gaussian opt+freq log or its feather
# notes/analysis/lucas_check: pdh2_extract.py (table pdh2_interim_fit30_2026-10-13.csv), filtered by pdh2_explorer_case.py build.
# The recipe runs those scripts' own functions (settings.scripts, held to the hashes in `code`).
def pdh2_log(path, settings, extra):
    import contextlib, io, os
    sys.path.insert(0, settings["scripts"])
    import verify_features as VF, extract_v6 as E, pdh2_extract as PE, pdh2_numbering_check as NC
    from M2_data_extractor.feather_extractor import gauss_file_handler, save_to_feather
    p = Path(path)
    if p.suffix.lower() == ".log":
        if p.read_text(errors="ignore").count("Normal termination") != 2: raise Refused("the log does not end two jobs normally (the set's logs are opt+freq)")
        df, msg = gauss_file_handler(str(p))
        if df is None or df.empty: raise Refused("the log could not be read: %s" % msg)
        cwd = os.getcwd(); os.chdir(p.parent)
        try:
            with contextlib.redirect_stdout(io.StringIO()): save_to_feather(df, p.stem)
        finally: os.chdir(cwd)
        p = p.with_suffix(".feather")
    PE.FE = p.parent                                                    # check() reads <PE.FE>/<name>.feather
    try: bad = NC.check(p.stem)["problems"]
    except Exception as e: bad = "%s: %s" % (type(e).__name__, e)
    if bad: raise Refused("the atoms are not numbered as the set's (the PdMeCl v7 numbering, pdh2_atom_map.json: C1 stereocentre, N5 oxazoline, N11 pyridine, Pd12, hydrides 13 and 14): %s" % bad)
    m = VF.load(p); X = m.xyz_df[["x", "y", "z"]].to_numpy(float)
    if VF.ang(X, 11, 12, 13) < VF.ang(X, 11, 12, 14):                  # as pdh2_extract.one: H13 is the hydride trans to N11
        with contextlib.redirect_stdout(io.StringIO()): m.renumber_atoms({13: 14, 14: 13}, bonds_are_0_based=False, save_feather=False)
        X = m.xyz_df[["x", "y", "z"]].to_numpy(float)
    d = lambda a, b: float(np.linalg.norm(X[a - 1] - X[b - 1]))
    if not (d(12, 13) < 1.70 and d(12, 14) < 1.70 and d(13, 14) > 1.2 and float(m.info_df["Frequency"].min()) > 0): raise Refused("not a valid dihydride minimum (Pd-H under 1.70 A, H-H over 1.2 A, no negative frequency), the rule the set's shapes were chosen by")
    import pandas as pd
    old = [c for c in pd.read_csv(Path(settings["scripts"]) / "lucas_features_v7_2026-10-08.csv", nrows=1).columns if c.startswith(("hirshfeld_atom_", "nbo_atom_", "bond_length_", "angle_"))]
    with contextlib.redirect_stdout(io.StringIO()): return PE.old_columns(m, old) | E.featurise(m, r1=PE.R1, stretches=PE.STRETCHES) | PE.metal(m)


# ------------------------------------------------------------------ cs3_*: Corminboeuf chelate complexes, xyz + its GFN2 single point
# zenodo single_structure_xtb/<build>/<lig>_<build>_xtbopt.xyz and the cluster's single points on them (CORMINBOEUF/cheap, the
# compact Q/W/D/E text of cs3_review/cluster/collect.sh) -> cs3_review/arms/build_pools.py (geometric_features with phi,
# electronic_features) -> bassa/build_wide.py (chelate placements, sym = mean of the arms, asym = |a - b|) -> filter50.py (columns only)
BUILD = {"cp": "CuCl", "oa": "CuCl", "cc": "Ni", "da_f": "CuCl2"}; METAL = {"CuCl": "Cu", "CuCl2": "Cu", "Ni": "Ni"}; COLIG = {"CuCl": ["Cl"], "CuCl2": ["Cl", "Cl"], "Ni": ["H", "H"]}


def _sp(path):
    """the compact GFN2 single point: Q charges, W metal Wiberg orders (`j wbo`, j 1-based), D dipole 'full:' in a.u., E HOMO LUMO in eV"""
    from M2_data_extractor.xtb_singlepoint import XtbSinglePoint
    d = {l[0]: [float(v) for v in l.split()[1:]] for l in Path(path).read_text(errors="replace").splitlines() if l[:2] in ("Q ", "W ", "D ", "E ")}
    if set("QWDE") - set(d): raise Refused("the single-point file lacks its %s line: it is the Q/W/D/E text cs3_review/cluster/collect.sh writes from xtb --sp" % ", ".join(sorted(set("QWDE") - set(d))))
    return XtbSinglePoint(charges=d["Q"], wbo={(0, int(d["W"][i]) - 1): d["W"][i + 1] for i in range(0, len(d["W"]) - 1, 2)}, dipole=d["D"], homo=d["E"][0], lumo=d["E"][1], dipole_unit="au")


def _topology(m):
    """tf_*: the package's LigandTopology indices per heavy atom, on the free ligand's graph read from the complex (metal and
    co-ligands left out, Kier-Hall delta_v from the H count); equal to the deposit's SMILES route on every complex the roles accept"""
    from M2_data_extractor.ligand_topology import _chi, _paths_of_length, _dist_matrix, ZV, Z, TF_FROM_T
    lig = [i for i in range(min(m.donor_1, m.donor_2), len(m.symbols)) if m.symbols[i] != "H"]; pos = {a: k for k, a in enumerate(lig)}
    adj = [[pos[j] for j in m.adj[a] if j in pos] for a in lig]; bonds = [(i, j) for i in range(len(lig)) for j in adj[i] if i < j]; delta = [len(a) for a in adj]; dv = []
    for a in lig:
        s, h = m.symbols[a], sum(m.symbols[j] == "H" for j in m.adj[a]); den = Z.get(s, 6) - ZV.get(s, 4) - 1
        dv.append((ZV.get(s, 4) - h) / den if den > 0 else float(ZV.get(s, 4) - h))
    n = float(len(lig)); t = {}
    for k in range(5): t[f"chi{k}"], t[f"chi{k}v"] = _chi(adj, delta, k), _chi(adj, dv, k)
    p = {k: len(_paths_of_length(adj, k)) for k in (1, 2, 3)}
    t["kappa1"] = n * (n - 1) ** 2 / p[1] ** 2; t["kappa2"] = (n - 1) * (n - 2) ** 2 / p[2] ** 2
    t["kappa3"] = ((n - 1) * (n - 3) ** 2 if len(lig) % 2 else (n - 3) * (n - 2) ** 2) / p[3] ** 2
    t["zagreb1"] = float(sum(x * x for x in delta)); t["zagreb2"] = float(sum(delta[i] * delta[j] for i, j in bonds))
    D = _dist_matrix(adj); fin = D[np.triu_indices(len(lig), 1)]; fin = fin[np.isfinite(fin)]
    t["wiener"] = float(fin.sum()); t["global_eff"] = float((1.0 / fin[fin > 0]).mean())
    ds = np.where(np.isfinite(D), D, 0).sum(axis=1); mu = len(bonds) - len(lig) + 1
    t["balaban"] = len(bonds) / (mu + 1) * sum(1.0 / np.sqrt(ds[i] * ds[j]) for i, j in bonds)
    A = np.zeros((len(lig), len(lig)))
    for i, j in bonds: A[i, j] = A[j, i] = 1.0
    t["estrada"] = float(np.exp(np.linalg.eigvalsh(A)).sum())
    return {f"tf_{k}": float(t[k]) / n for k in TF_FROM_T}


def _chelate_roles(m, reaction):
    """why a complex cannot play the roles of the set's (the metal first, then its co-ligands, then the two donors), or None"""
    from M2_data_extractor import metal_complex as mc, placement_search as ps
    b = BUILD[reaction]; S = m.symbols
    if S[0] != METAL[b]: return "atom 1 is %s: the %s complexes have %s first" % (S[0], b, METAL[b])
    anc = [S[i] for i in m.ancillary_indices]
    if anc != COLIG[b]: return "the atoms after the metal are %s; the %s build has %s there" % (anc, b, COLIG[b])
    if any(i not in m.adj[0] for i in m.ancillary_indices): return "a co-ligand is not bonded to the metal"
    d1, d2 = m.donor_1, m.donor_2
    for d in (d1, d2):
        if S[d] not in ("N", "P", "O", "S"): return "atom %d (%s) should be a donor (N, P, O or S) right after the co-ligands" % (d + 1, S[d])
        if d not in m.adj[0]: return "donor %s%d is not bonded to the metal (%.2f A)" % (S[d], d + 1, np.linalg.norm(m.coords[d] - m.coords[0]))
    pm = ps.Mol("x", S, m.coords)                                       # the placements find the donors by distance: they must be the same two
    dn = sorted(sorted([j for j in pm.adj[0] if S[j] in ("N", "P", "O", "S")], key=lambda j: np.linalg.norm(m.coords[j] - m.coords[0]))[:2])
    if dn != [d1, d2]: return "the two donors nearest the metal are atoms %s, not %d and %d" % ([x + 1 for x in dn], d1 + 1, d2 + 1)
    for tag, d, o in (("first", d1, d2), ("second", d2, d1)):
        if mc.stereocentre(S, m.adj, d, o, metal=0, coords=m.coords) is None: return "the %s arm (donor %s%d) has no stereocentre C* with a substituent R" % (tag, S[d], d + 1)
    v1, v2 = m.coords[d1] - m.coords[0], m.coords[d2] - m.coords[0]; bite = np.degrees(np.arccos(np.clip(v1 @ v2 / np.linalg.norm(v1) / np.linalg.norm(v2), -1, 1))); lo, hi = mc.BITE_WINDOW
    if not lo <= bite <= hi: return "the bite angle is %.1f deg, outside %g-%g" % (bite, lo, hi)
    return None


def chelate_xyz(path, settings, extra):
    import warnings
    from M2_data_extractor import metal_complex as mc, placement_search as ps
    if not (mc.STERIMOL_SCAN_STEP == 1 and mc.STERIMOL_FRAME == "fragment" and mc.STERIMOL_THETA_RULE == "soft" and mc.SUB_FRAGMENT_BOUND == "donor_ring" and mc.STEREOCENTRE_RULE == "alpha" and mc.DIPOLE_ABS):
        raise Refused("the package's MetalComplex settings are not the ones the CS3 columns were made with (1 deg scan, fragment frame, soft theta, donor-ring bound, alpha rule, unsigned dipoles)")
    if len(extra) != 1: raise Refused("the CS3 columns need two files: the complex's .xyz and its GFN2 single point (the Q/W/D/E text); got %d beside the structure" % len(extra))
    m = mc.MetalComplex.from_xyz(path); why = _chelate_roles(m, settings["reaction"])
    if why: raise Refused(why)
    sp = _sp(extra[0])
    if len(sp.charges) != len(m.symbols): raise Refused("the single point has %d atoms and the structure %d: not the same molecule" % (len(sp.charges), len(m.symbols)))
    keys, mc.STERIMOL_KEYS = mc.STERIMOL_KEYS, ("B1", "B5", "L", "theta", "angle")        # build_pools.py: phi ("angle") as well
    try:
        with warnings.catch_warnings(): warnings.simplefilter("ignore"); g = m.geometric_features()
    finally: mc.STERIMOL_KEYS = keys
    out = {k.replace("_angle_", "_phi_"): float(v) for k, v in g.items()}
    out.update({k: float(v) for k, v in m.electronic_features(sp).items() if k != "mu_outofplane"}); out.update(_topology(m))
    F, _ = ps.compute(ps.Dataset.from_xyz([path], ["x"]), ps.chelate_placements())
    for key in ("CsR", "DCs", "MD", "MCs", "MR", "DR", "RRp"):
        for q in ("B1", "B5", "L", "loc_B5", "theta", "phi"):
            a, b = (float(F[f"{key}_{t}_{q}"].iloc[0]) if f"{key}_{t}_{q}" in F else np.nan for t in "ab"); out[f"{key}_{q}_sym"], out[f"{key}_{q}_asym"] = (a + b) / 2, abs(a - b)
    return out


# ------------------------------------------------------------------ cs2: Ackermann oxazoline ligands, Gaussian feather or log
# Two routes, both the package's: the deposited pool (ackermann_features_sterimol_regenerated.csv: data_extractor.Molecule's
# charges, lengths, angles, and Sterimol 6->4 on the 1.82 A bond rule, regenerate_cs12_sterimol.py), and the build_wide.py
# columns (CS2 block of `organic`: payload charges, the dipole in the frame of each scaffold bond, placement_search Sterimol
# with the soft tie). The set fixed which columns exist; no value reads the rest of the set.
ACK_ROLE = {1: {"O"}, 2: {"C"}, 3: {"C", "N"}, 4: {"C"}, 5: {"N"}, 6: {"C"}, 7: {"C"}, 12: {"C"}, 15: {"O"}, 17: {"C"}, 21: {"H"}}
ACK_BONDS = [(1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (6, 12), (12, 15), (4, 15), (2, 17), (1, 21)]


def _ack_roles(S, X, de, ps):
    """why atoms 1-7, 12, 15, 17 and 21 cannot carry the roles of the set's (H21-O1-C2(-C17)=X3-C4(=N5-C6(H)(-C7)-C12-O15-C4)), or None"""
    import pandas as pd
    if len(S) < max(ACK_ROLE): return "only %d atoms" % len(S)
    for i, el in ACK_ROLE.items():
        if S[i - 1] not in el: return "atom %d is %s; the set has %s there" % (i, S[i - 1], " or ".join(sorted(el)))
    m = ps.Mol("x", S, X); b182 = {tuple(sorted(r)) for r in de.extract_connectivity(pd.DataFrame({"atom": S, "x": X[:, 0], "y": X[:, 1], "z": X[:, 2]}), threshold_distance=1.82).to_numpy().astype(int).tolist()}
    for a, b in ACK_BONDS:
        if (b - 1) not in m.adj[a - 1] or (min(a, b), max(a, b)) not in b182: return "atoms %d and %d are not bonded, as they are in every ligand of the set" % (a, b)
    if len(m.adj[5]) != 4: return "C6 has %d neighbours; the set's have 4 (N5, C7, C12, H)" % len(m.adj[5])
    if np.linalg.det(np.array([X[4] - X[5], X[6] - X[5], X[11] - X[5]])) > 0: return "C6 has the other configuration than every ligand of the set (the mirror image turns every mu_x and mu_z column round)"
    for o, a in ((6, 7), (2, 17)):
        hit = {1, 2, 3, 4, 5, 6, 12, 15} & {i + 1 for i in ps.fragment(m, o, a)}
        if hit: return "the substituent %d->%d reaches the scaffold atoms %s" % (o, a, sorted(hit))
    return None


def ackermann(path, settings, extra):
    import contextlib, importlib.util, io, warnings
    import pandas as pd
    warnings.filterwarnings("ignore")
    sys.path.insert(0, str(Path(S["descripytor"]) / "M2_data_extractor"))
    from M2_data_extractor import placement_search as ps
    spec = importlib.util.spec_from_file_location("fe_data", str(Path(S["descripytor"]) / "M2_data_extractor/feature_explorer/data.py")); TD = importlib.util.module_from_spec(spec); spec.loader.exec_module(TD)
    with contextlib.redirect_stdout(io.StringIO()):
        import data_extractor as de
        m = TD._molecule(path)
    p = TD.payload(path); L = p["xyz"].splitlines()[2:]; sym = [l.split()[0] for l in L]; X = np.array([[float(v) for v in l.split()[1:4]] for l in L])
    if p["dipole"] is None or not {"nbo", "hirshfeld", "cm5"} <= set(p["charges"] or {}): raise Refused("the file carries no dipole, or not all of the NBO, Hirshfeld and CM5 charges the columns read")
    why = _ack_roles(sym, X, de, ps)
    if why: raise Refused(why)
    mu = np.array(p["dipole"][:3]); mol = ps.Mol("x", sym, X, mu); xyz = m.xyz_df[["atom", "x", "y", "z"]].copy(); b182 = None; ster, out = {}, {}
    num = lambda s: [int(v) for v in re.findall(r"\d+", s)]
    for c in settings["columns"]:
        if c == "aniso":                                                 # not the anisotropy: the reader's slot for Gaussian's isotropic polarizability (bohr^3), as in the table
            if path.lower().endswith(".log"): v = float(m.polarizability_df["iso"].iloc[0])
            else: v = float(m.polarizability_df["iso" if "iso" in pd.read_feather(path).columns else "aniso"].iloc[0])
        elif c == "mu_total": v = np.linalg.norm(mu)
        elif r := re.fullmatch(r"(nbo|hirshfeld|cm5)_atom_(\d+)", c): v = m.charge_dict[r[1]]["charge"].to_numpy(float)[int(r[2]) - 1]
        elif r := re.fullmatch(r"0_bond_length_(\d+)-(\d+)", c): v = m.get_bond_length([[int(r[1]), int(r[2])]]).iloc[0, 0]
        elif r := re.fullmatch(r"0_(angle|dihedral)_\[(.*)\]", c): v = m.get_bond_angle(num(r[2])).iloc[0, 0]
        elif r := re.fullmatch(r"(B1|B5|L|loc_B5|B1_B5_angle)_(\d+)-(\d+)", c):
            if b182 is None: b182 = de.extract_connectivity(xyz, threshold_distance=1.82)
            key = (int(r[2]), int(r[3]))
            if key not in ster:
                with contextlib.redirect_stdout(io.StringIO()): ster[key] = de.get_sterimol_df(xyz, b182, list(key), None, radii="CPK").iloc[0]
            v = ster[key][r[1]]
        elif r := re.fullmatch(r"q_(nbo|hirshfeld|cm5)\[[A-Z][a-z]?(\d+)\]", c): v = p["charges"][r[1]][int(r[2]) - 1]
        elif r := re.fullmatch(r"mu_([xyz])\[[A-Z][a-z]?(\d+)->[A-Z][a-z]?(\d+)\|[A-Z][a-z]?(\d+)\]", c):
            o, a, t = (int(r[k]) - 1 for k in (2, 3, 4)); u = X[a] - X[o]; u /= np.linalg.norm(u); w = np.cross(u, X[t] - X[o]); w /= np.linalg.norm(w)
            v = (np.array([np.cross(w, u), u, -w]) @ mu)["xyz".index(r[1])]
        elif r := re.fullmatch(r"ster_(B1|B5|L|loc_B5|theta|phi)\[[A-Z][a-z]?(\d+)->[A-Z][a-z]?(\d+)\]", c):
            o, a = int(r[2]), int(r[3])
            if ("ps", o, a) not in ster: ster[("ps", o, a)] = ps.sterimol_tie(mol, o, a, ps.fragment(mol, o, a), tie="soft")
            v = ster[("ps", o, a)][r[1]]
        else: raise KeyError("no rule for column " + c)
        out[c] = float(v)
    return out


KINDS = {"acid_log": acid_log, "pdh2_log": pdh2_log, "chelate_xyz": chelate_xyz, "ackermann": ackermann}


def featurize(extraction, paths, raw=False):
    """the columns of one structure, rounded as the case's own X was (`round`: each rounding in turn, decimals or "6g"
    for 6 significant digits, or {column: those, "*": the rest}; default ["6g"], as build_case.compute writes X), so a
    measured compound added as new is
    priced exactly at its fitted value. paths: the structure file, and the file the recipe needs with it (`companion`)"""
    paths = [str(p) for p in ([paths] if isinstance(paths, (str, Path)) else paths)]; main = [p for p in paths if Path(p).suffix.lower() in extraction["accepts"]]
    if len(main) != 1: raise Refused("this case reads one %s file per compound%s; got %s" % (" or ".join(extraction["accepts"]), ", with " + extraction["companion"] if extraction.get("companion") else "", ", ".join(Path(p).name for p in paths) or "none"))
    extra = [p for p in paths if p != main[0]]; path = main[0]
    if extra and not extraction.get("companion"): raise Refused("this case reads the structure file alone; drop only the %s file" % " or ".join(extraction["accepts"]))
    import hashlib
    for f, h in extraction.get("code", {}).items():                    # the scripts a recipe runs must be the ones it was checked with
        if not Path(f).is_file() or hashlib.sha256(Path(f).read_bytes()).hexdigest() != h: raise Refused("%s has changed (or is gone) since this case's recipe was checked: re-check it (featurize.py check) before featurizing" % Path(f).name)
    x = {k: float(v) for k, v in KINDS[extraction["kind"]](path, extraction.get("settings", {}), extra).items()}
    r = extraction.get("round", ["6g"]); r = r if isinstance(r, dict) else {"*": r}
    def rnd(c, v):
        for k in r.get(c, r["*"]): v = float("%.6g" % v) if k == "6g" else float(np.round(v, k))     # numpy, as pandas rounded the tables: 1.20075 goes to 1.2008
        return v
    return x if raw else {k: rnd(k, v) for k, v in x.items()}


def case_file(case):
    from settings import CASES_DIR
    t = (CASES_DIR / f"{case}.js").read_text(encoding="utf-8"); return json.loads(t[t.index("] = ") + 4:].rstrip().rstrip(";"))


def check(case, folder, names=()):
    """re-extract compounds of a case (measured and external) from their own structures in a folder (<name>.<suffix>, and a
    companion as <name>.<another suffix>): per compound the
    largest difference from the case's source table (`extraction.table`, full precision; relative above 1), and that
    from the case's X of the rounded values, which is 0 when the recipe holds"""
    import pandas as pd
    d = case_file(case); e = d["extraction"]; T = pd.read_csv(e["table"], index_col=0) if e.get("table") else None; out = {}
    rows = {n: (d["X"], d["ligands"].index(n)) for n in d["ligands"]}; rows.update({n: (d["ext"]["X"], d["ext"]["names"].index(n)) for n in (d.get("ext") or {}).get("names", [])})
    files = {p.stem: p for p in Path(folder).iterdir() if p.suffix.lower() in e["accepts"]}
    together = lambda n: [files[n]] + [p for p in Path(folder).iterdir() if p.stem == n and p != files[n]]          # a companion is <name>.<another suffix>
    for n in names or [n for n in rows if n in files]:
        x = featurize(e, together(n), raw=True); y = featurize(e, together(n)); X, i = rows[n]; rel = lambda a, b: abs(a - b) / max(1.0, abs(b))
        out[n] = dict(table=None if T is None or n not in T.index else max(rel(x[c], float(T.loc[n, c])) for c in d["names"]), case_X=max(abs(y[c] - X[c][i]) for c in d["names"]))
    return out


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    if sys.argv[1:2] == ["check"]: print(json.dumps(check(sys.argv[2], sys.argv[3], sys.argv[4:]))); sys.exit()
    import contextlib
    q = json.loads(sys.stdin.read())
    try:
        with contextlib.redirect_stdout(sys.stderr): out = dict(x=featurize(q["extraction"], q["paths"]))      # what the readers print would spoil the answer
    except Refused as e: out = dict(error=str(e))
    print(json.dumps(out))
