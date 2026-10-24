"""Fold conformers that are the same shape after identical atoms trade labels.

A nitro's two oxygens, a CF3's fluorines, or the two faces of a symmetric para
ring can be relabeled and still be one geometry. A flip that puts a different
substituent on the other side of the scaffold does not overlay, and stays.

Auto-fold only pairs at or below 0.05 Angstrom. Between 0.05 and 0.15 Angstrom
the core torsions must agree within 15 degrees, and that pair is offered for
acceptance rather than hidden. Above 0.15 Angstrom the frames stay apart.
"""
from __future__ import annotations

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds, rdMolAlign

DUP_A = 0.15
TIGHT_A = 0.05
TORSION_DEG = 15.0
MAX_MATCHES = 128

_RADII = {
    "H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "P": 1.07,
    "S": 1.05, "Cl": 1.02, "Br": 1.20, "I": 1.39, "B": 0.84, "Si": 1.11,
}


def _symbol(el) -> str:
    s = str(el).strip()
    if not s:
        return "C"
    if len(s) == 1:
        return s.upper()
    return s[0].upper() + s[1:].lower()


def _coord(value):
    if value is None:
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if not np.isfinite(number):
        return None
    return number


def _row_xyz(row):
    if not row or len(row) < 4:
        return None
    coords = [_coord(row[1]), _coord(row[2]), _coord(row[3])]
    if any(v is None for v in coords):
        return None
    return coords


def heavy_rows(atoms):
    return [a for a in atoms if a and _row_xyz(a) is not None and _symbol(a[0]) != "H"]


def heavy_rows_idx(atoms):
    rows, idx = [], []
    for i, row in enumerate(atoms or []):
        if not row or _symbol(row[0]) == "H" or _row_xyz(row) is None:
            continue
        rows.append(row)
        idx.append(i)
    return rows, idx


def _fingerprint(rows):
    """Element key and sorted distances from the centre. Invariant to rotation and to atom order."""
    coords = [_row_xyz(r) for r in rows]
    if len(coords) < 3 or any(c is None for c in coords):
        return None
    xyz = np.array(coords, float)
    elem = [_symbol(r[0]) for r in rows]
    centre = xyz.mean(0)
    radial = np.linalg.norm(xyz - centre, axis=1)
    parts = []
    for el in sorted(set(elem)):
        parts.append(np.sort(radial[np.array(elem) == el]))
    return (tuple(sorted(elem)), np.concatenate(parts))


def _radial_lb(fa, fb):
    if fa is None or fb is None or fa[0] != fb[0]:
        return None
    return float(np.sqrt(np.mean((fa[1] - fb[1]) ** 2)))


def mol_from_heavy(heavy):
    """A heavy-atom mol with one conformer. Bonds from DetermineBonds, else distances."""
    rw = Chem.RWMol()
    for row in heavy:
        rw.AddAtom(Chem.Atom(_symbol(row[0])))
    base = rw.GetMol()
    conf = Chem.Conformer(len(heavy))
    for i, row in enumerate(heavy):
        conf.SetAtomPosition(i, (float(row[1]), float(row[2]), float(row[3])))
    base.AddConformer(conf, assignId=True)
    try:
        trial = Chem.Mol(base)
        rdDetermineBonds.DetermineBonds(trial, charge=0)
        if trial.GetNumBonds() >= max(1, len(heavy) - 3):
            return trial
    except Exception:
        pass
    em = Chem.RWMol(base)
    for i in range(len(heavy)):
        for j in range(i + 1, len(heavy)):
            ri = _RADII.get(_symbol(heavy[i][0]), 0.77)
            rj = _RADII.get(_symbol(heavy[j][0]), 0.77)
            d = float(np.linalg.norm(
                np.array(heavy[i][1:4], float) - np.array(heavy[j][1:4], float)
            ))
            if 0.4 < d < (ri + rj) * 1.3:
                em.AddBond(i, j, Chem.BondType.SINGLE)
    mol = em.GetMol()
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        mol.UpdatePropertyCache(strict=False)
        Chem.GetSymmSSSR(mol)
    return mol


def _xyz(mol) -> np.ndarray:
    conf = mol.GetConformer()
    return np.array(
        [[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y, conf.GetAtomPosition(i).z]
         for i in range(mol.GetNumAtoms())],
        float,
    )


def core_torsions(mol):
    """Rotatable heavy-atom bonds outside rings. Each end uses its lowest-index heavy neighbour."""
    mol.UpdatePropertyCache(strict=False)
    Chem.GetSymmSSSR(mol)
    out = []
    for bond in mol.GetBonds():
        if bond.IsInRing():
            continue
        a, b = bond.GetBeginAtom(), bond.GetEndAtom()
        if a.GetAtomicNum() == 1 or b.GetAtomicNum() == 1:
            continue
        left = [n.GetIdx() for n in a.GetNeighbors() if n.GetAtomicNum() > 1 and n.GetIdx() != b.GetIdx()]
        right = [n.GetIdx() for n in b.GetNeighbors() if n.GetAtomicNum() > 1 and n.GetIdx() != a.GetIdx()]
        if not left or not right:
            continue
        out.append((min(left), a.GetIdx(), b.GetIdx(), min(right)))
    return out


def _dihedral(xyz, i, j, k, l) -> float:
    p = xyz[[i, j, k, l]]
    b0, b1, b2 = p[0] - p[1], p[2] - p[1], p[3] - p[2]
    n1 = np.linalg.norm(b1)
    if n1 < 1e-12:
        return 0.0
    u = b1 / n1
    v = b0 - np.dot(b0, u) * u
    w = b2 - np.dot(b2, u) * u
    return float(np.degrees(np.arctan2(np.dot(np.cross(u, v), w), np.dot(v, w))))


def _circ(a, b) -> float:
    d = b - a
    while d > 180.0:
        d -= 360.0
    while d < -180.0:
        d += 360.0
    return abs(d)


def _best_rms(probe, ref):
    try:
        if probe.GetNumAtoms() != ref.GetNumAtoms() or probe.GetNumAtoms() < 3:
            return None
        return float(rdMolAlign.GetBestRMS(probe, ref))
    except Exception:
        return None


def _bare(mol):
    """Same atoms and bonds, every bond single, so a Kekule ring still has its flip."""
    rw = Chem.RWMol(mol)
    for bond in rw.GetBonds():
        bond.SetBondType(Chem.BondType.SINGLE)
    out = rw.GetMol()
    out.UpdatePropertyCache(strict=False)
    return out


def _maps(probe, ref):
    try:
        found = _bare(ref).GetSubstructMatches(_bare(probe), uniquify=False, maxMatches=MAX_MATCHES)
    except Exception:
        found = ()
    ident = tuple(range(probe.GetNumAtoms()))
    maps = [ident]
    for m in found:
        if m != ident and m not in maps:
            maps.append(m)
    return maps


def _kabsch_rmsd(P, Q) -> float:
    """RMSD after the best rotation of P onto Q. Rows already correspond."""
    P = np.asarray(P, float)
    Q = np.asarray(Q, float)
    Pc, Qc = P.mean(0), Q.mean(0)
    A = (P - Pc).T @ (Q - Qc)
    V, _S, Wt = np.linalg.svd(A)
    if np.linalg.det(V @ Wt) < 0:
        Wt = Wt.copy()
        Wt[-1] *= -1
    diff = (P - Pc) @ (V @ Wt) - (Q - Qc)
    return float(np.sqrt(np.mean(np.sum(diff * diff, axis=1))))


def _mapped_rms(probe, ref, amap, xyz_probe=None, xyz_ref=None) -> float:
    order = [int(j) for j in amap]
    xp = _xyz(probe) if xyz_probe is None else xyz_probe
    xr = _xyz(ref) if xyz_ref is None else xyz_ref
    return _kabsch_rmsd(xp, xr[order])


def _with_coords(template, rows):
    mol = Chem.Mol(template)
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i, row in enumerate(rows):
        conf.SetAtomPosition(i, (float(row[1]), float(row[2]), float(row[3])))
    mol.RemoveAllConformers()
    mol.AddConformer(conf, assignId=True)
    return mol


def _torsion_span(xyz_a, xyz_b, torsions, amap, da=None) -> float:
    if not torsions:
        return 0.0
    worst = 0.0
    for t, (i, j, k, l) in enumerate(torsions):
        a = _dihedral(xyz_a, i, j, k, l) if da is None else da[t]
        db = _dihedral(xyz_b, amap[i], amap[j], amap[k], amap[l])
        worst = max(worst, _circ(a, db))
    return worst


def _swaps(amap, orig):
    """1-based atom numbers that are not left in place, each pair once."""
    if not amap or not orig:
        return []
    seen, out = set(), []
    for i, j in enumerate(amap):
        if int(j) == i:
            continue
        a, b = orig[i] + 1, orig[int(j)] + 1
        key = (min(a, b), max(a, b))
        if key in seen:
            continue
        seen.add(key)
        out.append([key[0], key[1]])
    return out[:8]


def judge_pair(mol_a, mol_b, rmsd_cut=DUP_A, tight=TIGHT_A, torsion_cut=TORSION_DEG, maps=None):
    """kind is tight, torsion (close, torsions agree), rejected, apart, or unmatched.

    The atom map is the numbering of mol_b that best overlays mol_a. Identical
    numbering is checked first, so an unpermuted copy never enumerates symmetries.
    """
    n = mol_a.GetNumAtoms()
    if mol_b.GetNumAtoms() != n or n < 3:
        return {"kind": "unmatched", "rmsd": None, "amap": None, "torsion": None}
    ident = tuple(range(n))
    xyz_a, xyz_b = _xyz(mol_a), _xyz(mol_b)
    try:
        ident_rms = _mapped_rms(mol_a, mol_b, ident, xyz_a, xyz_b)
    except Exception:
        ident_rms = None
    if ident_rms is not None and ident_rms <= tight:
        return {"kind": "tight", "rmsd": ident_rms, "amap": ident, "torsion": 0.0}
    maps = maps if maps is not None else _maps(mol_a, mol_b)
    # every relabelling's RMSD once: their minimum is the symmetry-best RMSD, the value RDKit's GetBestRMS would
    # enumerate again; RDKit stays the judge only when the list of relabellings was cut at MAX_MATCHES
    scored = []
    for amap in maps:
        try:
            scored.append((_mapped_rms(mol_a, mol_b, amap, xyz_a, xyz_b), amap))
        except Exception:
            continue
    rms = min((r for r, _m in scored), default=None) if len(maps) < MAX_MATCHES else _best_rms(mol_a, mol_b)
    if rms is None:
        return {"kind": "unmatched", "rmsd": None, "amap": None, "torsion": None}
    if rms > rmsd_cut:
        return {"kind": "apart", "rmsd": rms, "amap": None, "torsion": None}
    tors = core_torsions(mol_a)
    da = [_dihedral(xyz_a, *t) for t in tors]       # the first molecule's torsions once, not once per relabelling
    best_tight, best_grey = None, None
    for mapped, amap in scored:
        span = _torsion_span(xyz_a, xyz_b, tors, amap, da)
        if mapped <= tight and (best_tight is None or mapped < best_tight[0]):
            best_tight = (mapped, amap, span)
            break
        if mapped <= rmsd_cut and span <= torsion_cut and (best_grey is None or mapped < best_grey[0]):
            best_grey = (mapped, amap, span)
    if best_tight is not None or rms <= tight:
        chosen = best_tight or (rms, None, None)
        return {"kind": "tight", "rmsd": chosen[0], "amap": chosen[1], "torsion": chosen[2]}
    if best_grey is not None:
        return {"kind": "torsion", "rmsd": best_grey[0], "amap": best_grey[1], "torsion": best_grey[2]}
    return {"kind": "rejected", "rmsd": rms, "amap": None, "torsion": None}


def same_conformer(mol_a, mol_b, rmsd_cut=DUP_A, tight=TIGHT_A, torsion_cut=TORSION_DEG):
    """(is_duplicate, rmsd or None, reason). reason is tight, torsion, or rejected."""
    judged = judge_pair(mol_a, mol_b, rmsd_cut, tight, torsion_cut)
    return judged["kind"] in ("tight", "torsion"), judged["rmsd"], judged["kind"]


def _record(child, parent, judged, orig):
    return {
        "child": child,
        "parent": parent,
        "rmsd": judged["rmsd"],
        "reason": judged["kind"],
        "swaps": _swaps(judged["amap"], orig[child]),
        "torsion": None if judged["torsion"] is None else round(float(judged["torsion"]), 1),
    }


def fold_duplicates(frames, rmsd_cut=DUP_A, tight=TIGHT_A, torsion_cut=TORSION_DEG, should_stop=None):
    """frames: [{atoms: [[el,x,y,z], ...], energy: float|None}, ...].

    pairs are auto-folded (at or below 0.05 Angstrom). suggestions are the
    0.05–0.15 Angstrom pairs whose core torsions agree; they stay in the list
    until accepted. A frame is compared with RDKit only when the
    rotation-invariant distance fingerprint is inside the cutoff of a cluster.
    """
    order = sorted(
        range(len(frames)),
        key=lambda i: (frames[i].get("energy") is None, float(frames[i].get("energy") or 0.0)),
    )
    prepared = [heavy_rows_idx(fr.get("atoms") or []) for fr in frames]
    prints = [_fingerprint(rows) for rows, _idx in prepared]
    graphs, built = {}, {}

    def mol(i):
        if i in built:
            return built[i]
        key = prints[i][0] if prints[i] else None
        if key is None:
            built[i] = (None, None)
            return built[i]
        if key not in graphs:
            template = mol_from_heavy(prepared[i][0])
            graphs[key] = (template, _maps(template, template))
        template, maps = graphs[key]
        built[i] = (_with_coords(template, prepared[i][0]), maps)
        return built[i]

    # the radial lower bound of frame i to every frame, in one numpy step (frames with another element key get None)
    keys = [p[0] if p else None for p in prints]
    stacks = {}
    for k in set(x for x in keys if x is not None):
        ids = [j for j, x in enumerate(keys) if x == k]
        stacks[k] = (ids, np.vstack([prints[j][1] for j in ids]))

    def lb_row(i):
        out = {}
        if keys[i] is None:
            return out
        ids, mat = stacks[keys[i]]
        vals = np.sqrt(np.mean((mat - prints[i][1]) ** 2, axis=1))
        for j, v in zip(ids, vals):
            out[j] = float(v)
        return out

    clusters, pairs, suggestions = [], [], []
    for i in order:
        if should_stop is not None and should_stop():     # a newer request from the page: stop, the caller drops this
            break
        if prints[i] is None:
            clusters.append([i])
            continue
        placed = False
        offer = None
        lbs = lb_row(i)
        for cl in clusters:
            rep = cl[0]
            bound = lbs.get(rep)
            if bound is None:
                continue
            if bound > 2.5 * rmsd_cut:
                continue
            members = [m for m in cl if lbs.get(m) is not None and lbs[m] <= rmsd_cut]
            probe, maps = mol(i)
            for member in members:
                other, _other_maps = mol(member)
                if probe is None or other is None:
                    continue
                judged = judge_pair(probe, other, rmsd_cut, tight, torsion_cut, maps=maps)
                if judged["kind"] == "tight":
                    # parent = the member judged against, so its rmsd and swaps describe that pair (the page follows
                    # the chain to the cluster's lowest)
                    pairs.append(_record(i, member, judged, [p[1] for p in prepared]))
                    cl.append(i)
                    placed = True
                    break
                if judged["kind"] == "torsion" and (offer is None or judged["rmsd"] < offer["rmsd"]):
                    offer = _record(i, member, judged, [p[1] for p in prepared])
            if placed:
                break
        if not placed:
            clusters.append([i])
            if offer is not None:
                suggestions.append(offer)
    return {"pairs": pairs, "suggestions": suggestions, "kept": len(clusters), "folded": len(pairs)}


def _source_tag(name):
    text = str(name or "").upper()
    if "CREST" in text:
        return "CREST"
    if "GOAT" in text:
        return "GOAT"
    return str(name or "file")


def cross_tags(files, rmsd_cut=DUP_A, tight=TIGHT_A, torsion_cut=TORSION_DEG):
    """Two ensembles of one molecule, on one energy zero.

    A frame is BOTH when the other file has the same heavy-atom shape (the
    0.15 Angstrom rule, torsions included). Otherwise it keeps its own file's
    tag. joint_kcal is the gap to the lower of the two minima.
    """
    if len(files) < 2:
        return {"tags": [], "minimum": None}
    left, right = files[0], files[1]
    energies = [
        float(fr["energy"]) for group in (left, right) for fr in group.get("frames") or []
        if fr.get("energy") is not None
    ]
    emin = min(energies) if energies else None
    prepared, prints = {}, {}
    templates, automaps, cross_maps = {}, {}, {}

    def prep(fr):
        key = id(fr)
        if key not in prepared:
            prepared[key] = heavy_rows_idx(fr.get("atoms") or [])
            prints[key] = _fingerprint(prepared[key][0])
        return key

    def sequence(fr):
        key = prep(fr)
        rows = prepared[key][0]
        seq = tuple(_symbol(row[0]) for row in rows)
        if seq not in templates and rows:
            templates[seq] = mol_from_heavy(rows)
            automaps[seq] = _maps(templates[seq], templates[seq])
        return seq

    def mol(fr):
        seq = sequence(fr)
        if seq not in templates:
            return None
        return _with_coords(templates[seq], prepared[prep(fr)][0])

    def maps_between(seq_a, seq_b):
        if seq_a == seq_b:
            return automaps.get(seq_a)
        key = (seq_a, seq_b)
        if key not in cross_maps and seq_a in templates and seq_b in templates:
            cross_maps[key] = _maps(templates[seq_a], templates[seq_b])
        return cross_maps.get(key)

    # every pair that is the same conformer, then one-to-one by an optimal assignment: as many matched as possible,
    # then the smallest total RMSD. (Greedy in file order let an early frame take a partner a later one needed.)
    from scipy.optimize import linear_sum_assignment

    lf, rf = list(left.get("frames") or []), list(right.get("frames") or [])
    found = {}
    for i, fr in enumerate(lf):
        for j, other in enumerate(rf):
            bound = _radial_lb(prints[prep(fr)], prints[prep(other)])
            if bound is None or bound > rmsd_cut:
                continue
            probe, other_mol = mol(fr), mol(other)
            if probe is None or other_mol is None:
                continue
            judged = judge_pair(
                probe, other_mol, rmsd_cut, tight, torsion_cut,
                maps=maps_between(sequence(fr), sequence(other)),
            )
            if judged["kind"] in ("tight", "torsion"):
                found[(i, j)] = judged
    matched_left, matched_right = {}, set()
    if found:
        miss = 1e6                                     # no pair: costs more than any set of real matches
        cost = np.full((len(lf), len(rf)), miss)
        for (i, j), judged in found.items():
            cost[i, j] = float(judged["rmsd"] or 0.0)
        rows, cols = linear_sum_assignment(cost)
        for i, j in zip(rows, cols):
            if (i, j) in found:
                matched_left[i] = j
                matched_right.add(j)
    tags = []
    for i, fr in enumerate(lf):
        if i not in matched_left:
            tags.append({"id": fr.get("id"), "tag": _source_tag(left.get("name")), "rmsd": None, "with": None, "swaps": []})
            continue
        j = matched_left[i]
        other, judged = rf[j], found[(i, j)]
        swaps = _swaps(judged["amap"], prepared[prep(fr)][1])
        tags.append({"id": fr.get("id"), "tag": "BOTH", "rmsd": judged["rmsd"], "with": other.get("id"), "swaps": swaps})
        tags.append({"id": other.get("id"), "tag": "BOTH", "rmsd": judged["rmsd"], "with": fr.get("id"), "swaps": swaps})
    for j, other in enumerate(rf):
        if j in matched_right:
            continue
        tags.append({"id": other.get("id"), "tag": _source_tag(right.get("name")), "rmsd": None, "with": None, "swaps": []})

    def rel(fr):
        if emin is None or fr.get("energy") is None:
            return None
        return round((float(fr["energy"]) - emin) * 627.509474, 4)

    by_id = {fr.get("id"): fr for group in (left, right) for fr in group.get("frames") or []}
    for row in tags:
        row["joint_kcal"] = rel(by_id.get(row["id"]) or {})
        if row["rmsd"] is not None:
            row["rmsd"] = round(float(row["rmsd"]), 4)
    return {"tags": tags, "minimum": emin}
