"""SMILES to starting geometries with a fixed atom numbering.

Three build types:

- ``organic``: the molecule alone. With ``core_smarts`` the atoms matching it come first, in
  SMARTS order, so the same descriptor atoms carry the same numbers in every molecule.
- ``metal_mono``: one metal on one donor (``donor_smarts``, default ``[PX3]``), on the donor's
  lone-pair direction. Order: metal, donor, then the rest.
- ``metal_chelate``: one metal on a bidentate ligand, plus ``n_ancillary`` ancillary ligands
  (H or Cl). Order: metal, ancillaries, donor 1, donor 2, then the rest: the order
  ``MetalComplex`` expects. The placement is CS3's ``build_general.py`` (validated against
  DFT-optimised Ni complexes); donors come from ``donor_smarts`` (atom maps 1 and 2) or from its
  graph heuristic.

Each molecule: ``n_conformers`` ETKDG embeds (``seed``), MMFF (UFF where MMFF has no parameters),
tried in force-field energy order; the first whose metal clears every atom is used. A chelate's
ligand is then relaxed around the fixed metal (see ``_relax_around_metal``), as is a monodentate
ligand when no conformer leaves room otherwise.
Writes ``<workdir>/build/<id>.xyz``, ``<workdir>/elements/<id>.elements``, ``ids.txt`` and
``manifest.json``.
"""
from __future__ import annotations

import csv
import json
import re
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

from .protocol import Protocol

# M-donor and M-ancillary distances, as in build_general.py
MD_BY_METAL = {
    "Ni": {"N": 2.09, "P": 2.20, "S": 2.25, "O": 2.05},
    "Cu": {"N": 2.00, "P": 2.20, "S": 2.30, "O": 2.00},
    "Pd": {"N": 2.10, "P": 2.30, "S": 2.30, "O": 2.10},
}
M_ANC = {("Ni", "H"): 1.48, ("Cu", "Cl"): 2.20, ("Ni", "Cl"): 2.24, ("Pd", "Cl"): 2.30}
CLEAR_HEAVY, CLEAR_H = 2.2, 1.8        # A, metal to any ligand atom other than the donors
ANC_CLEAR = 1.45                         # A, ancillary to ligand, as build_general


class BuildError(RuntimeError):
    pass


def read_molecules(p: Protocol) -> list[dict]:
    with open(p.molecules, newline="", encoding="utf-8-sig") as fh:
        rows = list(csv.DictReader(fh))
    if not rows or p.smiles_column not in rows[0]:
        raise BuildError(f"{p.molecules.name}: no column {p.smiles_column!r}")
    out = []
    for k, r in enumerate(rows, 1):
        mid = r[p.id_column].strip() if p.id_column else f"m{k:03d}"
        if not re.fullmatch(r"[A-Za-z0-9_-]+", mid):
            raise BuildError(f"id {mid!r} has characters other than letters, digits, _ and -")
        out.append(dict(id=mid, name=r.get(p.name_column, mid) if p.name_column else mid,
                        smiles=r[p.smiles_column].strip()))
    ids = [m["id"] for m in out]
    if len(set(ids)) != len(ids):
        raise BuildError("molecule ids are not unique")
    return out


def _embed(smiles: str, n_conf: int, seed: int, pair=None, pair_range=(2.5, 3.0)):
    """ETKDG conformers relaxed by MMFF (UFF where MMFF lacks parameters), in force-field energy
    order. With ``pair`` the two donor atoms are held at a chelating distance while embedding and
    relaxing: a free ligand like bipyridine otherwise comes out anti, lone pairs apart."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise BuildError(f"RDKit cannot read {smiles!r}")
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    if pair:
        from rdkit import DistanceGeometry
        from rdkit.Chem import rdDistGeom
        i, j = sorted(pair)
        bm = rdDistGeom.GetMoleculeBoundsMatrix(mol)
        bm[i][j], bm[j][i] = pair_range[1], pair_range[0]          # upper triangle = max, lower = min
        DistanceGeometry.DoTriangleSmoothing(bm)
        params.SetBoundsMat(bm)
        params.useExpTorsionAnglePrefs = False     # the torsion library holds bipyridine anti over the bound
    cids = list(AllChem.EmbedMultipleConfs(mol, n_conf, params))
    if not cids:
        raise BuildError("embedding failed" + (" with the donors at a chelating distance" if pair else ""))
    mmff = AllChem.MMFFHasAllMoleculeParams(mol)
    props = AllChem.MMFFGetMoleculeProperties(mol) if mmff else None
    energy = {}
    for cid in cids:
        ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=cid) if mmff else AllChem.UFFGetMoleculeForceField(mol, confId=cid)
        if pair and mmff:
            ff.MMFFAddDistanceConstraint(pair[0], pair[1], False, pair_range[0], pair_range[1] - 0.1, 100.0)
        elif pair:
            ff.UFFAddDistanceConstraint(pair[0], pair[1], False, pair_range[0], pair_range[1] - 0.1, 100.0)
        ff.Minimize(maxIts=2000)
        energy[cid] = ff.CalcEnergy()
    return mol, sorted(cids, key=energy.get)


def _lone_pair(X, i, neighbours):
    v = -sum((X[j] - X[i]) / np.linalg.norm(X[j] - X[i]) for j in neighbours)
    n = np.linalg.norm(v)
    return v / n if n > 1e-6 else None


def _clear(sym, X, point, skip, heavy=CLEAR_HEAVY, h=CLEAR_H):
    d = np.linalg.norm(X - point, axis=1)
    return all(d[i] >= (h if sym[i] == "H" else heavy) for i in range(len(sym)) if i not in skip), \
        float(min(d[i] for i in range(len(sym)) if i not in skip))


def _relax_around_metal(mol, cid, m, donors, dists):
    """Relax the ligand (UFF) around a metal fixed at ``m``, donors held at ``dists``.

    The force field has no metal types, so a silicon atom stands in at ``m``: unbonded, fixed, and
    bonded to the donors (so they feel no repulsion from it) and only there for its van der Waals
    repulsion on everything else, which turns a methyl or tBu hydrogen out of the
    metal's pocket. The stand-in is not written anywhere. Returns the relaxed ligand coordinates."""
    from rdkit.Geometry import Point3D

    from rdkit import RDLogger

    rw = Chem.RWMol(Chem.Mol(mol, confId=cid))
    k = rw.AddAtom(Chem.Atom("Si"))
    for d in donors:                       # bonded, so the donors feel no repulsion from the stand-in
        rw.AddBond(int(d), k, Chem.BondType.SINGLE)
    cx = rw.GetMol()
    cx.UpdatePropertyCache(strict=False)
    atom = cx.GetAtomWithIdx(k)
    atom.SetHybridization(Chem.HybridizationType.SP3)
    atom.SetNoImplicit(True)
    cx.GetConformer().SetAtomPosition(k, Point3D(*m))
    RDLogger.DisableLog("rdApp.*")        # the typer complains about the stand-in's charge state
    try:
        ff = AllChem.UFFGetMoleculeForceField(cx, ignoreInterfragInteractions=False)
    finally:
        RDLogger.EnableLog("rdApp.*")
    ff.AddFixedPoint(k)
    for d, r in zip(donors, dists):
        ff.UFFAddDistanceConstraint(int(d), k, False, r - 0.02, r + 0.02, 2000.0)
    ff.Minimize(maxIts=5000)
    return cx.GetConformer().GetPositions()[:-1]


def _match(mol, smarts, what):
    q = Chem.MolFromSmarts(smarts)
    if q is None:
        raise BuildError(f"cannot parse {what} {smarts!r}")
    hits = mol.GetSubstructMatches(q)
    if not hits:
        raise BuildError(f"{what} {smarts!r} does not match")
    return q, hits


# ---------------------------------------------------------------- chelate placement (build_general)
def _place_chelate_metal(X, i, j, li, lj, ri, rj):
    from scipy.optimize import minimize

    def cost(m):
        di, dj = m - X[i], m - X[j]
        ni, nj = np.linalg.norm(di), np.linalg.norm(dj)
        return 3.0 * ((ni - ri) ** 2 + (nj - rj) ** 2) + 0.5 * ((1 - di @ li / ni) + (1 - dj @ lj / nj))

    start = (X[i] + X[j]) / 2 + (li + lj) / np.linalg.norm(li + lj) * 1.6
    return minimize(cost, start, method="Nelder-Mead", options=dict(xatol=1e-4, fatol=1e-8, maxiter=4000)).x


def _ancillaries(m, d1, d2, ligand, dist, n, target):
    """n = 2: each trans to a donor (square planar); n = 1: on the external bisector (trigonal).
    Tilted off the ideal direction, up to 30 deg, when a ligand atom is in the way."""
    u1, u2 = (d1 - m) / np.linalg.norm(d1 - m), (d2 - m) / np.linalg.norm(d2 - m)
    nrm = np.cross(u1, u2)
    nrm = nrm / np.linalg.norm(nrm) if np.linalg.norm(nrm) > 1e-6 else np.array([0.0, 0.0, 1.0])
    if n == 0:
        return []
    ideal = [-u1, -u2] if n == 2 else [-(u1 + u2) / np.linalg.norm(u1 + u2)]
    out = []
    for u in ideal:
        best, bestd = u, -1.0
        perp = np.cross(nrm, u)
        perp /= np.linalg.norm(perp)
        for tilt in np.radians([0, 10, 20, 30]):
            for phi in np.linspace(0, 2 * np.pi, 12, endpoint=False):
                axis = np.cos(phi) * nrm + np.sin(phi) * perp
                v = u * np.cos(tilt) + np.cross(axis, u) * np.sin(tilt)
                v /= np.linalg.norm(v)
                score = min(np.linalg.norm(ligand - (m + v * dist), axis=1).min(), target) - 0.15 * tilt
                if score > bestd:
                    best, bestd = v, score
            if bestd > 0.9 * target:
                break
        out.append(m + best * dist)
    return out


def _chelate_donor_pairs(mol, b):
    if b.get("donor_smarts"):
        q, hits = _match(mol, b["donor_smarts"], "donor_smarts")
        maps = {a.GetAtomMapNum(): a.GetIdx() for a in q.GetAtoms() if a.GetAtomMapNum()}
        if 1 not in maps or 2 not in maps:
            raise BuildError("donor_smarts must mark the two donors with atom maps :1 and :2")
        return [(h[maps[1]], h[maps[2]]) for h in hits]
    sym = [a.GetSymbol() for a in mol.GetAtoms()]
    cand = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() in ("N", "P", "S") and a.GetDegree() <= 3]
    pairs = []
    for x in range(len(cand)):
        for y in range(x + 1, len(cand)):
            i, j = cand[x], cand[y]
            if len(Chem.GetShortestPath(mol, i, j)) - 1 in (3, 4):
                heavy = sum(n.GetSymbol() != "H" for n in mol.GetAtomWithIdx(i).GetNeighbors()) + \
                    sum(n.GetSymbol() != "H" for n in mol.GetAtomWithIdx(j).GetNeighbors())
                pairs.append((heavy, i, j))
    return [(i, j) for _, i, j in sorted(pairs)] or []


# ---------------------------------------------------------------- one molecule
def build_one(smiles: str, b: dict) -> dict:
    kind = b["type"]
    n_conf, seed = int(b.get("n_conformers", 10)), int(b.get("seed", 0))
    if kind == "metal_chelate":
        return _build_chelate(smiles, b, n_conf, seed)
    mol, confs = _embed(smiles, n_conf, seed)
    sym = [a.GetSymbol() for a in mol.GetAtoms()]
    n = len(sym)
    key = {}

    if kind == "organic":
        order = list(range(n))
        if b.get("core_smarts"):
            _, hits = _match(mol, b["core_smarts"], "core_smarts")
            core = list(hits[0])
            order = core + [i for i in range(n) if i not in core]
            key["core_1based"] = list(range(1, len(core) + 1))
        X = mol.GetConformer(confs[0]).GetPositions()
        return dict(symbols=[sym[i] for i in order], coords=X[order], conformer=confs[0], key=key, mol=mol)

    metal = b["metal"]
    md = dict(MD_BY_METAL.get(metal, {}), **b.get("metal_distance", {})) if isinstance(b.get("metal_distance", {}), dict) \
        else {e: float(b["metal_distance"]) for e in ("N", "P", "S", "O")}

    if kind == "metal_mono":
        q, hits = _match(mol, b.get("donor_smarts", "[PX3]"), "donor_smarts")
        donors = {h[0] for h in hits}
        if len(donors) != 1:
            raise BuildError(f"donor_smarts matches {len(donors)} atoms; it must pick exactly one donor")
        d = donors.pop()
        nb = [a.GetIdx() for a in mol.GetAtomWithIdx(d).GetNeighbors()]
        for cid in confs:
            X = mol.GetConformer(cid).GetPositions()
            lp = _lone_pair(X, d, nb)
            if lp is None:
                continue
            m = X[d] + md.get(sym[d], 2.2) * lp
            ok, closest = _clear(sym, X, m, {d})
            if ok:
                order = [d] + [i for i in range(n) if i != d]
                key.update(metal_1based=1, donor_1based=[2], closest_to_metal=round(closest, 3))
                return dict(symbols=[metal] + [sym[i] for i in order], coords=np.vstack([m, X[order]]),
                            conformer=cid, key=key, mol=mol)
        for cid in confs:                        # no conformer leaves room: relax the ligand around the metal
            X = mol.GetConformer(cid).GetPositions()
            lp = _lone_pair(X, d, nb)
            if lp is None:
                continue
            m = X[d] + md.get(sym[d], 2.2) * lp
            X = _relax_around_metal(mol, cid, m, [d], [md.get(sym[d], 2.2)])
            ok, closest = _clear(sym, X, m, {d})
            if ok:
                order = [d] + [i for i in range(n) if i != d]
                key.update(metal_1based=1, donor_1based=[2], closest_to_metal=round(closest, 3), relaxed_around_metal=True)
                return dict(symbols=[metal] + [sym[i] for i in order], coords=np.vstack([m, X[order]]),
                            conformer=cid, key=key, mol=mol)
        raise BuildError("no conformer leaves room for the metal, even relaxed around it")



def _metal_distances(b):
    md = b.get("metal_distance", {})
    if not isinstance(md, dict):
        return {e: float(md) for e in ("N", "P", "S", "O")}
    return dict(MD_BY_METAL.get(b["metal"], {}), **md)


def _build_chelate(smiles, b, n_conf, seed):
    metal, md = b["metal"], _metal_distances(b)
    anc, n_anc = b.get("ancillary", "H"), int(b.get("n_ancillary", 2))
    dist = M_ANC.get((metal, anc), 2.2)
    target = 2.0 if anc == "H" else 3.0
    graph = Chem.AddHs(Chem.MolFromSmiles(smiles))
    pairs = _chelate_donor_pairs(graph, b)
    if not pairs:
        raise BuildError("no bidentate donor pair found")
    tried = []
    for i, j in pairs:
        mol, confs = _embed(smiles, n_conf, seed, pair=(i, j))
        sym = [a.GetSymbol() for a in mol.GetAtoms()]
        n = len(sym)
        key = {}
        for cid in confs:
            X = mol.GetConformer(cid).GetPositions()
            li = _lone_pair(X, i, [a.GetIdx() for a in mol.GetAtomWithIdx(i).GetNeighbors()])
            lj = _lone_pair(X, j, [a.GetIdx() for a in mol.GetAtomWithIdx(j).GetNeighbors()])
            if li is None or lj is None or np.linalg.norm(li + lj) < 1e-6:
                continue
            ri, rj = md.get(sym[i], 2.1), md.get(sym[j], 2.1)
            m = _place_chelate_metal(X, i, j, li, lj, ri, rj)
            X = _relax_around_metal(mol, cid, m, [i, j], [ri, rj])
            extra = _ancillaries(m, X[i], X[j], X, dist, n_anc, target)
            ui, uj = (X[i] - m) / np.linalg.norm(X[i] - m), (X[j] - m) / np.linalg.norm(X[j] - m)
            bite = float(np.degrees(np.arccos(np.clip(ui @ uj, -1, 1))))
            ok, closest = _clear(sym, X, m, {i, j})
            anc_ok = all(np.linalg.norm(X - e, axis=1).min() >= ANC_CLEAR for e in extra)
            if ok and anc_ok and 65 < bite < 105:
                order = [i, j] + [k for k in range(n) if k not in (i, j)]
                na = len(extra)
                key.update(metal_1based=1, ancillary_1based=list(range(2, 2 + na)),
                           donor_1based=[2 + na, 3 + na], bite=round(bite, 1), closest_to_metal=round(closest, 3))
                return dict(symbols=[metal] + [anc] * na + [sym[k] for k in order],
                            coords=np.vstack([m] + extra + [X[order]]), conformer=cid, key=key, mol=mol)
            tried.append((cid, sym[i], sym[j], round(bite, 1), round(closest, 2)))
    raise BuildError(f"every conformer / donor pair rejected, e.g. {tried[:3]}")


def write_xyz(path: Path, symbols, coords, comment: str):
    lines = [str(len(symbols)), comment] + [f"{s:<2} {x:16.8f} {y:16.8f} {z:16.8f}" for s, (x, y, z) in zip(symbols, coords)]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")


def build_all(p: Protocol, strict: bool = True) -> list[dict]:
    """Build every molecule into p.workdir. A failed molecule stops the build when ``strict``;
    otherwise it is reported in the manifest and left out of ids.txt."""
    wd = p.workdir
    (wd / "build").mkdir(parents=True, exist_ok=True)
    (wd / "elements").mkdir(exist_ok=True)
    rows = []
    for m in read_molecules(p):
        try:
            r = build_one(m["smiles"], p.build)
        except BuildError as exc:
            if strict:
                raise BuildError(f"{m['id']} ({m['name']}): {exc}") from None
            rows.append(dict(m, status=f"failed: {exc}"))
            continue
        write_xyz(wd / "build" / f"{m['id']}.xyz", r["symbols"], r["coords"], f"{m['id']} {m['name']}")
        (wd / "elements" / f"{m['id']}.elements").write_text(" ".join(r["symbols"]) + "\n", encoding="utf-8", newline="\n")
        formula = rdMolDescriptors.CalcMolFormula(r["mol"]) + ("" if p.build["type"] == "organic" else f" + {p.build['metal']}")
        rows.append(dict(m, status="built", n_atoms=len(r["symbols"]), formula=formula, conformer=r["conformer"], **r["key"]))
    built = [r["id"] for r in rows if r["status"] == "built"]
    (wd / "ids.txt").write_text("\n".join(built) + "\n", encoding="utf-8", newline="\n")
    (wd / "manifest.json").write_text(json.dumps(dict(protocol=p.name, build=p.build, molecules=rows), indent=2) + "\n",
                                      encoding="utf-8", newline="\n")
    return rows
