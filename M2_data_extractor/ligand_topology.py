"""Graph-index descriptors of a free ligand from SMILES.

These are the Case Study 3 ``t_*`` / ``tf_*`` family: connectivity (χ), Kier
kappa, Zagreb, Wiener, Balaban J, Estrada, and related counts. They need no
complex, no conformer ensemble, and no QM.

``tf_*`` columns are the raw ``t_*`` values divided by the heavy-atom count.
That size-normalized block is what the reported models searched over.

RDKit is used only to parse SMILES into a graph. The indices themselves are
computed from the adjacency list so the formulae stay explicit (and so a
NumPy 2 / RDKit GraphDescriptors mismatch cannot silently change numbers).
"""

from __future__ import annotations

from collections import deque

import numpy as np
import pandas as pd

try:
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
except ImportError as exc:  # pragma: no cover
    Chem = None
    _RDKIT_IMPORT_ERROR = exc
else:
    _RDKIT_IMPORT_ERROR = None

# Valence electrons for Kier–Hall delta_v.
ZV = {"C": 4, "N": 5, "O": 6, "P": 5, "S": 6, "F": 7, "Cl": 7, "Br": 7, "I": 7}
Z = {"C": 6, "N": 7, "O": 8, "P": 15, "S": 16, "F": 9, "Cl": 17, "Br": 35, "I": 53}

TF_FROM_T = (
    "chi0", "chi0v", "chi1", "chi1v", "chi2", "chi2v", "chi3", "chi3v",
    "chi4", "chi4v", "kappa1", "kappa2", "kappa3", "zagreb1", "zagreb2",
    "wiener", "global_eff", "balaban", "estrada",
)


def _require_rdkit():
    if Chem is None:
        raise ImportError(
            "LigandTopology requires RDKit to parse SMILES"
        ) from _RDKIT_IMPORT_ERROR


def _graph(smiles: str):
    _require_rdkit()
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.RemoveHs(mol)
    n = mol.GetNumAtoms()
    adj = [[] for _ in range(n)]
    bonds = []
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        adj[i].append(j)
        adj[j].append(i)
        bonds.append((i, j))
    delta = [len(neighbours) for neighbours in adj]
    dv = []
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        zv, z = ZV.get(symbol, 4), Z.get(symbol, 6)
        h = atom.GetTotalNumHs()
        den = z - zv - 1
        dv.append((zv - h) / den if den > 0 else float(zv - h))
    return n, adj, bonds, delta, dv, mol


def _paths_of_length(adj, k):
    out = []

    def walk(path):
        if len(path) - 1 == k:
            if path[0] < path[-1]:
                out.append(tuple(path))
            return
        for neighbour in adj[path[-1]]:
            if neighbour not in path:
                walk(path + [neighbour])

    for start in range(len(adj)):
        walk([start])
    return out


def _chi(adj, weights, k):
    if k == 0:
        return sum(1.0 / np.sqrt(x) for x in weights if x > 0)
    total = 0.0
    for path in _paths_of_length(adj, k):
        product = 1.0
        ok = True
        for idx in path:
            if weights[idx] <= 0:
                ok = False
                break
            product *= weights[idx]
        if ok:
            total += 1.0 / np.sqrt(product)
    return total


def _dist_matrix(adj):
    n = len(adj)
    dist = np.full((n, n), np.inf)
    for start in range(n):
        dist[start, start] = 0
        queue = deque([start])
        while queue:
            u = queue.popleft()
            for v in adj[u]:
                if dist[start, v] == np.inf:
                    dist[start, v] = dist[start, u] + 1
                    queue.append(v)
    return dist


class LigandTopology:
    """Graph indices of one free-ligand SMILES string."""

    def __init__(self, smiles: str, name: str | None = None):
        self.smiles = smiles
        self.name = name
        self._graph = _graph(smiles)

    @classmethod
    def from_smiles(cls, smiles: str, name: str | None = None) -> "LigandTopology":
        return cls(smiles, name=name)

    @property
    def is_valid(self) -> bool:
        return self._graph is not None

    def raw_features(self) -> dict:
        """Un-normalized ``t_*`` indices, including atom/ring/rotatable counts."""
        if self._graph is None:
            return {}
        n, adj, bonds, delta, dv, mol = self._graph
        out = {}
        for k in range(5):
            out[f"t_chi{k}"] = float(_chi(adj, delta, k))
            out[f"t_chi{k}v"] = float(_chi(adj, dv, k))

        n_atoms = float(n)
        path_counts = {k: len(_paths_of_length(adj, k)) for k in (1, 2, 3)}
        out["t_kappa1"] = float(n_atoms * (n_atoms - 1) ** 2 / path_counts[1] ** 2) if path_counts[1] else 0.0
        out["t_kappa2"] = float((n_atoms - 1) * (n_atoms - 2) ** 2 / path_counts[2] ** 2) if path_counts[2] else 0.0
        if path_counts[3]:
            if n % 2:
                out["t_kappa3"] = float((n_atoms - 1) * (n_atoms - 3) ** 2 / path_counts[3] ** 2)
            else:
                out["t_kappa3"] = float((n_atoms - 3) * (n_atoms - 2) ** 2 / path_counts[3] ** 2)
        else:
            out["t_kappa3"] = 0.0

        out["t_zagreb1"] = float(sum(x * x for x in delta))
        out["t_zagreb2"] = float(sum(delta[i] * delta[j] for i, j in bonds))

        dist = _dist_matrix(adj)
        iu = np.triu_indices(n, 1)
        finite = dist[iu][np.isfinite(dist[iu])]
        out["t_wiener"] = float(finite.sum())
        positive = finite[finite > 0]
        out["t_global_eff"] = float((1.0 / positive).mean()) if len(positive) else 0.0

        dist_sums = np.where(np.isfinite(dist), dist, 0).sum(axis=1)
        mu = len(bonds) - n + 1
        jb = sum(
            1.0 / np.sqrt(dist_sums[i] * dist_sums[j])
            for i, j in bonds
            if dist_sums[i] > 0 and dist_sums[j] > 0
        )
        out["t_balaban"] = float(len(bonds) / (mu + 1) * jb) if mu + 1 else 0.0

        adj_mat = np.zeros((n, n))
        for i, j in bonds:
            adj_mat[i, j] = adj_mat[j, i] = 1.0
        out["t_estrada"] = float(np.exp(np.linalg.eigvalsh(adj_mat)).sum())

        out["t_heavy"] = float(n)
        out["t_nrings"] = float(mu)
        out["t_nrot"] = float(sum(
            1 for bond in mol.GetBonds()
            if bond.GetBondType() == Chem.BondType.SINGLE and not bond.IsInRing()
            and bond.GetBeginAtom().GetDegree() > 1
            and bond.GetEndAtom().GetDegree() > 1
        ))
        return out

    def size_normalized_features(self) -> dict:
        """``tf_*`` = ``t_*`` / heavy-atom count (the CS3 search-pool block)."""
        raw = self.raw_features()
        n = raw.get("t_heavy") or 0.0
        if n <= 0:
            return {}
        return {f"tf_{key}": raw[f"t_{key}"] / n for key in TF_FROM_T}

    def features(self, include_raw: bool = True, include_normalized: bool = True) -> dict:
        out = {}
        if include_raw:
            out.update(self.raw_features())
        if include_normalized:
            out.update(self.size_normalized_features())
        return out

    @classmethod
    def table_from_smiles(cls, smiles_by_name, normalized: bool = True) -> pd.DataFrame:
        """Build a ligand × descriptor table from ``{name: smiles}``."""
        rows = {}
        for name, smiles in smiles_by_name.items():
            topo = cls(smiles, name=name)
            if not topo.is_valid:
                continue
            rows[name] = (
                topo.size_normalized_features() if normalized
                else topo.features()
            )
        return pd.DataFrame.from_dict(rows, orient="index")
