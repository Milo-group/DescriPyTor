"""Regression tests for the feature-definition fixes of 2026-09-27 (docs/FEATURE_FIXES.md)."""
import contextlib
import io
from pathlib import Path

import networkx as nx
import numpy as np
import pytest

from M2_data_extractor.data_extractor import Molecule, Molecules
from M2_data_extractor.extractor_utils.vibrations_utils import get_benzene_ring_indices

EX = Path(__file__).resolve().parents[1] / "Getting_started_with_examples" / "feather_example"


def _quiet(fn, *a, **k):
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*a, **k)


@pytest.fixture(scope="module")
def mol():
    return _quiet(Molecule, str(EX / "m_Br.feather"))


def _graph(m):
    return nx.Graph([(int(a), int(b)) for a, b in m.bonds_df.values])


def test_carbon_bromine_bond_is_found(mol):
    br = [i + 1 for i, e in enumerate(mol.xyz_df.atom) if e == "Br"][0]
    assert any(br in pair for pair in map(tuple, mol.bonds_df.values.tolist()))


def test_ring_positions_are_measured_from_the_given_atom(mol):
    G = _graph(mol)
    ring = [c for c in nx.cycle_basis(G) if len(c) == 6][0]
    for atom in ring:
        z, x, c, v, b, n = get_benzene_ring_indices(mol.bonds_df, [atom])
        d = nx.single_source_shortest_path_length(G.subgraph(ring), atom)
        assert (d[x], d[z]) == (0, 3)                 # primary, para
        assert (d[c], d[v]) == (1, 1)                 # ortho pair
        assert (d[b], d[n]) == (2, 2)                 # meta pair


def test_ring_angles_are_not_a_constant(mol):
    ring = [c for c in nx.cycle_basis(_graph(mol)) if len(c) == 6][0]
    out = _quiet(mol.get_ring_vibrations, [ring[0]], verbose=False)
    assert out is not None
    assert not np.allclose(out[["cross_angle", "para_angle"]].to_numpy(), 90.0)


def _a_bond(m, heavy=True):
    el = list(m.xyz_df.atom)
    for a, b in m.bonds_df.values.tolist():
        if not heavy or (el[a - 1] != "H" and el[b - 1] != "H"):
            return [int(a), int(b)]


def test_tuples_give_the_same_answer_as_lists(mol):
    p = _a_bond(mol)
    q = [p[1], p[0]]
    assert np.allclose(mol.get_bond_length([p, q]).to_numpy(), mol.get_bond_length([tuple(p), tuple(q)]).to_numpy())
    trip = [p[0], p[1], [n for n in _graph(mol)[p[1]] if n != p[0]][0]]
    assert np.allclose(mol.get_bond_angle([trip]).to_numpy(), mol.get_bond_angle([tuple(trip)]).to_numpy())


def test_stretch_without_a_window_uses_the_default():
    ms = _quiet(Molecules, str(EX))
    ref = ms.molecules[0]
    pair = _a_bond(ref)
    got = _quiet(ms.get_molecules_features_set, {"Stretching": "%d,%d" % tuple(pair)})
    assert got is not None and any("Stretch" in str(c) for c in got.columns)


def test_a_bend_triplet_means_its_two_ends(mol):
    G = _graph(mol)
    centre = next(n for n in G if G.degree[n] >= 2)
    a, c = sorted(G[centre])[:2]
    by_pair = _quiet(mol.get_bend_vibration_single, [a, c])
    by_triplet = _quiet(mol.get_bend_vibration_single, [a, centre, c])
    assert np.allclose(by_pair.to_numpy(float), by_triplet.to_numpy(float))
