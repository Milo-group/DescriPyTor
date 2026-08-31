"""Bond-graph helpers used to replace python-igraph."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
M2 = ROOT / "M2_data_extractor"
for path in (str(ROOT), str(M2)):
    if path not in sys.path:
        sys.path.insert(0, path)

from extractor_utils.bond_graph import (
    all_simple_paths_from,
    diameter_vertex_path,
    graph_from_bonds_df,
    unique_atoms_on_paths,
)


def _chain_bonds():
    # 1-2-3-4 plus side branch 2-5
    return pd.DataFrame({0: [1, 2, 3, 2], 1: [2, 3, 4, 5]})


def _connections(bonds_df, source, direction, mode="all"):
    source, direction = int(source), int(direction)
    graph = graph_from_bonds_df(bonds_df, directed=False)
    if not graph.has_edge(source, direction):
        graph.add_edge(source, direction)
    paths = all_simple_paths_from(graph, source)
    started = [
        path for path in paths
        if len(path) >= 2 and path[0] == source and path[1] == direction
    ]
    if mode == "shortest":
        return min(started, key=len) if started else []
    return unique_atoms_on_paths(started)


def test_graph_from_bonds_df_is_undirected():
    G = graph_from_bonds_df(_chain_bonds())
    assert G.has_edge(1, 2) and G.has_edge(2, 1)
    assert set(G.nodes()) == {1, 2, 3, 4, 5}


def test_connections_follow_attached_bond():
    atoms = set(int(x) for x in _connections(_chain_bonds(), 1, 2, mode="all"))
    assert atoms == {1, 2, 3, 4, 5}


def test_connections_exclude_other_branch_when_direction_is_side():
    atoms = set(int(x) for x in _connections(_chain_bonds(), 2, 5, mode="all"))
    assert atoms == {2, 5}


def test_shortest_mode_starts_with_forced_edge():
    path = _connections(_chain_bonds(), 1, 2, mode="shortest")
    assert list(path)[:2] == [1, 2]


def test_diameter_path_is_longest_shortest_path():
    G = graph_from_bonds_df(_chain_bonds())
    path = diameter_vertex_path(G)
    assert len(path) == 4
    assert {path[0], path[-1]} <= {1, 4, 5}


def test_importing_data_extractor_does_not_need_igraph():
    import pytest

    pytest.importorskip("morfeus")
    sys.modules.pop("igraph", None)
    sys.modules.pop("data_extractor", None)
    sys.modules.pop("M2_data_extractor.data_extractor", None)
    import data_extractor  # noqa: F401

    assert "igraph" not in sys.modules
