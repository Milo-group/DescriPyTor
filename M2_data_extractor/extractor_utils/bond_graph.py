"""Bond graphs via networkx (replaces python-igraph)."""
from __future__ import annotations

import numpy as np
import pandas as pd
import networkx as nx


def graph_from_bonds_df(bonds_df, directed=False):
    """Undirected (default) or directed graph from a two-column bond table."""
    G = nx.DiGraph() if directed else nx.Graph()
    if bonds_df is None or len(bonds_df) == 0:
        return G
    u = pd.to_numeric(bonds_df.iloc[:, 0], errors="raise").astype(int)
    v = pd.to_numeric(bonds_df.iloc[:, 1], errors="raise").astype(int)
    G.add_edges_from(zip(u.tolist(), v.tolist()))
    return G


def all_simple_paths_from(G, source):
    """All simple paths that start at ``source`` (igraph ``get_all_simple_paths``)."""
    source = int(source)
    if source not in G:
        return []
    paths = []
    for target in list(G.nodes()):
        if int(target) == source:
            continue
        try:
            paths.extend(nx.all_simple_paths(G, source, int(target)))
        except (nx.NetworkXError, nx.NodeNotFound):
            continue
    return paths


def unique_atoms_on_paths(paths):
    if not paths:
        return np.array([], dtype=int)
    return np.unique([int(n) for path in paths for n in path])


def diameter_vertex_path(G):
    """Vertices along a graph-diameter shortest path (igraph ``get_diameter``)."""
    if G.number_of_nodes() == 0:
        return []
    work = G
    if isinstance(G, nx.Graph) and not G.is_directed() and not nx.is_connected(G):
        giant = max(nx.connected_components(G), key=len)
        work = G.subgraph(giant).copy()
    best = []
    for _, dests in nx.all_pairs_shortest_path(work):
        for path in dests.values():
            if len(path) > len(best):
                best = path
    return [int(v) for v in best]
