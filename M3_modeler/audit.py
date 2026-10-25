"""
Checks on a descriptor pool before any model sees the response.

    flags = zaudit(X)                          # columns one compound dominates (degenerate or a real outlier)
    kept, clusters = reduce_pool(X, 50)        # one representative per correlation cluster, response never read
    icc1([values_of_mol1, values_of_mol2, ...])  # how much of a quantity's spread is between molecules, not conformers
    reliability_gate(columns, icc_of)          # drop the columns whose conformer ICC is poor

Ported from the DescriPyTor paper notes (zaudit.py, filter50.py, feature_reliability/reliability.py and
reliability_gate.py, Aug-Oct 2026); each reproduces the stored paper output.
"""
from __future__ import annotations

import re

import numpy as np
import pandas as pd

__all__ = ["zaudit", "reduce_pool", "icc1", "icc_bootstrap", "reliability_gate"]


def _numeric(X: pd.DataFrame) -> pd.DataFrame:
    X = X.loc[:, X.apply(pd.api.types.is_numeric_dtype)]
    return X.loc[:, X.notna().all() & (X.std() > 1e-12)]


def zaudit(X: pd.DataFrame, zcut: float = 6.0, degenerate_ratio: float = 0.25) -> pd.DataFrame:
    """Columns where one held-out compound sits more than ``zcut`` leave-one-out SDs from the others.

    For each column and row i: z_i = |x_i - mean(x_-i)| / sd(x_-i) and sd_ratio = sd(x_-i) / sd(x). A flagged column with
    sd_ratio < ``degenerate_ratio`` is degenerate (one compound carries all its variance, e.g. an asym column on a C2 set,
    the cc sub_B1_asym retraction); at or above it the compound is a genuine outlier. Non-numeric, constant and
    incomplete columns are skipped. Returns one row per flagged column, worst first.
    """
    X = _numeric(X)
    n = len(X)
    a = X.to_numpy(float)
    m = (a.sum(0) - a) / (n - 1)                                    # leave-one-out mean
    v = ((a ** 2).sum(0) - a ** 2) / (n - 1) - m ** 2                # leave-one-out variance (ddof 0)
    sd = np.sqrt(np.maximum(v, 0.0))
    z = np.abs(a - m) / np.maximum(sd, 1e-12)
    j = z.argmax(0)
    cols = np.arange(a.shape[1])
    zmax, ratio = z[j, cols], sd[j, cols] / a.std(0)
    rows = [dict(column=X.columns[c], max_z=float(zmax[c]), sd_ratio=float(ratio[c]), row=X.index[j[c]],
                 value=float(a[j[c], c]), median_rest=float(np.median(np.delete(a[:, c], j[c]))),
                 kind="degenerate" if ratio[c] < degenerate_ratio else "outlier")
            for c in cols if zmax[c] > zcut]
    out = pd.DataFrame(rows, columns=["column", "max_z", "sd_ratio", "row", "value", "median_rest", "kind"])
    return out.sort_values("max_z", ascending=False, ignore_index=True)


def reduce_pool(X: pd.DataFrame, n_keep: int = 50, tier=None, min_levels: int = 4, drop: str | None = None):
    """Cut a wide table to ``n_keep`` columns without reading the response.

    1. Clean: drop incomplete and constant columns, columns with fewer than ``min_levels`` distinct values (to 3
       decimals; a two-level column is a class label, not a dial), columns whose name matches the regex ``drop``, and
       exact duplicates after z-scoring.
    2. Cluster on 1 - |Pearson r|, average linkage, cut at ``n_keep`` clusters.
    3. Keep one member per cluster: the lowest ``tier(name)`` (default: every column tier 0), then the member with the
       highest mean |r| to its cluster (the medoid).

    Returns (the reduced table, a cluster table: every cleaned column, its cluster, whether kept, |r| with the kept one).
    """
    from scipy.cluster.hierarchy import fcluster, linkage
    from scipy.spatial.distance import squareform

    tier = tier or (lambda c: 0)
    x = X.select_dtypes("number")
    bad = [c for c in x.columns if x[c].isna().any() or x[c].std() < 1e-9 or x[c].round(3).nunique() < min_levels
           or (drop and re.search(drop, str(c)))]
    x = x.drop(columns=bad)
    z = ((x - x.mean()) / x.std()).round(6)
    x = x.loc[:, ~z.T.duplicated().to_numpy()]
    r = x.corr().abs().fillna(0.0).to_numpy()
    np.fill_diagonal(r, 1.0)
    k = min(n_keep, x.shape[1])
    lab = fcluster(linkage(squareform(1 - r, checks=False), "average"), k, "maxclust")
    keep, rows = [], []
    for cl in sorted(set(lab)):
        idx = np.where(lab == cl)[0]
        best = min(tier(x.columns[i]) for i in idx)
        rep = max((i for i in idx if tier(x.columns[i]) == best), key=lambda i: r[i, idx].mean())
        keep.append(x.columns[rep])
        rows += [dict(column=x.columns[i], cluster=int(cl), kept=int(i == rep), tier=tier(x.columns[i]),
                      r_with_kept=round(float(r[i, rep]), 3)) for i in idx]
    return x[keep], pd.DataFrame(rows)


def icc1(groups) -> float:
    """One-way random-effects ICC(1) with unequal group sizes; ``groups`` holds one array of values per molecule
    (e.g. a quantity over that molecule's conformers). 1: conformers agree; about 0: conformer choice is the spread."""
    groups = [np.asarray(g, float) for g in groups]
    k = len(groups)
    n = np.array([len(v) for v in groups])
    N = n.sum()
    gm = np.concatenate(groups).mean()
    msb = sum(len(v) * (v.mean() - gm) ** 2 for v in groups) / (k - 1)
    msw = sum(((v - v.mean()) ** 2).sum() for v in groups) / (N - k)
    n0 = (N - (n ** 2).sum() / N) / (k - 1)
    return float((msb - msw) / (msb + (n0 - 1) * msw))


def icc_bootstrap(groups, n_boot: int = 1000, seed: int = 0) -> dict:
    """ICC(1) with a 95 % bootstrap interval over molecules (resampled with replacement)."""
    groups = [np.asarray(g, float) for g in groups]
    rng = np.random.default_rng(seed)
    bs = [icc1([groups[i] for i in rng.integers(0, len(groups), len(groups))]) for _ in range(n_boot)]
    return dict(icc=icc1(groups), lo=float(np.nanpercentile(bs, 2.5)), hi=float(np.nanpercentile(bs, 97.5)))


def reliability_gate(columns, icc_of, threshold: float = 0.50) -> pd.DataFrame:
    """Which columns to drop because their conformer ICC is poor (< ``threshold``, Koo & Li 2016).

    ``icc_of`` maps a column name to its ICC (a dict or a function); a column it does not know (None / missing) was not
    measured and is kept. Returns column, icc, drop."""
    get = icc_of.get if isinstance(icc_of, dict) else icc_of
    rows = []
    for c in columns:
        v = get(c)
        v = None if v is None or (isinstance(v, float) and np.isnan(v)) else float(v)
        rows.append(dict(column=c, icc=v, drop=v is not None and v < threshold))
    return pd.DataFrame(rows)
