"""
Controls for a descriptor search: is the selected model better than a search on noise, and is
it better than another model on the same samples.

    null = selection_null(X, y, combos, lam=0)          # the search itself, on permuted y
    null["p"], null["p95"]                               # p of the real best, floor of the null
    cmp = paired_sign_flip(y, pred_ours, pred_theirs)    # which model is closer, per sample

``selection_null`` builds each candidate's hat matrix once and scores every permutation through
it, so 1000 permutations cost about as much as one search.
"""
from __future__ import annotations

import numpy as np

from .ridge_search import hat_matrices, q2 as _q2

__all__ = ["selection_null", "paired_sign_flip"]


def selection_null(X, y, combos, lam: float = 1.0, n_perm: int = 1000, seed: int = 0,
                   batch: int = 500) -> dict:
    """Best Q2_LOO over ``combos`` for the real ``y`` and for ``n_perm`` permutations of it.

    The permutations are ``np.random.default_rng(seed).permutation(y)``, drawn in order, so a
    script that draws them the same way gets the same null. Restrict ``combos`` to test a
    constrained search (e.g. every triple that contains one hypothesised descriptor).

    Returns ``real`` (best Q2 on y), ``null`` (n_perm best Q2s), ``median``, ``p95``, ``max``,
    ``n_geq`` (permutations at or above ``real``) and ``p`` = (1 + n_geq) / (1 + n_perm).
    """
    y = np.asarray(y, float)
    rng = np.random.default_rng(seed)
    Y = np.vstack([y] + [rng.permutation(y) for _ in range(n_perm)])       # (1 + n_perm, n)
    sst = ((Y - Y.mean(1, keepdims=True)) ** 2).sum(1)
    best = np.full(len(Y), -np.inf)
    for s in range(0, len(combos), batch):
        H, h = hat_matrices(X, combos[s:s + batch], lam)                   # (b, n, n), (b, n)
        resid = (Y[None] - np.einsum("bij,pj->bpi", H, Y)) / h[:, None, :]
        best = np.maximum(best, (1 - (resid ** 2).sum(-1) / sst).max(0))
    real, null = float(best[0]), best[1:]
    n_geq = int((null >= real).sum())
    return dict(real=real, null=null, median=float(np.median(null)), p95=float(np.percentile(null, 95)),
                max=float(null.max()), n_geq=n_geq, p=(1 + n_geq) / (1 + n_perm))


def paired_sign_flip(y, pred_a, pred_b, n_draws: int = 100000, seed: int = 0) -> dict:
    """Is model A closer to the data than model B, sample by sample?

    The statistic is the mean of |y - pred_b| - |y - pred_a| (positive when A is closer). Its
    null flips the sign of each sample's difference at random. Returns ``closer`` (samples
    where A is strictly closer), ``mean_gain``, ``p`` (two-sided) and ``residual_r``, the
    correlation of the two models' residuals.
    """
    y, a, b = (np.asarray(v, float) for v in (y, pred_a, pred_b))
    diff = np.abs(y - b) - np.abs(y - a)
    rng = np.random.default_rng(seed)
    flips = (rng.choice([-1, 1], (n_draws, len(y))) * diff).mean(1)
    return dict(closer=int((np.abs(y - a) < np.abs(y - b)).sum()), mean_gain=float(diff.mean()),
                p=float((np.abs(flips) >= abs(diff.mean())).mean()),
                residual_r=float(np.corrcoef(y - a, y - b)[0, 1]))
