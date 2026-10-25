"""
Controls for a descriptor search: is the selected model better than a search on noise, and is
it better than another model on the same samples.

    null = selection_null(X, y, combos, lam=0)          # the search itself, on permuted y
    null["p"], null["p95"]                               # p of the real best, floor of the null
    cmp = paired_sign_flip(y, pred_ours, pred_theirs)    # which model is closer, per sample

``selection_null`` builds each candidate's hat matrix once and scores every permutation through
it, so 1000 permutations cost about as much as one search.

Further controls (ported from the DescriPyTor paper notes, Aug-Oct 2026, each reproducing its stored output):

    family_null(X, y, ["theta_a", "theta_b"])         # does this family beat random same-size column sets?
    slot_null(X, y, slot=metal_cols, family=angles)   # in a search forced to use a slot, is the family's win chance?
    equivalence(y, pred_ours, pred_theirs)            # TOST on the squared-error difference
    calibration(protocol)                             # do the 90 % intervals of a BASSA protocol run hold?
    domain_flags(X, weights); rank_check(d, err)      # nearest-neighbour out-of-domain flag, and whether it ranks error
"""
from __future__ import annotations

import itertools
import math

import numpy as np
from scipy import stats

from .ridge_search import hat_matrices, q2 as _q2

__all__ = ["selection_null", "paired_sign_flip", "family_null", "slot_null", "equivalence", "pooled_equivalence",
           "calibration", "nn_distance", "domain_flags", "rank_check"]


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


def _loo_q2_batches(X, Y, combos, lam, batch):
    """Q2_LOO (n_combos, n_responses) of every combination on every row of Y."""
    sst = ((Y - Y.mean(1, keepdims=True)) ** 2).sum(1)
    out = []
    for s in range(0, len(combos), batch):
        H, h = hat_matrices(X, combos[s:s + batch], lam)
        resid = (Y[None] - np.einsum("bij,pj->bpi", H, Y)) / h[:, None, :]
        out.append(1 - (resid ** 2).sum(-1) / sst)
    return np.concatenate(out)


def family_null(X, y, family, n: int = 1000, lam: float = 0.0, seed: int = 0, batch: int = 500) -> dict:
    """Q2_LOO of the model on ``family`` (column names of the DataFrame ``X``) against ``n`` random column sets of the
    same size drawn from the other columns. A family that beats few random families does not predict as a family; it
    predicts as any column set of that size does.

    Returns ``real``, ``null``, ``percentile`` (share of random families below the real one), ``p95`` and
    ``p`` = (1 + #null >= real) / (1 + n)."""
    cols = list(X.columns)
    fam = [cols.index(c) for c in family]
    rest = np.array([i for i in range(len(cols)) if i not in fam])
    rng = np.random.default_rng(seed)
    combos = np.vstack([fam] + [rng.choice(rest, len(fam), replace=False) for _ in range(n)])
    q = _loo_q2_batches(X.to_numpy(float), np.asarray(y, float)[None], combos, lam, batch)[:, 0]
    real, null = float(q[0]), q[1:]
    return dict(real=real, null=null, percentile=float((null < real).mean() * 100), p95=float(np.percentile(null, 95)),
                p=float((1 + (null >= real).sum()) / (1 + n)))


def slot_null(X, y, slot, family, k: int = 3, n_perm: int = 1000, lam: float = 0.0, seed: int = 0, batch: int = 500) -> dict:
    """In an exhaustive k-term search where every model must contain a ``slot`` column, how often does chance fill the
    slot with a ``family`` member? (e.g. slot = the metal-referenced columns, family = its angles.)

    The search runs on ``y`` and on ``n_perm`` permutations; the share of permutations whose best model contains a
    family member is the chance rate (often well above len(family) / len(slot)). Returns the real best ``combo`` and
    ``q2``, ``real_has_family``, ``rate``, ``naive_rate`` and ``p`` (only when the real winner has a family member).
    Cost: one n x n hat matrix per k-subset holding a slot column, computed in batches."""
    cols = list(X.columns)
    si, fi = {cols.index(c) for c in slot}, {cols.index(c) for c in family}
    combos = np.array([c for c in itertools.combinations(range(len(cols)), k) if si & set(c)])
    has = np.array([bool(fi & set(c)) for c in combos])
    y = np.asarray(y, float)
    rng = np.random.default_rng(seed)
    Y = np.vstack([y] + [rng.permutation(y) for _ in range(n_perm)])
    Xv = X.to_numpy(float)
    best_q, best_i = np.full(len(Y), -np.inf), np.zeros(len(Y), int)
    cols_r = np.arange(len(Y))
    for s in range(0, len(combos), batch):
        q = _loo_q2_batches(Xv, Y, combos[s:s + batch], lam, batch)
        j = q.argmax(0)
        better = q[j, cols_r] > best_q
        best_q[better], best_i[better] = q[j, cols_r][better], (s + j)[better]
    hit = has[best_i]
    return dict(combo=[cols[i] for i in combos[best_i[0]]], q2=float(best_q[0]), real_has_family=bool(hit[0]),
                rate=float(hit[1:].mean()), naive_rate=len(fi) / len(si),
                p=float((hit[1:].sum() + 1) / (n_perm + 1)) if hit[0] else None)


def _tost(z, se, df):
    t90, t95 = stats.t.ppf(0.95, df), stats.t.ppf(0.975, df)
    return z, z - t90 * se, z + t90 * se, z - t95 * se, z + t95 * se


def _tost_verdict(m, lo90, hi90, lo95, hi95, D):
    eq, diff = (-D < lo90 and hi90 < D), (lo95 > 0 or hi95 < 0)
    if eq:
        return f"equivalent within {D:g}" + (" (difference real but below the margin)" if diff else "")
    if diff:
        return "different, " + ("A better" if lo95 > 0 else "B better")
    return "underpowered"


def _tost_rows(q, n, margins):
    return [dict(n=n, margin=D, dq2=float(q[0]), ci90_lo=float(q[1]), ci90_hi=float(q[2]), ci95_lo=float(q[3]),
                 ci95_hi=float(q[4]), verdict=_tost_verdict(*q, D)) for D in margins]


def _gain(y, a, b):
    y, a, b = (np.asarray(v, float) for v in (y, a, b))
    return ((y - b) ** 2 - (y - a) ** 2) / y.var()


def equivalence(y, pred_a, pred_b, margins=(0.05, 0.01)) -> list:
    """Two one-sided tests on the per-compound gain z = ((y - b)² - (y - a)²) / var(y), whose mean is ΔQ² (A minus B).

    Per margin D: "equivalent within D" when the 90 % CI lies inside ±D, "different" when the 95 % CI excludes 0,
    otherwise "underpowered" -- which is what equal Q² on 20-30 compounds usually is. Equal accuracy is not equivalence."""
    z = _gain(y, pred_a, pred_b)
    return _tost_rows(_tost(z.mean(), z.std(ddof=1) / np.sqrt(len(z)), len(z) - 1), len(z), margins)


def pooled_equivalence(sets, margins=(0.05, 0.01)) -> list:
    """``equivalence`` pooled over data sets given as (y, pred_a, pred_b) triples; each set keeps its own var(y) and the
    SE is the block-pooled one (df = N - number of sets)."""
    zs = [_gain(*s) for s in sets]
    N = sum(map(len, zs))
    se = np.sqrt(sum(len(z) * z.var(ddof=1) for z in zs)) / N
    return _tost_rows(_tost(np.concatenate(zs).mean(), se, N - len(zs)), N, margins)


def _wilson(k, n, z=1.96):
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h


def calibration(protocol: dict, higher_is_better: bool = True) -> dict:
    """Coverage of the 90 % intervals of a BASSA protocol run (``bassa_explorer/protocol.py`` output or its json).

    5-fold coverage is the mean over deals (with the lowest and highest deal); the leave-one-out coverage with a Wilson
    95 % interval sets the verdict: calibrated (0.90 inside), overconfident (all below) or conservative (all above).
    ``best_quarter_*`` covers the ceil(n/4) best compounds: their coverage, and how much better they measured than
    predicted in response SDs (positive: the model pulls the best toward the mean)."""
    cp = protocol["compounds"]
    n = len(cp)
    y = np.array([q["y"] for q in cp])
    sd = y.std(ddof=1)
    deal = [q["cover"] for q in protocol["per_deal"]]
    k = sum(q["loo_lo"] <= q["y"] <= q["loo_hi"] for q in cp)
    lo, hi = _wilson(k, n)
    sign = 1 if higher_is_better else -1
    best = sorted(range(n), key=lambda i: y[i], reverse=higher_is_better)[:max(1, math.ceil(n / 4))]
    return dict(n=n, q2_5fold=protocol["k5"]["q2"], mae_5fold=protocol["k5"]["mae"],
                cover_5fold=sum(deal) / len(deal), cover_5fold_min=min(deal), cover_5fold_max=max(deal),
                cover_loo=k / n, cover_loo_ci_lo=lo, cover_loo_ci_hi=hi,
                verdict="calibrated" if lo <= 0.9 <= hi else "overconfident" if hi < 0.9 else "conservative",
                width_sd=sum(q["hi"] - q["lo"] for q in cp) / n / sd,
                cover_best_quarter=sum(cp[i]["inside"] for i in best) / len(best),
                best_quarter_under_sd=sum(sign * (cp[i]["y"] - cp[i]["pred"]) for i in best) / len(best) / sd,
                never_inside=[q["name"] for q in cp if q["inside"] == 0])


def nn_distance(Xtr, Xte, weights=None) -> np.ndarray:
    """Distance of each test row to its nearest training row, columns z-scored on the training rows (ddof 1) and
    weighted by ``weights`` (e.g. BASSA inclusion frequencies, normalised to sum 1). Constant columns count zero."""
    Xtr, Xte = np.asarray(Xtr, float), np.asarray(Xte, float)
    w = np.ones(Xtr.shape[1]) if weights is None else np.asarray(weights, float)
    m, s = Xtr.mean(0), Xtr.std(0, ddof=1)
    s[s == 0] = np.inf
    f = np.sqrt(w / w.sum())
    a, b = (Xtr - m) / s * f, (Xte - m) / s * f
    return np.sqrt(((b[:, None, :] - a[None, :, :]) ** 2).sum(-1)).min(1)


def _loo_nn(X, w):
    return np.array([nn_distance(np.delete(X, i, 0), X[i:i + 1], w)[0] for i in range(len(X))])


def domain_flags(X, weights=None, pct: float = 95) -> tuple:
    """Leave-one-out NN distance of each compound to the others, and whether it is out of domain: above the ``pct``-th
    percentile of the other compounds' own leave-one-out distances (computed without it). Returns (d, flag).
    For new compounds: nn_distance(X, X_new, w) > np.percentile(d, pct)."""
    X = np.asarray(X, float)
    d = _loo_nn(X, weights)
    flag = np.array([d[i] > np.percentile(_loo_nn(np.delete(X, i, 0), weights), pct) for i in range(len(X))])
    return d, flag


def rank_check(score, error, n_perm: int = 10000, seed: int = 0) -> dict:
    """Spearman rho between an out-of-domain score and the absolute error, with a one-sided permutation p."""
    rho = lambda a, b: float(np.corrcoef(stats.rankdata(a), stats.rankdata(b))[0, 1])
    r = rho(score, error)
    rng = np.random.default_rng(seed)
    null = np.array([rho(score, rng.permutation(error)) for _ in range(n_perm)])
    return dict(rho=r, p=float((1 + (null >= r).sum()) / (n_perm + 1)))
