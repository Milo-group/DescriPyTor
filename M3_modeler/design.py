"""
Design: does a model find the best compounds, and which compound to make next.

    top_k_recovery(y, pred)               of the k best measured, how many are among the k best predicted (hypergeometric p)
    leave_best_out(X, y)                  fit without the best quarter: are they still predicted above the rest?
    suggest(X, y, candidates)             rank new compounds by expected improvement, with an interval and a domain flag

Read the predictions for rank, not size. Across the nine paper cases the best quarter is under-predicted every time (by 0.2 to
1.2 sd, interval calibration, 18 Oct 2026): a linear model shrinks toward the training mean, so a candidate predicted at the top
is likely to be the best but its value is likely to be an underestimate. ``top_k_recovery`` reports that shrinkage.

Validated in the paper notes (design tests, 15 Oct 2026): top-quarter recovery passes on CS2, cyclopropanation, Wu–Doyle and
the CS1 acids; leave-best-out passes on three and fails on Wu–Doyle, whose best ligands sit at the extreme of every size column.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

__all__ = ["top_k_recovery", "leave_best_out", "suggest", "ols_fit"]


def _k(n, k):
    return max(1, int(np.ceil(n / 4))) if k is None else int(k)        # a quarter, rounded up (17 → 5, 30 → 8, as in the paper tests)


def top_k_recovery(y, pred, k: int | None = None, higher_is_better: bool = True) -> dict:
    """Hits = how many of the k best measured compounds are among the k best predicted (k defaults to a quarter).

    p = P(X ≥ hits) for k draws from n with k successes (hypergeometric); chance = k²/n. Also Spearman ρ, and the shrinkage of
    the best: the mean of (measured − predicted)/sd over the k best measured, signed so that positive = under-predicted."""
    y, pred = np.asarray(y, float), np.asarray(pred, float)
    n, k = len(y), _k(len(y), k)
    s = 1 if higher_is_better else -1
    best_meas = set(np.argsort(-s * y, kind="stable")[:k])
    best_pred = set(np.argsort(-s * pred, kind="stable")[:k])
    hits = len(best_meas & best_pred)
    idx = np.array(sorted(best_meas))
    return dict(n=n, k=k, hits=hits, chance=k * k / n, p=float(stats.hypergeom.sf(hits - 1, n, k, k)),
                spearman=float(stats.spearmanr(y, pred)[0]),
                best_shrinkage_sd=float(s * (y[idx] - pred[idx]).mean() / y.std(ddof=1)))


def ols_fit(Xtr, ytr):
    """The default fit: OLS with an intercept; returns a predict function."""
    Z = np.c_[np.ones(len(Xtr)), Xtr]
    b = np.linalg.lstsq(Z, ytr, rcond=None)[0]
    return lambda Xte: np.c_[np.ones(len(Xte)), Xte] @ b


def leave_best_out(X, y, k: int | None = None, higher_is_better: bool = True, fit=ols_fit, n_null: int = 20, seed: int = 500) -> dict:
    """Fit on all but the k best (a quarter by default) and predict them.

    S = mean over the held-out best of s·(prediction − mean y_train)/sd y_train: how far above the training compounds the model
    puts compounds it never saw (s = −1 when lower is better). Null: the same fit with y shuffled among the training compounds,
    ``n_null`` times (``default_rng(seed + i)``, as in the paper test). Passes when S is above every null value
    (p = 1/(n_null + 1)). Also: how many held-out compounds are predicted beyond the best training value."""
    X, y = np.asarray(X, float).reshape(len(y), -1), np.asarray(y, float)
    n, k = len(y), _k(len(y), k)
    s = 1 if higher_is_better else -1
    held = np.argsort(-s * y, kind="stable")[:k]
    tr = np.setdiff1d(np.arange(n), held)
    ym, ys = y[tr].mean(), y[tr].std(ddof=1)
    pred = fit(X[tr], y[tr])(X[held])
    S = float(np.mean(s * (pred - ym) / ys))
    null = []
    for i in range(n_null):
        yp = np.random.default_rng(seed + i).permutation(y[tr])
        null.append(float(np.mean(s * (fit(X[tr], yp)(X[held]) - ym) / ys)))
    null = np.array(null)
    best_train = (y[tr].max() if higher_is_better else y[tr].min())
    beyond = int(((pred > best_train) if higher_is_better else (pred < best_train)).sum())
    return dict(n=n, k=k, S=S, null_max=float(null.max()), p=float((1 + (null >= S).sum()) / (1 + n_null)),
                passes=bool(S > null.max()), held_out=held.tolist(), predicted=pred.tolist(), beyond_best_training=beyond)


def _nn_distances(Z):
    D = np.sqrt(((Z[:, None, :] - Z[None, :, :]) ** 2).sum(-1))
    np.fill_diagonal(D, np.inf)
    return D.min(1)


def suggest(X, y, candidates, higher_is_better: bool = True, batch: int = 1, xi: float = 0.0, alpha: float = 0.10,
            domain: str = "first") -> pd.DataFrame:
    """Rank candidate compounds for the next round from an OLS model on the same columns.

    X, y        : the measured compounds (DataFrame or array) and their response
    candidates  : the proposed compounds, same columns (a DataFrame keeps their names)
    batch       : how many to pick; after each pick the model treats it as measured at its prediction, which shrinks the
                  uncertainty around it, so the next pick goes elsewhere ("kriging believer")
    xi          : improvement margin for expected improvement
    domain      : "first" picks candidates inside the measured domain before any outside it (expected improvement rewards
                  extrapolation, and outside the domain the predictions were not reliable in the paper cases); "ignore" ranks
                  by expected improvement alone

    Columns of the result: prediction, its (1 − alpha) prediction interval (t distribution), sd, expected improvement over the
    best measured value, the nearest-neighbour distance to the measured compounds in standardised columns, ``outside_domain``
    when that distance exceeds the 95th percentile of the measured compounds' own nearest-neighbour distances (the flag that
    tracked error in the paper cases; leverage did not), the columns outside the measured range, and ``pick`` (1..batch)."""
    if domain not in ("first", "ignore"):
        raise ValueError('domain must be "first" or "ignore"')
    names = list(candidates.index) if isinstance(candidates, pd.DataFrame) else list(range(len(candidates)))
    cols = list(X.columns) if isinstance(X, pd.DataFrame) else [f"x{j}" for j in range(np.asarray(X).shape[1])]
    Xa, ya, Ca = np.asarray(X, float), np.asarray(y, float), np.asarray(candidates, float)
    s = 1 if higher_is_better else -1
    mu_x, sd_x = Xa.mean(0), Xa.std(0) + 1e-12
    Zt, Zc = (Xa - mu_x) / sd_x, (Ca - mu_x) / sd_x
    nn_train = _nn_distances(Zt)
    nn_cand = np.sqrt(((Zc[:, None, :] - Zt[None, :, :]) ** 2).sum(-1)).min(1)
    cut = float(np.quantile(nn_train, 0.95))
    lo, hi = Xa.min(0), Xa.max(0)
    outside = [", ".join(c for c, v, a, b in zip(cols, row, lo, hi) if v < a or v > b) for row in Ca]
    A = np.c_[np.ones(len(Xa)), Xa]
    b = np.linalg.lstsq(A, ya, rcond=None)[0]
    dof = max(len(ya) - A.shape[1], 1)
    s2 = float(((ya - A @ b) ** 2).sum() / dof)
    tq = stats.t.ppf(1 - alpha / 2, dof)
    Cc = np.c_[np.ones(len(Ca)), Ca]
    mean = Cc @ b
    best = s * (ya.max() if higher_is_better else ya.min())
    G = A.T @ A
    picks = np.zeros(len(Ca), int)
    first_sd = first_ei = None
    for r in range(1, min(batch, len(Ca)) + 1):
        Gi = np.linalg.pinv(G)
        sd = np.sqrt(s2 * (1 + np.einsum("ij,jk,ik->i", Cc, Gi, Cc)))
        z = (s * mean - best - xi) / sd
        ei = (s * mean - best - xi) * stats.norm.cdf(z) + sd * stats.norm.pdf(z)
        if r == 1:
            first_sd, first_ei = sd, ei
        ei_open = np.where(picks == 0, ei, -np.inf)
        if domain == "first" and ((picks == 0) & (nn_cand <= cut)).any():
            ei_open = np.where(nn_cand <= cut, ei_open, -np.inf)
        j = int(np.argmax(ei_open))
        picks[j] = r
        G = G + np.outer(Cc[j], Cc[j])                     # believed measured: the design gains this row, the fit is unchanged
    out = pd.DataFrame({"prediction": mean, "lower": mean - tq * first_sd, "upper": mean + tq * first_sd, "sd": first_sd,
                        "expected_improvement": first_ei, "nn_distance": nn_cand, "outside_domain": nn_cand > cut,
                        "outside_range": outside, "pick": picks}, index=names)
    out.attrs.update(best_measured=float(s * best), domain_cut=cut, note="rank, not size: the best are usually under-predicted")
    order = np.lexsort((-out.expected_improvement.to_numpy(), out.outside_domain.to_numpy() if domain == "first" else np.zeros(len(out)),
                        np.where(out.pick.to_numpy() > 0, out.pick.to_numpy(), np.iinfo(int).max)))
    return out.iloc[order]
