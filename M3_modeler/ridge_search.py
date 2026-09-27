"""
Exhaustive small-model ridge search scored by closed-form leave-one-out, with nested selection
and the size-matched constraint null -- the Case Study 3 protocol of the DescriPyTor paper.

Every k-descriptor combination is fitted as ridge regression on standardized columns
(penalty ``lam`` on the slopes, free intercept). Its LOO predictions come from the hat matrix,
so no model is refitted per fold, and combinations are processed in batches. This makes an
exhaustive search over all C(58, 3) = 30,856 triples take seconds.

    combos = all_combos(X.shape[1], 3)
    pred = loo_predictions(X, y, combos)          # (n_combos, n)
    score = adj_q2(pred, y, 3)
    best = combos[score.argmax()]
    nested = nested_loo(X, y, combos)             # selection repeated inside every LOO fold
    null = constraint_null(X, y, combos)          # every 2-column family as a required term

Unlike ``nested_cv.nested_cv_search`` (repeated k-fold, OLS, champion by count), selection here
is by adjusted Q2_LOO and the outer loop is leave-one-out.
"""
from __future__ import annotations

import itertools

import numpy as np

__all__ = ["all_combos", "loo_predictions", "adj_q2", "nested_loo", "constraint_null"]


def all_combos(n_features: int, k: int) -> np.ndarray:
    """Every k-subset of the columns, as an (n_combos, k) index array."""
    return np.array(list(itertools.combinations(range(n_features), k)))


def _standardize(X):
    return (X - X.mean(0)) / X.std(0, ddof=1)


def _penalty(k, lam):
    R = lam * np.eye(k + 1)
    R[0, 0] = 0.0                                  # the intercept is not penalized
    return R


def loo_predictions(X, y, combos, lam: float = 1.0, batch: int = 4000) -> np.ndarray:
    """LOO predictions of every combination, from the hat matrix. Returns (n_combos, n)."""
    X, y = np.asarray(X, float), np.asarray(y, float)
    m, k = len(y), combos.shape[1]
    Xs, R = _standardize(X), _penalty(k, lam)
    out = np.empty((len(combos), m))
    for s in range(0, len(combos), batch):
        idx = combos[s:s + batch]
        Z = np.ones((len(idx), m, k + 1))
        Z[:, :, 1:] = Xs[:, idx].transpose(1, 0, 2)
        H = Z @ np.linalg.inv(np.einsum("bmi,bmj->bij", Z, Z) + R) @ Z.transpose(0, 2, 1)
        h = 1 - np.einsum("bmm->bm", H)
        out[s:s + batch] = y - (y - H @ y) / h
    return out


def adj_q2(pred, y, k: int):
    """Adjusted Q2 of LOO predictions (last axis = samples) for k descriptors."""
    y = np.asarray(y, float)
    n = len(y)
    q2 = 1 - ((y - pred) ** 2).sum(-1) / ((y - y.mean()) ** 2).sum()
    return 1 - (1 - q2) * (n - 1) / (n - k - 1)


def _fit_predict(X_train, y_train, x_test, lam):
    mu, sd = X_train.mean(0), X_train.std(0, ddof=1)
    Z = np.c_[np.ones(len(y_train)), (X_train - mu) / sd]
    w = np.linalg.solve(Z.T @ Z + _penalty(X_train.shape[1], lam), Z.T @ y_train)
    return np.r_[1.0, (x_test - mu) / sd] @ w


def nested_loo(X, y, combos, lam: float = 1.0) -> dict:
    """Hold each sample out, select the best combination on the rest by adjusted Q2_LOO, then
    predict the held-out sample with it. Returns the per-fold selections, the held-out
    predictions and the nested R2."""
    X, y = np.asarray(X, float), np.asarray(y, float)
    n, k = len(y), combos.shape[1]
    selected, pred = [], np.empty(n)
    for i in range(n):
        keep = np.arange(n) != i
        best = combos[adj_q2(loo_predictions(X[keep], y[keep], combos, lam), y[keep], k).argmax()]
        selected.append(tuple(int(j) for j in best))
        pred[i] = _fit_predict(X[keep][:, best], y[keep], X[i, best], lam)
    r2 = 1 - ((y - pred) ** 2).sum() / ((y - y.mean()) ** 2).sum()
    return dict(selected=selected, pred=pred, nested_r2=float(r2))


def constraint_null(X, y, combos, lam: float = 1.0, batch: int = 4000) -> dict:
    """Nested R2 when the search is restricted to combinations containing a given 2-column
    family, for every one of the C(p, 2) families.

    Requiring a family shrinks the search, and a smaller search alone changes selection
    variance. Ranking a chosen family among all families of the same size controls for that,
    which a y-scramble does not. Returns ``pairs`` (C(p,2), 2) and ``nested_r2`` aligned to it;
    use ``nested_r2_for(pair)`` for one family.
    """
    X, y = np.asarray(X, float), np.asarray(y, float)
    n, k = len(y), combos.shape[1]
    Q = np.empty((n, len(combos)))                 # inner adj Q2 of every combo, per fold
    P = np.empty((n, len(combos)))                 # its prediction of the held-out sample
    R = _penalty(k, lam)
    for i in range(n):
        keep = np.arange(n) != i
        Xk = X[keep]
        mu, sd = Xk.mean(0), Xk.std(0, ddof=1)
        Xs, xi = (Xk - mu) / sd, (X[i] - mu) / sd
        Q[i] = adj_q2(loo_predictions(Xk, y[keep], combos, lam, batch), y[keep], k)
        for s in range(0, len(combos), batch):
            idx = combos[s:s + batch]
            Z = np.ones((len(idx), n - 1, k + 1))
            Z[:, :, 1:] = Xs[:, idx].transpose(1, 0, 2)
            w = np.linalg.solve(np.einsum("bmi,bmj->bij", Z, Z) + R,
                                np.einsum("bmi,m->bi", Z, y[keep])[..., None])[..., 0]
            P[i, s:s + batch] = w[:, 0] + np.einsum("bj,bj->b", w[:, 1:], xi[idx])
    sst = ((y - y.mean()) ** 2).sum()

    def nested_r2_for(pair):
        idx = np.where(np.isin(combos, pair).any(1))[0]
        sel = idx[np.argmax(Q[:, idx], 1)]
        return float(1 - ((y - P[np.arange(n), sel]) ** 2).sum() / sst)

    pairs = all_combos(X.shape[1], 2)
    return dict(pairs=pairs, nested_r2=np.array([nested_r2_for(p) for p in pairs]),
                nested_r2_for=nested_r2_for)
