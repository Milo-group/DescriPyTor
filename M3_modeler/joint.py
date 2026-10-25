"""
Joint descriptor selection across related reactions: the same columns for every reaction, each reaction with its own weights.

When a catalyst (ligand) set is run in several related reactions, choosing descriptors on all of them at once can beat choosing
per reaction, because the choice is no longer made on one small, noisy response. The weights stay per reaction, so a descriptor
can change sign with the electrophile. Responses are never averaged over reactions.

    res = joint_selection({"I": y_I, "II": y_II, "III": y_III}, F, max_terms=2)
    res["passes"], res["gain"], res["q2"]["joint"], res["columns"]       # the gate, per-reaction Q2, the shared pick
    fit_joint({"I": y_I, ...}, F, res["columns"])                        # per-reaction coefficients on all compounds

Validated on BUNNY (four Ni reactions, ΔΔG‡): joint selection picks the Ni charge and mean donor charge, mean leave-one-ligand-out
gain +0.28 Q², p 0.024, above the 95th percentile of random-partner gains (paper notes, ootb_models, 19 Oct 2026).
"""
from __future__ import annotations

from itertools import combinations

import numpy as np
import pandas as pd

__all__ = ["joint_selection", "fit_joint", "press_operators"]


def _ols_pred(Xtr, ytr, Xte):
    Z = np.c_[np.ones(len(Xtr)), Xtr]
    b = np.linalg.lstsq(Z, ytr, rcond=None)[0]
    return np.c_[np.ones(len(Xte)), Xte] @ b


def press_operators(X, designs):
    """M (P, m, m) with PRESS(y) = ||M y||² for OLS with an intercept on each design (a tuple of column indices).

    The leave-one-out residual of a linear smoother is (I - H) y / (1 - h), so one matrix per design scores any response."""
    out = []
    for k in sorted({len(d) for d in designs}):
        D = np.array([d for d in designs if len(d) == k], dtype=int)
        Z = np.concatenate([np.ones((len(D), len(X), 1)), X[:, D].transpose(1, 0, 2)], 2)
        G = np.linalg.pinv(np.einsum("pmi,pmj->pij", Z, Z))
        Hm = np.einsum("pmi,pij,pnj->pmn", Z, G, Z)
        h = np.einsum("pmm->pm", Hm)
        out.append((np.eye(len(X))[None] - Hm) / np.clip(1 - h, 1e-9, None)[..., None])
    order = [d for k in sorted({len(d) for d in designs}) for d in designs if len(d) == k]
    if order != list(designs):
        raise ValueError("designs must be listed by size (all singles, then all pairs, ...)")
    return np.concatenate(out)


def _q2(y, p):
    return float(1 - ((y - p) ** 2).sum() / ((y - y.mean()) ** 2).sum())


def _signflip(d, n=10_000, seed=0):
    """one-sided: the mean of d (separate loss minus joint loss) above zero"""
    S = np.random.default_rng(seed).choice((-1.0, 1.0), size=(n, len(d)))
    return float((1 + ((S * d).mean(1) >= d.mean()).sum()) / (1 + n))


def _align(Y, F):
    R = list(Y)
    ligs = sorted(set().union(*[set(Y[r].index) for r in R]) & set(F.index), key=str)
    if not ligs:
        raise ValueError("no compound is both in the feature table and in a response")
    X = F.loc[ligs].to_numpy(float)
    if np.isnan(X).any():
        raise ValueError("the feature table has missing values for compounds with a response")
    yv = {r: pd.Series(Y[r]).astype(float).reindex(ligs).to_numpy(float) for r in R}
    return R, ligs, X, yv, {r: ~np.isnan(yv[r]) for r in R}


def joint_selection(Y: dict, F: pd.DataFrame, max_terms: int = 2, n_null: int = 200, seed: int = 0, min_gain: float = 0.05) -> dict:
    """Separate against joint selection, leave-one-compound-out, the selection redone in every fold.

    Y : {reaction: Series of the response over compounds} (compounds may be missing from some reactions)
    F : feature table over compounds (one row per compound, numeric columns; the same columns for every reaction)

    In each fold the compound is removed from every reaction. *Separate*: each reaction takes the design (1 to ``max_terms``
    columns, OLS with an intercept) with the lowest PRESS on its own training compounds. *Joint*: every reaction takes the one
    design with the lowest sum over reactions of PRESS / SS. Each reaction is then fitted with its own coefficients and predicts
    the held-out compound. The null keeps the target reaction's response and shuffles the others' (``n_null`` draws), so it
    measures what any partner reactions would give.

    Gate (as in the paper's test): gain = mean over reactions of Q2(joint) - Q2(separate); passes when gain >= ``min_gain``,
    the one-sided sign-flip p on the per-compound loss differences < 0.05, and gain is above the null's 95th percentile.
    """
    R, ligs, X, yv, has = _align(Y, F)
    if len(R) < 2:
        raise ValueError("joint selection needs at least two reactions")
    p = X.shape[1]
    if not 1 <= max_terms <= min(3, p):
        raise ValueError("max_terms must be 1, 2 or 3 and at most the number of columns")
    designs = [d for k in range(1, max_terms + 1) for d in combinations(range(p), k)]
    rng = np.random.default_rng(seed)
    shuf = {r: np.array([rng.permutation(np.flatnonzero(has[r])) for _ in range(n_null)]) for r in R}
    n = len(ligs)
    pred = {a: {r: np.full(n, np.nan) for r in R} for a in ("separate", "joint")}
    predN = {r: np.full((n_null, n), np.nan) for r in R}
    picks = []
    for i in range(n):
        tr = np.arange(n) != i
        crit, critN = {}, {}
        for r in R:
            rows = np.flatnonzero(tr & has[r])
            M = press_operators(X[rows], designs)
            yt = yv[r][rows]
            crit[r] = ((M @ yt) ** 2).sum(1) / ((yt - yt.mean()) ** 2).sum()
            Ys = np.empty((len(rows), n_null))
            for s in range(n_null):
                full = np.full(n, np.nan)
                full[np.flatnonzero(has[r])] = yv[r][shuf[r][s]]
                Ys[:, s] = full[rows]
            critN[r] = (np.einsum("pmn,ns->pms", M, Ys) ** 2).sum(1) / ((Ys - Ys.mean(0)) ** 2).sum(0)
        dJ = designs[int(np.argmin(sum(crit.values())))]
        picks.append(dJ)
        for t in R:
            if not has[t][i]:
                continue
            rows = np.flatnonzero(tr & has[t])
            dS = designs[int(np.argmin(crit[t]))]
            pred["separate"][t][i] = _ols_pred(X[rows][:, dS], yv[t][rows], X[i:i + 1, dS])[0]
            pred["joint"][t][i] = _ols_pred(X[rows][:, dJ], yv[t][rows], X[i:i + 1, dJ])[0]
            totN = crit[t][:, None] + sum(critN[r] for r in R if r != t)            # target true, partners shuffled
            for s, dn in enumerate(np.argmin(totN, 0)):
                predN[t][s, i] = _ols_pred(X[rows][:, designs[dn]], yv[t][rows], X[i:i + 1, designs[dn]])[0]
    q = {a: {r: _q2(yv[r][has[r]], pred[a][r][has[r]]) for r in R} for a in pred}
    gain = float(np.mean([q["joint"][r] - q["separate"][r] for r in R]))
    d = np.concatenate([((yv[r] - pred["separate"][r]) ** 2 - (yv[r] - pred["joint"][r]) ** 2)[has[r]] / np.nanvar(yv[r], ddof=1) for r in R])
    null = np.array([np.mean([_q2(yv[r][has[r]], predN[r][s][has[r]]) - q["separate"][r] for r in R]) for s in range(n_null)])
    p_flip = _signflip(d)
    null95 = float(np.quantile(null, 0.95))
    common = max(set(picks), key=picks.count)
    cols = [F.columns[j] for j in common]
    return dict(reactions={r: int(has[r].sum()) for r in R}, n_compounds=n, n_designs=len(designs), q2=q, gain=gain, p=p_flip,
                null95=null95, passes=bool(gain >= min_gain and p_flip < 0.05 and gain > null95), columns=cols,
                pick_share=picks.count(common) / len(picks), compounds=ligs,
                predictions={a: {r: pd.Series(pred[a][r], index=ligs) for r in R} for a in pred},
                coefficients=fit_joint({r: pd.Series(yv[r], index=ligs) for r in R}, F.loc[ligs], cols))


def fit_joint(Y: dict, F: pd.DataFrame, columns) -> pd.DataFrame:
    """Each reaction's OLS on the shared ``columns``, on all its compounds: intercept, raw coefficients, standardised coefficients
    (coefficient × the column's sd over that reaction's compounds), n and in-sample R². One row per reaction."""
    rows = {}
    for r, y in Y.items():
        y = pd.Series(y).dropna().astype(float)
        idx = [c for c in y.index if c in F.index]
        X = F.loc[idx, list(columns)].to_numpy(float)
        yy = y.loc[idx].to_numpy(float)
        b = np.linalg.lstsq(np.c_[np.ones(len(yy)), X], yy, rcond=None)[0]
        fit = np.c_[np.ones(len(yy)), X] @ b
        rows[r] = {"n": len(yy), "intercept": b[0], **{c: b[k + 1] for k, c in enumerate(columns)},
                   **{f"z_{c}": b[k + 1] * X[:, k].std() for k, c in enumerate(columns)}, "r2_in_sample": _q2(yy, fit)}
    return pd.DataFrame(rows).T
