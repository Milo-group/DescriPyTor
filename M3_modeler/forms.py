"""
Model forms beyond a straight line, each judged against the straight line on the same compounds.

    logistic_link(X, y)            yield-type responses bounded at 0 and 100: binomial GLM against OLS, leave-one-out
    flexible_ceiling(X, y, ref)    a Gaussian process (ARD) and a random forest on the whole pool: does linear miss structure?
    threshold_analysis(F, labels)  one column and one cut that separate good from bad, chosen inside every fold

Every result carries a gate: a gain, a sign-flip p on the per-compound losses and a permutation null, and passes only when the
gain is ≥ 0.05, p < 0.05 and the gain is above the null's 95th percentile. Validated in the paper notes (readable_nonlinear,
ootb_models, 19 Oct 2026): the logistic link lifts Kraken yields (cone + %Vbur, LOO Q² 0.27 → 0.45) and nothing else; the GP
and forest beat the linear model only where the linear model is weak.
"""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd

__all__ = ["logistic_link", "flexible_ceiling", "threshold_analysis"]


def _q2(y, p):
    return float(1 - ((y - p) ** 2).sum() / ((y - y.mean()) ** 2).sum())


def _signflip(d, n=10_000, seed=0):
    """one-sided p that the mean of d (reference loss minus form loss) is above zero"""
    S = np.random.default_rng(seed).choice((-1.0, 1.0), size=(n, len(d)))
    return float((1 + ((S * d).mean(1) >= d.mean()).sum()) / (1 + n))


def _gate(gain, d, null, min_gain=0.05):
    p = _signflip(d)
    n95 = float(np.quantile(null, 0.95)) if len(null) else float("nan")
    return dict(gain=float(gain), p=p, null95=n95, passes=bool(gain >= min_gain and p < 0.05 and (not len(null) or gain > n95)))


def _ols_pred(Xtr, ytr, Xte):
    Z = np.c_[np.ones(len(Xtr)), Xtr]
    b = np.linalg.lstsq(Z, ytr, rcond=None)[0]
    return np.c_[np.ones(len(Xte)), Xte] @ b


def _glm_pred(Xtr, ytr, Xte, scale):
    import statsmodels.api as sm
    A = sm.add_constant(Xtr, has_constant="add")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = sm.GLM(np.clip(ytr / scale, 0, 1), A, family=sm.families.Binomial()).fit()
    return scale * m.predict(sm.add_constant(Xte, has_constant="add"))


def _loo_pair(X, y, scale):
    n = len(y)
    po, pg = np.zeros(n), np.zeros(n)
    for i in range(n):
        k = np.arange(n) != i
        po[i] = _ols_pred(X[k], y[k], X[i:i + 1])[0]
        pg[i] = _glm_pred(X[k], y[k], X[i:i + 1], scale)[0]
    return po, pg


def logistic_link(X, y, scale: float = 100.0, n_null: int = 200, seed: int = 0) -> dict:
    """Binomial GLM (logit link) on y / scale against OLS on the same columns, leave-one-out.

    For yields and conversions, where OLS predicts below 0 or above 100 and a floor-heavy response bends the relation. The gate
    is against OLS and, separately, against OLS clipped to [0, scale], so a pass is not only the bounds. Null: the descriptor
    rows shuffled against y, the same comparison repeated ``n_null`` times."""
    X = np.asarray(X, float).reshape(len(y), -1)
    y = np.asarray(y, float)
    po, pg = _loo_pair(X, y, scale)
    pc = np.clip(po, 0, scale)
    rng = np.random.default_rng(seed)
    null = []
    for _ in range(n_null):
        j = rng.permutation(len(y))
        a, b = _loo_pair(X[j], y, scale)
        null.append(_q2(y, b) - _q2(y, a))
    null = np.array(null)
    return dict(n=len(y), q2_ols=_q2(y, po), q2_logit=_q2(y, pg), q2_clipped=_q2(y, pc),
                mae_ols=float(np.abs(y - po).mean()), mae_logit=float(np.abs(y - pg).mean()),
                vs_ols=_gate(_q2(y, pg) - _q2(y, po), (y - po) ** 2 - (y - pg) ** 2, null),
                vs_clipped=_gate(_q2(y, pg) - _q2(y, pc), (y - pc) ** 2 - (y - pg) ** 2, null),
                pred_ols=po, pred_logit=pg)


def flexible_ceiling(X, y, reference_pred=None, seed: int = 0) -> dict:
    """Leave-one-out Q² of a Gaussian process with one length scale per column (ARD) and of a random forest, on the whole pool.

    ``reference_pred``: the linear model's own out-of-sample predictions (for example a BASSA or nested-CV result); when it is
    omitted the reference is OLS on every column, which overfits a wide pool and makes the ceiling look better than it is.
    The verdict is a reading, not a test: if neither flexible model beats the linear reference, linear is enough on this set;
    if one does, the linear model is missing something, which in the paper sets was one compound or a class more often than
    a curve (check those before reaching for a non-linear form)."""
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.gaussian_process import GaussianProcessRegressor
    from sklearn.gaussian_process.kernels import RBF, ConstantKernel, WhiteKernel
    X = np.asarray(X, float)
    y = np.asarray(y, float)
    n = len(y)
    pg, pr = np.zeros(n), np.zeros(n)
    for i in range(n):
        k = np.arange(n) != i
        Xt = X[k]
        sd = Xt.std(0)
        use = sd > 1e-12
        mu = Xt[:, use].mean(0)
        Z, zt = (Xt[:, use] - mu) / sd[use], (X[i, use] - mu) / sd[use]
        ym, ys = y[k].mean(), y[k].std()
        ker = ConstantKernel(1.0, (1e-3, 1e3)) * RBF(np.ones(use.sum()), (1e-2, 1e3)) + WhiteKernel(0.1, (1e-5, 1e1))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            gp = GaussianProcessRegressor(ker, n_restarts_optimizer=5, random_state=seed).fit(Z, (y[k] - ym) / ys)
        pg[i] = gp.predict(zt[None])[0] * ys + ym
        pr[i] = RandomForestRegressor(500, max_features=1 / 3, min_samples_leaf=2, random_state=seed, n_jobs=-1).fit(Xt, y[k]).predict(X[i][None])[0]
    ref = np.asarray(reference_pred, float) if reference_pred is not None else np.array(
        [_ols_pred(np.delete(X, i, 0), np.delete(y, i), X[i:i + 1])[0] for i in range(n)])
    q = dict(linear=_q2(y, ref), gp=_q2(y, pg), forest=_q2(y, pr))
    gain = max(q["gp"], q["forest"]) - q["linear"]
    reading = ("linear is enough here: neither flexible model beats it" if gain <= 0 else
               f"a flexible model beats the linear one by {gain:.2f}: look for one compound or a class before a curve")
    return dict(q2=q, ceiling_gain=float(gain), reading=reading, reference="given" if reference_pred is not None else "OLS on all columns",
                pred_gp=pg, pred_forest=pr)


def _best_cut(F, lab):
    """For every column, the cut and side with the highest balanced accuracy on (F, lab); returns (score, col, cut, side)."""
    best = (-1.0, 0, 0.0, 1)
    pos, neg = lab.sum(), (~lab).sum()
    for j in range(F.shape[1]):
        o = np.argsort(F[:, j], kind="stable")
        v, l = F[o, j], lab[o]
        tp_left = np.cumsum(l)[:-1]                       # positives at or below each gap
        tn_left = np.cumsum(~l)[:-1]
        ok = v[1:] > v[:-1]                               # a cut only between distinct values
        if not ok.any():
            continue
        ba_hi = 0.5 * ((pos - tp_left) / pos + tn_left / neg)    # "good above the cut"
        ba_lo = 0.5 * (tp_left / pos + (neg - tn_left) / neg)    # "good below the cut"
        for side, ba in ((1, ba_hi), (-1, ba_lo)):
            ba = np.where(ok, ba, -1)
            g = int(np.argmax(ba))
            if ba[g] > best[0]:
                best = (float(ba[g]), j, float((v[g] + v[g + 1]) / 2), side)
    return best


def threshold_analysis(F: pd.DataFrame, labels, n_perm: int = 200, seed: int = 0) -> dict:
    """Threshold analysis (one column, one cut) for a good/bad outcome, e.g. ``labels = ee >= 80``.

    In-sample: the column and cut with the best balanced accuracy over all compounds. Honest score: leave-one-out with the
    column and the cut chosen again in every fold. Null: the whole procedure on shuffled labels (``n_perm`` times); p is the
    share of shuffles whose leave-one-out balanced accuracy reaches the real one."""
    X = F.to_numpy(float)
    lab = np.asarray(labels, bool)
    if lab.all() or not lab.any():
        raise ValueError("labels need both classes")

    def loo_ba(lab_):
        hit = np.zeros(len(lab_), bool)
        cols = []
        for i in range(len(lab_)):
            k = np.arange(len(lab_)) != i
            if lab_[k].all() or not lab_[k].any():
                continue
            _, j, cut, side = _best_cut(X[k], lab_[k])
            cols.append(j)
            hit[i] = ((X[i, j] > cut) if side == 1 else (X[i, j] <= cut)) == lab_[i]
        ba = 0.5 * (hit[lab_].mean() + hit[~lab_].mean())
        return float(ba), cols

    score, j, cut, side = _best_cut(X, lab)
    loo, cols = loo_ba(lab)
    rng = np.random.default_rng(seed)
    null = np.array([loo_ba(rng.permutation(lab))[0] for _ in range(n_perm)])
    picked = pd.Series([F.columns[c] for c in cols]).value_counts(normalize=True)
    return dict(column=F.columns[j], cut=cut, good_side="above" if side == 1 else "below", balanced_accuracy_in_sample=score,
                balanced_accuracy_loo=loo, chance=0.5, null95=float(np.quantile(null, 0.95)),
                p=float((1 + (null >= loo).sum()) / (1 + n_perm)), column_share_in_folds=picked.round(3).to_dict())
