"""The sampler set-up every run of the explorer uses: BASSA's spike-and-slab sampler with the corrected core.

`fit` and `usable` are taken unchanged from the DescriPyTor paper's engine.py (notes/analysis/bassa); the comparison
methods of that file (forward selection, Lasso, exhaustive search) are not needed here and are left out.
Needs `bassa_reg` (pip install bassa-reg).
"""
import os, tempfile
from collections import Counter

import numpy as np
import pandas as pd

from bassa_reg.spike_and_slab.spike_and_slab import SpikeAndSlabConfigurations, SpikeAndSlabRegression
from bassa_reg.spike_and_slab.spike_and_slab_util_models import SpikeAndSlabPriors, TestSet

ITER, BURN, TMP = 20000, 0.2, os.path.join(tempfile.gettempdir(), "bassa_explorer_engine")
# The package sampler (github.com/Milo-group/BASSA, commit 661a164) is biased when descriptors are correlated:
# reconsidering a descriptor that is in the model, it reads the determinant ratio off the first row of the precision
# matrix, which is the candidate's row only when the candidate is the first active column. The default here is the
# corrected copy, core_fixed.py. BASSA_CORE=package runs the released core instead.
CORE = os.environ.get("BASSA_CORE", "fixed"); assert CORE in ("fixed", "package"), CORE
if CORE == "fixed":
    import bassa_reg.spike_and_slab.spike_and_slab as _SS, core_fixed
    _SS.SpikeAndSlabCore = core_fixed.SpikeAndSlabCore

FLOOR = 0.25


def usable(x, full=None):
    """columns a fit can use: not constant and, inside a leave-one-out fold (given the table with all rows), not carried
    by the row left out: fold sd at least FLOOR x the sd on all rows. A column below it is noise without the held-out
    row, which then lands tens of sd away. Reads X only."""
    s = x.std(); keep = s > 1e-9
    if full is not None: keep &= s >= FLOOR * full[x.columns].std()
    return x.loc[:, keep]


def held(X, y, tr, te, iters=ITER, seed=7, k=3, sig=(3, 0.4), tau=(1, 1)):
    """mean and 90 % interval of BASSA's predictive draws for the rows te, the sampler run on the rows tr (fit's set-up). A draw:
    at a kept step, that step's model applied to the compound plus noise with that step's sigma (bassa_cases/ensemble_protocol.py)."""
    x = usable(X[tr], X); dd = x.shape[1]; yv = np.asarray(y, float)
    pri = SpikeAndSlabPriors(a_alpha_prior=1, a_beta_prior=max(dd - k, 1) / k, sigma_alpha_prior=sig[0], sigma_beta_prior=sig[1], tau_alpha_prior=tau[0], tau_beta_prior=tau[1])
    cfg = SpikeAndSlabConfigurations(sampler_iterations=iters, save_meta_data=False, save_samples=False, seed=seed, turn_on_convergence_tracking=False, dont_save_anything=True, verbose=False)
    np.random.seed(seed); reg = SpikeAndSlabRegression(x=x, y=y[tr], priors=pri, config=cfg, test_set=None, project_path=TMP, experiment_name="h"); reg.run()
    names = list(reg.feature_names); b = int(BURN * iters); S = (np.asarray(reg.all_s)[:, b:] > 0.5)[:, ::4]; W = np.asarray(reg.all_w, float)[:, b:][:, ::4]; sg = np.asarray(reg.all_sigma, float).reshape(-1)[b:][::4]
    m = min(S.shape[1], len(sg)); S, W, sg = S[:, :m], W[:, :m], sg[:m]
    xs = (X.loc[te, names].to_numpy() - x[names].mean().to_numpy()) / x[names].std(ddof=0).to_numpy(); ym, ys = yv[tr].mean(), yv[tr].std()
    M = xs @ (W * S); D = (M + sg * np.random.default_rng(seed).standard_normal(M.shape)) * ys + ym
    return M.mean(1) * ys + ym, np.quantile(D, 0.05, axis=1), np.quantile(D, 0.95, axis=1)


def fit(x, y, xt=None, seed=7, k=3, sig=(3, 0.4), iters=ITER, full=None, tau=(1, 1)):
    x = usable(x, full); d = x.shape[1]
    pri = SpikeAndSlabPriors(a_alpha_prior=1, a_beta_prior=max(d - k, 1) / k, sigma_alpha_prior=sig[0], sigma_beta_prior=sig[1], tau_alpha_prior=tau[0], tau_beta_prior=tau[1])
    cfg = SpikeAndSlabConfigurations(sampler_iterations=iters, save_meta_data=False, save_samples=False, seed=seed, turn_on_convergence_tracking=False,
                                     dont_save_anything=True, verbose=False)
    np.random.seed(seed)
    ts = TestSet(x_test=xt[x.columns], samples_per_y=100, iterations=2000) if xt is not None else None
    reg = SpikeAndSlabRegression(x=x, y=y, priors=pri, config=cfg, test_set=ts, project_path=TMP, experiment_name="e")
    reg.run()
    names = list(reg.feature_names); b = int(BURN * iters)
    s = np.asarray(reg.all_s)[:, b:] > 0.5; w = np.asarray(reg.all_w)[:, b:]
    n_in = np.maximum(s.sum(1), 1)
    unit = y.std(ddof=0) / x[names].std(ddof=0).to_numpy()                     # scaled weight -> response units per descriptor unit
    sets = Counter(tuple(np.flatnonzero(c)) for c in s.T)
    return dict(freq=pd.Series(s.mean(1), names), pos=pd.Series(((w > 0) & s).sum(1) / n_in, names), coef=pd.Series((w * s).sum(1) / n_in * unit, names),
                size=float(s.sum(0).mean()), top=[(" + ".join(names[i] for i in t) or "(empty)", n / s.shape[1]) for t, n in sets.most_common(5)],
                pred=None if ts is None else np.asarray(ts.y_pred, float), s=pd.DataFrame(s, index=names))
