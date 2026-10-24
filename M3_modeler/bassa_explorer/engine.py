"""The sampler set-up every run of the explorer uses: BASSA's spike-and-slab sampler with the corrected core.

`fit` and `usable` are taken unchanged from the DescriPyTor paper's engine.py (notes/analysis/bassa); the comparison
methods of that file (forward selection, Lasso, exhaustive search) are not needed here and are left out.
Needs `bassa_reg` (pip install "bassa-reg>=0.1.16,<0.2"; see TESTED below).
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
# The copy was made from that commit and checked with bassa-reg 0.1.15 and 0.1.16. Another release may change the core it
# replaces: a release outside TESTED is reported, and one whose core lacks a method the copy provides is refused, rather
# than mixing a new sampler with an old core and quietly giving other numbers than the paper's.
TESTED = ("0.1.15", "0.1.16")
if CORE == "fixed":
    import sys
    from importlib.metadata import PackageNotFoundError, version
    import bassa_reg.spike_and_slab.spike_and_slab as _SS, core_fixed
    try: BASSA_VERSION = version("bassa-reg")
    except PackageNotFoundError: BASSA_VERSION = "unknown"
    _own = {m for m in dir(core_fixed.SpikeAndSlabCore) if not m.startswith("_")}
    _theirs = {m for m in dir(_SS.SpikeAndSlabCore) if not m.startswith("_")}
    if not _theirs <= _own:
        raise RuntimeError("bassa-reg %s has a sampler core this explorer does not know (new: %s). Install a tested release: "
                           "pip install \"bassa-reg>=0.1.16,<0.2\", or run the released core with BASSA_CORE=package"
                           % (BASSA_VERSION, ", ".join(sorted(_theirs - _own))))
    if BASSA_VERSION not in TESTED:
        print("bassa-reg %s is not one of the releases the corrected core was checked with (%s); results may differ from the paper's"
              % (BASSA_VERSION, ", ".join(TESTED)), file=sys.stderr, flush=True)
    _SS.SpikeAndSlabCore = core_fixed.SpikeAndSlabCore

FLOOR = 0.25
# The prior on the inclusion rate is Beta(a, b), mean k/d (k terms expected out of d). "k3", the recipe and the default since 6 Oct 2026:
# a = 1, b = (d - k)/k, a weak prior, so small sets still average 4-5 terms. "strong": a = k, b = d - k, the same mean held more firmly.
# "n5": the recipe with k = n/5, one term per five compounds of the fit. BASSA_PRIOR picks one for every run; held() and fit() also
# take prior=. Compared in notes/analysis/prior_test (PREREG_prior_test_2026-10-21.md) before any default changes.
PRIOR = os.environ.get("BASSA_PRIOR", "k3")


def inclusion_prior(d, n, k=3, prior=None):
    p = prior or PRIOR
    if p == "k3": return 1.0, max(d - k, 1) / k
    if p == "strong": return float(k), float(max(d - k, 1))
    if p == "n5": kn = max(1.0, n / 5); return 1.0, max(d - kn, 1) / kn
    raise ValueError("unknown BASSA prior %r (k3, strong, n5)" % p)


def usable(x, full=None):
    """columns a fit can use: not constant and, inside a leave-one-out fold (given the table with all rows), not carried
    by the row left out: fold sd at least FLOOR x the sd on all rows. A column below it is noise without the held-out
    row, which then lands tens of sd away. Reads X only."""
    s = x.std(); keep = s > 1e-9
    if full is not None: keep &= s >= FLOOR * full[x.columns].std()
    return x.loc[:, keep]


def held(X, y, tr, te, iters=ITER, seed=7, k=3, sig=(3, 0.4), tau=(1, 1), prior=None):
    """mean and 90 % interval of BASSA's predictive draws for the rows te, the sampler run on the rows tr (fit's set-up). A draw:
    at a kept step, that step's model applied to the compound plus noise with that step's sigma (bassa_cases/ensemble_protocol.py)."""
    x = usable(X[tr], X); dd = x.shape[1]; yv = np.asarray(y, float); a, b = inclusion_prior(dd, int(np.sum(tr)), k, prior)
    pri = SpikeAndSlabPriors(a_alpha_prior=a, a_beta_prior=b, sigma_alpha_prior=sig[0], sigma_beta_prior=sig[1], tau_alpha_prior=tau[0], tau_beta_prior=tau[1])
    cfg = SpikeAndSlabConfigurations(sampler_iterations=iters, save_meta_data=False, save_samples=False, seed=seed, turn_on_convergence_tracking=False, dont_save_anything=True, verbose=False)
    np.random.seed(seed); reg = SpikeAndSlabRegression(x=x, y=y[tr], priors=pri, config=cfg, test_set=None, project_path=TMP, experiment_name="h"); reg.run()
    names = list(reg.feature_names); b = int(BURN * iters); S = (np.asarray(reg.all_s)[:, b:] > 0.5)[:, ::4]; W = np.asarray(reg.all_w, float)[:, b:][:, ::4]; sg = np.asarray(reg.all_sigma, float).reshape(-1)[b:][::4]
    m = min(S.shape[1], len(sg)); S, W, sg = S[:, :m], W[:, :m], sg[:m]
    xs = (X.loc[te, names].to_numpy() - x[names].mean().to_numpy()) / x[names].std(ddof=0).to_numpy(); ym, ys = yv[tr].mean(), yv[tr].std()
    M = xs @ (W * S); D = (M + sg * np.random.default_rng(seed).standard_normal(M.shape)) * ys + ym
    return M.mean(1) * ys + ym, np.quantile(D, 0.05, axis=1), np.quantile(D, 0.95, axis=1)


def fit(x, y, xt=None, seed=7, k=3, sig=(3, 0.4), iters=ITER, full=None, tau=(1, 1), prior=None):
    x = usable(x, full); d = x.shape[1]; a, b = inclusion_prior(d, len(y), k, prior)
    pri = SpikeAndSlabPriors(a_alpha_prior=a, a_beta_prior=b, sigma_alpha_prior=sig[0], sigma_beta_prior=sig[1], tau_alpha_prior=tau[0], tau_beta_prior=tau[1])
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
