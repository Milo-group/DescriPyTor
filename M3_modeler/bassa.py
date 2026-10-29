"""
BASSA (spike-and-slab) from this Python, with the explorer's sampler set-up.

bassa-reg pins NumPy 2 and DescriPyTor needs NumPy < 2, so the sampler lives in its own environment (made once by
`descripytor bassa_setup`). These functions send the table to that Python and bring the answer back:

    X, y = load_table("features.csv")       # names first; the response column found by name (or the last column)
    r = fit(X, y)                           # inclusion frequency, sign, coefficient per descriptor; the most visited models
    pngs = figures(X, y, "bassa_figures")   # BASSA's own plots: model arrangement, Markov chain, survival, convergence
    cv = cross_validate(X, y)               # stratified 5-fold with the sampler rerun in every fold: Q2, MAE, coverage

Set-up as engine.fit (corrected core, recipe prior with 3 expected terms, 20,000 steps, a fifth burned, seed 7).
The full protocol of a kept case (20 deals + LOO, chance line, drop-one) is `descripytor bassa_protocol`.
"""
from __future__ import annotations

import json
import os
import subprocess
import tempfile
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

__all__ = ["load_table", "fit", "figures", "cross_validate", "sampler_python"]

HERE = Path(__file__).resolve().parent / "bassa_explorer"

_WORKER = r"""
import json, sys
sys.path.insert(0, sys.argv[1])
import numpy as np, pandas as pd
import engine as E, protocol as P
job = json.load(open(sys.argv[2], encoding="utf-8"))
X = pd.DataFrame(job["X"], index=job["names"], columns=job["cols"], dtype=float)
y = pd.Series(job["y"], index=job["names"], dtype=float)
if job["what"] == "fit":
    r = E.fit(X, y, seed=job["seed"], iters=job["iters"])
    out = dict(freq=r["freq"].to_dict(), pos=r["pos"].to_dict(), coef=r["coef"].to_dict(), size=r["size"], top=r["top"])
elif job["what"] == "figures":
    from bassa_reg import Bassa
    from bassa_reg.spike_and_slab.spike_and_slab import SpikeAndSlabConfigurations, SpikeAndSlabRegression
    from bassa_reg.spike_and_slab.spike_and_slab_util_models import SpikeAndSlabPriors
    x = E.usable(X); a, b = E.inclusion_prior(x.shape[1], len(y))
    pri = SpikeAndSlabPriors(a_alpha_prior=a, a_beta_prior=b, sigma_alpha_prior=3, sigma_beta_prior=0.4, tau_alpha_prior=1, tau_beta_prior=1)
    cfg = SpikeAndSlabConfigurations(sampler_iterations=job["iters"], seed=job["seed"], turn_on_convergence_tracking=True, verbose=False)
    import os; os.chdir(job["out"])   # bassa_reg makes its run folder relative to the working folder
    np.random.seed(job["seed"])
    reg = SpikeAndSlabRegression(x=x, y=y, priors=pri, config=cfg, test_set=None, project_path=job["out"], experiment_name="bassa")
    reg.run(); Bassa(model=reg).run()
    out = dict(folder=job["out"] + "/" + reg.full_experiment_name)
else:
    yv = y.to_numpy(float); out = dict(deals=[])
    for deal in range(job["deals"]):
        fold = P.folds(yv, deal); p, lo, hi = np.zeros(len(yv)), np.zeros(len(yv)), np.zeros(len(yv))
        for f in range(5):
            te = fold == f
            p[te], lo[te], hi[te] = P.held(X, y, ~te, te, seed=job["seed"])
        out["deals"].append(dict(pred=p.tolist(), lo=lo.tolist(), hi=hi.tolist()))
json.dump(out, open(sys.argv[3], "w", encoding="utf-8"))
"""


def sampler_python(bassa_python=None) -> str:
    """The Python with bassa_reg, found as bassa_explorer finds it (argument, DESCRIPYTOR_BASSA_PYTHON, bassa_setup's)."""
    from descripytor.cli import _sampler_python
    return _sampler_python(bassa_python, remember=False)


def _run(job, X, y, bassa_python):
    X = X.astype(float)
    y = pd.Series(y, index=X.index, dtype=float) if not isinstance(y, pd.Series) else y.astype(float).loc[X.index]
    if X.isna().any().any() or y.isna().any():
        raise ValueError("missing values in X or y: drop or fill those rows first")
    job.update(X=X.to_numpy().tolist(), y=y.loc[X.index].tolist(), names=[str(i) for i in X.index], cols=[str(c) for c in X.columns])
    with tempfile.TemporaryDirectory() as d:
        src, dst = Path(d) / "job.json", Path(d) / "out.json"
        src.write_text(json.dumps(job), encoding="utf-8")
        r = subprocess.run([sampler_python(bassa_python), "-c", _WORKER, str(HERE), str(src), str(dst)], capture_output=True, text=True)
        if r.returncode:
            raise RuntimeError("the sampler failed (is bassa_reg installed? run  descripytor bassa_setup):\n" + r.stderr[-2000:])
        return json.loads(dst.read_text(encoding="utf-8")), y


RESPONSE_NAMES = ("output", "y", "ddg", "response", "target", "yield", "ee", "er")
FIGURES = ("bassa_plot", "markov_chain_visualization", "survival_plot", "sparsity_convergence", "log_likelihood_convergence")


def load_table(path, response: str | None = None):
    """A CSV of names (first column), descriptors and one response -> (X, y).

    The response is ``response`` if given, else the first column named output, y, ddG, response, target, yield, ee or er
    (any case), else the last column. Text columns and columns with missing values are dropped, and so are rows without a
    response; each drop is printed."""
    # sep=None sniffs , ; or tab (Excel writes ; in many locales); utf-8-sig drops a byte-order mark
    t = pd.read_csv(path, index_col=0, keep_default_na=False, na_values=[""], sep=None, engine="python", encoding="utf-8-sig")
    t.index = t.index.astype(str)
    if t.index.duplicated().any():
        raise ValueError(f"these names appear twice: {sorted(set(t.index[t.index.duplicated()]))}")
    if response is None:
        low = {str(c).lower(): c for c in t.columns}
        response = next((low[n] for n in RESPONSE_NAMES if n in low), t.columns[-1])
    y = pd.to_numeric(t[response], errors="coerce")
    X = t.drop(columns=response).apply(pd.to_numeric, errors="coerce")
    text = [c for c in X.columns if X[c].isna().all()]
    keep = y.notna()
    X, y = X.loc[keep].drop(columns=text), y[keep]
    gaps = [c for c in X.columns if X[c].isna().any()]
    X = X.drop(columns=gaps)
    print(f"response: {response!r}; {len(y)} compounds, {X.shape[1]} descriptors")
    for what, items in (("text columns dropped", text), ("columns with missing values dropped", gaps),
                        ("compounds without a response dropped", list(keep.index[~keep]))):
        if items:
            print(f"  {what}: {items}")
    return X, y


def figures(X: pd.DataFrame, y, out_dir, seed: int = 7, iters: int = 20000, bassa_python=None) -> dict:
    """BASSA's own figures (as in its README), from one sampler run with this set-up and convergence tracking: the model
    arrangement plot, the Markov chain of the top features, the survival plot, and the sparsity and log-likelihood
    convergence plots. Written as PNG (and HTML) under ``out_dir``; returns {figure name: PNG path} for those made, plus
    "folder". BASSA keeps only models with training R2 >= 0.5 in its plot; with none, that plot is not made.
    shortcut: Windows refuses BASSA's file paths beyond 260 characters; keep ``out_dir`` short (no long-path handling)."""
    out_dir = Path(out_dir).resolve()
    if os.name == "nt" and len(str(out_dir)) > 150:
        warnings.warn(f"{out_dir} is {len(str(out_dir))} characters long; BASSA's files beyond 260 fail on Windows. Use a shorter folder.")
    out_dir.mkdir(parents=True, exist_ok=True)
    out, _ = _run(dict(what="figures", seed=seed, iters=iters, out=str(out_dir)), X, y, bassa_python)
    folder = Path(out["folder"])
    found = {f: folder / f"{f}.png" for f in FIGURES if (folder / f"{f}.png").exists()}
    return dict(found, folder=folder)


def fit(X: pd.DataFrame, y, seed: int = 7, iters: int = 20000, bassa_python=None) -> dict:
    """One sampler run on all compounds. Returns ``table`` (per descriptor: inclusion frequency, share of positive weight
    when included, mean coefficient in response units per descriptor unit; sorted by frequency), ``models`` (the five most
    visited descriptor sets with their share of steps) and ``size`` (mean number of terms). In-sample: not a validation."""
    out, _ = _run(dict(what="fit", seed=seed, iters=iters), X, y, bassa_python)
    table = pd.DataFrame(dict(inclusion=out["freq"], positive=out["pos"], coef=out["coef"])).sort_values("inclusion", ascending=False)
    return dict(table=table, models=pd.DataFrame(out["top"], columns=["model", "share"]), size=out["size"])


def cross_validate(X: pd.DataFrame, y, deals: int = 1, seed: int = 7, bassa_python=None) -> dict:
    """Stratified 5-fold (compounds sorted by y, blocks of 5, folds shuffled with default_rng(100 + deal), as the
    protocol), the sampler rerun on every training fold. Returns ``q2`` and ``mae`` (mean over deals), ``coverage`` of the
    90 % intervals, and ``predictions`` (y, pred, lo, hi of the first deal). Each fold is one sampler run (~10-30 s)."""
    if len(X) < 10:
        raise ValueError(f"5-fold cross-validation needs at least 10 compounds, this table has {len(X)}")
    out, y = _run(dict(what="cv", seed=seed, deals=deals), X, y, bassa_python)
    yv = y.to_numpy(float)
    q2 = [1 - ((yv - np.array(d["pred"])) ** 2).sum() / ((yv - yv.mean()) ** 2).sum() for d in out["deals"]]
    mae = [np.abs(yv - np.array(d["pred"])).mean() for d in out["deals"]]
    cov = [((np.array(d["lo"]) <= yv) & (yv <= np.array(d["hi"]))).mean() for d in out["deals"]]
    first = out["deals"][0]
    pred = pd.DataFrame(dict(y=yv, pred=first["pred"], lo=first["lo"], hi=first["hi"]), index=y.index)
    return dict(q2=float(np.mean(q2)), mae=float(np.mean(mae)), coverage=float(np.mean(cov)), q2_per_deal=q2, predictions=pred)
