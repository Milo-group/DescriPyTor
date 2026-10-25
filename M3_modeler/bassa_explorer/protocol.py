"""The validation protocol of a kept case, and the checks on it (from the paper's bassa_cases scripts, 13–16 October 2026).

    <a Python with bassa_reg> protocol.py protocol <case> [--deals 20] [--no-loo]  ->  protocol_<case>.json, and the case's `cv`
    <a Python with bassa_reg> protocol.py chance   <case> [--shuffles 20]          ->  chance_<case>.json, and cv.chance
    <a Python with bassa_reg> protocol.py seed     <case> [--seed 11]              ->  seed_<case>.json
    <a Python with bassa_reg> protocol.py drop_one <case> [<case> ...]             ->  dropone_<case>.json
    descripytor bassa_protocol ...   runs the same with the sampler's Python (the one bassa_explorer uses)

<case> is a case name in the cases folder (BASSA_EXPLORER_CASES, as the explorer uses) or a path to its .js file.

protocol: the sampler is rerun in every fold; nothing chosen on all compounds is scored. Stratified 5-fold over ``deals`` deals
  (compounds sorted by the response, cut into blocks of 5, fold labels shuffled inside each block with default_rng(100 + deal)),
  and leave-one-out as one more row. Q2 = 1 - SSE/SST. MAE beside Q2, the worst compound named and both numbers again without it.
  For every held-out compound the 90 % interval of the predictive draws (a kept step's model on the compound plus noise with that
  step's sigma) and the share of compounds it holds. The prediction is the mean over kept steps, checked on the first fold against
  BASSA's own test-set route (``first_fold_check``; above 5 % of the response sd, do not quote the interval numbers).
chance: the response permuted (default_rng(500 + s)), one deal each, the sampler rerun in every fold; p = (1 + shuffles at or above
  the real Q2) / (shuffles + 1), and how many real deals beat every shuffle. seed: the deals again with another sampler seed.
drop_one: BASSA refitted without each compound. A compound CARRIES the model when, without it, (a) the lead term's inclusion falls
  below 0.5, or (b) the most visited model changes and the lead term leaves it, or (c) the averaged model's predictions for the others
  shift by more than 0.25 sd and more than 3 x the sampler's own scatter (the full fit with seed 11). Cook's distance beside.
Sampler set-up everywhere: engine.fit's (recipe prior with 3 expected terms, corrected core, 20,000 steps, a fifth dropped, seed 7,
the fold screen engine.usable with the full table).
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import engine as E  # noqa: E402
import settings  # noqa: E402

K = 5


def case_path(case):
    p = Path(case)
    if p.suffix == ".js" and p.exists():
        return p
    q = Path(settings.CASES_DIR) / f"{case}.js"
    if not q.exists():
        raise SystemExit(f"no case {case!r}: neither a .js file nor {q}")
    return q


def load(case):
    """(path, text, start of the JSON, the case dict, X, y)"""
    f = case_path(case)
    t = f.read_text(encoding="utf-8")
    k0 = t.index("] = ") + 4
    d = json.loads(t[k0:len(t.rstrip().rstrip(";"))])
    X = pd.DataFrame({c: d["X"][c] for c in d["names"]}, index=d["ligands"]).astype(float)
    y = pd.Series([d["y"][l] for l in d["ligands"]], index=d["ligands"], dtype=float)
    return f, t, k0, d, X, y


def _save_case(f, t, k0, d):
    f.write_text(t[:k0] + json.dumps(d) + ";\n", encoding="utf-8")


def held(X, y, tr, te, seed=7, k=3, sig=(3, 0.4), tau=(1, 1), draws=True):
    """Mean (and, with ``draws``, the 90 % interval of the predictive draws) for the rows ``te``, from the sampler run on ``tr``."""
    from bassa_reg.spike_and_slab.spike_and_slab import SpikeAndSlabConfigurations, SpikeAndSlabRegression
    from bassa_reg.spike_and_slab.spike_and_slab_util_models import SpikeAndSlabPriors
    yv = y.to_numpy(float)
    x = E.usable(X[tr], X)
    dd = x.shape[1]
    pri = SpikeAndSlabPriors(a_alpha_prior=1, a_beta_prior=(dd - k) / k, sigma_alpha_prior=sig[0], sigma_beta_prior=sig[1], tau_alpha_prior=tau[0], tau_beta_prior=tau[1])
    cfg = SpikeAndSlabConfigurations(sampler_iterations=E.ITER, save_meta_data=False, save_samples=False, seed=seed, turn_on_convergence_tracking=False,
                                     dont_save_anything=True, verbose=False)
    np.random.seed(seed)
    reg = SpikeAndSlabRegression(x=x, y=y[tr], priors=pri, config=cfg, test_set=None, project_path=E.TMP, experiment_name="p")
    reg.run()
    names = list(reg.feature_names)
    b = int(E.BURN * E.ITER)
    S = (np.asarray(reg.all_s)[:, b:] > 0.5)[:, ::4]
    W = np.asarray(reg.all_w, float)[:, b:][:, ::4]
    sg = np.asarray(reg.all_sigma, float).reshape(-1)[b:][::4]
    m = min(S.shape[1], W.shape[1], len(sg) if draws else W.shape[1])
    S, W = S[:, :m], W[:, :m]
    xs = (X.loc[te, names].to_numpy() - x[names].mean().to_numpy()) / x[names].std(ddof=0).to_numpy()
    ym, ys = yv[tr].mean(), yv[tr].std()
    M = xs @ (W * S)
    if not draws:
        return M.mean(1) * ys + ym
    rng = np.random.default_rng(seed)
    D = (M + sg[:m] * rng.standard_normal(M.shape)) * ys + ym
    return M.mean(1) * ys + ym, np.quantile(D, 0.05, axis=1), np.quantile(D, 0.95, axis=1)


def folds(yv, deal):
    order = np.argsort(yv)
    rng = np.random.default_rng(100 + deal)
    fold = np.zeros(len(yv), int)
    for a in range(0, len(yv), K):
        blk = order[a:a + K]
        fold[blk] = rng.permutation(K)[:len(blk)]
    return fold


def _q2(yv, p, keep=slice(None)):
    return float(1 - ((yv[keep] - p[keep]) ** 2).sum() / ((yv[keep] - yv[keep].mean()) ** 2).sum())


def _mae(yv, p, keep=slice(None)):
    return float(np.abs(yv[keep] - p[keep]).mean())


def protocol(case, deals=20, loo=True, write=True, log=print):
    f, t, k0, d, X, y = load(case)
    name = f.stem
    yv = y.to_numpy(float)
    n = len(yv)
    t0 = time.time()
    rows, P, LO, HI = [], np.zeros((deals, n)), np.zeros((deals, n)), np.zeros((deals, n))
    FF = None
    for deal in range(deals):
        fold = folds(yv, deal)
        for q in range(K):
            te = fold == q
            if not te.any():
                continue
            P[deal, te], LO[deal, te], HI[deal, te] = held(X, y, ~te, te)
            if deal == 0 and q == 0:                     # the draws built here against BASSA's own test-set prediction
                own = np.asarray(E.fit(X[~te], y[~te], xt=X[te], full=X)["pred"], float).reshape(-1)
                FF = dict(difference=float(np.abs(own - P[0, te]).max()), sd=float(yv.std()))
                FF["ok"] = bool(FF["difference"] < 0.05 * FF["sd"])
        rows.append(dict(deal=deal, q2=_q2(yv, P[deal]), mae=_mae(yv, P[deal]), cover=float(((yv >= LO[deal]) & (yv <= HI[deal])).mean())))
        log("deal %2d  Q2 %.2f  MAE %.2f  inside the 90 %% interval %.0f %%  (%.0f s)" % (deal, rows[-1]["q2"], rows[-1]["mae"], 100 * rows[-1]["cover"], time.time() - t0))
    R = pd.DataFrame(rows)
    pm = P.mean(0)
    err = np.abs(yv - P).mean(0)
    w = int(np.argmax(err))
    keep = np.arange(n) != w
    out = dict(case=name, n=n, deals=deals, K=K, first_fold_check=FF, per_deal=rows, pred_per_deal=[[round(float(v), 4) for v in P[k]] for k in range(deals)],
               k5=dict(q2=float(R.q2.mean()), lo=float(R.q2.min()), hi=float(R.q2.max()), sd=float(R.q2.std(ddof=0)), mae=float(R.mae.mean()),
                       mae_lo=float(R.mae.min()), mae_hi=float(R.mae.max()), cover=float(R.cover.mean()), worst=d["ligands"][w], worst_err=float(err[w]),
                       q2_without=float(np.mean([_q2(yv, P[k], keep) for k in range(deals)])), mae_without=float(np.mean([_mae(yv, P[k], keep) for k in range(deals)]))))
    comp = dict(name=d["ligands"], y=yv, pred=pm, err=err, inside=((yv >= LO) & (yv <= HI)).mean(0), lo=LO.mean(0), hi=HI.mean(0))
    if loo:
        L, Ll, Lh = np.zeros(n), np.zeros(n), np.zeros(n)
        for i in range(n):
            te = np.arange(n) == i
            L[i], Ll[i], Lh[i] = held(X, y, ~te, te)
        wl = int(np.argmax(np.abs(yv - L)))
        kl = np.arange(n) != wl
        out["loo"] = dict(q2=_q2(yv, L), mae=_mae(yv, L), cover=float(((yv >= Ll) & (yv <= Lh)).mean()), worst=d["ligands"][wl], worst_err=float(abs(yv[wl] - L[wl])),
                          q2_without=_q2(yv, L, kl), mae_without=_mae(yv, L, kl))
        comp.update(loo=L, loo_lo=Ll, loo_hi=Lh)
    out["compounds"] = [{k: (v[i] if k == "name" else float(v[i])) for k, v in comp.items()} for i in range(n)]
    out["seconds"] = round(time.time() - t0)
    if write:
        (f.parent / f"protocol_{name}.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
        cv = d.get("cv", {})
        cv["k5"] = dict(cv.get("k5", {}), **out["k5"])
        cv.update(pred5=[round(float(v), 4) for v in pm], lo5=[round(float(v), 4) for v in LO.mean(0)], hi5=[round(float(v), 4) for v in HI.mean(0)],
                  deals=[round(float(r["q2"]), 3) for r in rows], split=f"stratified by the response, {deals} deals", d_ok=bool(FF["ok"]) if FF else None)
        if loo:
            cv["loo"] = out["loo"]
        d["cv"] = cv
        _save_case(f, t, k0, d)
    log("\n5-fold over %d deals: Q2 %.2f (%.2f to %.2f), MAE %.2f, inside the 90 %% interval %.0f %%; worst %s (%.2f), without it Q2 %.2f"
        % (deals, out["k5"]["q2"], out["k5"]["lo"], out["k5"]["hi"], out["k5"]["mae"], 100 * out["k5"]["cover"], out["k5"]["worst"], out["k5"]["worst_err"], out["k5"]["q2_without"]))
    if loo:
        log("leave-one-out: Q2 %.2f, MAE %.2f; worst %s (%.2f), without it Q2 %.2f" % (out["loo"]["q2"], out["loo"]["mae"], out["loo"]["worst"], out["loo"]["worst_err"], out["loo"]["q2_without"]))
    if FF and not FF["ok"]:
        log("first-fold check: the draws differ from BASSA's own prediction by %.3f (> 5 %% of the sd %.2f); rerun with another seed and do not quote the interval numbers" % (FF["difference"], FF["sd"]))
    return out


def _deal_q2(X, y, yv, k, seed):
    fold = folds(yv, k)
    p = np.zeros(len(yv))
    ys = pd.Series(yv, index=X.index)
    for q in range(K):
        te = fold == q
        if te.any():
            p[te] = held(X, ys, ~te, te, seed=seed, draws=False)
    return _q2(yv, p), _mae(yv, p)


def chance(case, shuffles=20, write=True, log=print):
    f, t, k0, d, X, y = load(case)
    name = f.stem
    real = json.loads((f.parent / f"protocol_{name}.json").read_text(encoding="utf-8"))
    y0 = y.to_numpy(float)
    rows, t0 = [], time.time()
    for s in range(shuffles):
        q, m = _deal_q2(X, y, y0[np.random.default_rng(500 + s).permutation(len(y0))], 0, 7)
        rows.append(dict(shuffle=s, q2=q, mae=m))
        log("shuffle %2d  Q2 %.2f  MAE %.2f  (%.0f s)" % (s, q, m, time.time() - t0))
    Q, M = np.array([r["q2"] for r in rows]), np.array([r["mae"] for r in rows])
    deals = [r["q2"] for r in real["per_deal"]]
    rq, rm = real["k5"]["q2"], real["k5"]["mae"]
    out = dict(case=name, shuffles=shuffles, real_q2=rq, real_mae=rm, q2_median=float(np.median(Q)), q2_max=float(Q.max()), mae_median=float(np.median(M)),
               mae_min=float(M.min()), p_q2=float((1 + (Q >= rq).sum()) / (shuffles + 1)), p_mae=float((1 + (M <= rm).sum()) / (shuffles + 1)),
               real_deals_above_max=int(sum(v > Q.max() for v in deals)), rows=rows, seconds=round(time.time() - t0))
    if write:
        (f.parent / f"chance_{name}.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
        d.setdefault("cv", {})["chance"] = dict(q2=[round(r["q2"], 3) for r in rows], p=out["p_q2"], above=out["real_deals_above_max"], deals=len(deals))
        _save_case(f, t, k0, d)
    log("\n%s: real Q2 %.2f | shuffled median %.2f, highest %.2f | p %.3f | real deals above every shuffle: %d of %d"
        % (name, rq, out["q2_median"], out["q2_max"], out["p_q2"], out["real_deals_above_max"], len(deals)))
    return out


def seed_check(case, seed=11, write=True, log=print):
    f, t, k0, d, X, y = load(case)
    name = f.stem
    real = json.loads((f.parent / f"protocol_{name}.json").read_text(encoding="utf-8"))
    y0 = y.to_numpy(float)
    rows, t0 = [], time.time()
    for k in range(real["deals"]):
        q, m = _deal_q2(X, y, y0, k, seed)
        rows.append(dict(deal=k, q2=q, mae=m))
        log("deal %2d  Q2 %.2f (seed 7: %.2f)  MAE %.2f  (%.0f s)" % (k, q, real["per_deal"][k]["q2"], m, time.time() - t0))
    Q, Q7 = np.array([r["q2"] for r in rows]), np.array([r["q2"] for r in real["per_deal"]])
    out = dict(case=name, seed=seed, q2=float(Q.mean()), q2_seed7=float(Q7.mean()), mae=float(np.mean([r["mae"] for r in rows])), mae_seed7=real["k5"]["mae"],
               largest_deal_difference=float(np.abs(Q - Q7).max()), rows=rows, seconds=round(time.time() - t0))
    if write:
        (f.parent / f"seed_{name}.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    log("\n%s: mean Q2 %.3f with seed %d, %.3f with seed 7; largest difference in one deal %.3f" % (name, out["q2"], seed, out["q2_seed7"], out["largest_deal_difference"]))
    return out


def _fit_rows(X, y, rows, seed=7):
    r = E.fit(X.loc[rows], y.loc[rows], seed=seed, full=X)
    fr = r["freq"]
    b = r["coef"] * fr                                   # marginal slope per descriptor unit
    mu = X.loc[rows, fr.index].mean()
    return fr, r["top"][0], y.loc[rows].mean() + (X[fr.index] - mu) @ b


def _cook(X, y, terms):
    A = np.c_[np.ones(len(y)), X[terms].to_numpy()]
    yv = y.to_numpy()
    Hm = A @ np.linalg.pinv(A.T @ A) @ A.T
    h = np.diag(Hm)
    e = yv - Hm @ yv
    p = A.shape[1]
    s2 = (e ** 2).sum() / (len(y) - p)
    return pd.Series(e ** 2 / (p * s2) * h / (1 - h) ** 2, y.index)


def drop_one(case, write=True, log=print):
    f, t, k0, d, X, y = load(case)
    name = f.stem
    n, sd, t0, all_ = len(y), y.std(ddof=0), time.time(), list(y.index)
    f0, top0, p0 = _fit_rows(X, y, all_)
    lead = f0.idxmax()
    _, _, pn = _fit_rows(X, y, all_, seed=11)
    noise = float(np.sqrt(((pn - p0) ** 2).mean()) / sd)
    terms = top0[0].split(" + ") if top0[0] != "(empty)" else []
    D = _cook(X, y, terms) if terms else pd.Series(0.0, y.index)
    rows = []
    for c in all_:
        rest = [k for k in all_ if k != c]
        fr, top, p = _fit_rows(X, y, rest)
        shift = float(np.sqrt(((p[rest] - p0[rest]) ** 2).mean()) / sd)
        lin = float(fr.get(lead, 0.0))
        why = dict(a=lin < 0.5, b=(top[0] != top0[0]) and lead not in top[0].split(" + "), c=shift > 0.25 and shift > 3 * noise)
        rows.append(dict(compound=c, y=float(y[c]), lead_incl=lin, top=top[0], top_share=float(top[1]), same_top=top[0] == top0[0], shift_sd=shift,
                         cook=float(D[c]), carries=any(why.values()), why="".join(k for k, v in why.items() if v),
                         top3=", ".join(f"{k} {v:.2f}" for k, v in fr.sort_values(ascending=False).head(3).items())))
        log(f"  {c:>14s}  lead {lin:.2f}  top {'same' if top[0] == top0[0] else top[0]:40.40s}  shift {shift:.2f} sd  Cook {D[c]:.2f}  "
            f"{'CARRIES ' + rows[-1]['why'] if rows[-1]['carries'] else ''}")
    R = pd.DataFrame(rows)
    out = dict(case=name, n=n, lead=lead, lead_incl=float(f0[lead]), top=top0[0], top_share=float(top0[1]), noise_sd=noise,
               top3_all=", ".join(f"{k} {v:.2f}" for k, v in f0.sort_values(ascending=False).head(3).items()), carriers=R.loc[R.carries, "compound"].tolist(),
               most_shift=R.loc[R.shift_sd.idxmax(), "compound"], max_shift=float(R.shift_sd.max()), lead_min=float(R.lead_incl.min()),
               lead_min_without=R.loc[R.lead_incl.idxmin(), "compound"], same_top_share=float(R.same_top.mean()), cook_gt1=R.loc[R.cook > 1, "compound"].tolist(),
               cook_max=R.loc[R.cook.idxmax(), "compound"], rows=rows, seconds=round(time.time() - t0))
    if write:
        (f.parent / f"dropone_{name}.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    log(f"{name}: lead {lead} ({f0[lead]:.2f}); noise {noise:.3f} sd; lead lowest {out['lead_min']:.2f} without {out['lead_min_without']}; "
        f"top unchanged in {out['same_top_share']:.0%}; largest shift {out['max_shift']:.2f} sd ({out['most_shift']}); carriers {out['carriers']}")
    return out


def main(argv=None):
    a = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    a.add_argument("check", choices=["protocol", "chance", "seed", "drop_one"])
    a.add_argument("cases", nargs="+", help="case names in the cases folder, or .js paths")
    a.add_argument("--deals", type=int, default=20)
    a.add_argument("--no-loo", action="store_true", help="protocol: skip leave-one-out (n more sampler runs)")
    a.add_argument("--shuffles", type=int, default=20)
    a.add_argument("--seed", type=int, default=11)
    a.add_argument("--dry", action="store_true", help="write nothing: print only")
    g = a.parse_args(argv)
    for c in g.cases:
        if g.check == "protocol":
            protocol(c, deals=g.deals, loo=not g.no_loo, write=not g.dry)
        elif g.check == "chance":
            chance(c, shuffles=g.shuffles, write=not g.dry)
        elif g.check == "seed":
            seed_check(c, seed=g.seed, write=not g.dry)
        else:
            drop_one(c, write=not g.dry)


if __name__ == "__main__":
    main()
