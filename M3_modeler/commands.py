"""
The `descripytor design`, `joint` and `forms` commands: read the CSVs, call design.py / joint.py / forms.py, print a short
report and write the tables. Tables are keyed by compound name (the first column, read as text: a compound called "NA"
stays "NA").
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from . import design, forms, joint


def _table(path):
    t = pd.read_csv(path, index_col=0, keep_default_na=False, na_values=[""])
    t.index = t.index.astype(str)
    if t.index.duplicated().any():
        raise SystemExit(f"{path}: these names appear twice: {sorted(set(t.index[t.index.duplicated()]))}")
    return t


def _columns(t, cols, y):
    cols = cols or [c for c in t.columns if c != y]
    missing = [c for c in cols if c not in t.columns]
    if missing:
        raise SystemExit(f"not in the table: {missing}. The columns are: {list(t.columns)}")
    X = t[cols].apply(pd.to_numeric, errors="coerce")
    bad = X.index[X.isna().any(axis=1)].tolist()
    if bad:
        raise SystemExit(f"non-numeric or empty cells for: {bad}")
    return X


def _ols_loo(X, y):
    n = len(y)
    return np.array([design.ols_fit(np.delete(X, i, 0), np.delete(y, i))(X[i:i + 1])[0] for i in range(n)])


def run_design(args):
    t = _table(args.table)
    if args.y not in t.columns:
        raise SystemExit(f"no column {args.y!r} in {args.table}")
    y = pd.to_numeric(t[args.y], errors="coerce")
    keep = y.notna()
    X = _columns(t.loc[keep], args.columns, args.y)
    y = y[keep].to_numpy(float)
    higher = not args.lower
    loo = _ols_loo(X.to_numpy(float), y)
    rec = design.top_k_recovery(y, loo, k=args.k, higher_is_better=higher)
    lbo = design.leave_best_out(X.to_numpy(float), y, k=args.k, higher_is_better=higher, n_null=args.null)
    q2 = float(1 - ((y - loo) ** 2).sum() / ((y - y.mean()) ** 2).sum())
    print(f"\n  OLS on {list(X.columns)}, n = {len(y)}, leave-one-out Q2 {q2:.3f}, {'higher' if higher else 'lower'} is better")
    print(f"  Top {rec['k']}: {rec['hits']} of the best {rec['k']} measured are among the best {rec['k']} predicted "
          f"(chance {rec['chance']:.1f}, p {rec['p']:.3f}); Spearman {rec['spearman']:.2f}; "
          f"the best are under-predicted by {rec['best_shrinkage_sd']:.2f} sd")
    print(f"  Leave the best {lbo['k']} out: S {lbo['S']:+.2f} against a null maximum of {lbo['null_max']:+.2f} "
          f"({'passes' if lbo['passes'] else 'does not pass'}, p {lbo['p']:.3f}); "
          f"{lbo['beyond_best_training']} predicted beyond the best training value")
    report = dict(columns=list(X.columns), n=len(y), q2_loo=q2, top_k=rec, leave_best_out=lbo)
    if args.candidates:
        C = _columns(_table(args.candidates), list(X.columns), None)
        out = design.suggest(X, y, C, higher_is_better=higher, batch=args.batch, domain="ignore" if args.allow_outside else "first")
        path = Path(args.output or Path(args.candidates).with_name(Path(args.candidates).stem + "_ranked.csv"))
        out.to_csv(path)
        print(f"\n  {len(C)} candidates ranked by expected improvement over the best measured value "
              f"({out.attrs['best_measured']:.3g}); {int(out.outside_domain.sum())} outside the measured domain. Picks:")
        print(out[out.pick > 0].round(3).to_string())
        print(f"  Wrote {path}. Read the predictions for rank, not size.")
        report["ranked"] = str(path)
    if args.json:
        Path(args.json).write_text(json.dumps(report, indent=1, default=float), encoding="utf-8")


def run_joint(args):
    F = _table(args.features)
    O = pd.read_csv(args.outcomes, keep_default_na=False, na_values=[""], dtype={args.name_col: str})
    for c in (args.name_col, args.reaction_col, args.y_col):
        if c not in O.columns:
            raise SystemExit(f"no column {c!r} in {args.outcomes}; it has {list(O.columns)}")
    O[args.y_col] = pd.to_numeric(O[args.y_col], errors="coerce")
    O = O.dropna(subset=[args.y_col])
    if O.duplicated([args.name_col, args.reaction_col]).any():
        raise SystemExit("a compound appears twice in one reaction; give one value per compound and reaction")
    Y = {str(r): g.set_index(args.name_col)[args.y_col].astype(float) for r, g in O.groupby(args.reaction_col)}
    cols = args.columns or [c for c in F.columns]
    F = _columns(F, cols, None)
    F = F.loc[:, F.std() > 1e-12]
    res = joint.joint_selection(Y, F, max_terms=args.max_terms, n_null=args.null)
    print(f"\n  {len(Y)} reactions, {res['n_compounds']} compounds, {res['n_designs']} candidate designs (1 to {args.max_terms} columns)")
    for r in Y:
        print(f"  {r}: separate selection Q2 {res['q2']['separate'][r]:+.3f}  ->  joint {res['q2']['joint'][r]:+.3f}")
    print(f"  Mean gain {res['gain']:+.3f}, p {res['p']:.4f}, null 95th {res['null95']:+.3f}: "
          f"{'joint selection helps' if res['passes'] else 'no detectable gain from selecting together'}")
    print(f"  Shared pick {res['columns']} in {res['pick_share']:.0%} of folds. Coefficients on all compounds:")
    print(res["coefficients"].round(3).to_string())
    if args.output:
        res["coefficients"].to_csv(args.output)
        print(f"  Wrote {args.output}")


def run_forms(args):
    t = _table(args.table)
    if args.y not in t.columns:
        raise SystemExit(f"no column {args.y!r} in {args.table}")
    y = pd.to_numeric(t[args.y], errors="coerce")
    keep = y.notna()
    X = _columns(t.loc[keep], args.columns, args.y)
    y = y[keep].to_numpy(float)
    line = lambda g: f"gain {g['gain']:+.3f}, p {g['p']:.3f}, null 95th {g['null95']:+.3f}: {'passes' if g['passes'] else 'does not pass'}"
    if args.logit:
        r = forms.logistic_link(X.to_numpy(float), y, scale=args.scale, n_null=args.null)
        print(f"\n  Logistic link on {list(X.columns)}: LOO Q2 OLS {r['q2_ols']:.3f} -> logit {r['q2_logit']:.3f} "
              f"(OLS clipped {r['q2_clipped']:.3f}); MAE {r['mae_ols']:.2f} -> {r['mae_logit']:.2f}")
        print(f"    against OLS: {line(r['vs_ols'])}\n    against clipped OLS: {line(r['vs_clipped'])}")
    if args.ceiling:
        r = forms.flexible_ceiling(X.to_numpy(float), y)
        print(f"\n  Ceiling on {X.shape[1]} columns: LOO Q2 linear ({r['reference']}) {r['q2']['linear']:.3f}, "
              f"Gaussian process {r['q2']['gp']:.3f}, random forest {r['q2']['forest']:.3f}\n    {r['reading']}")
    if args.threshold is not None:
        lab = (y >= args.threshold) if not args.lower else (y <= args.threshold)
        r = forms.threshold_analysis(X, lab, n_perm=args.null)
        print(f"\n  Threshold analysis ({'y >= ' if not args.lower else 'y <= '}{args.threshold:g} is good, {int(lab.sum())} of {len(lab)}): "
              f"{r['column']} {'>' if r['good_side'] == 'above' else '<='} {r['cut']:.4g}\n"
              f"    balanced accuracy in sample {r['balanced_accuracy_in_sample']:.2f}, leave-one-out {r['balanced_accuracy_loo']:.2f} "
              f"(chance 0.50, null 95th {r['null95']:.2f}, p {r['p']:.3f}); columns chosen in the folds: {r['column_share_in_folds']}")
    if not (args.logit or args.ceiling or args.threshold is not None):
        raise SystemExit("choose at least one of --logit, --ceiling, --threshold VALUE")
