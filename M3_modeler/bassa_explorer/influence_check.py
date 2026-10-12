"""Does a descriptor's effect rest on single compounds, or only the sampler's preference for it?

For a case of the explorer (its .js file) and one descriptor (default: the most sampled), from what the case already holds:
  - the share of models with the descriptor in every leave-one-compound-out rerun, and the compound that moves it most;
  - that compound among the others: its response and how extreme it is in the descriptors;
  - the descriptor alone against the response, and the most visited model that contains it, with and without that
    compound: slope, R2, leave-one-out Q2, Cook's distance.
A share that moves while slope and model do not is the sampler changing its choice among near-equivalent models.
Exploratory: it reads training compounds only and selects nothing.
    python influence_check.py cs3_cc [descriptor column]
"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent)); from settings import CASES_DIR as H, HERE
case = sys.argv[1]; t = (H / f"{case}.js" if (H / f"{case}.js").exists() else HERE / f"{case}.js").read_text(encoding="utf-8"); d = json.loads(t[t.index("] = ") + 4:-2])
X = pd.DataFrame(d["X"], index=d["ligands"]); y = pd.Series(d["y"])[d["ligands"]]; freq = pd.Series(d["freq"], d["names"]); lab = dict(zip(d["names"], d["label"]))
if not d.get("fold_freq"): sys.exit(f"{case} has no leave-one-out reruns; run the checks for it first")
F = pd.DataFrame(d["fold_freq"], index=d["ligands"], columns=d["names"]) if isinstance(d["fold_freq"], list) else pd.DataFrame(d["fold_freq"], index=d["ligands"])
col = sys.argv[2] if len(sys.argv) > 2 else freq.idxmax(); f = F[col]; who = (f - freq[col]).abs().idxmax(); rest = [i for i in y.index if i != who]; z = (X - X.mean()) / X.std(ddof=0)


def ols(cols, rows):
    A = np.c_[np.ones(len(rows)), z.loc[rows, cols]]; b = np.linalg.lstsq(A, y[rows], rcond=None)[0]; r = (y[rows] - A @ b).to_numpy(); h = np.einsum("ij,jk,ik->i", A, np.linalg.pinv(A.T @ A), A); sst = ((y[rows] - y[rows].mean()) ** 2).sum()
    cook = (r ** 2 / len(b) / (r @ r / (len(rows) - len(b)))) * h / (1 - h) ** 2
    return dict(b=b, r2=1 - r @ r / sst, q2=1 - ((r / (1 - h)) ** 2).sum() / sst, cook=pd.Series(cook, rows))


print(f"{d['title']}\n{lab[col]}  [{col}]: in {100 * freq[col]:.0f} % of models with all {len(y)} compounds")
print(f"1. Leaving one compound out: {100 * f.min():.0f} % to {100 * f.max():.0f} %. Most moved by {who}: {100 * f[who]:.0f} %. The five lowest: " + ", ".join(f"{k} {100 * v:.0f}%" for k, v in f.sort_values().head(5).items()))
print(f"2. {who}: response {y[who]:.2f} (rank {int(y.rank()[who])} of {len(y)}); z on this descriptor {z.loc[who, col]:+.2f}; most extreme compound in {int((z.abs().idxmax() == who).sum())} of {X.shape[1]} columns; largest |z| {z.loc[who].abs().max():.1f}")
a, b = ols([col], list(y.index)), ols([col], rest)
print(f"3. The descriptor alone (per sd): slope {a['b'][1]:+.2f}, R2 {a['r2']:.2f} with all; slope {b['b'][1]:+.2f}, R2 {b['r2']:.2f} without {who}. Cook's distance of {who}: {a['cook'][who]:.2f} (largest {a['cook'].max():.2f}, {a['cook'].idxmax()})")
m = next((M for M in d["models"] if d["names"].index(col) in M["cols"] and len(M["cols"]) > 1), None)
if m:
    cols = [d["names"][i] for i in m["cols"]]; k = cols.index(col) + 1; a, b = ols(cols, list(y.index)), ols(cols, rest)
    print(f"4. Most visited model with it ({' + '.join(lab[c] for c in cols)}; {100 * m['share']:.1f} % of the run):")
    print(f"   all: R2 {a['r2']:.2f}, Q2 {a['q2']:.2f}, its coefficient {a['b'][k]:+.2f}.  Without {who}: R2 {b['r2']:.2f}, Q2 {b['q2']:.2f}, coefficient {b['b'][k]:+.2f}.  Cook's of {who}: {a['cook'][who]:.2f} (largest {a['cook'].max():.2f}, {a['cook'].idxmax()})")
moved = abs(b["b"][k if m else 1] - a["b"][k if m else 1]) > 0.25 * abs(a["b"][k if m else 1])
print("Reading: the " + ("coefficient changes by more than a quarter without that compound: the effect itself leans on it." if moved else "coefficient holds without that compound; what moves is how often the sampler prefers this descriptor."))
if a["cook"][who] > 1: print(f"But {who} is a strongly influential point in that model (Cook's distance above 1): look at the compound itself before quoting the model.")
