"""Many tables at once: each CSV in a folder becomes a kept case of the explorer.

    <a Python with bassa_reg> batch.py <folder of csv> [--response NAME] [--structures FOLDER] [--checks]

Per table: a full run (three chains of 20,000), as from the page's New case. --response names the response column
(default: the last column). --structures is one folder for all tables, or a folder holding one sub-folder per table
(named like the CSV). --checks also computes the chance line, the nested Q2 and the fold frequencies, which takes a
hundred more runs per table: start it and leave it.
The cases are written beside the page and listed in manifest.js, exactly as Keep this case does.
"""
import argparse, sys, time
from pathlib import Path

import server as S

a = argparse.ArgumentParser(); a.add_argument("folder"); a.add_argument("--response"); a.add_argument("--structures", default=""); a.add_argument("--checks", action="store_true"); a = a.parse_args()
tables = sorted(Path(a.folder).glob("*.csv"))
if not tables: sys.exit("no .csv in " + a.folder)
for f in tables:
    t0 = time.time(); csv = f.read_text(encoding="utf-8-sig"); head = csv.splitlines()[0]; sub = Path(a.structures) / f.stem if a.structures else None
    try:
        import pandas as pd, io
        q = dict(csv=csv, response=a.response or pd.read_csv(io.StringIO(head)).columns[-1], iters=20000, name=f.stem, structures=str(sub if sub and sub.is_dir() else a.structures))
        d = S.run(q)
        if a.checks:
            X, y, _, _ = S.table(q); job = dict(done=0, total=0, chance=None, nested_q2=None, error=None); S.checks(job, X, y)
            if job["error"]: raise RuntimeError(job["error"])
            d.update({k: job[k] for k in ("chance", "nested_q2", "fold_top", "n_folds", "fold_freq", "p_top", "chance_lo", "chance_hi", "n_shuffles")})
        r = S.keep(dict(name=f.stem, case=dict(d, title=f.stem)))
        print("%-28s kept as %-20s %2d compounds, %3d descriptors, %d structures, %.0f s" % (f.name, r["name"], d["n"], len(d["names"]), len(d["mols"]), time.time() - t0), flush=True)
    except Exception as e: print("%-28s FAILED: %s" % (f.name, e), flush=True)
