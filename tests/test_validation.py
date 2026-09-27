"""M3_modeler.validation against the controls the paper deposits.

Set DESCRIPYTOR_ARCHIVE to the Zenodo archive root; the tests skip without it.
- CS1: models/null_mazet_cs1_dipole.json (scripts/cs1/selection_null_cs1_dipole.py)
- CS3: the comparison with the published model asserted in scripts/cs3/cs3_search.py
"""
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from M3_modeler.ridge_search import adj_q2, all_combos, loo_predictions
from M3_modeler.validation import paired_sign_flip, selection_null

ARCHIVE = Path(os.environ["DESCRIPYTOR_ARCHIVE"]) if os.environ.get("DESCRIPYTOR_ARCHIVE") else None
pytestmark = pytest.mark.skipif(ARCHIVE is None or not ARCHIVE.is_dir(), reason="set DESCRIPYTOR_ARCHIVE")


def test_cs1_constrained_selection_null():
    d = pd.read_csv(ARCHIVE / "feature_matrices" / "mazet_features.csv")
    feats = [c for c in d.columns if c not in ("name", "output", "Amplitude_Stretch_3_8")]
    X, y = d[feats].to_numpy(float), d["output"].to_numpy(float)
    dip = feats.index("dipole_y_{23,24,25,26,27,28}-23-1")
    others = [i for i in range(len(feats)) if i != dip]
    combos = np.array([(dip, a, b) for a, b in all_combos(len(others), 2).tolist()])
    combos[:, 1:] = np.array(others)[combos[:, 1:]]
    out = selection_null(X, y, combos, lam=0, n_perm=1000, seed=0)
    want = json.loads((ARCHIVE / "models" / "null_mazet_cs1_dipole.json").read_text())
    assert len(combos) == want["n_constrained"]
    assert round(out["real"], 4) == want["reported_q2"]
    assert round(out["median"], 4) == want["null_median"]
    assert round(out["p95"], 4) == want["null_p95"]
    assert round(out["max"], 2) == want["null_max"]
    assert round(out["p"], 3) == want["p_value"]


def test_cs3_against_the_published_model():
    d = pd.read_csv(ARCHIVE / "feature_matrices" / "single_structure" / "cp.csv")
    drop = {"name", "smiles", "ligand_class", "ddg", "published_loo", "n_conformers", "mu_outofplane"}
    cols = [c for c in d.columns if c not in drop]
    X, y = d[cols].to_numpy(float), d.ddg.to_numpy(float)
    combos = all_combos(len(cols), 3)
    pred = loo_predictions(X, y, combos)
    ours = pred[adj_q2(pred, y, 3).argmax()]
    out = paired_sign_flip(y, ours, d.published_loo.to_numpy(float))
    assert out["closer"] == 14
    assert round(out["p"], 2) == 0.62
    assert round(out["residual_r"], 2) == 0.64
