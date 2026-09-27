"""M3_modeler.ridge_search must reproduce the paper's Case Study 3 numbers from the deposited pool.

Set DESCRIPYTOR_ARCHIVE to the Zenodo archive root (the folder holding feature_matrices/); the
test skips without it. The asserted values are the ones scripts/cs3/cs3_search.py and
cs3_constraint_null.py assert in the deposit.
"""
import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from M3_modeler.ridge_search import adj_q2, all_combos, constraint_null, loo_predictions, nested_loo

ARCHIVE = os.environ.get("DESCRIPYTOR_ARCHIVE")
CP = Path(ARCHIVE) / "feature_matrices" / "single_structure" / "cp.csv" if ARCHIVE else None
DROP = {"name", "smiles", "ligand_class", "ddg", "published_loo", "n_conformers", "mu_outofplane"}

pytestmark = pytest.mark.skipif(CP is None or not CP.is_file(), reason="set DESCRIPYTOR_ARCHIVE")


@pytest.fixture(scope="module")
def pool():
    d = pd.read_csv(CP)
    cols = [c for c in d.columns if c not in DROP]
    return cols, d[cols].to_numpy(float), d.ddg.to_numpy(float), all_combos(len(cols), 3)


def test_exhaustive_search_finds_the_reported_equation(pool):
    cols, X, y, combos = pool
    score = adj_q2(loo_predictions(X, y, combos), y, 3)
    assert len(cols) == 58 and len(combos) == 30856
    assert sorted(cols[j] for j in combos[score.argmax()]) == ["homo", "sub_B1_sym", "sub_theta_sym"]
    assert round(score.max(), 3) == 0.813


def test_nested_selection_is_stable(pool):
    cols, X, y, combos = pool
    out = nested_loo(X, y, combos)
    terms = ["homo", "sub_B1_sym", "sub_theta_sym"]
    assert sum(sorted(cols[j] for j in s) == terms for s in out["selected"]) == 30
    assert round(out["nested_r2"], 3) == 0.832


@pytest.mark.slow
def test_metal_angle_constraint_rank(pool):
    cols, X, y, combos = pool
    null = constraint_null(X, y, combos)
    v = null["nested_r2_for"]([cols.index("fromM_theta_sym"), cols.index("fromM_theta_asym")])
    assert int((null["nested_r2"] > v + 1e-9).sum()) == 384
    assert abs(null["nested_r2"].max() - 0.8321) < 1e-4
