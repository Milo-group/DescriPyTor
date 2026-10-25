"""Joint selection, model forms and design tools (M3_modeler.joint / forms / design) on data where the answer is known."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from M3_modeler.design import leave_best_out, suggest, top_k_recovery  # noqa: E402
from M3_modeler.forms import flexible_ceiling, logistic_link, threshold_analysis  # noqa: E402
from M3_modeler.joint import fit_joint, joint_selection, press_operators  # noqa: E402


def test_press_operator_equals_explicit_refits():
    rng = np.random.default_rng(0)
    X, y = rng.normal(size=(15, 4)), rng.normal(size=15)
    designs = [(0,), (1,), (2,), (3,), (0, 1), (2, 3)]
    M = press_operators(X, designs)
    for p, d in enumerate(designs):
        loo = []
        for i in range(15):
            k = np.arange(15) != i
            Z = np.c_[np.ones(14), X[k][:, d]]
            b = np.linalg.lstsq(Z, y[k], rcond=None)[0]
            loo.append(y[i] - np.r_[1, X[i, d]] @ b)
        assert np.allclose(M[p] @ y, loo)


def _reactions(seed=1, n=24, noise=0.6):
    """four reactions on the same compounds: two shared columns, weights (and a sign) that change by reaction"""
    rng = np.random.default_rng(seed)
    F = pd.DataFrame(rng.normal(size=(n, 6)), columns=[f"c{j}" for j in range(6)], index=[f"L{i}" for i in range(n)])
    W = {"I": (1.0, 0.8), "II": (1.2, 0.6), "III": (-0.9, 0.7), "IV": (-1.1, 0.9)}
    Y = {r: pd.Series(w[0] * F.c0 + w[1] * F.c3 + noise * rng.normal(size=n), index=F.index) for r, w in W.items()}
    Y["II"] = Y["II"].drop(["L0", "L5"])                          # a compound missing from one reaction
    return Y, F


def test_joint_selection_finds_the_shared_columns_and_their_signs():
    Y, F = _reactions()
    res = joint_selection(Y, F, max_terms=2, n_null=30)
    assert set(res["columns"]) == {"c0", "c3"} and res["pick_share"] > 0.8
    assert res["reactions"]["II"] == 22
    c = res["coefficients"]
    assert (np.sign(c.loc[["I", "II"], "c0"]) > 0).all() and (np.sign(c.loc[["III", "IV"], "c0"]) < 0).all()
    assert res["gain"] >= 0                                       # sharing the choice never hurts when the columns are shared


def test_joint_selection_refuses_one_reaction():
    Y, F = _reactions()
    with pytest.raises(ValueError):
        joint_selection({"I": Y["I"]}, F)


def test_fit_joint_on_the_true_columns():
    Y, F = _reactions(noise=0.01)
    c = fit_joint(Y, F, ["c0", "c3"])
    assert np.allclose(c.loc["III", "c0"], -0.9, atol=0.05) and (c.r2_in_sample > 0.99).all()


def test_logistic_link_helps_a_floor_heavy_yield():
    rng = np.random.default_rng(3)
    x = rng.normal(size=(40, 1))
    y = 100 / (1 + np.exp(-(2.5 * x[:, 0] - 2.0)))                  # most yields near the floor, a few high
    y = np.clip(y + rng.normal(scale=3, size=40), 0, 100)
    r = logistic_link(x, y, n_null=20)
    assert r["q2_logit"] > r["q2_ols"] and r["vs_ols"]["gain"] > 0.05


def test_flexible_ceiling_reads_linear_as_enough_on_a_line():
    rng = np.random.default_rng(4)
    X = rng.normal(size=(20, 2))
    y = 2 * X[:, 0] - X[:, 1] + 0.1 * rng.normal(size=20)
    r = flexible_ceiling(X, y)
    assert r["q2"]["linear"] > 0.95 and r["ceiling_gain"] <= 0 and r["reading"].startswith("linear is enough")


def test_threshold_analysis_finds_the_column_and_the_cut():
    rng = np.random.default_rng(5)
    F = pd.DataFrame(rng.normal(size=(30, 4)), columns=list("abcd"))
    lab = F["c"] > 0.2
    r = threshold_analysis(F, lab, n_perm=30)
    assert r["column"] == "c" and r["good_side"] == "above" and abs(r["cut"] - 0.2) < 0.4
    assert r["balanced_accuracy_loo"] > 0.85 and r["p"] < 0.05


def test_top_k_recovery_perfect_and_direction():
    y = np.arange(12, dtype=float)
    r = top_k_recovery(y, y)
    assert r["k"] == 3 and r["hits"] == 3 and np.isclose(r["p"], 1 / 220)          # 1 / C(12, 3)
    assert top_k_recovery(y, -y, higher_is_better=False)["hits"] == 0
    assert top_k_recovery(y, y * 0.5 + 3)["best_shrinkage_sd"] > 0                  # a shrunk model under-predicts the best


def test_leave_best_out_signal_and_noise():
    rng = np.random.default_rng(6)
    x = rng.normal(size=(24, 1))
    assert leave_best_out(x, 3 * x[:, 0] + 0.2 * rng.normal(size=24))["passes"]
    assert not leave_best_out(x, rng.normal(size=24))["passes"]


def test_suggest_ranks_up_the_slope_spreads_a_batch_and_flags_the_far():
    rng = np.random.default_rng(7)
    X = pd.DataFrame(rng.normal(size=(20, 2)), columns=["a", "b"])
    y = 2 * X.a + 0.1 * rng.normal(size=20)
    C = pd.DataFrame({"a": [0.0, 1.5, 1.6, 9.0], "b": [0.0, 0.0, 0.05, 0.0]}, index=["mid", "good", "good_twin", "far"])
    out = suggest(X, y, C, batch=3)
    assert out.loc["good", "expected_improvement"] > out.loc["mid", "expected_improvement"]
    assert bool(out.loc["far", "outside_domain"]) and "a" in out.loc["far", "outside_range"]
    assert sorted(out.pick[out.pick > 0]) == [1, 2, 3]
    inside = out.index[~out.outside_domain]
    assert out.index[0] in inside or not len(inside)                                # inside the domain is picked first
    assert suggest(X, y, C, domain="ignore").index[0] == "far"                     # expected improvement alone extrapolates
    assert (out.lower < out.prediction).all() and (out.prediction < out.upper).all()
