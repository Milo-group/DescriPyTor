"""Pool audit (M3_modeler.audit) and the further controls of M3_modeler.validation on data where the answer is known."""
import itertools
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from M3_modeler.audit import icc1, icc_bootstrap, reduce_pool, reliability_gate, zaudit  # noqa: E402
from M3_modeler.ridge_search import loo_predictions, q2  # noqa: E402
from M3_modeler.validation import (calibration, domain_flags, equivalence, family_null, nn_distance,  # noqa: E402
                                   pooled_equivalence, rank_check, slot_null)


def test_zaudit_tells_a_degenerate_column_from_an_outlier():
    rng = np.random.default_rng(0)
    X = pd.DataFrame(rng.normal(size=(25, 3)), columns=list("abc"))
    X["asym"] = 1e-4 * rng.normal(size=25)
    X.loc[7, "asym"] = 2.0                                   # one compound breaks the symmetry: all the variance
    X.loc[3, "b"] = 12.0                                     # a real outlier on a column with spread of its own
    Z = zaudit(X).set_index("column")
    assert set(Z.index) == {"asym", "b"}
    assert Z.loc["asym", "kind"] == "degenerate" and Z.loc["asym", "row"] == 7
    assert Z.loc["b", "kind"] == "outlier" and Z.loc["b", "row"] == 3


def test_reduce_pool_keeps_one_per_cluster_and_respects_tier():
    rng = np.random.default_rng(1)
    base = rng.normal(size=(30, 3))
    X = pd.DataFrame({"a": base[:, 0], "a_copy": 2 * base[:, 0] + 1, "a_new": base[:, 0] + 0.01 * rng.normal(size=30),
                      "b": base[:, 1], "c": base[:, 2], "flag": np.r_[np.zeros(15), np.ones(15)], "vib_Stretch": rng.normal(size=30)})
    kept, cl = reduce_pool(X, 3, tier=lambda c: 0 if c.endswith("_new") else 1, drop="Stretch")
    assert "a_copy" not in cl.column.values                  # exact duplicate after z-scoring
    assert "flag" not in cl.column.values and "vib_Stretch" not in cl.column.values
    assert "a_new" in kept.columns and "a" not in kept.columns   # the tier-0 member wins its cluster
    assert len(kept.columns) <= 3 and cl[cl.kept == 1].column.tolist() == list(kept.columns)
    assert cl.groupby("cluster").kept.sum().eq(1).all()


def test_icc1_limits_and_gate():
    r = np.random.default_rng(1)
    assert abs(icc1([np.full(3, float(i)) + 1e-12 * r.standard_normal(3) for i in range(10)]) - 1) < 1e-6
    assert abs(icc1([r.standard_normal(5) for _ in range(400)])) < 0.05
    b = icc_bootstrap([np.full(4, float(i)) + 0.1 * r.standard_normal(4) for i in range(12)], n_boot=200)
    assert b["lo"] <= b["icc"] <= b["hi"] and b["icc"] > 0.9
    g = reliability_gate(["theta", "B1", "dipole"], {"theta": 0.3, "B1": 0.9})
    assert list(g["drop"]) == [True, False, False]


def test_equivalence_verdicts():
    rng = np.random.default_rng(0)
    y = rng.normal(size=40)
    good, bad = y + 0.1 * rng.normal(size=40), y + 1.0 * rng.normal(size=40)
    assert equivalence(y, good, bad)[0]["verdict"] == "different, A better"
    near = good + 0.001 * rng.normal(size=40)
    assert equivalence(y, good, near, margins=(0.05,))[0]["verdict"].startswith("equivalent within 0.05")
    small = rng.normal(size=8)
    assert equivalence(small, small + rng.normal(size=8), small + rng.normal(size=8))[0]["verdict"] == "underpowered"
    one, pooled = equivalence(y, good, bad)[0], pooled_equivalence([(y, good, bad)])[0]
    assert np.isclose(one["dq2"], pooled["dq2"]) and np.isclose(one["ci90_lo"], pooled["ci90_lo"])


def test_calibration_counts_coverage():
    y = np.arange(8, dtype=float)
    comp = [dict(name=f"m{i}", y=v, pred=v - 0.5, lo=v - 1, hi=v + 1, loo_lo=v - 1 if i else v + 1, loo_hi=v + 2,
                 inside=0.0 if i == 0 else 1.0) for i, v in enumerate(y)]
    P = dict(compounds=comp, per_deal=[dict(cover=0.875), dict(cover=1.0)], k5=dict(q2=0.9, mae=0.5))
    r = calibration(P)
    assert r["cover_loo"] == 7 / 8 and r["verdict"] == "calibrated" and r["never_inside"] == ["m0"]
    assert np.isclose(r["width_sd"], 2 / y.std(ddof=1)) and np.isclose(r["best_quarter_under_sd"], 0.5 / y.std(ddof=1))


def test_domain_flag_finds_the_far_compound_and_ranks_error():
    rng = np.random.default_rng(3)
    X = rng.normal(size=(20, 3))
    X[5] = [8, 8, 8]
    d, flag = domain_flags(X)
    assert flag[5] and flag.sum() <= 3
    assert np.isclose(d[5], nn_distance(np.delete(X, 5, 0), X[5:6])[0])
    err = d + 0.05 * rng.normal(size=20)
    rc = rank_check(d, err, n_perm=500)
    assert rc["rho"] > 0.8 and rc["p"] < 0.01


def test_family_null_signal_and_noise():
    rng = np.random.default_rng(4)
    X = pd.DataFrame(rng.normal(size=(25, 12)), columns=[f"c{i}" for i in range(12)])
    y = X.c0 - X.c1 + 0.2 * rng.normal(size=25)
    r = family_null(X, y, ["c0", "c1"], n=200)
    assert r["percentile"] > 99 and r["p"] < 0.01
    assert np.isclose(r["real"], q2(loo_predictions(X.to_numpy(), y, np.array([[0, 1]]), lam=0), y)[0])
    assert family_null(X, rng.normal(size=25), ["c0", "c1"], n=200)["p"] > 0.01


def test_slot_null_matches_brute_force():
    rng = np.random.default_rng(5)
    X = pd.DataFrame(rng.normal(size=(15, 6)), columns=["m_angle", "m_len", "m_q", "s1", "s2", "s3"])
    y = X.m_angle + X.s1 + 0.3 * rng.normal(size=15)
    r = slot_null(X, y, slot=["m_angle", "m_len", "m_q"], family=["m_angle"], k=2, n_perm=30)
    combos = np.array([c for c in itertools.combinations(range(6), 2) if {0, 1, 2} & set(c)])
    best = combos[q2(loo_predictions(X.to_numpy(), y, combos, lam=0), y).argmax()]
    assert r["combo"] == [X.columns[i] for i in best] and r["real_has_family"]
    assert np.isclose(r["naive_rate"], 1 / 3) and 0 <= r["rate"] <= 1 and r["p"] is not None
