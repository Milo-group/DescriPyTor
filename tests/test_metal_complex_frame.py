"""MetalComplex Sterimol must not depend on how the input file is oriented.

The B1 rotation scan starts at e1. It used to be np.cross(axis, lab z), so a
rigid rotation of the file moved theta by up to 25 deg and phi by up to 83 deg.
"""
from pathlib import Path

import numpy as np
import pytest

from M2_data_extractor.metal_complex import MetalComplex

XYZ = Path(__file__).resolve().parents[1] / "M2_data_extractor" / "theta_explorer" / "fixtures" / "026_box_CuCl.xyz"


def _rotations(n, seed=7):
    rng = np.random.default_rng(seed)
    for _ in range(n):
        q, r = np.linalg.qr(rng.normal(size=(3, 3)))
        q = q * np.sign(np.diag(r))
        yield q if np.linalg.det(q) > 0 else -q


def test_sterimol_is_invariant_to_rigid_rotation():
    mc = MetalComplex.from_xyz(str(XYZ))
    want = mc.geometric_features()
    cols = [c for c in want if c.split("_")[0] in ("fromM", "sub")
            and any(k in c for k in ("_B1_", "_B5_", "_L_", "_theta_"))]
    assert len(cols) == 16
    for R in _rotations(5):
        got = MetalComplex(mc.symbols, mc.coords @ R.T + [0.3, -1.2, 2.0]).geometric_features()
        worst = max(abs(got[c] - want[c]) for c in cols)
        assert worst < 1e-6, {c: (got[c], want[c]) for c in cols if abs(got[c] - want[c]) >= 1e-6}
