"""The B1 plane normal must be the direction B1 was measured along.

The scan rotates the projected cloud, t = R p, and reads extents along the
rotated axes, so a rotated axis is the original direction R^T applied to it.
Storing R x̂ instead gives its mirror image across the frame x-axis. Nothing
errors: B1_B5_angle is still a number in [0, 90], it just depends on how the
frame happens to be turned about the Sterimol axis. On the Case Study 3 arms
the mirrored angle moved by a median 46 deg under such turns.

Two properties pin it, neither of which the mirrored normal has:
  - the cloud's extent along the stored normal IS the row's B1;
  - B1_B5_angle does not change when the substituent is turned about y.

The functions are sliced out of sterimol_utils rather than imported, for the
reason given in test_sterimol_loc_b1: the module pulls utils.help_functions ->
rdkit, built against NumPy 1.x in this interpreter. The slice is the shipped
function text. Runs under pytest or as a script.
"""
import ast
import os
import re

import numpy as np
import pandas as pd

PKG = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _load(path, names, ns):
    src = open(os.path.join(PKG, path), encoding="utf-8").read()
    for node in ast.parse(src).body:
        if isinstance(node, ast.FunctionDef) and node.name in names:
            exec(compile(ast.get_source_segment(src, node), path, "exec"), ns)


STERIMOL = "M2_data_extractor/extractor_utils/sterimol_utils.py"
SRC = open(os.path.join(PKG, STERIMOL), encoding="utf-8").read()
ns = {"np": np, "pd": pd,             # constants first: calc_loc_b1 reads one as a default
      "STERIMOL_CIRCLE_POINTS": int(re.search(r"STERIMOL_CIRCLE_POINTS = (\d+)", SRC).group(1)),
      "B1_TANGENT_TOL": float(re.search(r"B1_TANGENT_TOL = ([\d.]+)", SRC).group(1))}
_load(STERIMOL, {"generate_circle", "get_transfomed_plane_for_sterimol", "b1s_for_loop_function",
                 "scan_b1_over_angles", "get_b1s_list", "calc_loc_b1",
                 "filter_atoms_for_sterimol", "calc_sterimol"}, ns)
_load("utils/help_functions.py", {"adjust_indices"}, ns)

# A substituent in the Sterimol frame (y along the axis): x, y, z, radius.
# It wraps the axis, so no single atom sets the width over a range of
# directions: an atom centred on the axis would make every direction it alone
# touches equally narrow, and B1's direction -- hence the angle -- undefined.
# Here the narrowest plane beats any plane more than 20 deg away by 0.5 A, and
# the widest atom (the fifth) is 0.4 A clear of the next.
ATOMS = np.array([(0.0, 0.0, 0.0, 1.0), (2.0, 1.2, 0.2, 1.6), (-0.9, 2.0, 2.1, 1.6),
                  (-1.0, 2.6, -0.2, 1.1), (0.9, 3.5, -2.3, 1.7), (-2.2, 3.0, 1.0, 1.2)])


def frame(turn_deg):
    """extended_df / bonded_atoms_df for ATOMS turned about y by turn_deg."""
    g = np.radians(turn_deg)
    x = ATOMS[:, 0] * np.cos(g) - ATOMS[:, 2] * np.sin(g)
    z = ATOMS[:, 0] * np.sin(g) + ATOMS[:, 2] * np.cos(g)
    y, r = ATOMS[:, 1], ATOMS[:, 3]
    ext = pd.DataFrame({"x": x, "y": y, "z": z, "radius": r,
                        "B5": r + np.hypot(x, z), "L": y + r, "loc_B5": y})
    bonds = pd.DataFrame({"index_1": [1, 2, 3, 2, 4], "index_2": [2, 3, 4, 5, 6]})
    return bonds, ext


def test_stored_normal_carries_the_b1_extent():
    th = np.linspace(0, 2 * np.pi, ns["STERIMOL_CIRCLE_POINTS"])
    rng = np.random.default_rng(1)
    for scan in (ns["b1s_for_loop_function"], ns["scan_b1_over_angles"]):
        for _ in range(20):
            n = rng.integers(3, 10)
            P, R = rng.normal(0, 1.8, (n, 2)), rng.uniform(1.1, 1.9, n)
            cloud = np.vstack([np.c_[px + rr * np.cos(th), pz + rr * np.sin(th)]
                               for (px, pz), rr in zip(P, R)])
            args = (list(range(90)), cloud) if scan is ns["b1s_for_loop_function"] \
                else (cloud, list(range(90)))
            for _, row in scan(*args).iterrows():
                along = cloud @ row["b1_normal"]
                assert abs(min(abs(along.max()), abs(along.min())) - row["B1"]) < 1e-9, \
                    (scan.__name__, row["degree"])


def test_b1_b5_angle_ignores_how_the_frame_is_turned():
    ref = ns["calc_sterimol"](*frame(0)).iloc[0]
    assert ref["B1_B5_angle"] > 5          # a non-trivial tilt, or the test proves nothing
    for turn in (17, 45, 73, 120, 250):
        got = ns["calc_sterimol"](*frame(turn)).iloc[0]
        assert abs(got["B1"] - ref["B1"]) < 1e-2, (turn, got["B1"], ref["B1"])
        assert abs(got["B1_B5_angle"] - ref["B1_B5_angle"]) < 1.0, \
            (turn, got["B1_B5_angle"], ref["B1_B5_angle"])


if __name__ == "__main__":
    ok = True
    for t in (test_stored_normal_carries_the_b1_extent,
              test_b1_b5_angle_ignores_how_the_frame_is_turned):
        try:
            t()
            print("  ok    %s" % t.__name__)
        except AssertionError as e:
            ok = False
            print("  FAIL  %s  %s" % (t.__name__, e))
    print("\n%s" % ("ALL PASS" if ok else "FAILURES ABOVE"))
