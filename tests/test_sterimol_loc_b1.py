"""calc_loc_b1 on cases whose answer is known by hand.

The function is sliced out of sterimol_utils rather than imported, because the
module pulls utils.help_functions -> rdkit, which is built against NumPy 1.x in
this interpreter. The slice is the function text itself, so what is tested is
what ships.

Each case places four atoms so that exactly one of the four extents is smallest
in absolute value, which fixes which edge B1 sits on, and so that the atoms
reaching that edge are known. Two flanking atoms at large |z| (or |x|) keep the
perpendicular extents out of the way.
"""
import os
import re

import numpy as np
import pandas as pd

PKG = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = open(os.path.join(PKG, "M2_data_extractor/extractor_utils/"
                        "sterimol_utils.py"), encoding="utf-8").read()

ns = {"np": np, "pd": pd}
ns["STERIMOL_CIRCLE_POINTS"] = int(
    re.search(r"STERIMOL_CIRCLE_POINTS = (\d+)", SRC).group(1))
ns["B1_TANGENT_TOL"] = float(
    re.search(r"B1_TANGENT_TOL = ([\d.]+)", SRC).group(1))
exec(compile(SRC[SRC.index("def calc_loc_b1"):SRC.index("def calc_sterimol")],
             "calc_loc_b1", "exec"), ns)
calc_loc_b1 = ns["calc_loc_b1"]
NPT = ns["STERIMOL_CIRCLE_POINTS"]
ok = True


def build(atoms):
    """atoms: (x, z, radius, y) -> the plane and frame the B1 scan works on."""
    th = np.linspace(0, 2 * np.pi, NPT)
    plane = np.vstack([np.column_stack((x + r * np.cos(th), z + r * np.sin(th)))
                       for x, z, r, _ in atoms])
    df = pd.DataFrame({"x": [a[0] for a in atoms], "z": [a[1] for a in atoms],
                       "radius": [a[2] for a in atoms],
                       "y": [a[3] for a in atoms]})
    return plane, df


def run(label, atoms, want_loc, want_n, want_edge):
    global ok
    plane, df = build(atoms)
    ev = [plane[:, 0].max(), plane[:, 0].min(),
          plane[:, 1].max(), plane[:, 1].min()]
    edge = ["+x", "-x", "+z", "-z"][int(np.argmin(np.abs(ev)))]
    loc, nt = calc_loc_b1(plane, df)
    good = (abs(loc - want_loc) < 1e-9 and nt == want_n and edge == want_edge)
    ok &= good
    print("  %-44s edge %-2s  loc_B1 %6.3f (want %5.2f)  atoms %d (want %d)  %s"
          % (label, edge, loc, want_loc, nt, want_n, "ok" if good else "FAIL"))


# Flankers hold the extents we are not testing away from zero. FZ sits far
# enough left that min_x stays outside the +x edge under test; FZ4 sits at x=0
# so that a test atom can own min_x; FX does the same job for the z edges.
FZ = [(-2.5, 3.0, 0.5, 2.0), (-2.5, -3.0, 0.5, 3.0)]
FZ4 = [(0.0, 3.0, 0.5, 2.0), (0.0, -3.0, 0.5, 3.0)]
FX = [(3.0, 0.0, 0.5, 1.0), (-3.0, 0.0, 0.5, 3.0)]

print("calc_loc_b1")
run("single contact on the +x edge",
    [(1.0, 0.0, 0.5, 5.0)] + FZ, 5.0, 1, "+x")

run("two atoms tied on the edge -> mean of their y",
    [(1.0, -0.6, 0.5, 4.0), (1.0, 0.6, 0.5, 8.0)] + FZ, 6.0, 2, "+x")

run("tie broken by radius, not by centre position",
    [(1.0, -0.6, 0.7, 4.0), (1.0, 0.6, 0.5, 8.0)] + FZ, 4.0, 1, "+x")

run("the -x edge when that is the narrow one",
    [(-1.0, 0.0, 0.5, 7.0), (5.0, 0.0, 0.5, 2.0)] + FZ4, 7.0, 1, "-x")

run("the +z edge exercises the other column",
    [(0.0, 1.0, 0.5, 9.0), (0.0, -4.0, 0.5, 2.0)] + FX, 9.0, 1, "+z")

run("three-way contact averages all three",
    [(1.0, -0.8, 0.5, 3.0), (1.0, 0.0, 0.5, 6.0),
     (1.0, 0.8, 0.5, 9.0)] + FZ, 6.0, 3, "+x")

print("\ndegenerate inputs return nan rather than raising")
for lab, pl, fr in (("empty frame", np.zeros((0, 2)), pd.DataFrame({"y": []})),
                    ("plane length mismatch", np.zeros((7, 2)),
                     pd.DataFrame({"y": [1.0, 2.0]}))):
    loc, nt = calc_loc_b1(pl, fr)
    good = np.isnan(loc) and nt == 0
    ok &= good
    print("  %-30s -> (%s, %d)  %s" % (lab, loc, nt, "ok" if good else "FAIL"))

print("\n%s" % ("ALL PASS" if ok else "FAILURES ABOVE"))
