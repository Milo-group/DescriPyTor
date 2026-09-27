"""Reference Sterimol values for the fixtures, from the live Python, for gate.js.

Each fixture on its case-study axes, in both directions.

    python make_refs.py   -> refs.json
"""
import json
import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)                       # M2_data_extractor
for p in (PKG, os.path.dirname(PKG)):
    if p not in sys.path:
        sys.path.insert(0, p)
import data_extractor as de                        # noqa: E402
from extractor_utils import sterimol_utils as su    # noqa: E402

# calc_sterimol replaces the scan's B1_B5_angle (phi, in the plane) with theta before returning,
# so phi is only reachable on the way past.
_phi = {}
_b1s = su.get_b1s_list


def _spy(df, *a, **k):
    r = _b1s(df, *a, **k)
    _phi['v'] = float(pd.Series(r[1]).iloc[int(pd.Series(r[0]).idxmin())])
    return r


su.get_b1s_list = _spy

AXES = {'p-Otf.xyz': [(1, 23), (23, 1), (1, 3), (3, 1)],
        'FL_lig_13_str_target.xyz': [(6, 7), (7, 6), (6, 4), (4, 6)],
        '041_lig.xyz': [(7, 6), (6, 7)],
        'm-CN.xyz': [(1, 23), (23, 1), (1, 3)],
        '045_lig.xyz': [(8, 12), (12, 8), (19, 21), (21, 19)]}
KEYS = ('B1', 'B5', 'L', 'loc_B5', 'B1_B5_angle')

cases = []
for name, axes in AXES.items():
    xyz = pd.read_csv(os.path.join(HERE, 'fixtures', name), sep=r'\s+', skiprows=2, names=['atom', 'x', 'y', 'z'])
    bonds = de.extract_connectivity(xyz, threshold_distance=1.82)
    for a, b in axes:
        r = de.get_sterimol_df(xyz, bonds, [a, b], None, radii='CPK').iloc[0]
        cases.append(dict(file=name, a=a, b=b, ref={k: float(r[k]) for k in KEYS}, phi=round(_phi['v'], 4)))

with open(os.path.join(HERE, 'refs.json'), 'w', newline='\n') as f:
    json.dump(cases, f, indent=1)
print('%d reference axes -> refs.json' % len(cases))
