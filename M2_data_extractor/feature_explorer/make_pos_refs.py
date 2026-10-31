"""Reference %Vbur and exact cone angles for the fixtures, for gate_pos.js.

    <mol_ml python> make_pos_refs.py   -> pos_refs.json
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))
from M2_data_extractor import metal_complex as mc            # noqa: E402
import morfeus                                                # noqa: E402
import warnings
warnings.filterwarnings("ignore")

# fixture -> (vbur centres as 1-based atom indices with radii, cone apices)
CASES = {
    '009_phos_CuCl.xyz': dict(vbur=[(1, 3.0), (1, 3.5), (1, 5.0), (3, 3.5)], cone=[1, 3]),
    '009_phos.xyz': dict(vbur=[(15, 3.5), (15, 5.0)], cone=[15]),
    '055_phos.xyz': dict(vbur=[(10, 3.5), (10, 5.0)], cone=[10]),
    '045_lig.xyz': dict(vbur=[(8, 3.5), (12, 3.5)], cone=[]),   # a carbon apex has its neighbours inside
}


def xyz(name):
    L = open(os.path.join(HERE, 'fixtures', name)).read().split('\n')
    n = int(L[0])
    el = [r.split()[0] for r in L[2:2 + n]]
    X = np.array([[float(v) for v in r.split()[1:4]] for r in L[2:2 + n]])
    return el, X


out = []
for name, job in CASES.items():
    el, X = xyz(name)
    for i, r in job['vbur']:
        out.append(dict(file=name, kind='vbur', atom=i, radius=r,
                        ref=round(float(mc.buried_volume(el, X, X[i - 1], radius=r)), 6)))
    for i in job['cone']:
        ca = morfeus.ConeAngle(el, X, i)
        out.append(dict(file=name, kind='cone', atom=i,
                        ref=round(float(ca.cone_angle), 6),
                        tangent=sorted(int(t) for t in getattr(ca, 'tangent_atoms', []))))
with open(os.path.join(HERE, 'pos_refs.json'), 'w', newline='\n') as f:
    json.dump(out, f, indent=1)
print('%d reference values -> pos_refs.json' % len(out))
