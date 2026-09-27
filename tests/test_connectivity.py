"""extract_connectivity: covalent radii x BOND_SCALE by default, the old flat cutoff on request."""
from pathlib import Path

import pandas as pd

from utils.help_functions import extract_connectivity

FIX = Path(__file__).resolve().parents[1] / "M2_data_extractor" / "theta_explorer" / "fixtures"


def _xyz(name):
    lines = (FIX / name).read_text().splitlines()
    n = int(lines[0].split()[0])
    df = pd.DataFrame([l.split()[:4] for l in lines[2:2 + n]], columns=["atom", "x", "y", "z"])
    df[["x", "y", "z"]] = df[["x", "y", "z"]].astype(float)
    return df


def _bonds(df, **kw):
    return {tuple(sorted(r)) for r in extract_connectivity(df, **kw).values.tolist()}


def test_long_single_bonds_are_kept():
    """S-CF3 of the CS1 triflate is 1.841 A: dropped by the flat 1.82 A rule, kept now."""
    x = _xyz("p-Otf.xyz")
    assert (34, 37) in _bonds(x)
    assert (34, 37) not in _bonds(x, threshold_distance=1.82)


def test_phosphine_p_c_bonds_are_kept():
    x = _xyz("009_phos.xyz")
    p = [i + 1 for i, e in enumerate(x.atom) if e == "P"]
    new = sum(1 for b in _bonds(x) if p[0] in b)
    old = sum(1 for b in _bonds(x, threshold_distance=1.82) if p[0] in b)
    assert new == 3 and old < 3


def test_agostic_metal_h_contact_is_not_a_bond():
    """Cu...H at 2.68 A in 045_lig is inside the 2.8 A metal cutoff but is not a hydride."""
    x = _xyz("045_lig.xyz")
    cu = [i + 1 for i, e in enumerate(x.atom) if e == "Cu"][0]
    h = {i + 1 for i, e in enumerate(x.atom) if e == "H"}
    assert not any(cu in b and (set(b) & h) for b in _bonds(x))


def test_the_rule_changes_nothing_else_on_the_fixtures():
    """On every fixture the covalent rule only adds S-C / P-C bonds the flat rule dropped."""
    for path in sorted(FIX.glob("*.xyz")):
        x = _xyz(path.name)
        new, old = _bonds(x), _bonds(x, threshold_distance=1.82)
        assert old <= new, (path.name, old - new)
        for i, j in new - old:
            assert {x.atom[i - 1], x.atom[j - 1]} in ({"S", "C"}, {"P", "C"}), (path.name, i, j)
