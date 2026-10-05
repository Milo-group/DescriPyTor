"""0.2.1 arm fixes: the C*->R fragment stops at the donor ring, theta is the soft average over tied
B1 directions, and the dipole projections are unsigned (docs/STERIMOL_FIXES.md, FEATURE_FIXES.md)."""
import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from M2_data_extractor import metal_complex as mc

HERE = Path(__file__).resolve().parent
LEAKY = HERE / "data" / "cs3_reference" / "024_lig_CuCl_xtbopt.xyz"          # fused substituent, leaked on both arms
TBU = HERE.parent / "M2_data_extractor" / "theta_explorer" / "fixtures" / "045_box_CuCl.xyz"
V5 = os.environ.get("DESCRIPYTOR_V5")


def _sub_atoms(m, bound, arm=0):
    donor, other = (m.donor_1, m.donor_2)[arm], (m.donor_2, m.donor_1)[arm]
    stereo, subst = mc.stereocentre(m.symbols, m.adj, donor, other, metal=m.metal)
    if bound == "donor_ring":
        block = tuple(sorted({x + 1 for x in mc.donor_ring(m.adj, donor, m.metal) | {m.metal}} - {stereo + 1}))
    else:
        block = (donor + 1, m.metal + 1)
    return set(mc.fragment(m.symbols, m.coords, stereo + 1, subst + 1, block)), donor, other


def test_donor_ring_is_the_smallest_ring_without_the_metal():
    # metal 0 - donor 1; 5-ring 1-2-3-4-5; 6-ring 1-2-6-7-8-9-5 also passes through the donor
    adj = {0: {1}, 1: {0, 2, 5}, 2: {1, 3, 6}, 3: {2, 4}, 4: {3, 5}, 5: {4, 1, 9}, 6: {2, 7}, 7: {6, 8}, 8: {7, 9}, 9: {8, 5}}
    assert mc.donor_ring(adj, 1, 0) == {1, 2, 3, 4, 5}
    assert mc.donor_ring({0: {1}, 1: {0, 2}, 2: {1}}, 1, 0) == {1}


def test_fused_substituent_no_longer_reaches_the_backbone():
    m = mc.MetalComplex.from_xyz(str(LEAKY))
    old, donor, other = _sub_atoms(m, "donor")
    new, _, _ = _sub_atoms(m, "donor_ring")
    assert other in old and other not in new                   # the old walk ran round to the other donor
    assert new < old and donor not in new and m.metal not in new


def test_soft_theta_ignores_orientation_on_a_tied_tbu():
    m = mc.MetalComplex.from_xyz(str(TBU))
    want = m.geometric_features()["sub_theta_sym"]
    rng = np.random.default_rng(3)
    for _ in range(4):
        q, r = np.linalg.qr(rng.normal(size=(3, 3)))
        q = q * np.sign(np.diag(r))
        q = q if np.linalg.det(q) > 0 else -q
        m2 = mc.MetalComplex.from_xyz(str(TBU))
        m2.coords = m2.coords @ q.T
        assert abs(m2.geometric_features()["sub_theta_sym"] - want) < 0.05


@pytest.mark.skipif(not V5, reason="set DESCRIPYTOR_V5 to the v5 paper repository")
def test_defaults_reproduce_v5_cp_fix_soft():
    """feature_matrices/single_structure/cp_fix_soft.csv, 16 Sterimol/theta columns (written at %.6f)."""
    A = Path(V5) / "zenodo_archive"
    ref = pd.read_csv(A / "feature_matrices" / "single_structure" / "cp_fix_soft.csv").set_index("name")
    geo = A / "geometries" / "case_study_3_corminboeuf" / "single_structure_xtb" / "CuCl"
    cols = [f"{fr}_{k}_{p}" for fr in ("sub", "fromM") for k in ("B1", "B5", "L", "theta") for p in ("sym", "asym")]
    got = pd.DataFrame({n: mc.MetalComplex.from_xyz(str(geo / f"{n}_CuCl_xtbopt.xyz")).geometric_features() for n in ref.index}).T
    assert float((got[cols].astype(float) - ref[cols]).abs().max().max()) < 1e-6


def test_arm_values_do_not_depend_on_atom_order():
    """The (C*, R) pick and everything read from it survive a shuffle of the ligand atoms."""
    m = mc.MetalComplex.from_xyz(str(LEAKY))
    first = m.donor_2 + 1
    order = list(range(first)) + list(np.random.default_rng(0).permutation(np.arange(first, len(m.symbols))))
    shuffled = mc.MetalComplex([m.symbols[i] for i in order], m.coords[order])
    a, b = m.geometric_features(), shuffled.geometric_features()
    assert a.keys() == b.keys()
    for k in a:
        assert abs(a[k] - b[k]) < 1e-6, k


def test_p_aryl_arm_has_no_substituent_axis():
    """P bonded to two phenyls and a backbone carbon: ipso -> ortho is a ring bond, not a substituent."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    mol = Chem.AddHs(Chem.MolFromSmiles("[Cu](Cl)(N(C)(C)CC1)P1(c1ccccc1)c1ccccc1"))
    AllChem.EmbedMolecule(mol, randomSeed=1)
    sym = [a.GetSymbol() for a in mol.GetAtoms()]
    adj = [[n.GetIdx() for n in a.GetNeighbors()] for a in mol.GetAtoms()]
    p, n = sym.index("P"), sym.index("N")
    assert mc.stereocentre(sym, adj, p, n, metal=0) is None
