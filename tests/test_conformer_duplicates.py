"""Symmetry-aware conformer duplicates: a ring flip is one conformer, a real torsion change is not."""
from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("rdkit")
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdMolTransforms import SetDihedralDeg

from M2_data_extractor.conformer_duplicates import fold_duplicates, same_conformer, mol_from_heavy, heavy_rows


def _rows(mol):
    conf = mol.GetConformer()
    return [
        (a.GetSymbol(), conf.GetAtomPosition(a.GetIdx()).x, conf.GetAtomPosition(a.GetIdx()).y, conf.GetAtomPosition(a.GetIdx()).z)
        for a in mol.GetAtoms()
    ]


def _rotate(rows, angle=0.7):
    c, s = np.cos(angle), np.sin(angle)
    r = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
    out = []
    for el, x, y, z in rows:
        v = r @ np.array([x, y, z]) + np.array([1.5, -2.0, 0.4])
        out.append((el, float(v[0]), float(v[1]), float(v[2])))
    return out


def _embed(smiles, seed):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(mol, randomSeed=seed) == 0
    return mol


def test_a_rotated_copy_is_the_same_conformer():
    rows = _rows(_embed("CCO", 1))
    got = fold_duplicates([{"atoms": rows, "energy": -1.0}, {"atoms": _rotate(rows), "energy": -0.9}])
    assert got["folded"] == 1
    assert got["pairs"][0]["parent"] == 0
    assert got["pairs"][0]["child"] == 1
    assert got["pairs"][0]["rmsd"] < 0.05
    assert got["pairs"][0]["reason"] == "tight"


def test_a_para_ring_flip_is_one_conformer_and_names_the_atoms():
    mol = _embed("Fc1ccc(F)cc1", 2)
    rows = _rows(mol)
    heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    flip = next(
        m for m in mol.GetSubstructMatches(mol, uniquify=False, maxMatches=32)
        if any(m[i] != i for i in heavy)
    )
    flipped = [None] * len(rows)
    for i, j in enumerate(flip):
        src = rows[j]
        flipped[i] = (rows[i][0], src[1], src[2], src[3])
    got = fold_duplicates([{"atoms": rows, "energy": -2.0}, {"atoms": flipped, "energy": -1.9}])
    assert got["folded"] == 1
    assert got["pairs"][0]["rmsd"] < 0.05
    assert got["pairs"][0]["swaps"]


def test_a_close_shift_is_offered_not_hidden():
    rows = _rows(_embed("CCO", 6))
    found = None
    for step in (0.15, 0.2, 0.28, 0.35):
        nudged = [list(row) for row in rows]
        nudged[1][1] = float(nudged[1][1]) + step
        got = fold_duplicates([{"atoms": rows, "energy": -1.0}, {"atoms": nudged, "energy": -0.9}])
        rms = None
        if got["suggestions"]:
            rms = got["suggestions"][0]["rmsd"]
        elif got["pairs"]:
            rms = got["pairs"][0]["rmsd"]
        if rms is not None and 0.05 < rms <= 0.15:
            found = got
            break
    assert found is not None
    assert found["pairs"] == []
    assert found["suggestions"][0]["reason"] == "torsion"


def test_a_real_torsion_change_stays_two_conformers():
    mol = _embed("CCCC", 3)
    carbons = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == "C"]
    anti, gauche = Chem.Mol(mol), Chem.Mol(mol)
    SetDihedralDeg(anti.GetConformer(), carbons[0], carbons[1], carbons[2], carbons[3], 180.0)
    SetDihedralDeg(gauche.GetConformer(), carbons[0], carbons[1], carbons[2], carbons[3], 60.0)
    got = fold_duplicates([
        {"atoms": _rows(anti), "energy": -3.0},
        {"atoms": _rows(gauche), "energy": -2.5},
    ])
    assert got["folded"] == 0
    assert got["kept"] == 2


def test_a_frame_with_missing_coordinates_does_not_crash():
    rows = _rows(_embed("CCO", 3))
    broken = [list(r) for r in rows]
    broken[1] = [broken[1][0], None, None, None]
    blank = {"atoms": [["C", None, None, None], ["C", None, None, None], ["H", None, None, None]], "energy": 0.0}
    got = fold_duplicates([
        {"atoms": rows, "energy": -2.0},
        {"atoms": _rotate(rows), "energy": -1.0},
        {"atoms": broken, "energy": -1.5},
        blank,
    ])
    assert got["folded"] == 1
    assert got["pairs"][0]["child"] == 1
    assert got["pairs"][0]["parent"] == 0


def test_the_lower_energy_frame_is_the_one_kept():
    rows = _rows(_embed("CCO", 4))
    got = fold_duplicates([
        {"atoms": _rotate(rows, 0.3), "energy": -1.0},
        {"atoms": rows, "energy": -2.0},
    ])
    assert got["pairs"][0]["parent"] == 1
    assert got["pairs"][0]["child"] == 0


def test_grey_zone_rejects_a_large_core_torsion():
    base = _embed("CCCC", 5)
    carbons = [a.GetIdx() for a in base.GetAtoms() if a.GetSymbol() == "C"]
    anti, gauche = Chem.Mol(base), Chem.Mol(base)
    SetDihedralDeg(anti.GetConformer(), *carbons, 180.0)
    SetDihedralDeg(gauche.GetConformer(), *carbons, 60.0)
    mol_a = mol_from_heavy(heavy_rows(_rows(anti)))
    mol_b = mol_from_heavy(heavy_rows(_rows(gauche)))
    ok, rms, why = same_conformer(mol_a, mol_b, rmsd_cut=5.0, tight=0.0, torsion_cut=15.0)
    assert rms is not None and rms > 0.0
    assert ok is False
    assert why == "rejected"
