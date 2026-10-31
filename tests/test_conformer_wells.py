"""The well vote on a torus with four known centres, and a CREST/GOAT tag."""
from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("rdkit")

from M2_data_extractor.conformer_duplicates import cross_tags
from M2_data_extractor.conformer_wells import vote_wells


def test_four_separated_wells_are_the_vote():
    rng = np.random.default_rng(0)
    centres = [(0.0, 0.0), (0.0, 180.0), (180.0, 0.0), (180.0, 180.0)]
    angles, groups = [], []
    for ligand in range(3):
        for a, b in centres:
            for _copy in range(5):
                angles.append([a + float(rng.normal(0, 3)), b + float(rng.normal(0, 3))])
                groups.append(f"L{ligand}")
    got = vote_wells(angles, groups)
    assert got["criteria"]["density"] == 4
    assert got["suggested"] == 4
    assert set(got["labels"]["4"]) == {0, 1, 2, 3}
    # copies of every ligand sit in every well, so coverage stays at least at the four real wells
    assert got["criteria"]["coverage"] >= 4


def test_coverage_stops_when_a_file_misses_a_well():
    rng = np.random.default_rng(2)
    angles, groups = [], []
    for ligand in ("A", "B"):
        for centre in (0.0, 180.0):
            for _copy in range(8):
                angles.append([centre + float(rng.normal(0, 2))])
                groups.append(ligand)
    for _copy in range(8):
        angles.append([90.0 + float(rng.normal(0, 2))])
        groups.append("A")
    got = vote_wells(angles, groups, k_max=4)
    assert got["criteria"]["coverage"] == 2


def test_one_file_has_no_coverage_vote():
    rng = np.random.default_rng(1)
    angles = [[float(rng.normal(0, 2)), float(rng.normal(90, 2))] for _i in range(12)]
    angles += [[float(rng.normal(180, 2)), float(rng.normal(-90, 2))] for _i in range(12)]
    got = vote_wells(angles, ["only"] * len(angles))
    assert got["criteria"]["coverage"] is None
    assert got["criteria"]["density"] == 2


def test_crest_and_goat_share_one_zero():
    from rdkit import Chem
    from rdkit.Chem import AllChem
    mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    assert AllChem.EmbedMolecule(mol, randomSeed=1) == 0
    conf = mol.GetConformer()
    rows = [(a.GetSymbol(), conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y, conf.GetAtomPosition(i).z) for i, a in enumerate(mol.GetAtoms())]
    other = [(el, z + 1.0, x - 2.0, y) for el, x, y, z in rows]
    lone = [(el, x + 5.0, y, z) for el, x, y, z in rows]
    lone[0] = (lone[0][0], lone[0][1], lone[0][2] + 2.0, lone[0][3])
    got = cross_tags([
        {"name": "basic · CREST", "frames": [
            {"id": 1, "atoms": rows, "energy": -10.0},
            {"id": 2, "atoms": lone, "energy": -9.9},
        ]},
        {"name": "basic · GOAT", "frames": [
            {"id": 7, "atoms": other, "energy": -9.99},
        ]},
    ])
    by = {row["id"]: row for row in got["tags"]}
    assert by[1]["tag"] == "BOTH" and by[7]["tag"] == "BOTH"
    assert by[1]["with"] == 7
    assert by[2]["tag"] == "CREST"
    assert by[1]["joint_kcal"] == 0.0
    assert got["minimum"] == -10.0
