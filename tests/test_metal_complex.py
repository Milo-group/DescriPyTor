"""Tests for GOAT / xTB metal-complex extractors.

Fast tests cover energy parsing, donor indexing, the B1 scan, theta and the xTB parsers.
"""

from __future__ import annotations

import importlib.util
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import M2_data_extractor.metal_complex as metal_complex
from M2_data_extractor.ligand_topology import LigandTopology
from M2_data_extractor.metal_complex import (
    MetalComplex,
    MetalComplexEnsemble,
    donor_indices,
)
from M2_data_extractor.xtb_singlepoint import XtbSinglePoint, parse_props_file
from M2_data_extractor.xyz_io import XYZEnsemble, energy_from_comment

HERE = Path(__file__).resolve().parent
FIXTURES = HERE / "fixtures"


@pytest.fixture(autouse=True)
def earlier_scan_grid(monkeypatch):
    """These tests run on the earlier settings: the 18° coarse B1 grid, with phi, in the lab frame,
    with the C*->R fragment bounded only at the donor and signed dipole projections."""
    monkeypatch.setattr(metal_complex, "STERIMOL_SCAN_STEP", 18)
    monkeypatch.setattr(metal_complex, "SUB_FRAGMENT_BOUND", "donor")
    monkeypatch.setattr(metal_complex, "DIPOLE_ABS", False)
    monkeypatch.setattr(metal_complex, "STERIMOL_FRAME", "lab")
    monkeypatch.setattr(metal_complex, "STERIMOL_KEYS", ("B1", "B5", "L", "angle"))
    monkeypatch.setattr(metal_complex, "STERIMOL_THETA_RULE", "scan")
    monkeypatch.setattr(metal_complex, "STEREOCENTRE_RULE", "first")




class TestEnergyParser:
    def test_last_float_skips_trailing_text(self):
        comment = "-72.9194827032 converged=true"
        assert energy_from_comment(comment, "last") == pytest.approx(-72.9194827032)
        assert energy_from_comment(comment, "first") == pytest.approx(-72.9194827032)

    def test_goat_rmsd_then_energy(self):
        comment = "0.123456 -72.91"
        assert energy_from_comment(comment, "first") == pytest.approx(0.123456)
        assert energy_from_comment(comment, "last") == pytest.approx(-72.91)
        assert energy_from_comment(comment, "auto") == pytest.approx(-72.91)

    def test_crest_relative_zero_auto(self):
        assert energy_from_comment("0.00000000", "auto") == pytest.approx(0.0)


class TestDonorIndices:
    def test_nickel_dihydride(self):
        symbols = ["Ni", "H", "H", "N", "N", "C"]
        assert donor_indices(symbols) == (0, 3, 4)

    def test_copper_monochloride(self):
        symbols = ["Cu", "Cl", "N", "N", "C"]
        assert donor_indices(symbols) == (0, 2, 3)

    def test_copper_dichloride(self):
        symbols = ["Cu", "Cl", "Cl", "N", "N"]
        assert donor_indices(symbols) == (0, 3, 4)




def test_default_scan_finds_the_plane_the_18_degree_grid_misses(monkeypatch):
    path = FIXTURES / "072_lig.finalensemble.xyz"
    if path is None:
        pytest.skip("072_lig GOAT ensemble not available")
    conformer = MetalComplexEnsemble.from_xyz(path).complexes[0]
    coarse = conformer.geometric_features()["sub_B1_sym"]
    monkeypatch.undo()                          # back to the module default
    assert conformer.geometric_features()["sub_B1_sym"] < coarse - 0.01


def _toy_chelate(bite_deg):
    half = np.radians(bite_deg / 2)
    coords = np.array([[0.0, 0.0, 0.0], [-2.1, 0.0, 0.0],
                       [2.0 * np.cos(half), 2.0 * np.sin(half), 0.0],
                       [2.0 * np.cos(half), -2.0 * np.sin(half), 0.0]])
    return MetalComplex(["Cu", "Cl", "N", "N"], coords, name="toy")


def test_bite_outside_the_builder_window_warns():
    """The builder checks the bite at placement only; the relaxed structure is re-checked."""
    with pytest.warns(UserWarning, match="toy: bite 50.0 deg"):
        _toy_chelate(50.0).geometric_features()
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        _toy_chelate(90.0).geometric_features()


def test_theta_of_a_lone_hydrogen_is_zero():
    """A lone H on the axis has no B1 direction, so phi is a tie-break; theta is 0."""
    symbols = ["C", "H", "C"]
    coords = np.array([[0.0, 0.0, 0.0], [0.3, 0.4, 1.0], [1.5, 0.0, 0.0]])
    s = metal_complex.sterimol(symbols, coords, 1, 2, [1.70, 1.00, 1.70])
    assert s["n"] == 1
    assert s["theta"] == pytest.approx(0.0, abs=1e-9)


def test_theta_does_not_depend_on_how_the_molecule_is_oriented(monkeypatch):
    """theta must use the B1 normal the scan measured along (R^T, not R)."""
    path = FIXTURES / "081_lig.finalensemble.xyz"
    if path is None:
        pytest.skip("081_lig GOAT ensemble not available")
    monkeypatch.undo()                          # module defaults: 1° scan, theta
    conformer = MetalComplexEnsemble.from_xyz(path).complexes[0]
    axis = np.array([1.0, 2.0, 3.0]) / np.sqrt(14.0)
    K = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    rot = np.eye(3) + np.sin(1.1) * K + (1 - np.cos(1.1)) * K @ K
    turned = MetalComplex(conformer.symbols, conformer.coords @ rot.T + [1.0, -2.0, 0.5])
    got, want = turned.geometric_features(), conformer.geometric_features()
    cols = [c for c in want if "_theta_" in c]
    assert len(cols) == 4
    assert max(abs(got[c] - want[c]) for c in cols) < 2.0, {c: (got[c], want[c]) for c in cols}




class TestXtbParsers:
    def test_parse_props_file(self, tmp_path):
        path = tmp_path / "props.txt"
        path.write_text("081_lig 0.10 -0.20 0.30 -8.9532 -8.1422\n", encoding="utf-8")
        out = parse_props_file(path)
        assert out["081_lig"]["homo"] == pytest.approx(-8.9532)
        assert out["081_lig"]["lumo"] == pytest.approx(-8.1422)
        assert out["081_lig"]["dipole"][0] == pytest.approx(0.10)

    def test_parse_dump_and_map_onto_tiny_ni(self, tmp_path):
        dump = tmp_path / "ni.txt"
        dump.write_text(
            "F cheap 0\n"
            "Q -0.50 0.10 0.10 -0.05 -0.05\n"
            "W 2 0.40 3 0.40 4 0.90 5 0.90\n"
            "D 0.10 0.00 0.00\n"
            "E -8.00 -7.00\n",
            encoding="utf-8",
        )
        sp = XtbSinglePoint.parse_dump(dump)["cheap"][0]
        assert sp.gap == pytest.approx(1.0)
        xyz = ROOT / "tests" / "data" / "small_set" / "xyz" / "tiny_nih2n2.xyz"
        if not xyz.is_file():
            pytest.skip("tiny_nih2n2.xyz missing")
        elec = MetalComplex.from_xyz(xyz).electronic_features(sp)
        assert elec["q_metal"] == pytest.approx(-0.50)
        assert elec["q_anc_sum"] == pytest.approx(0.20)
        assert elec["wbo_MD_sym"] == pytest.approx(0.90)
        assert elec["gap"] == pytest.approx(1.0)


