"""Tests for GOAT / xTB metal-complex extractors.

Fast tests cover energy parsing, donor indexing, topology, and the two-conformer
081_lig GOAT file (copied into tests/fixtures when available).

Recreation of the published CS3 tables runs when ``CS3_SCRATCHPAD`` points at
the original scratchpad: GOAT ensembles in ``cbens`` / ``gcu``,
xTB charges in ``chg_*``, ``props*.txt``, cheap XYZ, and SMILES in ``mc/``.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

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
CORNMINBUF = (
    ROOT / "Getting_started_with_examples" / "case_study_notebooks" / "cornminbuf"
)

def scratchpad():
    import os

    env = os.environ.get("CS3_SCRATCHPAD")
    if not env:
        return None
    path = Path(env)
    if path.is_dir() and (path / "cbens").is_dir() and (path / "chg_ni").is_dir():
        return path
    return None

GEOM = CORNMINBUF / "geom"
GOAT_TABLES = CORNMINBUF / "data" / "goat_ensemble"
CHEAP_TABLES = CORNMINBUF / "data" / "single_structure"


def _max_abs(got: pd.Series, want: pd.Series) -> float:
    return float(np.nanmax(np.abs(got.astype(float) - want.astype(float))))


def _compare_row(got: dict, want: pd.Series, cols, atol=1e-8, rtol=1e-9):
    mismatches = []
    for col in cols:
        if col not in got or col not in want.index or pd.isna(want[col]):
            continue
        g, w = float(got[col]), float(want[col])
        if not np.isclose(g, w, atol=atol, rtol=rtol, equal_nan=True):
            mismatches.append((col, g, w, abs(g - w)))
    return mismatches


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


class TestTopology:
    @pytest.mark.skipif(
        importlib.util.find_spec("rdkit") is None,
        reason="rdkit is not installed in this interpreter",
    )
    def test_081_matches_cs3_tf_block(self):
        table_path = GOAT_TABLES / "cc.csv"
        if not table_path.is_file():
            pytest.skip("CS3 goat_ensemble tables not in this checkout")
        smiles = "c1ccc(C[C@H]2COC(C3=N[C@@H](Cc4ccccc4)CO3)=N2)cc1"
        table = pd.read_csv(table_path, index_col="name")
        got = LigandTopology.from_smiles(smiles, name="081_lig").size_normalized_features()
        want = table.loc["081_lig"]
        cols = [c for c in got if c in want.index]
        assert cols, "no overlapping tf_ columns"
        bad = _compare_row(got, want, cols, atol=1e-10)
        assert not bad, bad[:5]


EXPECTED_SMALL = FIXTURES / "expected_cc_goat_ens_small.csv"
SMALL_GOAT = ("081_lig", "072_lig", "080_lig", "083_lig")


def _ensemble_path(stem: str) -> Path | None:
    local = FIXTURES / f"{stem}.finalensemble.xyz"
    if local.is_file():
        return local
    pad = scratchpad()
    if pad is not None:
        candidate = pad / "cbens" / f"{stem}.finalensemble.xyz"
        if candidate.is_file():
            return candidate
    return None


def _081_ensemble_path() -> Path | None:
    return _ensemble_path("081_lig")


def _expected_geom_table() -> pd.DataFrame | None:
    if EXPECTED_SMALL.is_file():
        return pd.read_csv(EXPECTED_SMALL, index_col=0)
    table = GEOM / "desc_cc_goat_ens.csv"
    if table.is_file():
        return pd.read_csv(table, index_col=0)
    return None


@pytest.mark.parametrize("stem", SMALL_GOAT)
def test_small_goat_matches_published_geom(stem):
    path = _ensemble_path(stem)
    want_table = _expected_geom_table()
    if path is None or want_table is None or stem not in want_table.index:
        pytest.skip(f"fixture or expected row missing for {stem}")
    ens = MetalComplexEnsemble.from_xyz(path, name=stem)
    got = ens.geometric_features()
    want = want_table.loc[stem]
    cols = [c for c in want.index if c in got]
    bad = _compare_row(got, want, cols, atol=1e-8)
    assert not bad, bad
    assert got["n_conformers"] == want["n_conformers"]


@pytest.mark.skipif(_081_ensemble_path() is None, reason="081_lig GOAT ensemble not available")
class TestGoatGeometry081:
    def test_xyz_ensemble_last_float_energy(self):
        ens = XYZEnsemble(_081_ensemble_path(), energy_convention="last")
        assert ens.n_conformers == 2
        assert ens.energies_hartree[0] == pytest.approx(-72.9194827032)


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
        assert sp.dipole_unit == "au"
        xyz = ROOT / "tests" / "data" / "small_set" / "xyz" / "tiny_nih2n2.xyz"
        if not xyz.is_file():
            pytest.skip("tiny_nih2n2.xyz missing")
        mc = MetalComplex.from_xyz(xyz)
        elec = mc.electronic_features(sp)
        assert elec["q_metal"] == pytest.approx(-0.50)
        assert elec["q_anc_sum"] == pytest.approx(0.20)
        assert elec["wbo_MD_sym"] == pytest.approx(0.90)
        assert elec["wbo_MD_asym"] == pytest.approx(0.0)
        assert elec["gap"] == pytest.approx(1.0)

    def test_from_files_q_and_wbo(self, tmp_path):
        q_path = tmp_path / "x.q"
        w_path = tmp_path / "x.wbo"
        q_path.write_text("-0.50 0.10 0.10 -0.05 -0.05\n", encoding="utf-8")
        w_path.write_text("1 4 0.91\n1 5 0.89\n", encoding="utf-8")
        # electronic_features always rotates the dipole into the N–M–N frame
        sp = XtbSinglePoint.from_files(
            q_path, wbo_path=w_path, homo=-8.0, lumo=-7.5, dipole=[0.10, 0.0, 0.0]
        )
        xyz = ROOT / "tests" / "data" / "small_set" / "xyz" / "tiny_nih2n2.xyz"
        if not xyz.is_file():
            pytest.skip("tiny_nih2n2.xyz missing")
        elec = MetalComplex.from_xyz(xyz).electronic_features(sp)
        assert elec["wbo_MD_sym"] == pytest.approx(0.90)
        assert elec["homo"] == pytest.approx(-8.0)
        assert elec["lumo"] == pytest.approx(-7.5)


CHARGE_WBO_COLS = [
    "homo", "lumo", "gap", "q_metal", "q_donor_sym", "q_donor_asym",
    "q_anc_sum", "q_absmax", "q_spread", "wbo_MD_sym", "wbo_MD_asym",
    "wbo_metal_tot",
]


def _rebuild_elec(pad, elec_csv, chg_dir, geo_dir, props_file):
    want = pd.read_csv(pad / elec_csv, index_col=0)
    props = parse_props_file(pad / props_file)
    chg, geo = pad / chg_dir, pad / geo_dir
    got_rows = {}
    skipped = []
    for stem in want.index:
        q_path = chg / f"{stem}.q"
        xyz_path = geo / f"{stem}.xyz"
        if not q_path.is_file() or not xyz_path.is_file() or stem not in props:
            skipped.append(stem)
            continue
        mc = MetalComplex.from_xyz(xyz_path)
        sp = XtbSinglePoint.from_xtb_dir(str(chg), stem, props=props[stem])
        if len(sp.charges) != len(mc.symbols):
            skipped.append(stem)
            continue
        got_rows[stem] = mc.electronic_features(sp)
    got = pd.DataFrame.from_dict(got_rows, orient="index")
    return got, want, skipped


def _scratch():
    pad = scratchpad()
    if pad is None:
        pytest.skip("CS3 scratchpad not found (set CS3_SCRATCHPAD)")
    return pad


def _assert_charge_wbo(got, want, label):
    cols = [c for c in CHARGE_WBO_COLS if c in got.columns and c in want.columns]
    idx = got.index.intersection(want.index)
    assert len(idx) > 10, f"{label}: too few overlapping ligands ({len(idx)})"
    bad = {}
    for col in cols:
        g, w = got.loc[idx, col].astype(float), want.loc[idx, col].astype(float)
        if col.startswith("wbo_MD"):
            mask = ~((w.abs() < 1e-12) & (g.abs() > 1e-6))
            if mask.sum() == 0:
                continue
            d = float(np.nanmax(np.abs(g[mask] - w[mask])))
        else:
            d = float(np.nanmax(np.abs(g - w)))
        if d > 1e-9:
            bad[col] = d
    assert not bad, f"{label} charge/WBO/orbital mismatches: {bad}"
    return idx


class TestElectronicsFromScratchpad:
    def test_elec_ni_table(self):
        pad = _scratch()
        got, want, _ = _rebuild_elec(pad, "elec_ni.csv", "chg_ni", "cb_complexes", "props.txt")
        idx = _assert_charge_wbo(got, want, "elec_ni")
        assert len(idx) >= 90

    def test_elec_cu1_table(self):
        pad = _scratch()
        got, want, _ = _rebuild_elec(pad, "elec_cu1.csv", "chg_cu1", "cb_cu1", "props_cu1.txt")
        _assert_charge_wbo(got, want, "elec_cu1")

    def test_elec_cu2_table(self):
        pad = _scratch()
        got, want, _ = _rebuild_elec(pad, "elec_cu2.csv", "chg_cu2", "cb_cu2", "props_cu2.txt")
        _assert_charge_wbo(got, want, "elec_cu2")


class TestCheapArm:
    def test_cc_cheap_geometry_and_electronics(self):
        table = GEOM / "desc_cc_cheap.csv"
        if not table.is_file():
            pytest.skip("CS3 geom table not in this checkout")
        pad = _scratch()
        want_g = pd.read_csv(table, index_col=0)
        want_e = pd.read_csv(pad / "elec_cc_cheap.csv", index_col=0)
        geo_rows, elec_rows = {}, {}
        for stem in want_g.index:
            xyz = pad / "cheap_opt" / f"ni_{stem}.xyz"
            dump = pad / "cheap_out" / f"ni_{stem}.txt"
            if not xyz.is_file() or not dump.is_file():
                continue
            mc = MetalComplex.from_xyz(xyz)
            geo_rows[stem] = mc.geometric_features()
            parsed = XtbSinglePoint.parse_dump(dump)
            sp = parsed.get("cheap", {}).get(0)
            if sp is not None:
                elec_rows[stem] = mc.electronic_features(sp)
        assert geo_rows, "no cheap Ni geometries found"
        got_g = pd.DataFrame.from_dict(geo_rows, orient="index")
        cols = [c for c in want_g.columns if c != "n_conformers" and c in got_g.columns]
        idx = got_g.index.intersection(want_g.index)
        bad = {c: _max_abs(got_g.loc[idx, c], want_g.loc[idx, c]) for c in cols}
        bad = {c: d for c, d in bad.items() if d > 1e-8}
        assert not bad, f"desc_cc_cheap mismatches: {bad}"

        got_e = pd.DataFrame.from_dict(elec_rows, orient="index")
        ecols = [c for c in want_e.columns if c in got_e.columns]
        eidx = got_e.index.intersection(want_e.index)
        ebad = {c: _max_abs(got_e.loc[eidx, c], want_e.loc[eidx, c]) for c in ecols}
        ebad = {c: d for c, d in ebad.items() if d > 1e-8}
        assert not ebad, f"elec_cc_cheap mismatches: {ebad}"


class TestAgreesWithExtractCb:
    """Package kernel must match the scratchpad extractor, including H-arms."""

    def test_040_lig_matches_extract_cb(self):
        pad = _scratch()
        extract_cb_path = pad / "extract_cb.py"
        if not extract_cb_path.is_file():
            pytest.skip("scratchpad extract_cb.py not found")
        sys.path.insert(0, str(pad))
        import extract_cb  # noqa: E402

        path = pad / "gcu" / "cu1_040_lig.ens.xyz"
        orig, n = extract_cb.boltzmann(str(path))
        got = MetalComplexEnsemble.from_xyz(path, name="040_lig").geometric_features()
        bad = _compare_row(got, pd.Series(orig), list(orig), atol=1e-10)
        assert not bad, bad[:10]
        assert n == got["n_conformers"]


@pytest.mark.slow
class TestGoatTables:
    """Recreate Boltzmann-averaged geometric tables from GOAT XYZ files.

    Slow: thousands of conformers × Sterimol + %Vbur. Skip with
    ``pytest -m "not slow"``.
    """

    def _check_goat_dir(self, rxn, directory, pattern, name_from_path):
        pad = _scratch()
        want = pd.read_csv(GEOM / f"desc_{rxn}_goat_ens.csv", index_col=0)
        paths = []
        names = []
        for stem in want.index:
            path = name_from_path(pad / directory, stem)
            if path is not None and path.is_file():
                paths.append(path)
                names.append(stem)
        assert names, f"no GOAT files for {rxn}"
        print(f"{rxn}: {len(names)} GOAT ensembles", flush=True)
        rows = {}
        for i, (path, stem) in enumerate(zip(paths, names), 1):
            ens = MetalComplexEnsemble.from_xyz(path, name=stem)
            rows[stem] = ens.geometric_features()
            print(f"  {i}/{len(names)} {stem} n={ens.n_conformers}", flush=True)
        got = pd.DataFrame.from_dict(rows, orient="index")
        cols = [c for c in want.columns if c in got.columns]
        idx = got.index.intersection(want.index)
        bad = {}
        for col in cols:
            d = _max_abs(got.loc[idx, col], want.loc[idx, col])
            if d > 1e-8:
                bad[col] = d
        return idx, bad, got, want

    def _assert_goat_table(self, idx, bad, got, want, label):
        """All columns must match the CSV except H-substituent ``sub_angle``.

        A 1-atom (H) Sterimol fragment has B1 ≈ B5 and a ~180° B1–B5 angle.
        Current ``extract_cb.py`` (and this package) report that angle; the
        on-disk Cu CSVs do not. ``sub_B1`` / ``sub_B5`` / ``sub_L`` still
        match. Near-ties in the 1° B1 scan can also move ``sub_angle`` by
        < 0.1°; those are accepted.
        """
        angle_cols = {"sub_angle_sym", "sub_angle_asym"}
        rest = {c: d for c, d in bad.items() if c not in angle_cols}
        assert not rest, f"{label} mismatches: {rest}"
        for name in idx:
            diffs = []
            for col in angle_cols:
                if col in got.columns and col in want.columns:
                    diffs.append(abs(float(got.loc[name, col]) - float(want.loc[name, col])))
            if not diffs or max(diffs) <= 1e-8:
                continue
            if max(diffs) <= 0.1:
                continue
            b1_asym = float(got.loc[name, "sub_B1_asym"]) if "sub_B1_asym" in got.columns else 0.0
            assert b1_asym > 1.0, (
                f"{label} {name} sub_angle drifted {max(diffs):.4g}° but "
                f"sub_B1_asym={b1_asym} (expected an H vs heavy substituent)"
            )

    def test_cc_goat_ensemble_table(self):
        idx, bad, got, want = self._check_goat_dir(
            "cc", "cbens", "*.finalensemble.xyz",
            lambda d, stem: d / f"{stem}.finalensemble.xyz",
        )
        assert len(idx) == len(want), f"cc coverage {len(idx)}/{len(want)}"
        self._assert_goat_table(idx, bad, got, want, "desc_cc_goat_ens")

    def test_cp_goat_ensemble_table(self):
        idx, bad, got, want = self._check_goat_dir(
            "cp", "gcu", "*.ens.xyz",
            lambda d, stem: d / f"cu1_{stem}.ens.xyz",
        )
        assert len(idx) >= 20
        self._assert_goat_table(idx, bad, got, want, "desc_cp_goat_ens")

    def test_oa_goat_ensemble_table(self):
        idx, bad, got, want = self._check_goat_dir(
            "oa", "gcu", "*.ens.xyz",
            lambda d, stem: d / f"cu1_{stem}.ens.xyz",
        )
        assert len(idx) >= 10
        self._assert_goat_table(idx, bad, got, want, "desc_oa_goat_ens")

    def test_da_f_goat_ensemble_table(self):
        idx, bad, got, want = self._check_goat_dir(
            "da_f", "gcu", "*.ens.xyz",
            lambda d, stem: d / f"cu2_{stem}.ens.xyz",
        )
        assert len(idx) >= 20
        self._assert_goat_table(idx, bad, got, want, "desc_da_f_goat_ens")


class TestTopologyTables:
    @pytest.mark.skipif(
        importlib.util.find_spec("rdkit") is None,
        reason="rdkit is not installed in this interpreter",
    )
    def test_all_four_reactions(self):
        pad = _scratch()
        for rxn in ("cc", "cp", "oa", "da_f"):
            preds = pd.read_csv(pad / "mc" / f"mc_preds_{rxn}.csv")
            preds["name"] = preds["name"].str.replace(".xyz", "", regex=False)
            smiles = dict(zip(preds["name"], preds["c_smiles"]))
            got = LigandTopology.table_from_smiles(smiles, normalized=True)
            want = pd.read_csv(pad / f"topo_{rxn}.csv", index_col=0)
            cols = [c for c in want.columns if c.startswith("tf_") and c in got.columns]
            idx = got.index.intersection(want.index)
            bad = {c: _max_abs(got.loc[idx, c], want.loc[idx, c]) for c in cols}
            bad = {c: d for c, d in bad.items() if d > 1e-10}
            assert not bad, f"topo_{rxn} mismatches: {bad}"
            assert len(idx) == len(preds)
