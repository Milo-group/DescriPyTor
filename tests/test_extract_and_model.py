"""End-to-end extraction and modeling on the small test set."""
from __future__ import annotations

import importlib.util
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

SMALL = ROOT / "tests" / "data" / "small_set"
XYZ_DIR = SMALL / "xyz"


def _load_extract():
    spec = importlib.util.spec_from_file_location(
        "small_set_extract", SMALL / "extract.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


_extract = _load_extract()
extract_table = _extract.extract_table
extract_tiny_ni = _extract.extract_tiny_ni


@pytest.fixture(scope="module")
def alcohol_table():
    if not (XYZ_DIR / "methanol.xyz").is_file():
        pytest.skip("small_set XYZ files missing; run tests/data/small_set/build_xyz.py")
    return extract_table(with_topology=False)


class TestAlcoholExtraction:
    def test_all_eight_load(self, alcohol_table):
        assert len(alcohol_table) == 8
        assert set(alcohol_table.index) == {
            "methanol", "ethanol", "n_propanol", "isopropanol",
            "n_butanol", "isobutanol", "sec_butanol", "tert_butanol",
        }

    def test_atom_order_and_bonds(self, alcohol_table):
        for name, row in alcohol_table.iterrows():
            assert 1.35 < row["CO"] < 1.50
            assert 90 < row["HOC"] < 125
            assert np.isfinite(row["B1"]) and np.isfinite(row["B5"]) and np.isfinite(row["L"])
            assert row["B5"] + 1e-6 >= row["B1"]

    def test_sterimol_trends(self, alcohol_table):
        # Chain length stretches L; branching widens B5.
        assert alcohol_table.loc["methanol", "L"] < alcohol_table.loc["n_butanol", "L"]
        assert alcohol_table.loc["ethanol", "L"] < alcohol_table.loc["n_propanol", "L"]
        assert alcohol_table.loc["tert_butanol", "B5"] > alcohol_table.loc["methanol", "B5"]
        assert alcohol_table.loc["n_butanol", "n_heavy"] == 5
        assert alcohol_table.loc["methanol", "n_heavy"] == 2

    def test_committed_table_matches_live_extract(self, alcohol_table):
        committed = SMALL / "modeling_table.csv"
        if not committed.is_file():
            pytest.skip("modeling_table.csv not written yet")
        want = pd.read_csv(committed, index_col="name")
        cols = [c for c in ("B1", "B5", "L", "CO", "HOC", "vbur_3.5") if c in want.columns]
        idx = alcohol_table.index.intersection(want.index)
        for col in cols:
            delta = float(np.nanmax(np.abs(
                alcohol_table.loc[idx, col].astype(float) - want.loc[idx, col].astype(float)
            )))
            assert delta < 1e-6, f"{col} drifted by {delta}"


class TestTinyNickel:
    def test_geometric_features(self):
        if not (XYZ_DIR / "tiny_nih2n2.xyz").is_file():
            pytest.skip("tiny_nih2n2.xyz missing")
        geom = extract_tiny_ni()
        assert geom["n_atoms"] == 5
        assert geom["bite"] == pytest.approx(180.0, abs=0.5)
        assert geom["MD_mean"] == pytest.approx(1.90, abs=0.02)
        assert geom["MD_asym"] == pytest.approx(0.0, abs=1e-6)
        assert 0.0 < geom["vbur_3.5"] < 100.0
        assert np.isfinite(geom["fromM_L_sym"])


class TestTopology:
    @pytest.mark.skipif(
        importlib.util.find_spec("rdkit") is None,
        reason="rdkit is not installed in this interpreter",
    )
    def test_tf_columns_scale_with_size(self):
        table = extract_table(with_topology=True)
        assert "tf_zagreb1" in table.columns
        assert table.loc["n_butanol", "tf_zagreb1"] != table.loc["methanol", "tf_zagreb1"]


class TestRegressionModeling:
    def test_ridge_loo_on_extracted_sterimol(self, alcohol_table, tmp_path):
        try:
            from M3_modeler.modeling import LinearRegressionModel, _analytic_loo_linear
        except Exception as exc:
            pytest.skip(f"M3_modeler import failed: {exc}")

        df = alcohol_table.reset_index()[["name", "B1", "B5", "L", "n_heavy", "bp_c"]]
        csv_path = tmp_path / "alcohols.csv"
        df.to_csv(csv_path, index=False)

        # In-sample + hat-matrix LOO without writing the search DB.
        X = df[["L", "n_heavy"]].to_numpy(float)
        y = df["bp_c"].to_numpy(float)
        q2, mae, rmsd = _analytic_loo_linear(X, y, alpha=1.0)
        assert np.isfinite(q2) and np.isfinite(mae) and np.isfinite(rmsd)
        assert mae < 40.0  # °C; n=8 toy geometries, not a paper model

        prev = os.getcwd()
        os.chdir(tmp_path)
        try:
            model = LinearRegressionModel(
                {"features_csv_filepath": str(csv_path)},
                process_method="one csv",
                y_value="bp_c",
                names_column="name",
                min_features_num=1,
                max_features_num=2,
                db_path=str(tmp_path / "alcohols.db"),
                n_splits=8,
            )
            X_fit = model.features_df[["L", "n_heavy"]].to_numpy(float)
            y_fit = model.target_vector.to_numpy(float)
            model.fit(X_fit, y_fit, alpha=1.0)
            scores, pred = model.evaluate(X_fit, y_fit)
            assert scores["r2"] > 0.5
            assert pred.shape == y_fit.shape
        finally:
            os.chdir(prev)

    def test_search_two_features(self, alcohol_table, tmp_path):
        try:
            from M3_modeler.modeling import LinearRegressionModel
        except Exception as exc:
            pytest.skip(f"M3_modeler import failed: {exc}")

        df = alcohol_table.reset_index()[["name", "B1", "B5", "L", "n_heavy", "bp_c"]]
        csv_path = tmp_path / "alcohols.csv"
        df.to_csv(csv_path, index=False)
        prev = os.getcwd()
        os.chdir(tmp_path)
        try:
            model = LinearRegressionModel(
                {"features_csv_filepath": str(csv_path)},
                process_method="one csv",
                y_value="bp_c",
                names_column="name",
                min_features_num=1,
                max_features_num=2,
                db_path=str(tmp_path / "search.db"),
            )
            results = model.search_models(
                top_n=5,
                n_jobs=1,
                threshold=0.0,
                bool_parallel=False,
                min_models_to_keep=1,
                min_threshold=0.0,
            )
            assert results is not None and len(results) >= 1
        finally:
            os.chdir(prev)

    def test_linear_example_csv(self, tmp_path):
        example = ROOT / "Getting_started_with_examples" / "modeling_example" / "Linear_Dataset_Example.csv"
        if not example.is_file():
            pytest.skip("Linear_Dataset_Example.csv missing")
        try:
            from M3_modeler.modeling import LinearRegressionModel
        except Exception as exc:
            pytest.skip(f"M3_modeler import failed: {exc}")

        raw = pd.read_csv(example)
        # Keep a handful of columns so the search is tiny.
        name_col = raw.columns[0]
        keep = [name_col, "B1", "B5", "L", "Dist_1-2", "output"]
        slim = raw[keep].head(10)
        csv_path = tmp_path / "linear_example_slim.csv"
        slim.to_csv(csv_path, index=False)
        prev = os.getcwd()
        os.chdir(tmp_path)
        try:
            model = LinearRegressionModel(
                {"features_csv_filepath": str(csv_path)},
                process_method="one csv",
                y_value="output",
                min_features_num=1,
                max_features_num=2,
                db_path=str(tmp_path / "linear_example.db"),
            )
            X = model.features_df[["B1", "L"]].to_numpy(float)
            y = model.target_vector.to_numpy(float)
            model.fit(X, y, alpha=1.0)
            scores, _ = model.evaluate(X, y)
            assert np.isfinite(scores["r2"])
        finally:
            os.chdir(prev)


class TestClassificationModeling:
    def test_fit_toy_labels(self, tmp_path):
        toy = SMALL / "classification_toy.csv"
        try:
            from M3_modeler.modeling import ClassificationModel
        except Exception as exc:
            pytest.skip(f"M3_modeler import failed: {exc}")

        prev = os.getcwd()
        os.chdir(tmp_path)
        try:
            model = ClassificationModel(
                {"features_csv_filepath": str(toy)},
                process_method="one csv",
                y_value="class",
                names_column="name",
                min_features_num=1,
                max_features_num=1,
                db_path=str(tmp_path / "cls.db"),
            )
            X = model.features_df.to_numpy(float)
            y = model.target_vector.to_numpy()
            model.fit(X, y)
            scores, pred = model.evaluate(X, y)
            assert pred.shape[0] == len(y)
            acc = scores.get("accuracy", scores.get("Accuracy"))
            if acc is not None:
                assert 0.0 <= float(acc) <= 1.0
        finally:
            os.chdir(prev)


class TestGoatFixtureStillExtracts:
    def test_081_ensemble(self):
        path = ROOT / "tests" / "fixtures" / "081_lig.finalensemble.xyz"
        if not path.is_file():
            pytest.skip("081_lig fixture missing")
        from M2_data_extractor.metal_complex import MetalComplexEnsemble

        ens = MetalComplexEnsemble.from_xyz(path)
        geom = ens.geometric_features()
        assert ens.n_conformers == 2
        assert geom["n_conformers"] == 2
        assert 70 < geom["bite"] < 100
        assert np.isfinite(geom["fromM_L_sym"])
