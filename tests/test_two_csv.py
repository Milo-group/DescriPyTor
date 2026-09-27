"""`descripytor model -f features.csv -t target.csv`: features and target in separate files."""
import numpy as np
import pandas as pd
import pytest

from M3_modeler.modeling import merge_features_and_target


def _files(tmp_path):
    rng = np.random.default_rng(0)
    names = [f"mol{i}" for i in range(12)]
    X = rng.normal(size=(12, 3))
    feats = pd.DataFrame(X, columns=["a", "b", "c"], index=names)
    feats.to_csv(tmp_path / "features.csv")                     # names as an unnamed first column, like the extractor
    y = pd.DataFrame({"name": names[::-1] + ["extra"], "output": list(2 * X[::-1, 0] + 1) + [0.0]})
    y.to_csv(tmp_path / "outcomes.csv", index=False)             # other order, one molecule without features
    return tmp_path / "features.csv", tmp_path / "outcomes.csv", X


def test_merge_joins_on_the_name_not_the_row_order(tmp_path):
    f, t, X = _files(tmp_path)
    m = merge_features_and_target(f, t, "output").set_index("name")
    assert len(m) == 12
    assert np.allclose(m.loc["mol3", "output"], 2 * X[3, 0] + 1)


def test_missing_target_column_is_a_clear_error(tmp_path):
    f, t, _ = _files(tmp_path)
    with pytest.raises(KeyError, match="yield"):
        merge_features_and_target(f, t, "yield")


def test_linear_model_loads_two_files(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    from M3_modeler.modeling import LinearRegressionModel
    f, t, X = _files(tmp_path)
    model = LinearRegressionModel({"features_csv_filepath": str(f), "target_csv_filepath": str(t)},
                                  process_method="two csvs", y_value="output",
                                  min_features_num=1, max_features_num=1)
    assert len(model.target_vector) == 12
    assert set(model.features_df.columns) == {"a", "b", "c"}
