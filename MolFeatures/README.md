# MolFeatures

MolFeatures is a Python toolkit for molecular descriptor workflows. It helps convert quantum-chemistry outputs into tabular descriptors, assemble feature matrices, and search regression or classification models for structure-property analysis.

The project currently contains three main workflow layers:

- `M1_pre_calculations`: helpers for preparing and submitting molecular calculations.
- `M2_data_extractor`: feature extraction from Gaussian/log-derived feather files, cube files, and XYZ geometries.
- `M3_modeler`: regression and classification model search, validation, plotting, and reporting.

It also includes `MolAlign` utilities for molecular alignment and atom renumbering, plus practical notebooks under `Getting_started_with_examples`.

## Status

This repository is being prepared for publication. The source code is usable, but the package metadata, dependency boundaries, examples, and tests should be reviewed before a public release.

See `PUBLISHING_CHECKLIST.md` for the remaining release checklist.

Maintainer documentation is in `docs/`, especially `docs/ARCHITECTURE.md` and
`docs/CLI.md`.

## Installation

Create and activate a Python environment, then install the project in editable mode:

```bash
python -m pip install -U pip
python -m pip install -e .
```

For notebooks:

```bash
python -m pip install -e ".[notebooks]"
```

Some workflows also require external chemistry software that is not installed through pip, such as Gaussian and Open Babel.

The current dependency set pins `numpy<2` because common compiled chemistry/data packages in this project environment, including RDKit and PyArrow, may fail when installed against NumPy 2.x.

## Command-line usage

From the package directory:

```bash
python __main__.py --help
```

After installation:

```bash
molfeatures --help
```

Available commands:

```text
gui              Run the desktop GUI
model            Run regression or classification model search
extractor        Extract a feature set from a saved input JSON file
logs_to_feather  Convert Gaussian log files to feather files
cube             Calculate cube Sterimol descriptors
sterimol         Calculate Sterimol descriptors from XYZ files
```

## Feature extraction

Use the GUI to choose molecules, inspect structures, select atom indices, and save extraction inputs:

```bash
python __main__.py gui
```

If you already have a saved input JSON file, run the extractor directly:

```bash
python __main__.py extractor \
  --input Getting_started_with_examples/feather_example/input_example.json \
  --output feature_set \
  --feather_directory Getting_started_with_examples/feather_example
```

The extractor writes a CSV feature matrix and, when configured, correlation summaries.

## Modeling

Run model search from the command line:

```bash
python __main__.py model \
  --mode regression \
  --features_csv path/to/features.csv \
  --target_csv path/to/targets.csv \
  --y_value output \
  --min-features 1 \
  --max-features 4 \
  --top-n 20 \
  --threshold 0.70
```

Classification uses the same command with `--mode classification`.

## Python API

Feature extraction:

```python
from MolFeatures.M2_data_extractor.data_extractor import Molecules

molecules = Molecules("Getting_started_with_examples/feather_example")
features = molecules.get_molecules_features_set(
    entry_widgets={
        "Sterimol": "[[1, 6], [3, 4]]",
        "Bond-Angle": "[1, 2, 3]",
        "Dipole": "[[1, 2, 3], 5, 6]",
    },
    parameters={"Radii": "CPK", "Isotropic": True},
    save_as=True,
    csv_file_name="molecule_features",
)
```

Modeling:

```python
from MolFeatures.M3_modeler.modeling import LinearRegressionModel

csvs = {
    "features_csv_filepath": "features.csv",
    "target_csv_filepath": "targets.csv",
}

model = LinearRegressionModel(
    csvs,
    process_method="two csvs",
    y_value="yield",
    min_features_num=1,
    max_features_num=3,
)
model.search_models(top_n=10, threshold=0.6, bool_parallel=True)
```

## Examples

Start with:

- `Getting_started_with_examples/README.md`
- `Getting_started_with_examples/Practical_Notebook_Features.ipynb`
- `Getting_started_with_examples/Practical_Notebook_Modeling.ipynb`

Before publishing, review the examples and study cases for data size, licensing, and redistribution permissions.

## Development

Install development tools:

```bash
python -m pip install -e ".[dev]"
```

Run a lightweight syntax check:

```bash
python -m compileall .
```

Run CLI help:

```bash
python __main__.py --help
```

If CLI help fails before printing usage, check that the active environment satisfies `numpy<2`. The current entry point imports several heavy scientific packages before argument parsing; lazy-loading those imports is a recommended cleanup before release.

## License

This project currently includes an MIT license placeholder. Confirm the intended license before public release.
