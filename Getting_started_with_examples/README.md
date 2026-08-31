# Getting started

Example data ships with the repo — no Gaussian job, nothing to download.

```bash
pip install -e .
jupyter notebook Getting_started_with_examples/Practical_Notebook_Features.ipynb
```

---

## What's in here

| | |
|---|---|
| `Practical_Notebook_Features.ipynb` | Load a molecule set, pick atoms, extract descriptors |
| `Practical_Notebook_Modeling.ipynb` | Search models over a feature matrix, read the report |
| bundled benzene set | 26 substituted-benzene `.feather` files (`descripytor.examples.feather_example_dir()`) |
| bundled Baptiste set | 18 product `.feather` files plus outcomes (`descripytor.examples.baptiste_example_dir()`) |
| `modeling_example/` | Two ready-made datasets, regression and classification |
| `cube_example/` | Three electron-density cubes for cube Sterimol |
| `input_example.json` | A saved extraction input — replayable from the CLI |
| `descriptor_extraction_toolkit/` | The 3D atom picker, multi-engine extraction, web app |

---

## The example data

**Substituted benzenes** — one `.feather` per molecule (geometry, connectivity, charges,
dipoles, vibrations). 26 structures including `basic.feather`. Load with
`descripytor.examples.feather_example_dir()`. This is the GUI default (**Use example set**).

**Baptiste products** — 18 structures including `unsub.feather`, plus `outcomes.csv` for
modeling. Load with `descripytor.examples.baptiste_example_dir()`.

**`modeling_example/`** — feature matrices you can model immediately:

| File | Shape | Target |
|---|---|---|
| `Linear_Dataset_Example.csv` | 15 samples × 22 descriptors | `output`, continuous |
| `Logistic_Dataset_Example.csv` | 55 samples × 3 descriptors | `class`, categorical |

One CSV: first column is the sample name, last is the target. Fifteen samples is typical
here — that is why modeling uses Q² and Y-randomization, not R².

**`cube_example/`** — density cubes for `descripytor cube` (Sterimol from electron density).

---

## Three things to try

### 1. Extract descriptors without writing code

Start the 3D picker on the bundled benzene set:

```bash
descripytor visual
```

Open http://localhost:7432, click **Use example set**, pick atoms, then **Extract CSV**.
To replay a saved JSON on the same molecules:

```bash
descripytor extractor \
  --input Getting_started_with_examples/input_example.json \
  --output feature_set \
  --feather_directory descripytor/examples/feather_example
```

You get a CSV of vibrations, dipole, charges, Sterimol, bonds, and angles, plus a
correlation table.

Keys match extractors by name prefix (`Sterimol…`, `Bond_length…`, …). Values are 1-based
atom indices. An unknown key is printed, not ignored.

### 2. Model a dataset

```bash
descripytor model \
  --mode regression \
  --features_csv Getting_started_with_examples/modeling_example/Linear_Dataset_Example.csv \
  --y_value output \
  --min-features 1 --max-features 3 \
  --top-n 10 --threshold 0.6
```

Writes `runs/` (SQLite, figures, PDF with SHAP and sanity checks). For classification, use
`Logistic_Dataset_Example.csv` with `--mode classification` and `--y_value class`.

### 3. Pick atoms in 3D

```bash
cd Getting_started_with_examples/descriptor_extraction_toolkit
python make_picker.py --feather-dir ../../descripytor/examples/feather_example --index 0
```

Opens the picker on a molecule from the set. Click atoms, then **Extract CSV**. See
[descriptor_extraction_toolkit/README.md](descriptor_extraction_toolkit/README.md).

---

## From the Python side

```python
from M2_data_extractor.data_extractor import Molecules
from descripytor.examples import feather_example_dir

molset = Molecules(str(feather_example_dir()), threshold=1.82)
print(molset.success_molecules, molset.failed_molecules)

features = molset.get_molecules_features_set(
    entry_widgets={
        'Sterimol':    '[[1,4],[4,1]]',
        'Bond-Length': '[[1,2],[4,7]]',
        'Dipole':      '[[1,2,19], 20, 19]',
    },
    parameters={'Radii': 'CPK', 'Isotropic': True},
    save_as=True,
    csv_file_name='my_features',
)
```

Full API: [M2_data_extractor/README.md](../M2_data_extractor/README.md) for extraction,
[M3_modeler/README.md](../M3_modeler/README.md) for modelling.

---

## If something goes wrong

- **`failed_molecules` is not empty** — that file did not parse. Check the list before
  trusting the matrix.
- **Empty descriptor** — indices usually do not match what the extractor needs (a stretch
  pair must be bonded; a bend pair must share a centre). Indices are 1-based.
- **Wrong bonds** — raise or lower `threshold` in `Molecules` (default 1.82 Å).
- **Import errors** — `pip install -e .` from the repo root; keep `numpy<2`.

A guided walkthrough with more context:
[DescriPyTor_Tutorial.md](descriptor_extraction_toolkit/DescriPyTor_Tutorial.md).
