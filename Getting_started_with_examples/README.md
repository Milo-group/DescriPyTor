# Getting started

Everything here runs against data that ships with the repository — no Gaussian job, no
cluster, nothing to download. Start a notebook and you have output in a few minutes.

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
| `feather_example/` | 26 substituted benzenes, ready to extract from |
| `modeling_example/` | Two ready-made datasets, regression and classification |
| `cube_example/` | Three electron-density cubes for cube Sterimol |
| `input_example.json` | A saved extraction input — replayable from the CLI |
| `descriptor_extraction_toolkit/` | The 3D atom picker, multi-engine extraction, web app |

---

## The example data

**`feather_example/`** — one `.feather` per molecule (geometry, connectivity, charges,
dipoles, vibrations). 26 substituted benzenes: `basic` plus ortho/meta/para halides and a
spread of para substituents. The same files are installed with the package as
`descripytor.examples.feather_example_dir()`.

**`modeling_example/`** — feature matrices you can model immediately:

| File | Shape | Target |
|---|---|---|
| `Linear_Dataset_Example.csv` | 15 samples × 22 descriptors | `output`, continuous |
| `Logistic_Dataset_Example.csv` | 55 samples × 3 descriptors | `class`, categorical |

Both are in **one-CSV** form: first column is the sample name, last is the target. Fifteen
samples is realistic for this field, and it is exactly why the modelling side leans on
out-of-fold Q² and Y-randomization rather than R².

**`cube_example/`** — `Ad_1_a`, `Bn_1_a`, `Cy_1_a` density cubes (~6.5 MB each) for
`descripytor cube`, where Sterimol radii come from the electron density instead of a
lookup table.

---

## Three things to try

### 1. Extract descriptors without writing code

`input_example.json` records a full set of atom selections. Replay it against the example
molecules:

```bash
descripytor extractor \
  --input Getting_started_with_examples/input_example.json \
  --output feature_set \
  --feather_directory Getting_started_with_examples/feather_example
```

You get a CSV with ring and bending vibrations, dipole components, charges, charge
differences, Sterimol, bond lengths and angles — plus a correlation table.

The file is hand-editable. Each key selects an extractor by the descriptor name it starts
with — `Sterimol…`, `Bond_length…`, `Charge difference…` — and the value is the atom indices,
1-based. Change the numbers, rerun, compare. A key that matches no descriptor is reported on
the console rather than ignored, so a typo can't quietly cost you a column.

### 2. Model a dataset

```bash
descripytor model \
  --mode regression \
  --features_csv Getting_started_with_examples/modeling_example/Linear_Dataset_Example.csv \
  --y_value output \
  --min-features 1 --max-features 3 \
  --top-n 10 --threshold 0.6
```

Writes a `runs/` folder with the scored models in SQLite, the figures, and a PDF report
carrying the SHAP and sanity-check pages. Swap in `Logistic_Dataset_Example.csv` with
`--mode classification` and `--y_value class`.

### 3. Pick atoms in 3D

```bash
cd Getting_started_with_examples/descriptor_extraction_toolkit
python make_picker.py --feather-dir ../feather_example --index 0
```

Opens the browser picker with a real molecule from the set. Click atoms, watch the Sterimol
vectors and dipole update live, then export a run config or a `answers_dict` to paste into a
notebook. See [descriptor_extraction_toolkit/README.md](descriptor_extraction_toolkit/README.md).

---

## From the Python side

```python
from M2_data_extractor.data_extractor import Molecules

molset = Molecules("Getting_started_with_examples/feather_example", threshold=1.82)
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

- **`failed_molecules` is not empty** — that file didn't parse. It is skipped silently in
  batch extractors, so check the list before trusting a feature matrix.
- **A descriptor came back empty** — usually the atom indices don't describe what the
  extractor needs (a stretch pair must actually be bonded; a bend pair must share a centre).
  Visualise the molecule and re-read the indices; they are 1-based.
- **A bond looks wrong** — adjust `threshold` in the `Molecules` constructor. The default
  1.82 Å suits organic molecules; longer bonds need more.
- **Import errors on a fresh environment** — `pip install -e .` from the repository root,
  and check `numpy<2` was honoured.

A guided walkthrough with more context:
[DescriPyTor_Tutorial.md](descriptor_extraction_toolkit/DescriPyTor_Tutorial.md).
