# DescriPyTor

**Chemical-intuition-based molecular feature extraction and modeling, for computational chemists.**

DescriPyTor turns quantum-chemistry output into descriptors *you choose* — pick the atoms
that matter mechanistically, get a model-ready feature matrix, then search for the linear or
classification model that explains your measured outcome.

Modeled after the R package *MoleculaR* ([docs](https://barkais.github.io/)).

<p align="center">
  <img src="docs/images/pipeline-overview.png" width="620" alt="From SMILES or structure, through conformer search and DFT, to a feature matrix and a cross-validated model."><br>
  <em>The full pipeline: structure in, descriptors out, model at the end.</em>
</p>

---

## Contents

- [Why it looks like this](#why-it-looks-like-this)
- [Install](#install)
- [The three ways to use it](#the-three-ways-to-use-it)
- [Feature extraction](#feature-extraction)
  - [Descriptor families](#descriptor-families)
- [Modeling](#modeling)
- [Command-line reference](#command-line-reference)
- [Molecular alignment and renumbering](#molecular-alignment-and-renumbering)
- [Docker](#docker)
- [Examples and further reading](#examples-and-further-reading)
- [Repository layout](#repository-layout)

---

## Why it looks like this

Generic descriptor packages hand you thousands of columns and leave interpretation to the
model. DescriPyTor goes the other way: you cut the molecule where the chemistry happens,
put every structure in a common reference frame, and read off a handful of descriptors
whose meaning you already understand.

<p align="center">
  <img src="docs/images/feature-concept-blackboard.png" width="720" alt="Three panels: trimming a symmetric molecule to its unique fragment, aligning conformers, and labelling the chosen descriptors (Sterimol, dipole, charge difference, angle)."><br>
  <em>Trim by symmetry → align to a shared frame → pick the descriptors that mean something.</em>
</p>

Alignment is what makes descriptors comparable across a series. Every molecule is rotated
and translated onto the same origin and axes, so a "dipole along x" means the same thing in
every row of your table.

<p align="center">
  <img src="docs/images/aligned-common-frame.png" width="680" alt="Three substituted benzenes each with local axes, rotated and translated onto one common ring-centered frame."><br>
  <em>Ring center or nuclear-charge center as origin; substituent direction fixes the axes.</em>
</p>

You give it three selections and it builds the frame by Gram–Schmidt: the origin is a
centroid, `ŷ` points at your y-atom, and `x̂` is whatever is left of the plane atom once its
`ŷ` component is subtracted off.

<p align="center">
  <img src="docs/animations/frame.svg" width="760" alt="Animation: the origin is set to the centroid of the chosen atoms, y is normalized toward the y-atom, and x is obtained by removing the y-component from the plane-atom vector.">
</p>

That is `_build_basis` in [`dipole_utils.py`](M2_data_extractor/extractor_utils/dipole_utils.py),
and every descriptor with a direction in it — dipole components, NPA, ring vibrations — is
read off in that frame.

---

## Install

```bash
git clone https://github.com/Milo-group/DescriPyTor.git
cd DescriPyTor
```

Create an environment (Python 3.9–3.11) and install:

```bash
conda create -n descripytor python=3.10 && conda activate descripytor
pip install -r requirements.txt
```

`numpy<2` is pinned: RDKit and PyArrow builds commonly used here fail against NumPy 2.x.

Optional extras, only if you need the feature they unlock:

| Install | Unlocks |
|---|---|
| `torch` | MolAlign automatic atom renumbering (see note below) |
| `streamlit umap-learn` | the Streamlit extraction web app |
| `mordred deepchem ase` | extra descriptor engines in the toolkit |
| `aqme autoqchem` (+ xTB) | xTB / Gaussian-log descriptor engines |

For `torch`, pick the build that matches your machine:

```bash
# with an NVIDIA GPU
pip install torch --index-url https://download.pytorch.org/whl/cu121
# CPU only
pip install torch --index-url https://download.pytorch.org/whl/cpu
```

Gaussian and Open Babel are external programs and are not installed by pip.

---

## The three ways to use it

| | Best for | Start with |
|---|---|---|
| **Desktop GUI** | Loading a feather set, clicking through descriptor prompts, saving a reusable input file | `python __main__.py gui` |
| **3D atom picker** | Choosing atoms visually, live Sterimol/dipole/vibration overlays, driving extraction + modeling from one page | [`Getting_started_with_examples/descriptor_extraction_toolkit/`](Getting_started_with_examples/descriptor_extraction_toolkit/README.md) |
| **CLI / notebooks** | Reproducible batch runs, scripting, HPC | `python __main__.py --help` |

---

## Feature extraction

### Desktop GUI

```bash
python __main__.py gui
```

<p align="center">
  <img src="docs/images/gui_main.png" width="760" alt="Molecule Data Extractor main window: buttons for browsing a feather directory, filtering molecules, extracting features, modeling, and exporting DataFrames or xyz files."><br>
  <em>The main window. Every action logs to the right-hand pane, which you can save as text.</em>
</p>

**Load molecules.** *Browse for Feather Files Directory* reads every `.feather` in a folder
and reports what loaded:

```text
Molecules initialized : ['basic', 'm_Br', 'm_Cl', 'm_F', ... , 'p_tfm']
Failed to load Molecules: []
```

If you only have Gaussian logs, convert them first with **File Handler → Log to Feather**,
or from the command line:

<p align="center">
  <img src="docs/images/logs_to_feather.jpg" width="820" alt="Terminal running the logs_to_feather subcommand: it prompts for a log directory and reports the feather file it saved."><br>
  <em>One <code>.feather</code> per molecule — geometry, charges, dipoles and vibrations in a single file.</em>
</p>

**Look before you pick.** *Visualize Molecules* opens an interactive viewer with atom
indices, bond lengths and the dipole vector — this is how you find the index numbers the
extractor asks for.

<p align="center">
  <img src="docs/images/molvisualizer.png" width="700" alt="Interactive molecule viewer with numbered atoms and a side menu toggling atom indices, bond lengths, and the dipole vector."><br>
  <em>Atom indices are 1-based throughout, matching Gaussian.</em>
</p>

**Fill in the descriptors you want.** *Extract Features* opens one prompt per descriptor
family, each with an example of the expected input. Leave a field blank to skip it.

<p align="center">
  <img src="docs/images/gui_questions.png" width="700" alt="Feature extraction window listing every descriptor family with an input box and worked example: ring vibration, stretching, bending, dipole, NPA, charges, Sterimol, bond length, bond angle."><br>
  <em>Only the fields you fill get computed.</em>
</p>

*Choose Parameters* sets the radii system used for Sterimol — **Pyykkö** covalent radii
(defined for every element), **bondi**, or **CPK/VDW** (a subset of elements) — and whether
to append isotropic polarizability and energy.

<p align="center">
  <img src="docs/images/feature_extraction_choose_parameters.png" width="640" alt="The Parameters popup over the feature extraction window, with a radii dropdown offering bondi, CPK and Pyykko, and a True/False dropdown for isotropic values."><br>
  <em>Radii system and isotropic toggle, applied to every molecule in the set.</em>
</p>

**Save input** writes your indices to a JSON file; **Load input** replays it, so a run is
reproducible and a second dataset costs no re-typing. **Submit** prints the feature matrix
to the dashboard and writes a `.csv`.

<p align="center">
  <img src="docs/images/input_example.png" width="520" alt="A saved input JSON file mapping each descriptor prompt to the atom indices chosen."><br>
  <em>A saved input file — hand-editable. An example ships in <code>Getting_started_with_examples/feather_example/</code>.</em>
</p>

Already have an input file? Skip the GUI:

```bash
python __main__.py extractor \
  --input Getting_started_with_examples/feather_example/input_example.json \
  --output feature_set \
  --feather_directory Getting_started_with_examples/feather_example
```

### 3D atom picker

The browser-based picker is the richer path: click atoms in 3D, watch the descriptor
you're defining get drawn live, then export a ready-to-run config.

<p align="center">
  <img src="docs/images/atom-picker.png" width="900" alt="Browser atom picker: a 3D molecule with Sterimol L/B1/B5 arrows, dipole vector and transformation axes on the left; a scrollable list of descriptor fields with committed atom groups on the right."><br>
  <em>Pick a field, click atoms, see the vectors. Pairs and triplets auto-commit.</em>
</p>

It also overlays the molecular dipole, animates normal modes with a playable amplitude
slider, and stacks conformer ensembles with per-conformer RMSD:

<p align="center">
  <img src="docs/images/conformer-viewer.png" width="900" alt="Conformer viewer showing six overlaid conformers in different colors with energies and RMSD values against the reference conformer."><br>
  <em>Kabsch-aligned conformer overlay, aligned on a chosen substructure.</em>
</p>

### Data prep and pre-flight check

Before a long run, `data_gathering_validation.html` assembles the whole command sequence for
you. Pick a starting point — existing feathers, Gaussian logs, or a SMILES CSV — fill in the
paths, and it writes out the exact commands to copy, including a dependency check:

```bash
python _dx_check.py --config run_config.json --install-report
```

That reports which optional engines are actually importable in your environment, so a
six-hour extraction doesn't die on a missing package at hour five.

Full instructions: [descriptor_extraction_toolkit/README.md](Getting_started_with_examples/descriptor_extraction_toolkit/README.md).

### Beyond the built-in descriptors

The toolkit can run external descriptor packages in the same pass and merge everything into
one CSV, each engine namespaced by a `prefix` so columns never collide. Every engine is
optional and skips itself with a one-line message if its package is missing, so a run never
dies on a dependency:

| Engine | Gives you | Needs |
|---|---|---|
| `descripytor_full` / `_steric` | the built-in descriptor set above | — |
| `xyz_sterimol`, `xyz_geometry`, `xyz_buried_volume` | Sterimol, angles, bond lengths and buried volume from bare `.xyz` | — |
| `morfeus_suite` | Sterimol, buried volume, cone angle, SASA, dispersion, pyramidalization | `morfeus-ml` |
| `rdkit`, `rdkit_fp` | full 2D descriptor list; ECFP/FCFP/MACCS/AtomPair/Torsion + USRCAT | `rdkit` |
| `mordred` | Mordred 2D/3D descriptors | `mordred` |
| `deepchem` | CircularFingerprint + RDKitDescriptors | `deepchem` |
| `qm` | Gaussian-log descriptors via autoqchem | `autoqchem` |
| `aqme_qdescp` | xTB/MORFEUS table: IP, EA, HOMO–LUMO, FOD, dispersion, … | `aqme` + xTB |

```bash
python descriptor_extractor.py --config my_run.json
```

**Starting from SMILES.** AQME CSEARCH turns a SMILES CSV into a conformer ensemble, so the
pipeline can begin without any structures at all:

```bash
python descriptor_extractor.py --csearch molecules.csv --csearch-out xyz_out \
    --name-col name --smiles-col smiles --csearch-program rdkit \
    --csearch-sample 10 --keep lowest
```

### Python API

```python
from M2_data_extractor.data_extractor import Molecules

molset = Molecules("Getting_started_with_examples/feather_example", threshold=1.82)

features = molset.get_molecules_features_set(
    entry_widgets={
        "Sterimol":    "[[1, 6], [3, 4]]",
        "Stretching":  "[[1, 2], [3, 4]]",
        "Bond-Angle":  "[1, 2, 3]",
        "Dipole":      "[[1, 2, 3], 5, 6]",
    },
    parameters={"Radii": "CPK", "Isotropic": True},
    save_as=True,
    csv_file_name="molecule_features",
)
```

Full class reference: [M2_data_extractor/README.md](M2_data_extractor/README.md).

### Descriptor families

| Family | What you give it | What you get |
|---|---|---|
| **Sterimol** | `[origin, attached]` atom pair | `L`, `B1`, `B5` — steric length and the narrow/wide widths |
| **Cube Sterimol** | a density `.cube` file | Sterimol from the real electron density, so radii respond to stereoelectronics |
| **Charges** | single atoms | NBO / Hirshfeld / CM5 / NPA values per atom |
| **Charge difference** | atom pair `[a, b]` | `q_a − q_b`, per charge scheme |
| **Dipole** | `[origin(s), y-axis, xy-plane]` | dipole components in your transformed frame, plus total |
| **NPA** | base trio (+ optional sub-atoms) | natural population analysis in a local frame |
| **Bond length / angle** | pairs, triads, or quartets | distances, angles, dihedrals |
| **Stretching vibration** | a bonded pair | frequency and amplitude of the mode aligned with that bond |
| **Bending vibration** | two atoms sharing a center | frequency and amplitude of the strongest bending mode |
| **Ring vibration** | one ring atom | `cross` / `para` frequencies and angles for a benzene-like ring |
| **Buried volume, cone angle, SASA, dispersion** | metal/apex atom | via the Morfeus suite in the toolkit |

**Sterimol, in pictures.** `L` runs along the bond axis; `B1` and `B5` are the smallest and
largest perpendicular extents of the van der Waals envelope.

<p align="center">
  <img src="docs/images/sterimol-3d-view.png" width="620" alt="A molecule with the Sterimol L axis drawn as a black arrow along the primary bond, B1 in blue to the nearest atom and B5 in red to the furthest, with a swept envelope."><br>
  <em>Side view: L along the axis, B1 and B5 perpendicular to it.</em>
</p>

`B5` is one number — the furthest any atom reaches from the axis, plus its radius. `B1` is
a search: rotate a supporting line all the way around the cross-section and keep the angle
where it sits closest.

<p align="center">
  <img src="docs/animations/sterimol.svg" width="760" alt="Animation: a supporting line sweeps 360 degrees around the substituent cross-section while a plot traces its distance from the axis; the minimum of that curve is B1.">
</p>

The traced curve is exactly what `scan_b1_over_angles` tabulates in
[`sterimol_utils.py`](M2_data_extractor/extractor_utils/sterimol_utils.py) — B1 is its
minimum, which is why B1 and B5 are rarely perpendicular.

You can measure against the whole molecule's envelope, or against just the substructure
hanging off your chosen bond — the two answers differ, and the choice is yours
(`sub_structure`, `drop_atoms`).

<p align="center">
  <img src="docs/images/sterimol-end-on-global.png" width="400" alt="End-on Sterimol view using the global projection: the substructure atoms are highlighted green against the full molecule's grey envelope, B1 = 2.08 A."> <img src="docs/images/sterimol-end-on-local.png" width="400" alt="End-on Sterimol view using local slices: only the slices at the B1 and B5 heights are drawn, giving B1 = 1.68 A."><br>
  <em>Left: global projection over the full envelope (B1 = 2.08 Å). Right: local slices at the relevant heights (B1 = 1.68 Å).</em>
</p>

**Picking the right vibration.** A Gaussian frequency job gives you every normal mode; only
one of them is *your* bond stretching. Each mode is scored by how much of its displacement
lies along the bond, and the best scorer inside a frequency window wins.

<p align="center">
  <img src="docs/animations/vibration.svg" width="760" alt="Animation: each normal mode's displacement vectors are projected onto the bond axis, giving a score per mode; a cursor scans the modes and the highest-scoring one inside the frequency window is selected.">
</p>

The score is `|dᵃ·û| + |dᵇ·û|` — `calc_vibration_dot_product` in
[`vibrations_utils.py`](M2_data_extractor/extractor_utils/vibrations_utils.py). The window
matters: a C–H stretch at 3055 cm⁻¹ can out-score your carbonyl if you let it.

**Ring positions.** Ring vibration descriptors need one ring atom; the primary/ortho/meta/para
pattern is resolved for you.

<p align="center">
  <img src="docs/images/rings.png" width="640" alt="Two numbered aromatic rings with primary, ortho, meta and para positions color-coded and labelled."><br>
  <em>Give it the primary atom; it finds the rest.</em>
</p>

---

## Modeling

Model search is exhaustive over feature subsets, scored with leakage-free cross-validation:
scaling is fit inside each fold, never on the full set.

```bash
python __main__.py model \
  --mode regression \
  --features_csv path/to/features.csv \
  --target_csv path/to/targets.csv \
  --y_value output \
  --min-features 1 --max-features 4 \
  --top-n 20 --threshold 0.70 \
  --bool-parallel --n_jobs -1
```

<p align="center">
  <img src="docs/animations/crossval.svg" width="760" alt="Animation: five folds each hold out a different fifth of the samples; the held-out predictions accumulate into one out-of-fold vector, which Q-squared is then computed from.">
</p>

With a few dozen molecules and thousands of candidate subsets, a model can fit beautifully
and mean nothing. Q² is computed only from predictions made by models that never saw the
sample in question — which is why it is the number worth reading.

- **Regression** — R², Q² (from out-of-fold predictions), MAE, RMSD; prediction intervals;
  VIF multicollinearity checks.
- **Classification** — accuracy, F1, McFadden's R²; stratified K-fold; class balancing by
  stratified or similarity-based sampling.
- Results persist to SQLite per dataset, so runs are incremental and reproducible.
- The top 5 models are written to a PDF report automatically.

Candidate models are ranked so you can see which descriptors keep earning their place:

<p align="center">
  <img src="docs/images/bassa-plot.jpeg" width="520" alt="Bar chart of three candidate models with R-squared 0.925, 0.918 and 0.901, and a dot matrix below showing which of features s1-s4 each model uses."><br>
  <em>Model R² against the feature subset used — a two-feature model within noise of the best.</em>
</p>

Every top model gets a PDF report, and these pages are generated for it automatically — no
extra flags:

```text
runs/<dataset>_<target>_<date>/
  db/    results_<dataset>.csv, .db     every model scored, for incremental reruns
  figs/  coefficients, parity, SHAP, sanity checks, diagnostics
  pdf/   <dataset>_top_models_report.pdf
  logs/  regression_results.txt
```

**Sanity checks that try to break your model.** The report refits the model against
deliberately corrupted data — targets shuffled (Y-randomization), descriptors shuffled
globally and one feature at a time, and a one-hot baseline that knows only substituent
identity. A real model has to beat all of them.

<p align="center">
  <img src="docs/images/sanity-checks-yscramble.png" width="640" alt="Histogram of RMSD from models fitted to randomly shuffled targets, clustered near 0.50, with the real model's RMSD marked at 0.351 far to the left, plus one-hot and X-shuffle baselines."><br>
  <em>The real model (red, 0.351) sits clear of the Y-randomized distribution (~0.50). If it landed inside that histogram, the model is fitting noise.</em>
</p>

**Feature attribution.** SHAP values show which descriptor pushed which sample, and in which
direction:

<p align="center">
  <img src="docs/images/shap-beeswarm.png" width="760" alt="SHAP beeswarm for the top three features: total dipole, dipole z and stretch amplitude, with individual substituents labelled at the extremes of each row."><br>
  <em>Named outliers make it obvious which substituents drive the model.</em>
</p>

Each model can also be decomposed per sample, so a prediction is auditable rather than opaque:

<p align="center">
  <img src="docs/images/model-components-chart.png" width="900" alt="Stacked contribution chart across 18 substituents, showing how the dipole, CM5 charge and O-C bond length terms each push the prediction above or below the intercept, with measured values as open circles and predictions as diamonds."><br>
  <em>Per-substituent breakdown: which descriptor moved which prediction, and by how much.</em>
</p>

Bayesian model averaging (BASSA, spike-and-slab priors) is available through the toolkit GUI
for feature-inclusion probabilities rather than a single winning subset:

<p align="center">
  <img src="docs/images/bassa-markov-chain.png" width="760" alt="MCMC trace showing which of ten features are included at each of 8000 iterations; two features are included almost always, the rest sporadically."><br>
  <em>MCMC inclusion trace — s1 and s2 are in nearly every sampled model.</em>
</p>

Python API and full reference: [M3_modeler/README.md](M3_modeler/README.md).

---

## Command-line reference

```text
python __main__.py {gui,model,extractor,logs_to_feather,cube,sterimol}

  gui               Desktop GUI
  model             Regression or classification model search
  extractor         Extract a feature set from a saved input JSON
  logs_to_feather   Convert Gaussian .log files to .feather
  cube              Sterimol from density cube files
  sterimol          Sterimol from plain .xyz files
```

Every subcommand takes `-h`.

**Sterimol from bare XYZ** — no logs or feather files needed:

```bash
python __main__.py sterimol
```

<p align="center">
  <img src="docs/images/sterimol_cmd.jpg" width="820" alt="Terminal session running the sterimol subcommand and printing a table of L, B1 and B5 values per molecule."><br>
</p>

**Cube Sterimol** — same idea, but radii come from the electron density, so they respond to
the electronic environment instead of being read off a table:

<p align="center">
  <img src="docs/images/cube_sterimol.jpg" width="820" alt="Terminal session running the cube subcommand over density cube files and printing the resulting Sterimol values."><br>
</p>

**Model search arguments**

| Flag | Type | Required | Description |
|---|---|:--:|---|
| `-m`, `--mode` | `regression` \| `classification` | **yes** | Task; sets metrics and thresholds |
| `-f`, `--features_csv` | path | **yes** | Descriptor CSV |
| `-t`, `--target_csv` | path | **yes** | Target/label CSV |
| `-y`, `--y_value` | str | **yes** | Target column name, also used to name outputs |
| `-j`, `--n_jobs` | int | no | Cores; `-1` = all. Defaults to `$NSLOTS` or all cores |
| `--min-features` / `--max-features` | int | no | Bounds on subset size |
| `--top-n` | int | no | How many top models to keep |
| `--bool-parallel` | flag | no | Parallel evaluation |
| `--threshold` | float | no | Quality cutoff: min R² (regression) or McFadden's R² (classification) |
| `--leave-out` | list | no | Sample IDs held out entirely, for external validation |

---

## Molecular alignment and renumbering

`MolAlign/` matches atom numbering across a series so descriptors line up row to row. It
finds a maximum common substructure, then optimizes the atom mapping. Needs `torch`:

```bash
pip install torch
python MolAlign/renumbering.py --help
```

The renumbering hook is imported lazily, so nothing else in DescriPyTor requires torch.

---

## Docker

Reproducible Linux environment with RDKit, Morfeus, AQME, Mordred, xTB and the app, from
the repository root:

```bash
docker compose up --build
```

- Atom picker + modeling GUI: <http://localhost:7432/visual>
- Streamlit extraction app: <http://localhost:8503>

Put input files in `work/` — it is mounted at `/work` inside the container, so use paths
like `/work/features.csv` in the GUI. Anything written to `/work` shows up there too.

Setup from scratch (WSL2, Docker Desktop): [DOCKER_QUICKSTART.md](Getting_started_with_examples/descriptor_extraction_toolkit/webapp/DOCKER_QUICKSTART.md).

---

## Examples and further reading

Start in `Getting_started_with_examples/`:

| File | What it covers |
|---|---|
| `Practical_Notebook_Features.ipynb` | Feature extraction end to end |
| `Practical_Notebook_Modeling.ipynb` | Model search, validation, reporting |
| `feather_example/` | 26 substituted molecules, plus `input_example.json` |
| `modeling_example/` | Ready-made linear and logistic datasets |
| `cube_example/` | Density cubes for cube Sterimol |
| `descriptor_extraction_toolkit/` | The 3D picker, multi-engine extraction, BASSA |

Deeper documentation:

- [M2_data_extractor/README.md](M2_data_extractor/README.md) — `Molecules` / `Molecule` API
- [M3_modeler/README.md](M3_modeler/README.md) — modeling classes, CV, sampling, reports
- [descriptor_extraction_toolkit/README.md](Getting_started_with_examples/descriptor_extraction_toolkit/README.md) — picker, engines, config format
- [MolFeatures_Tutorial.md](Getting_started_with_examples/descriptor_extraction_toolkit/MolFeatures_Tutorial.md) — guided walkthrough

---

## Repository layout

```text
__main__.py                     CLI entry point
M1_pre_calculations/            Prepare and submit calculations (SMILES to xyz, .com files)
M2_data_extractor/              Descriptor extraction; Molecules / Molecule
M3_modeler/                     Regression and classification model search
MolAlign/                       Alignment and atom renumbering
utils/                          Shared file handling, geometry, visualization
Getting_started_with_examples/  Notebooks, example data, the 3D picker toolkit, webapp
docs/images/                    Screenshots and figures used by this README
docs/animations/                Animated SVGs, and the script that generates them
work/                           Docker-visible scratch folder
```

Atom indices are **1-based** everywhere, matching Gaussian.

The four animated figures are generated, not drawn — the geometry in them is computed by
`docs/animations/build_animations.py` from the same formulas the extractors use, so they
cannot silently drift away from the code:

```bash
python docs/animations/build_animations.py
```

---

## License

MIT.
