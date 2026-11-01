# Workflow and options

From Gaussian output to a feature table, three ways: the command line, Python, and the browser. Then every option
the steps take. Atom numbers are **1-based** everywhere, as in Gaussian. What to decide and check at each step,
before trusting a table: [EXTRACTION_PROTOCOL.md](EXTRACTION_PROTOCOL.md).

- [1. Install](#1-install)
- [2. The workflow on the command line](#2-the-workflow-on-the-command-line)
- [3. The same workflow in Python](#3-the-same-workflow-in-python)
- [4. In the browser](#4-in-the-browser)
- [5. Metal complexes](#5-metal-complexes)
- [6. Option reference](#6-option-reference)
- [7. Running the tests](#7-running-the-tests)

## 1. Install

```bash
conda create -n descripytor python=3.10 -y
conda activate descripytor
pip install descripytor            # or, from a clone: pip install -e ".[dev]"
```

If RDKit fails to install with pip: `conda install -c conda-forge rdkit`. `numpy<2` is required.

## 2. The workflow on the command line

The example uses the 18 aryl substrates and measured outcomes that ship with the package
(`descripytor/examples/baptiste_products`).

**Step 1: Gaussian logs to `.feather`.** One file per molecule, holding the geometry, charges
(NBO, Hirshfeld, CM5), dipole, polarizability and vibrational modes. The command asks for the
log folder and writes `feather_files/` inside it.

```bash
descripytor logs_to_feather
```

Skip this step for the example: its `.feather` files are already made.

The step is optional: a folder of Gaussian `.log`, ORCA `.out` or `.xyz` files is read as it is, by `extractor` and
the visual extractor alike. What each gives:

| files | geometry | charges | dipole | polarizability | vibrations |
|---|---|---|---|---|---|
| `.feather`, Gaussian `.log` | yes | NBO, Hirshfeld, CM5 | yes | yes | yes |
| ORCA `.out` | yes | Mulliken, Löwdin (Hirshfeld if printed) | yes | if computed | frequency jobs |
| structured `.json` | yes | as written | as written | as written | as written |
| `.xyz` | yes | no | no | no | no |

A descriptor the files cannot give (a charge from an `.xyz`, a vibration from a single point) is struck through in the
visual extractor and reported beside the count by the extraction.

Atom numbers in the input JSON mean the same atom in every molecule only if the folder is numbered alike: build or
number the structures from one template, and use **Check numbering** in the visual extractor (§4) to confirm it.

**Step 2: say which atoms each descriptor reads.** An input JSON maps a descriptor name to atom
numbers. Groups are separated by spaces and atoms within a group by commas:

```json
{
  "Dipole": "23,24,25,26,27,28,23,1",
  "Charges": "8,1",
  "Charge difference": "8,3 8,1 3,1",
  "Sterimol": "1,23 1,3",
  "Bond length": "1,2 1,23 1,3 3,8",
  "Bond angle": "23,1,3 23,1,2,3"
}
```

That asks for the dipole in the ring frame (origin at the centroid of atoms 23–28, y towards
23, the plane through 1), charges on atoms 8 and 1, three charge differences, Sterimol along
1→23 and 1→3, four bond lengths, and one angle plus one dihedral. Every key is listed in
[§6.1](#61-input-json-keys).

**Step 3: extract.**

```bash
descripytor extractor -i input.json -o features -f path/to/feather_files
```

This writes `features_<timestamp>.csv` (one row per molecule, one column per descriptor) and
`features_<timestamp>_correlation_table.csv`. It also prints diagnostics: columns with missing
values, and near-duplicate columns (r > 0.999). For the example: 18 molecules × 46 descriptors.

## 3. The same workflow in Python

```python
from M2_data_extractor.data_extractor import Molecules
from descripytor.examples import baptiste_example_dir

folder = baptiste_example_dir()
mols = Molecules(str(folder))                     # every .feather / .json in the folder
features = mols.get_molecules_features_set({
    "Dipole": "23,24,25,26,27,28,23,1",
    "Charges": "8,1",
    "Charge difference": "8,3 8,1 3,1",
    "Sterimol": "1,23 1,3",
    "Bond length": "1,2 1,23 1,3 3,8",
    "Bond angle": "23,1,3 23,1,2,3",
})
features.to_csv("features.csv")
```

Every descriptor also has its own method, for one molecule (`Molecule`) or the whole set
(`Molecules.*_dict`):

```python
mol = mols.molecules[0]
mol.get_sterimol([1, 23])                       # B1, B5, L, loc_B1, loc_B5, B1_B5_angle (theta), ...
mol.get_dipole_gaussian_df([[23, 24, 25, 26, 27, 28], 23, 1])
mol.get_charge_df([8, 1])                       # every charge type in the file
mol.get_charge_diff_df([[8, 3], [8, 1]])
mol.get_bond_length([[1, 2], [3, 8]])
mol.get_bond_angle([[23, 1, 3], [23, 1, 2, 3]])  # three atoms: angle; four: dihedral
mol.get_stretch_vibration([1, 2])               # the pair must be bonded
mol.get_bend_vibration([2, 3])                  # a pair sharing a centre (here atom 1), or a triplet 2-1-3
mol.get_ring_vibrations([23])                   # one ring atom; ortho/meta/para are found
mols.get_sterimol_dict([[1, 23], [1, 3]])       # the same for every molecule, as one table
```

## 4. In the browser

```bash
descripytor visual            # --port 7432 (default), --host 127.0.0.1, --no-browser
```

This opens the 3D atom picker at `http://127.0.0.1:7432/visual`. Load a folder, click atoms to build the selections
and extract, all on one page. The feature explorer is one click away in the header; it also opens directly with
`descripytor feature_explorer` and has its own manual with screenshots, [FEATURE_EXPLORER.md](FEATURE_EXPLORER.md).

## 5. Metal complexes

For metal–ligand complexes from xTB or GOAT, `MetalComplex` finds the donors, the metal and the
stereocentres itself and returns arm-averaged (`_sym`) and arm-difference (`_asym`) descriptors:

```python
from M2_data_extractor.metal_complex import MetalComplex, MetalComplexEnsemble

geo = MetalComplex.from_xyz("box_3_CuCl.xyz").geometric_features()
# fromM_* : Sterimol B1/B5/L/theta read from the metal; sub_* : from the stereocentre
# bite, M-D distances, %Vbur at 3 / 3.5 / 5 A, ...
ens = MetalComplexEnsemble.from_xyz("box_6.finalensemble.xyz").geometric_features()   # Boltzmann-averaged
```

| setting | default | meaning |
|---|---|---|
| `metal_complex.STERIMOL_SCAN_STEP` | `1` | B1 rotation scan step in degrees; `18` is the earlier coarse scan |
| `metal_complex.STERIMOL_FRAME` | `"fragment"` | where the scan starts; `"lab"` is the earlier lab-frame start |
| `metal_complex.STERIMOL_THETA_RULE` | `"soft"` | θ over tied B1 directions; `"scan"` is the 0.2.0 value |
| `metal_complex.SUB_FRAGMENT_BOUND` | `"donor_ring"` | where the C\*→R fragment stops; `"donor"` is the 0.2.0 walk |
| `metal_complex.DIPOLE_ABS` | `True` | unsigned `mu_desym` / `mu_outofplane`; `False` is the 0.2.0 value |
| `metal_complex.STERIMOL_KEYS` | `("B1", "B5", "L", "theta")` | add `"angle"` for the in-plane azimuth φ |
| `metal_complex.BITE_WINDOW` | `(65, 105)` | warn when the bite angle says a donor has come off the metal |

Electronic descriptors from xTB single points, and the ensembles: [METAL_COMPLEX.md](METAL_COMPLEX.md).

## 6. Option reference

### 6.1 Input JSON keys

A key is matched by its start, case-insensitively, so the long labels saved by the older GUIs
("Sterimol atoms - Primary axis along: …") work as well.

| key | value | gives |
|---|---|---|
| `Sterimol` | pairs `base,direction` | B1, B5, L, loc_B1, loc_B5, θ (`B1_B5_angle`), the in-plane azimuth and its acute form |
| `Drop atoms` | atoms | leave these atoms out of every Sterimol fragment |
| `Dipole` | `o,y,p` or `o1,o2,…,y,p` | dipole x/y/z and total in the frame: origin at `o` (or the centroid of `o1…`), y towards `y`, xy-plane through `p` |
| `Center atoms` | atoms | move the dipole frame's origin to the centroid of these atoms, keeping its axes |
| `Charges` | atoms | NBO, Hirshfeld and CM5 charge of each atom (whichever the file has) |
| `Charge difference` | pairs | charge of the first atom minus the second, per charge type |
| `Bond length` | pairs | distance |
| `Bond angle` | triplets or quadruplets | angle, or (unsigned) dihedral |
| `Stretching` | bonded pairs | frequency and amplitude of the mode that moves most along the bond |
| `Stretch` / `Upper stretch` | cm⁻¹ | the window searched: default 1400–3500. Widen the upper edge for O–H (about 3650) |
| `Bending` | pairs sharing a centre, or triplets | frequency and amplitude of the strongest bending mode |
| `Bend` | cm⁻¹ | lowest frequency considered, default 1300 |
| `Ring` | one ring atom (optionally a second, to choose between fused rings) | the two ring modes along and across the primary–para axis, with their angles to it |

`parameters={"Radii": "CPK", "Isotropic": True}` is the second argument of
`get_molecules_features_set`. `Radii` sets the Sterimol radii (`"CPK"`, `"bondi"`, `"Pyykko"`);
`Isotropic` adds polarizability and energy.

### 6.2 Bonds

Every descriptor that follows the molecular graph (Sterimol fragments, stretch/bend, rings,
atom typing) uses `utils.help_functions.extract_connectivity`:

| parameter | default | meaning |
|---|---|---|
| `threshold_distance` | `None` | `None`: a non-metal pair is bonded below `scale` × the sum of covalent radii (Pyykkö). A number: a flat cutoff in Å instead (`1.82` is the rule before September 2026). |
| `scale` | `1.15` | covalent-radius multiple |
| `metals` | all metals | elements given the metal rule; `[]` turns it off |
| `metal_threshold` | `2.8` | Å, metal–ligand; a metal–H pair must also pass the covalent test |
| `max_coordination` | `6` | bonds kept per metal, shortest first |

`Molecule(..., threshold=…)` and `Molecules(..., threshold=…)` pass `threshold_distance` on.
Check `mol.bonds_df` when a descriptor looks wrong.

### 6.3 Command line

```text
descripytor visual_extractor [--port N] [--host H] [--no-browser]
descripytor feature_explorer [--mol NAME] [--preset NAME] [--port N] [--host H] [--no-browser]
descripytor extractor       -i input.json -o out_name [-f feather_dir]
descripytor logs_to_feather (asks for the log folder)
descripytor sterimol        (asks for an .xyz folder, atom pairs and radii; writes xyz_sterimol.csv)
descripytor cube            (asks for a .cube folder and atom pairs; writes cube_sterimol.csv)
```

## 7. Running the tests

```bash
pip install -e ".[dev]"
python -m pytest -m "not slow"
```

A clean clone passes with no data files. The explorer's JavaScript is checked against the Python by three scripts in
`M2_data_extractor/feature_explorer/`:

```bash
python make_refs.py && node gate.js          # Sterimol, 102 values: 0 must differ
python make_pos_refs.py && node gate_pos.js  # %V_bur and cone angle
node smoke.js                                # every drawing
```
