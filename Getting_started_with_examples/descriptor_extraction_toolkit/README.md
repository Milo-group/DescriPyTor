# DescriPytor descriptor-extraction toolkit

A small, configurable harness that runs **many descriptor-extraction options**
over a molecule set in one go. You pick which **engines/packages** to use, give
each its own config, and choose atoms visually in a **3Dmol.js picker**.

Built on the same `Molecules` / `Molecules_xyz` API as
`Practical_Notebook_Features.ipynb`.

## Files

| File | Purpose |
|------|---------|
| `descriptor_extractor.py` | Main runner. Loads a set, runs the enabled engines, merges everything into one CSV. CLI + importable. |
| `feature_packages.py` | The external descriptor packages ported from `feature_extraction_and_regression.ipynb` (rdkit, mordred, deepchem, rafbl/moltop, qm, aqme_qdescp, morfeus suite). |
| `atom_picker.html` | Primary visual workflow: 3D atom selection, descriptor-package configuration, extraction export, CSV validation, and M3/BASSA modeling controls. |
| `make_picker.py` | Embeds a *real* molecule (from a feather set or `.xyz`) into the picker and opens it. |
| `config_example.json` | A starter run config with every engine block. |

## 1. Pick atoms visually

```bash
# load one molecule of a feather set into the picker
python make_picker.py --feather-dir ..\feather_example --index 0

# or a single xyz / a folder of xyz
python make_picker.py --xyz mol.xyz
python make_picker.py --xyz-dir path\to\xyz --name conf_2
```

In the picker: choose a **field** (Sterimol, Bond length, Dipole, Charges, …),
then **click atoms** in the 3D view. Pair/triplet fields auto-commit at the
right size; variable fields (bond angle/dihedral) use **Commit group**. Atoms
are 1-indexed, matching DescriPytor/Gaussian.

### Visualize panel (3Dmol.js)

Under the 3D view:

- **Sterimol** — tick it and choose one of your picked `[origin, attached]`
  pairs; the L / B1 / B5 vectors are computed live in JavaScript (from
  coordinates + van der Waals radii) and drawn as arrows with their values.
- **Dipole** — draws the molecular dipole vector (µ) as an arrow from the
  molecular center. Needs dipole data, so load the molecule from a feather set
  via `make_picker.py` (not a bare `.xyz`).
- **Vibration** — pick a normal mode (frequency / IR intensity) and press
  **Play**: the atoms oscillate along their displacement vectors. The amplitude
  slider scales the motion; **arrows** overlays the per-atom displacement
  vectors. Vibration data also comes from the feather molecule.

`make_picker.py --feather-dir … --index N` embeds the dipole vector and all
normal-mode displacements into the page automatically; a bare `.xyz` still gives
you the geometry-only Sterimol overlay.

Export either:
- **Download `answers.json`** → feed it to the runner (`atoms_file`), or
- **Copy the Python `answers_dict`** → paste into a notebook and call
  `mols.get_molecules_features_set(entry_widgets=answers_dict, ...)`.

`atom_picker.html` also opens on its own (with a demo molecule) — use the
*Load structure* button to drop in any `.xyz`.

## 2. Configure engines

`descriptor_extractor.py --print-template > my_run.json`, then edit. Each engine
has `"enabled"` plus its own settings:

| Engine | Package / call | Set | Config |
|--------|----------------|-----|--------|
| `descripytor_full` | `Molecules.get_molecules_features_set` | feather | `atoms` (the picker output: IR, dipole, charges, sterimol, bonds, angles) |
| `descripytor_steric` | `Molecules.get_sterimol_dict` | feather | `pairs` |
| `xyz_sterimol` | `Molecules_xyz.get_sterimol_df` | xyz | `pairs`, `radii` (`CPK`/`bondi`) |
| `xyz_geometry` | `Molecules_xyz.get_angles_df` / `get_bond_lengths_df` | xyz | `angles`, `bond_lengths` |
| `xyz_buried_volume` | `Molecules_xyz.get_buried_volume_df` | xyz | `metal_index`, `radius` |
| `morfeus_sterimol` | `morfeus.Sterimol` | xyz | `pairs` (independent cross-check) |

### Expanded packages (from `feature_extraction_and_regression.ipynb`)

These run on the `.xyz` set (SMILES is perceived once from the structures and
shared) or on `.log` files. Each is optional and skips itself with a one-line
message if its package isn't installed, so a run never dies on a missing dep.

| Engine | Package | Needs | Config |
|--------|---------|-------|--------|
| `rdkit` | RDKit full 2D descriptor list | `rdkit` | `embed` |
| `rdkit_fp` | ECFP/FCFP/MACCS/AtomPair/Torsion + USRCAT | `rdkit` | `n_bits`, `include_usrcat` |
| `mordred` | Mordred 2D/3D descriptors | `mordred` | `ignore_3d` |
| `deepchem` | CircularFingerprint + RDKitDescriptors | `deepchem` | `n_bits` |
| `rafbl` | moltop topological graph descriptors | `moltop` | `mult` |
| `qm` | autoqchem Gaussian-log descriptors | `autoqchem` | `log_dir` |
| `aqme_qdescp` | AQME QDESCP xTB/MORFEUS table | `aqme` + xTB | `level`, `atoms`, `charge`, `solvent`, … |
| `morfeus_suite` | Sterimol, BuriedVolume, ConeAngle, SASA, Dispersion, Pyramidalization | `morfeus-ml` | `descriptors`, `sterimol_pairs`, `metal_index`, `cone_atoms`, `pyramid_atoms`, `bv_radius` |

Install what you want to use:

```
pip install rdkit mordred deepchem morfeus-ml ase networkx
pip install autoqchem      # qm engine (Gaussian logs)
pip install aqme           # aqme_qdescp engine (also needs a working xTB)
```

`rafbl` additionally needs `moltop`. Install it only if you have a working
`moltop` source or wheel; the Docker image leaves it out because it is not
reliably available from public package indexes.

### AQME: end-to-end SMILES → descriptors

AQME is used two ways here:

1. **`aqme_qdescp` engine** — xTB/MORFEUS descriptors from your `.xyz` set
   (IP, EA, HOMO–LUMO, dipole, polarizability, FOD, SASA, dispersion, buried
   volume, …), molecular and atom-level. Needs `aqme` + a working **xTB**.
   `level` picks the table width (`denovo` / `interpret` / `full`); `atoms`
   selects atom-centered descriptors (`[]` auto-detects, or e.g. `['P']`,
   `['C=O']`); `solvent` adds ALPB solvation terms. Run
   `python descriptor_extractor.py --list-aqme-descriptors` to see the names.

2. **`--csearch` preprocessing** — AQME CSEARCH turns a SMILES CSV into a
   conformer `.xyz` ensemble, so the whole pipeline can start from SMILES:

   ```
   python descriptor_extractor.py --csearch molecules.csv --csearch-out xyz_out \
       --name-col name --smiles-col smiles --csearch-program rdkit \
       --csearch-sample 10 --keep lowest
   ```

   `--keep lowest` writes one representative xyz per molecule; `--keep all`
   writes `name_conf_1.xyz`, `name_conf_2.xyz`, …  `--csearch-program crest`
   uses CREST (must be installed). Then point your run config's `xyz_dir` at
   `xyz_out` and enable the xyz engines (`aqme_qdescp`, `morfeus_suite`,
   `rdkit`, `mordred`, …).

`prefix` namespaces each engine's columns so they never collide. In the picker,
the **Feature packages** checkboxes toggle these engines and the **morfeus_suite
atoms** rows (Sterimol pair, buried-volume metal, cone-angle apex,
pyramidalization center) feed `morfeus_suite` — all of it lands in the exported
`run_config.json`.

## 3. Run

```bash
python descriptor_extractor.py --config my_run.json
# -> writes output_csv with every enabled engine merged column-wise
```

Or from a notebook:

```python
from descriptor_extractor import run_from_config
df = run_from_config("my_run.json")
```

## 4. Model the extracted CSV

Expand **Modeling** in `atom_picker.html`, then click **Use extraction output**.
If the descriptor CSV does not yet contain the measured target, paste the
output values into **Output values to merge**. Match them either to CSV row
order or provide one `name,value` pair per line. **Merge outputs into CSV**
writes a new CSV and selects it for modeling; the descriptor CSV is not
overwritten.

The panel then validates that the selected target (normally `output`) is
numeric, reports the usable descriptors, estimates the M3 feature-combination
count, and generates a copyable command. With
`M2_data_extractor/gui_server.py` running, it can execute M3 regression,
classification, or BASSA.

### Where calculations run

- The molecule display, picking, and visual overlays run in the browser.
- Extraction, CSV merging, M3, and BASSA run in the Python backend.
- M3 uses CPU workers selected by **CPU jobs**. BASSA/PyMC also runs on CPU.
- In Docker, both tools use the CPUs and memory assigned to Docker Desktop.
  The GUI status reports the number of processors visible inside the container.

### Docker

From the `MolFeatures` directory:

```powershell
docker compose up --build
```

Open the combined atom-picker and modeling GUI at
`http://localhost:7432/visual`. The Streamlit extraction app is available at
`http://localhost:8503`. Put local input files in the visible
`MolFeatures/work` folder and use container paths such as `/work/features.csv`
in the GUI. Generated files written to `/work` appear in that folder too.

## Notes

- Set `root_dir` (or the `DESCRIPYTOR_ROOT` env var) to your `MolFeatures`
  folder so the DescriPytor packages import. Defaults to the path used in the
  notebooks.
- Engines fail independently: if one errors it's reported and the rest still
  run.
- Rows are natural-sorted (`conf_2` before `conf_10`) like the notebooks.
