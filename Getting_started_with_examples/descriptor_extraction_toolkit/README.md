# Descriptor-extraction toolkit

Run several descriptor engines on one molecule set. Pick atoms in 3D, extract a CSV in the
GUI, or save `run_config.json` for the CLI.

Uses the same `Molecules` / `Molecules_xyz` API as `Practical_Notebook_Features.ipynb`.
## Files

| File | Purpose |
|------|---------|
| `descriptor_extractor.py` | Main runner. Loads a set, runs the enabled engines, merges everything into one CSV. CLI + importable. |
| `feature_packages.py` | The external descriptor packages ported from `feature_extraction_and_regression.ipynb` (rdkit, mordred, deepchem, rafbl/moltop, qm, aqme_qdescp, morfeus suite). |
| `atom_picker.html` | 3D picks, visualizations, Extract CSV, optional `run_config.json` |
| `make_picker.py` | Embeds a *real* molecule (from a feather set or `.xyz`) into the picker and opens it. |
| `config_example.json` | A starter run config with every engine block. |

## 1. Pick atoms visually

Start the picker from the repo root with `descripytor visual`, or embed one molecule:

```bash
# load one molecule of a feather set into the picker
python make_picker.py --feather-dir ..\..\descripytor\examples\feather_example --index 0

# or a single xyz / a folder of xyz
python make_picker.py --xyz mol.xyz
python make_picker.py --xyz-dir path\to\xyz --name conf_2
```

In the picker: choose a **field**, **click atoms** (1-based, Gaussian style). Pairs and
triplets auto-commit; variable fields use **Commit group**. **Visualizations** (Sterimol,
dipole, vibration, feature preview) sit under the picks. Set the feather folder, then
**Extract CSV** — the table appears in the page; **Download CSV** saves it. **Save picks**
writes `run_config.json` for later. Extra engines are under **Advanced**.

### Visualizations

Next to the 3D view:

- **Sterimol** — pick an `[origin, attached]` pair; L / B1 / B5 draw live.
- **Dipole** — molecular dipole from a feather (not a bare `.xyz`).
- **Vibration** — play a normal mode; amplitude slider and optional arrows.
- **Feature preview** — overlay the atoms you have already picked.

`make_picker.py --feather-dir … --index N` embeds dipole and modes. A `.xyz` still gets
geometry-only Sterimol.

Or export:

- **Save picks** → `run_config.json` for `descriptor_extractor.py --config`
- **Python `answers_dict`** → paste into a notebook

`atom_picker.html` opens on bundled `basic.feather` (unsubstituted benzene) when you run
`descripytor visual`. Load any `.xyz` or `.feather`.

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

These run on `.xyz` (SMILES is perceived once) or `.log` files. A missing package is skipped
with a one-line message.

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

`rafbl` needs `moltop`. Skip it if you have no wheel; the Docker image omits it.

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

`prefix` namespaces columns. In the picker, **Advanced → Feature packages** toggles engines;
morfeus atom rows feed `morfeus_suite`. All of that goes into `run_config.json`.

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

The 3D picker **Model** tab adds an output column (choose one from the CSV or paste values) and runs ordinary least-squares with leave-one-out Q². You can still use the CLI:

```bash
descripytor model -m regression -f merged_features.csv -y output --min-features 1 --max-features 3
```

Or Streamlit tab **5 · Modeling** (`pip install -e ".[webapp]"`, then the webapp / Docker at
http://localhost:8503). Add a numeric target column (usually `output`) if the CSV is
descriptors only.

Picking and overlays run in the browser. Extraction and modeling run in Python (CPU).

### Docker

From the repository root:

```powershell
docker compose up --build
```

Open the picker at `http://localhost:7432` (forms: `/forms`). Streamlit is at
`http://localhost:8503`. Put files in `work/`; use `/work/yourfile.feather` in the GUI.

## Notes

- Set `root_dir` or `DESCRIPYTOR_ROOT` to the clone (Docker: `/workspace/descripytor`).
- Engines fail independently; the rest still run.
- Rows are natural-sorted (`conf_2` before `conf_10`).
