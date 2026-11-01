# Metal-complex descriptors from GOAT ensembles and xTB single points

Metal complexes from GOAT, CREST or xTB do **not** go through Gaussian feathers. Their descriptors come from:

1. metal-complex XYZ (GOAT ensemble, CREST ensemble, or one GFN2 relaxation)
2. xTB charges / Wiberg pairs / dipole / HOMO–LUMO
3. SMILES graph indices of the free ligand

Those extractors now live in `M2_data_extractor` as first-class classes.
They do **not** import morfeus.

## Classes

| Class | Module | Role |
|---|---|---|
| `XYZEnsemble` | `xyz_io.py` | Parse multi-XYZ, Boltzmann weights |
| `GoatEnsemble` | `xyz_io.py` | GOAT multi-XYZ with last-float energies |
| `XtbSinglePoint` | `xtb_singlepoint.py` | Parse `.q` / `.wbo` / `.dip`, `props.txt`, or compact dumps |
| `MetalComplex` | `metal_complex.py` | One geometry → `fromM_*`, `sub_*`, `%Vbur`, electronics |
| `MetalComplexEnsemble` | `metal_complex.py` | Boltzmann-average `MetalComplex` over a GOAT/CREST XYZ |
| `MetalComplexSet` | `metal_complex.py` | Batch → DataFrame |
| `LigandTopology` | `ligand_topology.py` | SMILES → `t_*` and size-normalized `tf_*` |

`Molecules` / `Molecules_xyz` are unchanged. Use them for Gaussian feathers and
for generic XYZ Sterimol when **you** pick the atom indices. Use `MetalComplex`
when atom 0 is the metal and the next atoms are ancillaries then the two donors.

## Energy on the comment line

- CREST: first float is the GFN2 energy. Pass `energy_convention='first'`.
- GOAT: take the **last** float. An RMSD written first would otherwise become
  the Boltzmann weight. `MetalComplexEnsemble`, `XYZEnsemble`, and `GoatEnsemble`
  default to `'last'`.

GOAT files that look like `-72.91 converged=true` work with either rule; last-float is the safe one.

## Canonical atom order

```
0 = metal
1 … = H / Cl / F / Br ancillaries
then the two chelate donors
```

`Ni(H)2L` donors are atoms 3 and 4 (0-based). `Cu(Cl)L` donors are 2 and 3.
Do not hard-code 3 and 4.

## Geometric block

Per frame, then Boltzmann-averaged at 298.15 K (`HARTREE_TO_KCAL = 627.5095`):

- `fromM_*_sym` / `_asym` — Sterimol of each chelate half seen from the metal
- `sub_*`, `a_R_Cstereo_N_*` — substituent Sterimol and R–C*–N angle
- `bite`, `MD_mean`, `MD_asym`
- `vbur_3`, `vbur_3.5`, `vbur_5`, `vbur_noH_3.5` — Cavallo %Vbur, Bondi
  table (Ni 1.63 Å, Cu 1.40 Å), not Alvarez

Sterimol uses Verloop CPK types from coordination number, not morfeus. Each arm reports B1, B5, L and θ (the tilt of the
B5 vector out of the B1 plane, 0° for a lone H). Module settings:

| setting | default | set it to |
|---|---|---|
| `STERIMOL_SCAN_STEP` | `1` (every degree) | `18`: the earlier coarse scan |
| `STERIMOL_FRAME` | `"fragment"` (scan starts from the molecule) | `"lab"`: the earlier lab-frame start |
| `STERIMOL_THETA_RULE` | `"soft"` (θ averaged over tied B1 directions) | `"scan"`: θ at the direction the scan lands on |
| `SUB_FRAGMENT_BOUND` | `"donor_ring"` (C\*→R walk stops at the donor's ring) | `"donor"`: blocked only at the donor and metal |
| `DIPOLE_ABS` | `True` (\|mu_desym\|, \|mu_outofplane\|) | `False`: signed, flips with donor order |
| `STERIMOL_KEYS` | `("B1", "B5", "L", "theta")` | add `"angle"` for the in-plane azimuth φ |

## Electronic block

From one xTB single point on a geometry in that same atom order:

- `q_anc_sum` — sum of charges on the co-ligands between metal and donors
- `q_metal`, `q_donor_sym` / `_asym`, `q_absmax`, `q_spread`
- `wbo_MD_sym` / `_asym`, `wbo_metal_tot` (missing Wiberg pairs count as 0)
- `mu_bisector`, `mu_outofplane`, `mu_desym` in the M–D1–D2 frame
- `homo`, `lumo`, `gap`

φ (`sub_angle_*`, only with `"angle"` in `STERIMOL_KEYS`) on a **hydrogen substituent**
(1-atom Sterimol fragment) is ill-defined, which is why θ is the default: B1 ≈ B5 ≈ 1.0 Å and the B1–B5 angle comes out near 180°.

Electronic descriptors are read from one single point, not Boltzmann-averaged. The dipole components (`mu_*`) are
only meaningful when the XYZ is the structure the xTB single point was run on.


## Topology block

`LigandTopology.from_smiles(smi).features()`:

- `t_*` — raw graph indices
- `tf_*` — `t_* / n_heavy` 

No geometry, no QM, no stereochemistry.

## Examples

GOAT ensemble geometry:

```python
from M2_data_extractor import MetalComplexEnsemble

ens = MetalComplexEnsemble.from_xyz("box_6.finalensemble.xyz")
row = ens.geometric_features()
# row['fromM_theta_sym'], row['n_conformers'], ...
```

xTB electronics:

```python
from M2_data_extractor import (
    MetalComplex, XtbSinglePoint, parse_props_file,
)

props = parse_props_file("props.txt")["box_6"]
sp = XtbSinglePoint.from_xtb_dir("chg_ni", "box_6", props=props)
mc = MetalComplex.from_xyz("box_6.xyz")
elec = mc.electronic_features(sp)
```

Cheap GFN2 structure + compact dump:

```python
from M2_data_extractor import MetalComplex, XtbSinglePoint

mc = MetalComplex.from_xyz("ni_box_6.xyz")
dump = XtbSinglePoint.parse_dump("ni_box_6.txt")
elec = mc.electronic_features(dump["cheap"][0])
geom = mc.geometric_features()
```

Free-ligand graph indices:

```python
from M2_data_extractor import LigandTopology

topo = LigandTopology.from_smiles("C[C@H]1COC(C2=N[C@@H](C)CO2)=N1")
tf = topo.size_normalized_features()  # tf_zagreb1, tf_kappa3, ...
```

Batch a directory of GOAT ensembles:

```python
from M2_data_extractor import MetalComplexSet

s = MetalComplexSet.from_xyz_dir("cbens", pattern="*.finalensemble.xyz")
geom = s.geometric_dataframe()
```

## What this package still does not do

- Run GOAT, CREST, or xtb. Point the classes at the files those programs write.
- Place the metal on a free ligand. That is a
  separate geometry step; extraction assumes the complex XYZ already exists.
- NBO / Gaussian charges. Those still go through `Molecules` and feathers.

## Tests

Fast checks only (this is what to run on a laptop):

```bash
pytest tests/test_metal_complex.py -m "not slow"
```
