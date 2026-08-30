# Metal-complex descriptors from GOAT ensembles and xTB single points

Case Study 3 (Corminboeuf / Waser four-reaction benchmark) does **not** go
through Gaussian feathers. The tables were built from:

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

The CS3 GOAT files look like `-72.91 converged=true`, so first-float also
works, but last-float is the rule the tables were built with.

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
- `vbur_3`, `vbur_3.5`, `vbur_5`, `vbur_noH_3.5` — Cavallo %Vbur, CS3 Bondi
  table (Ni 1.63 Å, Cu 1.40 Å), not Alvarez

Sterimol uses Verloop CPK types from coordination number, same kernel as the
CS3 scratchpad (not morfeus).

## Electronic block

From one xTB single point on a geometry in that same atom order:

- `q_anc_sum` — sum of charges on the co-ligands between metal and donors
- `q_metal`, `q_donor_sym` / `_asym`, `q_absmax`, `q_spread`
- `wbo_MD_sym` / `_asym`, `wbo_metal_tot` (missing Wiberg pairs count as 0)
- `mu_bisector`, `mu_outofplane`, `mu_desym` in the M–D1–D2 frame
- `homo`, `lumo`, `gap`

`sub_angle_*` on a **hydrogen substituent** (1-atom Sterimol fragment) is
ill-defined: B1 ≈ B5 ≈ 1.0 Å and the B1–B5 angle comes out near 180°. The
package matches the current scratchpad `extract_cb.py` on those ligands.
`geom/desc_cp_goat_ens.csv` does **not** — nine Cu cyclopropanation ligands
differ by 30–90° on `sub_angle` while `sub_B1` / `sub_B5` / `sub_L` still
match to 1e-15. Tests follow `extract_cb`, not the stale CSV angle.

The **reported** CS3 electronics are one SP on the GOAT minimum
(`elec_ni.csv` / `elec_cu1.csv` / `elec_cu2.csv`), not a Boltzmann average.
Note: dipole components (`mu_*`) are only reproduced when the XYZ is the
structure the xTB single point was run on. The on-disk `cb_complexes` /
`cb_lowest` files no longer match `props.txt` Cartesian dipoles (the original
`elec_rebuild.py` check also fails on those files). Charge, WBO, and HOMO/LUMO
columns do reproduce. Dipoles are checked on the cheap arm, where geometry and
the compact dump come from the same GFN2 job.


## Topology block

`LigandTopology.from_smiles(smi).features()`:

- `t_*` — raw graph indices
- `tf_*` — `t_* / n_heavy` (the block the models searched)

No geometry, no QM, no stereochemistry.

## Examples

GOAT ensemble geometry (recreates a row of `desc_cc_goat_ens.csv`):

```python
from M2_data_extractor import MetalComplexEnsemble

ens = MetalComplexEnsemble.from_xyz("081_lig.finalensemble.xyz")
row = ens.geometric_features()
# row['fromM_angle_sym'], row['n_conformers'], ...
```

xTB electronics (recreates a row of `elec_ni.csv`):

```python
from M2_data_extractor import (
    MetalComplex, XtbSinglePoint, parse_props_file,
)

props = parse_props_file("props.txt")["081_lig"]
sp = XtbSinglePoint.from_xtb_dir("chg_ni", "081_lig", props=props)
mc = MetalComplex.from_xyz("081_lig.xyz")
elec = mc.electronic_features(sp)
```

Cheap GFN2 structure + compact dump:

```python
from M2_data_extractor import MetalComplex, XtbSinglePoint

mc = MetalComplex.from_xyz("ni_081_lig.xyz")
dump = XtbSinglePoint.parse_dump("ni_081_lig.txt")
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
- Place the metal on a free ligand. That builder (`build_general.py`) is a
  separate geometry step; extraction assumes the complex XYZ already exists.
- NBO / Gaussian charges. Those still go through `Molecules` and feathers.

## Tests

Fast checks only (this is what to run on a laptop):

```bash
pytest tests/test_metal_complex.py -m "not slow"
```

Do **not** run `-m slow` on a workstation. That path Boltzmann-averages every
GOAT frame (~4000 structures × Sterimol + %Vbur) and will pin the CPU for
many minutes. Recreate the full CS3 geometric tables on the cluster, or skip
them. The 081_lig fixture plus electronics / cheap-arm tests already check
the same kernel.
