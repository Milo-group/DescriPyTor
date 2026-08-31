# M2 — feature extraction

Turns quantum-chemistry output into model-ready descriptors.

- **`Molecule`** — one molecule: geometry, connectivity, charges, dipoles, vibrations.
- **`Molecules`** — a whole set, with batch extractors returning wide DataFrames.
- **`Molecules_xyz`** — the same idea for bare `.xyz` files, no quantum data needed.
- **`MetalComplex` / `MetalComplexEnsemble`** — GOAT/xTB metal-complex descriptors (atom 0 = metal). See [docs/METAL_COMPLEX.md](../docs/METAL_COMPLEX.md).

> **Atom indices are 1-based**, matching Gaussian. `adjust_indices` normalises them
> internally; you never pass 0-based indices.

---

## Loading a set

```python
from M2_data_extractor.data_extractor import Molecules
from descripytor.examples import feather_example_dir

molset = Molecules(str(feather_example_dir()), threshold=1.82)

molset.success_molecules   # files that parsed
molset.failed_molecules    # files that did not — check these before trusting a run
```

`Molecules` loads `*.feather` (and `*.json`). `threshold` is the bond-distance cutoff.

Don't have feathers yet? Convert Gaussian logs first:

```python
from M2_data_extractor.feather_extractor import logs_to_feather
logs_to_feather("path/to/logs", "path/to/feathers")
```

### Selection and export

```python
molset.filter_molecules([0, 3, 5])   # keep a subset by position
molset.export_all_xyz()              # write per-molecule XYZ
molset.extract_all_dfs()             # dump every internal DataFrame as CSV
molset.visualize_molecules([0, 1])   # 3D preview
molset.show_ring_atoms(n=5)          # inspect detected ring systems
```

---

## The one call that does everything

`get_molecules_features_set` is what the picker and CLI call. Pass atom groups per family;
you get one feature matrix.

```python
features = molset.get_molecules_features_set(
    entry_widgets={
        'Ring':        '[[8],[14]]',
        'Stretching':  '[[1,2],[3,4]]',
        'Bending':     '[[1,2],[3,4]]',
        'Dipole':      '[[1,2,3], 5, 6]',
        'Charges':     '[3,5,7,9]',
        'Charge-Diff': '[[1,2],[3,4]]',
        'Sterimol':    '[[1,6],[3,4]]',
        'Bond-Angle':  '[1,2,3]',
        'Bond-Length': '[[1,2],[3,4]]',
    },
    parameters={'Radii': 'CPK', 'Isotropic': True},
    save_as=True,
    csv_file_name='features_csv',
)
```

**Keys match the longest descriptor name they start with** (case, punctuation, and help
text ignored). All of these hit `bond_length`:

```text
'Bond-Length'   'Bond length  [a, b]'   'Bond_length - Atom pairs to calculate difference: …'
```

That is why the picker can use long labels. An unmatched label is **reported**, not dropped:

```text
Warning: ignoring unrecognized feature-set entry 'Wibble atoms' - no descriptor matches it.
```

then that field added no columns. Two labels for the same descriptor: the populated one wins.

| Key | Expects | Produces |
|---|---|---|
| `Ring` | one ring atom per group | cross/para ring-mode frequencies and angles |
| `Stretching` | bonded pairs | frequency and amplitude along each bond |
| `Bending` | pairs sharing a centre | strongest bending mode |
| `Dipole` | `[origin(s), y-atom, plane-atom]` | dipole components in the transformed frame |
| `Charges` | single atoms | one column per charge scheme |
| `Charge-Diff` | pairs `[a,b]` | `q_a − q_b` per scheme |
| `Sterimol` | `[origin, attached]` pairs | `L`, `B1`, `B5` |
| `Bond-Angle` | triads or quartets | angles, dihedrals |
| `Bond-Length` | pairs | distances |

Auxiliary keys tune the above rather than adding columns: `Stretch` and `Upper-Stretch`
(frequency window), `Bend` (threshold), `Drop-Atoms` (exclude from Sterimol),
`Center_Atoms` (move the dipole frame's origin — see [Dipoles](#dipoles)).

`parameters`: `Radii` (`'CPK'`, `'bondi'`, `'Pyykko'`) and `Isotropic` (append polarizability
and energy). Skipped steps are skipped. One failed molecule is reported; the run continues.

`save_as=True` writes `features_csv_<timestamp>.csv` plus a correlation table.

---

## Batch extractors

Call these for one family. Each returns a wide DataFrame.

```python
molset.get_sterimol_dict([[1,6],[3,4]], radii='CPK', drop_atoms=None)
molset.get_ring_vibration_dict([[8],[14]], freq_min=1550, freq_max=1700)
molset.get_dipole_dict([[1,2,3], 5, 6])
molset.get_bond_angle_dict([1,2,3])
molset.get_bond_length_dict([[1,2],[3,4]])
molset.get_stretch_vibration_dict([[1,2],[3,4]], threshold=1400, upper_threshold=3500)
molset.get_bend_vibration_dict([[1,2],[3,4]], threshold=1300)
molset.get_charge_df_dict([3,5,7,9])
molset.get_charge_diff_df_dict([[1,2],[3,4]], type='all')
```

---

## One molecule at a time

```python
from M2_data_extractor.data_extractor import Molecule
mol = Molecule("LS1716_optimized.feather", threshold=1.82)
```

Populates geometry (`xyz_df`, `coordinates_array`), connectivity (`bonds_df`, `atype_list`)
and quantum output (`gauss_dipole_df`, `polarizability_df`, `energy_value`, `charge_dict`,
`vibration_dict`).

### Sterimol

```python
mol.get_sterimol(base_atoms=[1, 6], radii='CPK', sub_structure=True,
                 drop_atoms=None, mode='all')
```

`base_atoms=None` auto-detects a base. Accepts a single pair or a list of pairs.
`sub_structure=True` measures only the fragment hanging off the bond; `drop_atoms` excludes
specific atoms. `L = max(yᵢ + rᵢ)`, `B5 = max(‖(xᵢ,zᵢ)‖ + rᵢ)`, and `B1` is the smallest
half-width found by rotating a supporting line around the cross-section.

### Charges

```python
mol.get_charge_df([3,5,7,9], type='all')          # dict of DataFrames, one per scheme
mol.get_charge_df([3,5,7,9], type=['nbo','cm5'])
mol.get_charge_diff_df([[1,2],[3,4]], type='all') # labelled diff_i-j
```

`type='all'` covers every scheme present in `charge_dict` and silently skips absent ones.

### Dipoles

```python
# base spec is one of: [o, y, plane] · [o1, o2, …, y, plane] · [[o1, o2, …], y, plane]
mol.get_dipole_gaussian_df([[1,2,3], 5, 6], visualize_bool=True)

# anchor the frame somewhere else while keeping the same axis directions
mol.get_dipole_gaussian_df([[1,2,3], 5, 6], center_atoms=[1,2,3,4,5,6])
```

The origin is the centroid of the origin set; `ŷ` points at the y-atom; `x̂` is the
plane-atom direction with its `ŷ` component removed (Gram–Schmidt). Components are reported
in that frame, which is what makes them comparable across a series.

**`center_atoms`** moves the origin to the centroid of the atoms you name, while `ŷ` and the
plane atom still come from the base spec. Because both axes are measured *from* the origin,
moving it re-aims them — this is how you get a ring-centred frame while a substituent still
sets the y direction. Columns from a run with a centre are suffixed (`…_c{1,2,3}`) so they
don't collide with the uncentred ones. Omit it, or pass `[]`, for the original behaviour.

The same argument is on `get_dipole_gaussian_df_single` and on the batch
`Molecules.get_dipole_dict(...)`, where it applies to every molecule in the set. In the
feature-set call it is the **`Center_Atoms`** entry.

`visualize_bool=True` draws the molecule in that same frame, so the arrow you see and the
components you get are measured from the same origin.

### Vibrations

```python
mol.get_stretch_vibration([1,2], threshold=1600, upper_threshold=3000)
mol.get_bend_vibration([1,2], threshold=1300)
mol.get_ring_vibrations([[1,4],[2,5]], return_nan_on_empty=True)
```

Stretch validates the pair is bonded, then scores every mode by `|dᵃ·û| + |dᵇ·û|` and keeps
the best inside the frequency window. Bend finds the shared centre and uses cross-product
magnitudes. Ring vibrations resolve a benzene-like pattern and return `cross`, `cross_angle`,
`para`, `para_angle` — with `return_nan_on_empty=True` you get a NaN row instead of `None`
when filters eliminate everything, which keeps table shapes aligned.

### Geometry helpers

```python
mol.get_coordinates_mean_point([1,2,3])                 # centroid
mol.get_coordination_transformation_df([1,2,3])         # frame transform
mol.get_coordination_transformation_df([1,2,3], origin=[4,5,6])   # centred elsewhere
mol.renumber_atoms({old: new, ...})                     # rebuilds everything consistently
mol.swap_atom_pair((1, 2))
mol.write_xyz_file(); mol.write_csv_files()
mol.visualize_molecule()
```

`renumber_atoms` rebuilds geometry, connectivity, types **and** reindexes `charge_dict` and
`vibration_dict`, so descriptors stay consistent after a renumber.

---

## Without quantum data: `Molecules_xyz`

Geometry-only descriptors straight from `.xyz` — no Gaussian run required.

```python
from M2_data_extractor.sterimol_standalone import Molecules_xyz

xyzset = Molecules_xyz("path/to/xyz_dir")
xyzset.get_sterimol_df([[1,6],[3,4]], radii='CPK')
xyzset.get_angles_df([[1,2,3]])
xyzset.get_bond_lengths_df([[1,2],[3,4]])
xyzset.get_buried_volume_df(metal_index=1, radius=3.5)
```

## Sterimol from electron density

```python
from M2_data_extractor.cube_sterimol import cube_many
```

Radii come from the density isosurface rather than a lookup table, so they respond to the
electronic environment. Driven from the CLI with `descripytor cube`.

---

## Practical notes

- **Connectivity**: default cutoff is 1.82 Å. Check `bonds_df` if a descriptor looks wrong —
  a missing bond is the usual cause of a failed stretch lookup.
- **Failures are per-molecule**: batch methods report and skip, so one bad log doesn't lose
  the run. Always check `failed_molecules`.
- **Re-extract after renumbering** — descriptors are index-based.
- **Correlations**: `get_molecules_features_set` writes a correlation table alongside the
  matrix; scan it before modelling, since near-duplicate descriptors inflate subset counts.

Full pipeline context: [root README](../README.md).
