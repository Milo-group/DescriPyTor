# Extraction protocol: from structures to a feature table you can trust

[WORKFLOW.md](WORKFLOW.md) says how to run each step. This page says what to decide and what to
check at each step, so the feature table is correct. Atom numbers are **1-based** throughout.

- [0. Before you compute](#0-before-you-compute)
- [1. Pick the route](#1-pick-the-route)
- [2. Structures](#2-structures)
- [3. Bonds](#3-bonds)
- [4. Atom selection](#4-atom-selection)
- [5. Descriptor rules](#5-descriptor-rules)
- [6. Extract and record](#6-extract-and-record)
- [7. Check the table](#7-check-the-table)
- [8. Hand-off](#8-hand-off)

## 0. Before you compute

**Build structures from names or SMILES**, not from figures or image recognition.

**Check every SMILES with RDKit against the source:**

- the molecular formula;
- the exact mass, when HRMS is given. A published "calculated" HRMS value is often the neutral
  mass, 0.00055 away from the ion.

**Read the source publication's computational method before building anything:** which metal, which
complex, which level of theory.

## 1. Pick the route

| input | use | typical case |
|---|---|---|
| Gaussian logs (NBO / Hirshfeld / CM5 charges, dipole, polarizability, frequencies) | `descripytor logs_to_feather`, then `Molecules(folder)` | organic substrates and ligands |
| xyz files, atoms chosen by you | `Molecules_xyz`, `descripytor sterimol` | steric descriptors only |
| metal-complex xyz (atom 1 = metal, then ancillary ligands, then the two donors) | `MetalComplex`, `MetalComplexEnsemble`, `MetalComplexSet`, `XtbSinglePoint` | chelating ligands on a metal |
| free-ligand SMILES | `LigandTopology` (`tf_*`: graph indices divided by heavy-atom count) | topology indices |

## 2. Structures

Read every descriptor on the final geometry, and keep every intermediate file.

**Conformer ensembles:**

- Report which conformer sampler you used.
- The energy on the comment line differs by program:
  - GOAT writes it as the **last** number;
  - CREST writes it as the **first**.

  Set `energy_convention` to match. An RMSD read as an energy silently corrupts the weights.

**CREST needs a relaxed starting geometry.** A placeholder metal–donor distance can make it
dissociate the complex.

**Placing a metal on a free ligand:**

- put the metal along the donor lone-pair direction;
- make sure no atom sits closer than about 2.2 Å (heavy) or 1.8 Å (H) to it;
- keep the builder script and a manifest (formula, metal and donor indices) next to the inputs.

**Check before going on:**

- the element sequence of every file matches the manifest;
- no imaginary frequencies in DFT minima;
- for chelates, the bite angle is inside `BITE_WINDOW` (65–105°). A complex whose donor came
  off the metal is removed from the set.

## 3. Bonds

Every graph-based descriptor depends on the bond list: Sterimol fragments, stretch and bend
modes, rings, atom typing. The rule since 0.2.0:

- **non-metal pairs:** bonded below 1.15 × the sum of Pyykkö covalent radii;
- **metal–ligand pairs:** bonded below 2.8 Å, and a metal–H pair must also pass the covalent test.

The old flat 1.82 Å cutoff silently dropped S–CF3 and P–C bonds.

**Before extracting a whole set**, look at `mol.bonds_df` (or `MetalComplex.adj`) for one molecule
of each ligand class:

- a P,N-chelating ligand must come out as P,N, not N,N;
- a monodentate ligand is not a chelate.

## 4. Atom selection

**Use the same numbering across the set,** or renumber first. Check by listing the element at
every chosen index in every molecule, not by looking at one molecule.

**Input JSON** (full list in [WORKFLOW.md §8.1](WORKFLOW.md#81-input-json-keys)): groups are separated
by spaces, atoms within a group by commas.

- Keys: `Sterimol`, `Dipole`, `Charges`, `Charge difference`, `Bond length`, `Bond angle`,
  `Stretching`, `Bending`, `Ring`.
- Modifiers: `Drop atoms`, `Center atoms`, `Stretch` / `Upper stretch`, `Bend`.

**Ring frames:**

- take the origin from the ring's own atoms, never from an index range that shifts when the
  substituent changes;
- put the plane-defining atom on the intended side of the ring.

## 5. Descriptor rules

### Sterimol (B1, B5, L, loc_B1, loc_B5, θ)

**Measure the whole molecule, and restrict the measurement** with `block` or `Drop atoms`. Never
measure a fragment cut out of the molecule: the frame's third atom changes, and B1 moves by up
to 0.09 Å.

**Metal complexes:**

- `fromM_*` is read from the metal and `sub_*` from the stereocentre.
- The defaults since 0.2.1 are the corrected ones:
  - the stereocentre→substituent fragment stops at the donor's ring;
  - θ is averaged over tied B1 directions (τ = 0.02 Å).

**θ edge cases:**

- θ is 0 for a lone-H substituent.
- When the attachment atom lies on the axis, the B1 minimum is flat: report θ as a range.

**B1 resolution is about ±0.01 Å.** Smaller differences are not signal.

### Dipole

**State the frame:** origin, y-axis atom, plane atom.

**Components along an axis whose sign depends on atom order are reported as magnitudes.** In a
metal complex these are `mu_desym` and `mu_outofplane`, whose axes flip with the order of the two
donors (`DIPOLE_ABS`).

### Charges and xTB properties

- NBO, Hirshfeld and CM5 charges come from the Gaussian feathers.
- xTB charges, Wiberg bond orders, HOMO/LUMO and the dipole come from one single point on exactly
  the geometry the other descriptors used. The dipole only reproduces on that geometry.

### Vibrations

**Windows:**

- stretch 1400–3500 cm⁻¹; raise the upper edge to about 3650 for O–H;
- bend from 1300 cm⁻¹.

**Inputs:** a stretch pair must be bonded, and a bend triplet a-b-c means its two ends, a and c.

**Check that the selected mode is the bond you meant.** A C–H stretch can outscore the bond of
interest.

### %V_bur and cone angle

State all four settings:

- the sphere radius;
- the centre: the metal, or for phosphines 2.28 Å from P along P→M, as in Tolman's protocol;
- whether hydrogens are included;
- the radii set (for example Bondi × 1.17, as in SambVca).

### Topology

The `tf_*` columns are the graph indices divided by heavy-atom count.

## 6. Extract and record

**Extract:**

- Gaussian route: `descripytor extractor -i input.json -o name -f feathers/`, or
  `Molecules(folder).get_molecules_features_set(json_dict)`.
- Metal complexes: `MetalComplexSet.from_xyz_dir(...).geometric_dataframe()`, plus the electronic
  block per molecule.

**Record next to the table:**

- `descripytor.__version__` and the git SHA;
- every non-default module setting: `STERIMOL_*`, `SUB_FRAGMENT_BOUND`, `DIPOLE_ABS`, the bond
  threshold;
- where the structures came from.

**Numbers from an earlier release** are reproduced by passing the earlier settings, not by changing the defaults.

## 7. Check the table

1. **Shape.** One row per molecule, and names that match the response file exactly. Quote names
   that contain commas.
2. **Bad columns.**
   - No missing values and no constant columns.
   - Near-duplicate pairs (r > 0.999) are flagged; the extractor prints both.
3. **Symmetric sets.** For C2-symmetric ligands every `_asym` column is about 0 by construction.
   Drop those columns: they hold only noise.
4. **Outliers.** Standardise each column. For any |z| > 6, open that molecule and look at the
   structure, not just the number.
5. **Spot-checks.**
   - Compare three molecules against an independent measurement: the `/feature_explorer` explorer, the
     published values, or a second implementation.
   - Rotate one input structure at random and re-extract. No geometric value may change.
6. **Column map.** Write down which physical quantity each column is (for example, `iso`/`aniso`
   are polarizability, not NMR shielding).

## 8. Hand-off

1. Save the table with a README: source, structures, settings, and
   the checks it passed.
