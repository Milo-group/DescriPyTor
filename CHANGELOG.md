# Changelog

## 0.5.3 — 2026-11-03

- **The arm's stereocentre is read at the carbon next to the donor** (`metal_complex.STEREOCENTRE_RULE = "alpha_carbon"`).
  The previous rule took the largest substituent anywhere on the donor's ring, so a 5,5-diphenyl oxazoline was read
  through a phenyl about 5 Å from the metal. R is now the group on that carbon with the larger B1. `"alpha"` restores the
  previous pick.
- **The metal bonds only its co-ligands and two donors** (`metal_complex.METAL_BONDS = "donors"`). The covalent rule also
  bonded it to nearby carbons (an imine carbon at 2.70 Å), aryl chlorides and a third donor, and each such edge cut the arm
  fragments. A mer tridentate is now read at its outer two donors, with a warning. `"all"` restores the previous graph.
- Arm values change for metal complexes with such a contact or such a substituent.

## 0.5.2 — 2026-11-01

- Built-in structures renamed; the example run script and the mode viewer removed (ORCA reading stays, in `utils/orca_text.py`).

## 0.5.1 — 2026-11-01

- Package description and documentation wording.

## 0.5.0 — 2026-10-31

DescriPyTor is now a descriptor-extraction package. It reads Gaussian, ORCA, xTB and `.xyz` output and returns a
feature table defined atom by atom. Earlier releases are withdrawn.

- **Descriptors:** Sterimol (L, B1, B5, loc_B1, loc_B5, θ and the in-plane azimuth, from tabulated radii or a density
  cube), charges and charge differences, the dipole in a frame you define, bond lengths, angles and dihedrals,
  stretching, bending and ring vibrations, buried volume and cone angle, and metal-complex descriptors.
- **Pages:** the visual extractor (`descripytor visual_extractor`) picks atoms on the molecule and extracts the table
  for a folder; the feature explorer (`descripytor feature_explorer`) draws every descriptor on the molecule.
- **Command line:** `extractor`, `logs_to_feather`, `sterimol`, `cube`, and `pipeline` (SMILES to structures and
  chained cluster stages).
- Bonds follow covalent radii (Pyykkö, × 1.15) with a separate metal rule; the B1 scan runs every degree.
