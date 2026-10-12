# Changelog

## 0.3.3 — 2026-10-12

### Fixed

- Visual extractor: the dipole arrow was switched on whenever the loaded file carried a dipole, and its origin jumped to
  whatever atom was clicked last, in any field. It is now off until you tick **Dipole** or pick a Dipole frame, and it
  starts at the Dipole pick's origin, or at the centre of the molecule when there is none.
- Removing a Dipole pick (or **Clear picks**) left its three atoms in the basis boxes, so the frame stayed drawn and
  stayed in the exported configuration. The boxes now follow the picks.
- Every atom of a group that is still being built is ringed, not only the last one clicked, and nothing stays ringed
  after the group is committed.
- Clicking an atom that is already in the unfinished group takes it out again; the same atom could be picked twice
  into one group.
- Ctrl+Z with no unfinished group takes back the last committed group of the active field. Pairs and triplets commit by
  themselves, so before this it did nothing after a wrong pick.

### Removed

- The Streamlit web app, the Docker files (`docker-compose*.yml`, the image's `Dockerfile`), `make_picker.py`, the second
  copy of `atom_picker.html` in the examples folder, the `webapp` extra and the empty `work/` folder. The "Show 2D graph"
  button went with them: only a page built by `make_picker.py` ever had a graph to show.
- The form-based page at `/forms` (`feature_extraction_gui.html`) and the eight server routes only it used (`/sterimol`,
  `/charges`, `/dipole`, `/vibrations`, `/molecule_info`, `/steriplot`, `/batch_sterimol`, `/features_set`). The visual
  extractor and the feature explorer do what it did.
- The `/modes` route of the visual extractor's server. The vibration-mode viewer itself stays as the standalone tool
  `python utils/mode_viewer.py --serve`.

## 0.3.2 — 2026-10-12

### Fixed

- A visual extractor server left running from before an update was reused by the next `descripytor visual_extractor`
  or `feature_explorer`. Its pages come from disk and were new, its routes were old, so the links to the feature
  explorer and the BASSA model explorer gave "Not Found". A new start now recognises the older server and says how
  to stop it.
- The same for the BASSA model explorer: its server answers `/status` with its version, and the `bassa_explorer` command,
  the visual extractor's link and the one-click hand-off refuse a server on its port that is older or is something else,
  and say how to stop it. A second `descripytor bassa_explorer` on a current server opens it instead of failing.
  From the BASSA page, the links to the other two pages refuse an older visual extractor server in the same way.

## 0.3.1 — 2026-10-12

### Changed

- The extraction engine behind **Extract CSV** (`descriptor_extractor.py`, `feature_packages.py`) moved from the examples
  folder into `M2_data_extractor/`: an installed package did not have it, so Extract CSV failed outside a clone.
  The old files forward to the new ones.
- BASSA model explorer: the columns `descripytor extractor` and the visual extractor write are drawn on the molecule:
  `dipole_x/y/z_<origin>-<y>-<plane>` and `total_dipole_…`, `loc_B1_`, `B1_tangent_atoms_`, `B1_B5_azimuth_`,
  `B1_B5_acute_`, and the vibration columns (`Frequency_Stretch_a_b`, `Frequency_Bending_a-b`, `cross_`/`para_` and
  their angles). `iso`, `aniso` and `energy` get readable names; they have no atoms to be drawn on.

### Removed

- The Tk desktop app left in `descripytor/cli.py` after the `gui` command went (about 1,300 lines), the interactive
  prompts whose command was commented out, and the `gui` extra (`customtkinter`).
- Nine pictures and one recording under `docs/` that no page referenced.

## 0.3.0 — 2026-10-12

### Added

- Visual extractor: **Model in BASSA explorer** in the Model tab. One click sends the extracted table, its target
  and the feather folder to the BASSA model explorer, which starts the run (`/bassa_explorer/handoff`).
- The three pages link to each other; `/visual_extractor`, `/feature_explorer` and `/bassa_explorer` mean the same
  page on both servers, and a link starts the other server when it is not running.
- `pip install "descripytor[bassa]"` installs the sampler (`bassa-reg`). The corrected sampler core is the default.
- `descripytor --help` describes every command, with examples. README rewritten around the three pages.
- `descripytor bassa_explorer` (alias `bassa`): the BASSA model explorer, now in the package
  (`M3_modeler/bassa_explorer/`, [docs/BASSA_EXPLORER.md](docs/BASSA_EXPLORER.md)). It opens with a bundled example; kept
  cases go to `--cases DIR` (default `./bassa_cases`); new tables need `bassa_reg` in the Python given by `--bassa-python`.
  The corrected sampler core (`core_fixed.py`) comes with it.
- Feature explorer: a figure-only embed mode (`?embed=figure`): the framing page posts sessions and is told which atom
  was clicked (`feature:atom`). A mark can ring its atoms (`ring`), and `L.arrows` draws arrows between points.

### Changed

- `pyarrow<26`: pyarrow 26 needs NumPy 2, and the package needs `numpy<2`; without the pin a fresh install could not
  read a feather.
- The Sterimol theta explorer is renamed **feature explorer**: `descripytor feature_explorer`, `/feature_explorer`,
  `M2_data_extractor/feature_explorer/`, [docs/FEATURE_EXPLORER.md](docs/FEATURE_EXPLORER.md), sessions saved as
  `.featurescene.json`. The descriptor θ keeps its name. The old names still work: `descripytor theta`, `/theta`
  (redirects), `.thetascene.json` files, the `theta:` message prefix, and the folder `M2_data_extractor/theta_explorer/`
  (forwarding `data.py`, `build.py`, `sterimol.js` and a copy of the built page).
- `descripytor visual_extractor`: the atom picker draws its molecule with the feature explorer instead of 3Dmol. Every
  selection is drawn by the explorer's own layers, with the values the package computes; the picker's own approximate
  Sterimol and %V_bur previews are removed (they disagreed with the package). The exported answers are unchanged.
  Screenshot now saves an SVG. 3Dmol remains only in the conformer viewer. The picker's server address follows `--port`.

- `get_ring_vibrations`: a mode is considered only if at least half of its heavy-atom motion sits on the
  requested ring (`min_ring_share=0.5`; `None` restores the old pick). In a molecule with a second ring
  the old pick could be that ring's mode: on the 33 Lucas ligands 19 pyridine-ring picks change, on the
  CS1 aryl ring 2 of 17. Rings that are the only ring of their kind are unaffected.
- `get_ring_vibrations` accepts a five-membered ring: `cross` and `para` are then the higher and lower of
  its two ring C-C stretching modes (new `ring_stretch_pair`), with NaN angles.
- `get_bend_vibration`: a mode is considered only if at least 30% of its motion sits on the two end atoms
  and their shared centre (`min_share=0.3`; `None` restores the old pick), and the row is NaN when none
  qualifies. The old pick could jump between unrelated modes across a set of molecules.

- `MetalComplex`: the arm's (C\*, R) pair no longer depends on the atom order of the file
  (`STEREOCENTRE_RULE = "alpha"`, [docs/STERIMOL_FIXES.md](docs/STERIMOL_FIXES.md) #12). Values change for
  arms with several equally sized substituents (`a_R_Cstereo_N` on unsubstituted arms) and for P-aryl
  arms, which now have no C\*→R block. `descripytor.compat.paper_v3()` and `STEREOCENTRE_RULE = "first"`
  restore the old pick.
- `sterimol`: φ (`angle`) follows `STERIMOL_THETA_RULE`; under `"soft"` it is averaged over tied B1
  directions like θ.

### Added

- `descripytor theta`: opens the Sterimol θ explorer in the browser, as `descripytor visual` opens the atom
  picker (`--mol`, `--preset`, `--port`, `--no-browser`); it reuses a server that is already running.
- **Cluster pipeline** (`descripytor pipeline`, `M1_pre_calculations/pipeline/`, [docs/PIPELINE.md](docs/PIPELINE.md)):
  SMILES → numbered starting geometries (organic with a fixed core numbering, one metal on one
  donor, or a chelate with ancillaries in `MetalComplex` order) → chained SGE array stages
  (ORCA GOAT, UMA, GFN2-xTB, ORCA DFT, Gaussian), each molecule's stages chained with `-hold_jid`, status, retry,
  watch, fetch, `adopt` (finished outside results join a run as a completed stage) and `reset`. Every
  stage checks program errors, imaginary frequencies (xTB `--ohess` by default, with one restart along
  the mode) and the structure against the build (bonds, clashes, donors, bite angle); a failed
  structure never reaches the next stage. Driven by one protocol JSON.
- [docs/THETA_EXPLORER.md](docs/THETA_EXPLORER.md): a manual for the Sterimol θ explorer with 15 screenshots,
  regenerated by `docs/images/theta_explorer/make_images.py`.
- [docs/EXTRACTION_PROTOCOL.md](docs/EXTRACTION_PROTOCOL.md): what to decide and check from
  structures to a feature table.

## 0.2.1 — 2026-09-28

Metal-complex (`MetalComplex`) values change; nothing else does. `descripytor.compat.paper_v3()`
restores the old definitions, and each one also has its own switch.

- **Substituent-arm fragment** (`sub_*`): the C\*→R walk stops at the donor's own ring, so a
  substituent fused to that ring no longer reaches the backbone and the other donor.
  `SUB_FRAGMENT_BOUND = "donor"` for the old walk.
- **θ at tied B1 directions**: averaged over every direction with weight exp(−(width − B1)/0.02 Å),
  so chemically identical groups (tBu: 4.7 / 38.8 / 45.5° before) get one θ. B1, B5 and L are
  unchanged. `STERIMOL_THETA_RULE = "scan"` for the old value.
- **`mu_desym`, `mu_outofplane`** are magnitudes; the sign depended on which donor was listed
  first. `DIPOLE_ABS = False` for the signed values.

With the defaults, the 16 Sterimol/θ columns of the paper's corrected cyclopropanation table
(`cp_fix_soft.csv`) rebuild to within 5e-7 from the deposited geometries.

## 0.2.0 — 2026-09-27

Descriptor values change in this release; the list is below. `descripytor.compat.paper_v3()`
restores the 0.1-era definitions the DescriPyTor paper was computed with (git tag `paper-v3`).

### Changed definitions

- **Bonds**: a non-metal pair is bonded below 1.15 × the sum of covalent radii instead of a
  flat 1.82 Å. S–CF₃, P–C, C–Br, C–I, Si–C and S–S bonds are kept; an agostic metal···H
  contact is not a bond. `threshold_distance=1.82` gives the old rule.
- **Sterimol**:
  - the B1 plane normal is no longer mirrored, so θ (`B1_B5_angle`) no longer depends on how
    the frame is turned about the axis;
  - every copy scans every degree;
  - `MetalComplex` reports θ instead of φ, and starts its scan from the molecule rather than the
    lab frame, so its values no longer depend on how the input file is oriented.
- **Vibrations**: one window everywhere: stretch 1400–3500 cm⁻¹, bend from 1300 cm⁻¹. The 3D
  picker's old bend default (3000) only ever found X–H stretches.
- **Rings**: ortho/meta/para are measured from the atom you give (they were measured from another
  atom in most cases), and the two ring angles are real values (they were always 90°).

Details and the measured effect of each: [docs/STERIMOL_FIXES.md](docs/STERIMOL_FIXES.md),
[docs/FEATURE_FIXES.md](docs/FEATURE_FIXES.md).

### Added

- `M3_modeler.ridge_search`: exhaustive k-descriptor search by closed-form ridge / OLS LOO,
  nested LOO and the constraint null.
- `M3_modeler.validation`: `selection_null` (the search on permuted y, 1000 permutations in about
  a second) and `paired_sign_flip`.
- `descripytor.compat.paper_v3()`.
- The Sterimol / θ / cone-angle / %V_bur explorer (`/theta` in `descripytor visual`), gated
  against the Python.
- `Molecules({label: folder, ...})`: several folders pooled into one set.
- New columns: `loc_B1`, `B1_tangent_atoms`, `B1_B5_azimuth`, `B1_B5_acute`.
- [docs/WORKFLOW.md](docs/WORKFLOW.md): the workflow step by step, and every option.

### Fixed

- `descripytor extractor` crashed on Windows consoles before writing its CSV.
- `descripytor model -t target.csv` ignored the target file; features and target are now joined
  on the molecule name.
- The model search waited for input with no terminal attached (cluster jobs, pipes).
- A stretch or bend with no frequency window silently dropped its columns.
- Tuples of atom pairs gave one meaningless number.
- A bend given as a triplet a-b-c (as the forms GUI sends it) used the wrong atoms.
- Opening a feather file needed dash and torch.

## 0.1.0 — 2026-08-30

First release on PyPI: the pip package and the `descripytor visual` first-run path.
