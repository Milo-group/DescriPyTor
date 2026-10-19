# Changelog

## Unreleased

### Added: atom renumbering with Conformaty

- **`descripytor renumber FOLDER -r REFERENCE`** numbers every molecule of a folder like the reference with
  [Conformaty](https://github.com/Milo-group/Conformaty) (the lab's package, now installable:
  `pip install git+https://github.com/Milo-group/Conformaty`). The copies and `renumbering.json` (the `{old: new}` dict
  of every molecule) go to `<folder>/renumbered`; the originals are not changed. A progress bar with the time left, one
  line per molecule: what moved, which reference atoms are another element there, the RMSD of the overlay.
- **Align & renumber** at the bottom of the visual extractor does the same in the browser: molecules appear as they
  finish, with a time bar each and those to check first; a click shows the molecule renumbered on its reference (drawn
  faint), with renumbered atoms, other elements and the picks ringed; *Swap two atoms* edits a numbering by clicking;
  *Compare fits* fits a molecule every way side by side; *Write renumbered folder* and *Use it for extraction* switch
  the page to the result and check the picks' numbering there.
- The overlay's fit is a choice: Conformaty's losses (`rmsd`, the default; `msd`; `mixed`, its own default; `mad`;
  `sig`) or `kabsch`, on the heavy atoms as Conformaty's `align` fits. On the example set rmsd, msd and Kabsch agree on
  all nine; mixed and mad stall about 6 Å away on five (ortho and para), sig on all, also in Conformaty's own `align`,
  so a fit more than 0.5 Å worse than Kabsch is reported as stalled. The numbering itself keeps Conformaty's default.
- `Molecules.get_renumbering_dict(reference)` uses it too (no xyz export, no change of directory) and returns the dicts.

### Removed

- `MolAlign/`, the old copy of the alignment code that Conformaty replaces.
- A debugging `print` of the energy table in `save_to_feather`.

### Security

- **The visual extractor's server answers only pages opened from this computer.** It allowed every origin, so any
  website open in the browser could read the first lines of a text file on the machine (`/model/dataset`) or write a
  CSV anywhere (`/model/merge-outputs`). Pages on other sites (and `file://` or sandboxed pages, Origin `null`), a Host
  name that is not this computer's (DNS rebinding) and cross-site `<img>`/`<form>` requests now get 403; this computer's
  pages on any port still work. Started with `--host 0.0.0.0`, it answers only at 127.0.0.1/localhost: give the
  machine's own address to open it from another computer.
- A request can no longer name the folder put first on `sys.path` (`root`, `root_dir`).
- File names, SDF titles, xyz comment lines and the conformer-search folder are shown as text, not inserted as HTML.

### Fixed (visual extractor)

- **A descriptor that fails is reported.** A pick that failed on some or all molecules (Sterimol on an atom they do not
  have, say) left its columns out and the page said only "10 molecules × 6 features". The count is now followed, in
  amber, by each failure and on how many molecules it happened.
- **Two Ring vibration picks are two rings.** They left the page as `[13, 20]`, which the package reads as one ring
  through both atoms; they now leave as `[[13], [20]]`.
- **Groups are checked when committed.** A 1- or 5-atom bond angle, a 3-atom Sterimol axis and the like are refused
  with the reason; several atoms committed on a one-atom field (Charges) are several picks, not the first one only.
  Typed atoms read `1,4 4,1` as two groups (commas inside a group, spaces between groups, as in `WORKFLOW.md` §8.1);
  without a comma they are clicks, as before. An atom number the molecule does not have adds nothing.
- **Buried volume no longer falls back to atom 1.** With the engine on and no centre atom chosen, Extract says so.
- **Picks made on another molecule are named** when a molecule is loaded: the atoms past its end or of another element.
  **Check numbering** compares the folder with the molecule the picks were made on, not with `basic.feather` or the
  folder's first file.
- The Sterimol preview leaves out the **Drop atoms**, as the extraction does.
- **morfeus Sterimol pairs** stay with the morfeus suite; they were also sent to the DescriPyTor Sterimol engines.
- **outcomes.csv saved by Excel** ("CSV UTF-8", with a byte-order mark) joined no row; it is read now. A file in another
  encoding is reported instead of losing the extracted table, and one that matches no molecule is said to.
- The page is revalidated on every load (`Cache-Control: no-cache`): after an update the browser kept the old page.

### Removed

- The modeling panel's leftovers: about 190 lines of its JavaScript (its elements were gone), the routes only it called
  (`/model/dataset`, `/model/merge-outputs`, the `features_csv` form of `/model/run`) and the `/theta/check` alias.
- The layout cycler (its button was hidden), the `__XYZ_DATA__`/`__MOL_DATA__` placeholders nothing filled, unused CSS
  and helpers, and the in-process Tk fallback of the folder dialog, which hangs on a server thread.

## 0.3.4 — 2026-10-13

### BASSA model explorer

- **A plain page by default.** A case opens on two views, Models and Recipes, with Key models on, the rows grouped by
  family, and models of up to one descriptor per five compounds. The other views, the filters, the row orders, the
  links, pins, threshold analysis and the extra Sterimol drawings wait behind **advanced**. Older links open as before.
- **Families as ideas.** Every family, and every descriptor on its own, has a letter and a colour. A family is one
  block showing its four most used descriptors, with the rest behind "more"; lettered squares under each column and a
  line joining its dots show what a model combines; "with" names an idea's most frequent partners.
- **Exploring models.** A hollow dot swaps a family member into the model; recipe cards pick and step through their
  models; **Compare** keeps one model in the plot as hollow marks beside a second one.
- **Checks on a model.** A mark under every column and a list under the shown model: fit against leave-one-out, high
  or low only (a switch, not a scale), the leading descriptor against the chance line and across the leave-one-out
  reruns; the ends of the range and far-off compounds are shown and not counted. The thresholds are not calibrated.
- **On the molecule.** "on molecule" in a family heading draws its descriptors together; a slider sets the length of
  the dipole arrows.
- **First visit and keyboard.** Three steps above the table instead of a dialog; Enter and space work on recipes;
  explanations appear on keyboard focus; the family colours carry white letters at 4.9 to 1 or better.
- A family now counts as often as a model draws on it: one and two descriptors of the same family are two columns.
- **Less on the page.** A case opens with every family folded, one row per idea. The table's head is shorter with the
  key in its corner, the toolbar is one row, and the rerun column, the partner badges, the leave-one-out strips and the
  full list of checks wait behind advanced or a click.

## 0.3.3 — 2026-10-12

### Fixed

- **Measured outputs from another data set are never used.** The bundled Baptiste outcomes were pasted into the Model
  tab for the benzene example, and joined onto any extracted table whose compound names matched theirs (p-Cl, p-NO2,
  ...). After an extraction without outcomes the page also filled invented values. Now only an `outcomes.csv` in the
  extracted folder is joined, as the column `output`, and the paste box holds only what you paste.
- Pasted outputs are read strictly: a line that cannot be read stops the paste (it was skipped, moving every later value
  up a row); the count must match the rows; a pasted name that matches no row is an error; names are matched ignoring
  case, `.feather`, and `-` against `_`; tab, semicolon and space separators work. A target column chosen together with
  pasted values is refused instead of quietly preferring the paste.
- An engine that failed during extraction is named beside the result; the status used to stay green.
- Every picked Dipole frame is exported; with two picked, only one left the page.
- Buried volume: the atom and the radius that are drawn are the ones exported.
- Loading a molecule with fewer atoms than an existing Dipole pick no longer breaks the page; a vibration that was
  playing stops; an unfinished group does not carry over. Picks are restored with their molecule after a reload.
- The "labels" box of the feature preview works; the numbering check includes Dipole atoms given as groups.
- Feature explorer: a saved session or permalink kept an axis on its stacked structure by an id that changes on reload,
  so it was re-measured on the active structure with the other one's atom numbers. Sessions now name the structure.
  A layer naming a structure that is not loaded is reported, not moved onto another molecule. Sheet panels no longer
  share their layers and look. A dipole given as three components, a linear molecule, a charge site on a structure
  without charges and a session without layers no longer stop the page. The Python check is not run for axes it does
  not compute (stacked, blocked). A framed figure neither reads nor stores a look in the browser.
- BASSA model explorer: a column named `B1`, `L`, ... was labelled as the Wu–Doyle ligand's ("on Ni→P") in every table;
  atom names with two-letter elements (`q_cm5[Cl12]`) are read; a column whose whole spread is below 1e-4 was sent to
  the page as zeros; the Q² and R² sliders at zero hid models with negative values; header names with spaces and
  cells holding `inf` are handled; an atom number beyond a compound's atoms no longer stops the page.
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
