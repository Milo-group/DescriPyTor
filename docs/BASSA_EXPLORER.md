# BASSA model explorer

A page that puts BASSA's output beside the molecule. Columns are the models the sampler visited, rows are descriptors, and
the feature explorer on the right draws the selected descriptors on a compound of the set. It is a way of looking at a run,
not an analysis: nothing here changes a reported number.

![The BASSA model explorer: descriptors against the models the sampler visited, and the shown model drawn on a compound](images/readme/bassa_explorer.png)

## Start

    descripytor bassa_explorer

The page opens at http://127.0.0.1:8765/ with the example that comes with the package (17 ligands, 6 descriptors).
The molecule view is the package's [feature explorer](FEATURE_EXPLORER.md), framed in its figure mode.

| option | |
|---|---|
| `--cases DIR` | the folder where kept cases live; default `./bassa_cases`. Nothing is written inside the package |
| `--bassa-python PYTHON` | the Python that has `bassa_reg` (`pip install bassa-reg`), for running new tables. Default: the Python that runs the command, or the environment variable `DESCRIPYTOR_BASSA_PYTHON` |
| `--port`, `--host`, `--no-browser` | as for the other pages; the default port is 8765 |

Without `bassa_reg` the page still opens the example and every kept case; **New case…** then says what is missing.
Structures (`.feather`, `.log`, `.xyz`) are always read by the Python that runs the command.

The sampler set-up is in `M3_modeler/bassa_explorer/`: `engine.py` (priors, chains) and `core_fixed.py`, a corrected copy
of the sampler core. The released core is wrong when descriptors are correlated; `BASSA_CORE=package` runs it anyway.

The other tools of the folder are scripts, run with the Python that has `bassa_reg` and `BASSA_EXPLORER_CASES` set to
the cases folder:

    python M3_modeler/bassa_explorer/batch.py <folder of csv> --response y [--structures FOLDER] [--checks]
    python M3_modeler/bassa_explorer/influence_check.py <case> [descriptor]
    python M3_modeler/bassa_explorer/build_case.py attach <case> <folder> ["title"] ["response"] ["unit"]
    python M3_modeler/bassa_explorer/build_case.py export [cases]      writes explorer_standalone.html in the cases folder
    python M3_modeler/bassa_explorer/test_page.py                      the page's click-through test, in headless Chrome

## From the visual extractor, in one click

In `descripytor visual_extractor`, after **Extract CSV**, open the **Model** tab, choose the target column or paste the
measured values, and press **Model in BASSA explorer**. The table, the target and the feather folder are handed over,
this page opens in a new tab and the run starts.

![The visual extractor's Model tab with the button Model in BASSA explorer](images/readme/handoff.png)

## Your own table

**New case…** in the page:

- a CSV whose first column is the compound name, with one response column and the descriptors;
- optionally a folder with one `.feather`, `.log` or `.xyz` per compound. A compound takes the file with its exact name,
  else one that contains its name as a whole word, else the only file that carries its number (`L1` and `FL_lig_1_str`).
  Feathers and logs bring charges and the dipole, xyz the geometry only;
- columns named as the DescriPyTor platform names them (`B1_6-4`, `loc_B5_6-4`, `nbo_atom_5`, `mu_y[C4->N5|C6]`,
  `0_angle_[5, 6, 7]`) are drawn on the molecule. Other columns work everywhere else.

A full run is three chains of 20,000 iterations. **Add chance line and nested Q²** then runs, in the background, 100
shuffles of the response and one rerun per compound left out. **Keep this case** writes it beside the page and adds it
to the menu.

## Reading a case

- **Models**: the table. **Key models** merges models that differ only within a family of correlated descriptors and
  hides groups visited in less than 1 % of the run; "visited ≥" and the size, Q² and R² filters do the rest.
  Do not rank by Q²: that brings up large models the sampler hardly visited.
- **Recipes**: every model reduced to its families; inside a slot, the descriptors that fill it replace each other.
- **Network**, **Map**: which descriptors share models (green) or stand in for each other (red). Network circles can be dragged.
- **Compounds**: for the ten most used descriptors, how the frequency moves when each compound is left out.
  A frequency that moves is the sampler changing its choice among near-equivalent models; check the coefficient before
  reading it as "the effect rests on that compound" (`../bassa/cs2_l1_check.py` is the worked example).
- `influence_check.py <case>` prints, for the leading descriptor, the compound that moves its frequency most and whether
  its coefficient and the model move with it.
- Right-click a point in the plot to leave a compound out of the fit shown; **Rerun BASSA without them** runs the sampler
  again and shows every frequency before and after.
- **The plot**: above it, R², Q² and MAE of the shown model (Q² and MAE from fits that leave each compound out), then the
  external compounds and the ones you held out. The shaded band is one MAE either side of the diagonal.
  **Threshold analysis** asks how well one descriptor splits the compounds at a response value you set.
- **On the molecule**: Sterimol axes, charges and dipole components of the shown model, each with its value. A length, an
  angle or a dihedral is a dashed line through its atoms with the table's value (d, ∠, τ). A dipole feature shows its own
  arrow and label only. Untick **descriptors** to see the structures alone.
- **Compare**: shift-click compounds in the plot to stack them, superposed on the atoms the descriptors name; every
  compound carries its own values in its own colour.
- A column the page cannot place on the molecule: click its name under the plot and say what it measures.
- The browser's Back button undoes the last selection.
- **All cases** (top right) is one table of every case: leading descriptor, its share and leave-one-out range, p against
  chance, nested Q².
- The grey band around the chance line is how far the line could move with other shuffles. With 100 shuffles it is wide;
  `../bassa/null1000.py` computes it from 1,000 for a prepared case.

## Structures for a case kept without them

    <anaconda python> build_case.py attach <case> <folder> ["title"] ["response label"] ["unit"]

Reads the folder as for a new case. A folder of metal complexes with two donors also gets the atom roles, so the
metal-referenced descriptors are drawn on both arms. Where the report has a chance line and nested Q² for a case of
that name, they replace the ones from the batch run.

## Many tables

    <BASSA venv python> batch.py <folder of csv> --response y [--structures FOLDER] [--checks]

Each CSV becomes a kept case. `--checks` adds the chance line, the nested Q² and the fold frequencies (a hundred more
runs per table).

## Prepared cases of the DescriPyTor paper

They are not in the package. `build_case.py models <case>` and `build_case.py page <case>` rebuild them from the paper's
folder (`BASSA_EXPLORER_PAPER`, and `BASSA_EXPLORER_CS2_FEATHERS` for cs2); to look at cases that already exist, point
`--cases` at the folder that holds their `.js` files and `manifest.js`.
