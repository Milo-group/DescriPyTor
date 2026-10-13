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

A case opens on the plain page: the views **Models** and **Recipes**, one column per combination of families, the
rows grouped by family, and only models of up to one descriptor per five compounds (rounded down); **Filters** changes the sizes and sets a least Q² or R². The **advanced** box
adds the rest (listed at the end of this section). A link made before the plain page existed opens in advanced.

- **The table**: a column is a model, a row a descriptor. Over each bar: its 5-fold Q² (the grey bar behind is R²). Click a
  column to draw its model; the arrow keys move along. **A column is** switches between a combination of families
  (models that differ only inside a family counted together, combinations under 1 % of sampled models left out) and a
  single model. A family counts as often as a model draws on it. Columns stand strongest first: by the 5-fold Q² of
  the column's most visited model, ties by visits; the number under a column is its rank. **?** in the header opens
  with **How this case is split**: which models are shown, the order, how the folds are dealt, and every family with
  the reason each member is in it.
- **Families**: descriptors that stand in for one another. A descriptor joins the most used descriptor for which its
  column is correlated at 0.7 or more, or at 0.5 or more when the sampler clearly takes one or the other; it never joins
  one the sampler uses it together with more than chance. The sampler's evidence is the 2x2 table of the two
  descriptors over the kept steps (log odds ratio, a tenth of the steps counted as independent, judged against the
  unrelated pairs of the same run). The "?" at the end of a row gives the reason with its numbers.
- **Ideas**: every family, and every descriptor on its own, has a letter and a colour (A is the most used). A family is
  a coloured block named after the descriptor the model in hand takes from it. It opens when that model uses it, shows
  its four most used descriptors, and keeps the rest behind **more**. **Fold families** leaves one row per idea.
- **Dots**: one per descriptor of the model. A hollow dot is used by another model of the same column; click it to swap
  it in. A folded family carries its dot on the heading. Lettered squares under a column name its ideas; a line joins
  its dots.
- **The mark under a column** (! or ‼, none when nothing is flagged): how many checks flag the model. Counted: R² more than 0.15 above the
  leave-one-out Q²; a leave-one-out Q² within 0.05 when every descriptor is cut to high or low at its
  median (a switch, not a scale); the leading descriptor among the three most used in fewer than 80 % of the
  leave-one-out reruns. Shown and not counted: the ends of the range predicted from the rest, compounds beyond two
  standard deviations, and the leading descriptor against the chance line (the same for every model of a case, and
  said once as a pill beside the case name). The thresholds are choices, not calibrated.
- **The tiles above the plot**: the model in hand, kept fixed (5-fold Q², with R² as the grey part of the line);
  **BASSA's averaged model**, sampler rerun in every fold, which is the number to quote (each dot is one deal of the
  folds; a case with one deal says so); the typical miss (MAE against the range of the response); and, for a case run
  through `ensemble_protocol.py`, the share of held-out values inside the 90 % interval with the compound furthest off.
  These come from the case file (`cv`). Every small **?** opens a drawing of that number.
- **Weights and equation**: under the plot, a bar per descriptor (its weight per standard deviation, ▲ raises and ▼
  lowers the response). **as an equation** unfolds the model on normalized descriptors, and in their own units.
- **Checks**: three chips for the counted checks (fit against test, a scale and not a switch, holds without one
  compound); hover one for its drawing, **all checks** for the list.
- **Averaged model**: the average over every sampled model, in the plot and as its five largest weights (`ens` in the
  case file). **Compare models** keeps one model as hollow marks beside a second. **Design**: a slider per descriptor
  (not on the averaged model). **Export**: the picture as SVG, or the model as text.
- **Recipes**: every model reduced to its families; click a descriptor to draw the model of that card that contains it,
  the arrows step through the card's models, and "together" is the card's models averaged by visits.
- **On the molecule**: Sterimol axes, charges and dipole components of the shown model, each with its value. A length, an
  angle or a dihedral is a dashed line through its atoms with the table's value (d, ∠, τ). **dipole arrow ×** and
  **Sterimol drawing** appear when the model has such a term. Untick **descriptors** to see the structures alone.
- **Header**: **All cases**, **More** (new case, keep this case, add the chance line and leave-one-out reruns, save as
  one file, copy link), **Reaction**, and **?**.

What **advanced** adds:

- **Network**, **Map**: which descriptors share models (green) or stand in for each other (red). **Compounds**: for the
  ten most used descriptors, how the frequency moves when each compound is left out. A frequency that moves is the
  sampler changing its choice among near-equivalent models; check the coefficient before reading it as "the effect
  rests on that compound". `influence_check.py <case>` prints this for the leading descriptor.
- Sorting the models by visits, Q² or R², and the filters (these use the leave-one-out Q², not the 5-fold one over the
  bars). Do not rank by Q²: that brings up large models the sampler hardly visited.
- Rows by frequency or by what goes with a selection; the links as arcs; pins; the count of merged models.
- Beside each descriptor, the leave-one-out reruns (top three, range); in family headings, the partner ideas and
  **on molecule**, which draws a family's descriptors together.
- Under the plot: where each descriptor's share goes without single compounds; **Threshold analysis**; point names.
- Right-click a point to leave a compound out of the fit shown; **Rerun BASSA without them** runs the sampler again.
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
