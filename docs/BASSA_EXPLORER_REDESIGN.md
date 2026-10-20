# BASSA explorer: the simple page (plan, 20 October 2026)

The PI's request: an interface like Apple's, limited to exactly what the user needs to analyse and view the models; clean,
little text, everything extra behind a toggle or advanced settings; fewer options so the user finds the way.

Today the plain page of one case shows 71 controls and about 690 words at once. The analysis stays; what changes is what is
shown first. Nothing is deleted: every current tool stays reachable through the ⋯ menu or **Advanced**.

## Decisions (Eden, 20 Oct)

- Two simple modes behind one switch in the top bar: **Overview** (named Design until 20 Oct; a chemist choosing what to make) and **Models** (a modeller
  checking the models).
- The table of visited models becomes a **compact strip of the top models** (share of sampled models, Q²) under the drivers;
  the full table is Advanced.
- Besides *Why this prediction*, the simple page keeps only **How to make it better**. Design sliders and candidates, comparing
  two compounds and *Your families* move to Advanced.
- **"Averaged model" is renamed "ensemble model"** on the page and in the docs.

## The page answers three questions, in order

1. **Is the model any good?** One number and one verdict in the top bar: "Ensemble model · Q² 0.61 · better than chance". The
   method (deals, chance line, tests) opens from that number.
2. **What drives it?** The three to five families as bars with their direction ("lower B1 on 6→4 is better").
3. **This compound, and what to make next?** The plot and one compound card: the molecule, *Why*, *How to make it better*.

## Layout

| area | Design mode | Models mode |
|---|---|---|
| top bar | case picker · verdict · Design / Models · ⋯ | same |
| left | **Drivers** (family bars), then the top-model strip | the top-model strip first (share, Q², R²), then the families |
| centre | predicted vs measured, ensemble model; click a point to open the compound | the plot of the model clicked in the strip, with "← ensemble model" to return |
| right | compound card: molecule, *Why*, *How to make it better* | compound card: molecule with the model's descriptors, the model's weights |

- **Top bar:** case picker, verdict, the mode switch and one ⋯ menu: All cases, New case, Export (picture, model,
  predictions), Copy link, Print, the visual extractor and feature explorer, Help, **Advanced**.
- **Interaction instead of controls:** clicking a family colours the plot by it and draws it on the molecule; clicking a
  model shows it, and the chip returns to the ensemble model; the plot needs no "all models / the selected model" switch.
- **Molecule settings** (descriptors on/off, dipole arrow size, Sterimol drawing, reset view) sit behind a small gear.
- **Text:** the explanatory sentences ("New here?", "over each column…", the method paragraphs) are replaced by a one-time
  three-step tour and **?** popovers that open only when asked.

## Advanced (one switch, remembered)

The full model table and Recipes, *Your families*, Filters, View (column type, family rule, fold), the Network, Map and
Influence views, ranking and colour options, "More about this model" (every weight, what is left over, how it was checked),
Design sliders and candidates, comparing two compounds, the threshold analysis, the per-test evidence badges, the
descriptor search (also on Ctrl+K in both).

## Polish (Eden, 20 Oct)

Speed and motion
- **A case loads when it is opened**, not all of them at start (today every load fetches all 14 case files, Denmark's 6.8 MB
  among them, before anything is drawn).
- **No layout jumps**: panels keep their size when the mode, the compound or the model changes; content fades in place; the
  molecule keeps its orientation from one compound to the next.

Less on screen
- **One verdict chip** (green, amber or red) instead of the three badges ("#2 of 4", "! 2 of 3 tests passed", "↓ lower is
  better"); the details open from it.
- **Say each thing once**: the family bars stay, their repeat inside the table goes; error bars for the selected compound
  only (the rest on hover); plot resize handles on hover.
- **The molecule's switches** (the unlabelled "on" boxes, dipole size, Sterimol drawing, reset) go into one gear menu; only
  the atoms a descriptor uses are labelled.

Consistency
- **One visual system**: one accent colour for controls, family colours for families only; three type sizes; an 8 px
  spacing grid; numbers always as "0.91", tabular, Q² and kcal/mol to two decimals.
- **One vocabulary**: ensemble model, family, driver, share of models; no "deal" or "kept fixed" on the simple page; the help
  in one **?** drawer instead of a "?" after every heading.
- **Short links**: only the settings that differ from the default go into the address (today about 40 do, and old defaults
  travel with them, as `fb=r` did).
- **Plain messages for missing pieces**: a case without structures, without ensemble predictions, or with its server
  stopped says so in one line where the panel would be.

## Steps

1. **Layout** — the simple page as the default with the two modes; the current page is Advanced. The self-test covers
   both modes and checks that every tool of today is still reachable.
2. **Rename** — "averaged model" → "ensemble model" in the page, the tooltips, the self-test and `BASSA_EXPLORER.md`.
3. **Visual pass** — one accent colour, more white space, a fixed type scale, fewer borders, the tour.
4. **User test** — two or three lab members, three tasks (find the best compound, explain one prediction, propose a
   candidate); count clicks and time; trim again.

Done when: a new user can answer the three questions without opening Advanced, the default page shows at most ~20 controls
and ~150 words, and the self-test passes on every case in both modes.
