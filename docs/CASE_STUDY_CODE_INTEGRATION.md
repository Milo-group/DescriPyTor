# Integrating the case-study code into DescriPyTor

2026-09-16. What the paper's analysis scripts do that the package cannot, and the
order to bring it in.

Surveyed: `MolFeatures` (M1/M2/M3/utils), the paper's working scripts
(`descripytor_paper_v3_submission/notes/analysis`, 75 files) and the deposit copies
(`zenodo_archive/scripts`, 48 files).

## 1. What is already in the package

| Area | Module | State |
|---|---|---|
| Feature extraction | `M2_data_extractor` | Gaussian/feather parsing, Sterimol, %V<sub>bur</sub>, `MetalComplex`, `LigandTopology`, `XtbSinglePoint`, cube and morfeus paths, GUI |
| Model search | `M3_modeler/modeling.py` | `LinearRegressionModel`: exhaustive k-combination search, VIF, spike-and-slab, `_analytic_loo_linear`, SQLite persistence |
| Nested CV | `M3_modeler/nested_cv.py` | k-fold / grouped nested CV, RSS shortlist then validation pick, `compare_representations`, Dietterich 5×2cv |
| Reporting | `M3_modeler/plot.py` | per-model PDF reports, scatter/violin/SHAP/threshold pages, `y_randomization_check`, correlation heatmaps |
| Structure views | `utils/visualize.py` | plotly molecule views, `plot_b1_visualization`, `plot_L_B5_plane` |

So "nested CV" and "plotting" exist. What is missing is narrower and specific: the
**validation protocol the case studies actually report**, and a **Sterimol geometry
interface** the figures can be drawn from.

## 2. What the scripts have and the package does not

Counted across the 123 script files (a capability is counted once per filename):

| Capability | Scripts that re-implement it | Package equivalent |
|---|---|---|
| Closed-form LOO through the hat matrix | **41** | `_analytic_loo_linear` (one combination, private) |
| Exhaustive combination scan, batched | **17** | combination loop exists, not batched, no LOO ranking |
| Nested selection: outer LOO, inner exhaustive | **16** | k-fold with an RSS shortlist — a different protocol |
| Response-permutation null | **20** | `y_randomization_check` (fixed features, no re-selection) |
| Size-matched constraint null over all pairs | **12** | none |
| Group holdout / transfer with a random-holdout percentile | **15** | none |
| Paired sign-flip or bootstrap model comparison | **5** | Dietterich 5×2cv only |
| Ridge with an unpenalised intercept on standardised slopes | **7** | scattered inside model classes |
| Sterimol rotation-scan kernel | **5** | 4 copies inside the package as well (see §4) |
| Trivial baselines: mean, best single, class identity | 2 | none |

Duplication is not only between the package and the scripts. 45 filenames appear in
both script folders, and **17 of those differ in content** — including
`cs2_verify.py`, `holdout.py`, `constraint_null.py`, `nested.py` and `sterimol.py`.
The deposit and the working copies have drifted apart.

## 3. Proposed additions

Two new modules carry almost everything. Both are small because the maths is already
written; the work is choosing one signature and deleting 40 copies.

### `M3_modeler/search.py` — the selection protocol

```python
loo_q2(X, y, combos, alpha=1.0, chunk=400)      # batched closed-form ridge LOO -> Q2 per combo
best_combos(X, y, k=3, alpha=1.0, top=1)        # exhaustive scan, adjusted Q2 ranking
nested_loo(X, y, k=3, alpha=1.0, mask=None)     # outer LOO, inner exhaustive, per-fold scaler
```

`nested_loo` returns the nested R², the per-fold selections and the stability
summary the SI prints (distinct equations, most frequent count, in-fold Q² spread).
`mask` is what makes the constraint null one line.

### `M3_modeler/validation.py` — the controls

```python
selection_null(X, y, n_perm=1000, seed=...)          # floor (95th pct) and p for a search maximum
family_null(X, y, family, n_perm=1000, require=None) # is this family in the winner beyond chance
constraint_null(X, y, pairs="all")                   # nested R2 of every size-matched pair constraint
paired_sign_flip(y, pred_a, pred_b, n_draws=20000)   # closer-on count, MAE difference, p
group_holdout(X, y, group, n_random=300)             # fit outside, predict the group, percentile
baselines(X, y, classes=None)                        # mean, best single descriptor (OLS), class identity
```

Two implementation notes that decide whether these are usable interactively:
`selection_null` and `family_null` permute only the response, so the hat matrix of
each candidate model is permutation-invariant — build it once per triple and
contract all 1000 permuted responses through it in one einsum. That is 20 seconds
on a 58-descriptor pool, against the ~35 minutes a naive re-fit implies, and it is
the difference between a control you run every time and one you schedule.
`group_holdout` must take the random-holdout draws as an argument, so two
equations can be scored against the same splits; percentiles from separate draw
sets are not comparable and have already caused one wrong comparison.

### `M2_data_extractor/extractor_utils/sterimol_geometry.py` — what a figure needs

```python
sterimol_frame(coords, bonds, base_atoms, radii="CPK")
# -> atoms in the Sterimol frame, the B1 normal, the B5 atom, B1/B5/L/loc_B5/theta/phi/rho
```

plus `get_sterimol_plot_data` (already written in the `gaussian-fairshare-jobs`
worktree, missing from main) and two drawing helpers: the down-axis projection
`utils/visualize.plot_b1_visualization` already covers, and the **slice view**
(B1 plane edge-on, θ true, the axis carrying loc B₅ and L), which nothing covers.

### `M1_pre_calculations/complex_builder.py` and `xtb_runner.py`

`build_general.py` — the metal placement CS3 depends on — lives in the paper v2 notes
and is **not in the package at all**, so `MetalComplex` can featurise a complex the
package cannot build. The xTB drivers (optimise, single point, constrained
optimisation, parse energies/dipole/orbitals) are likewise re-written in three
scripts. Both belong in M1.

## 4. Consolidation the survey exposed

- **Four Sterimol kernels ship in the package**: `extractor_utils/sterimol_utils.py`,
  `sterimol_standalone.py`, `cube_sterimol.py` and `metal_complex.sterimol`, plus
  `sterimol_utils_old.py`. Their defaults were only aligned by hand on 2026-09-15
  (1° scan). One kernel with options, and thin wrappers, or they will drift again.
- **`main` and the worktree disagree** on `sterimol_utils.py`: main lacks
  `STERIMOL_CIRCLE_POINTS`, `calc_loc_b1`, `prepare_sterimol_inputs` and
  `get_sterimol_plot_data`. Reconcile before porting anything else.
- Dead weight to drop: `M2_data_extractor/_tmp_*.py`, `M3_modeler/plot_backup.py`,
  `extractor_utils/sterimol_utils_old.py`, `tmp/goat_table_cluster/pkg/`.

## 5. Order of work

**Phase 1 — search and validation (highest value).** Port `search.py` and
`validation.py`. Every function already has published numbers to gate against, so
each port lands with a regression test rather than a hope:

| Test | Expected |
|---|---|
| `nested_loo` on the deposited cp pool | 0.602, 9 equations, most frequent in 15 of 30 |
| in-fold Q² spread, same run | 0.802 ± 0.012 (0.773–0.827) |
| `constraint_null` on the 1° cp pool | metal pair 243rd of 1653, p 0.174, median 0.692 |
| `paired_sign_flip`, cc and cp against the published models | p 0.85 and 0.87 |
| `baselines` on cp / cc / oa | 0.648 (0.330), 0.680 (0.168), 0.585 (0.147) |
| `group_holdout`, cp BOX → PyOx | MAE 0.390 (deposited φ pool, 99th pct); published DFT triple 0.248, 54th |
| `selection_null` floors | cc 0.39, cp 0.36 |

**Phase 2 — Sterimol geometry and the slice view.** Reconcile the two
`sterimol_utils.py` copies, add `sterimol_frame`, move the figure code
(`figure_three_cases.py`, `build_figure_page.py`) behind it. Gate: the three figure
molecules reproduce their matrix rows, which they already do.

**Phase 3 — builder and xTB runner.** Move `build_general` in, wrap the xTB calls,
and expose `complex_from_ligand(xyz, metal, ancillary)` so the CS3 route runs from
the package alone. The post-relaxation bite check added on 2026-09-16 belongs with it.

**Phase 4 — kernel consolidation and cleanup.** One Sterimol kernel, delete the
`_tmp_`/`_old`/`_backup` files, and regenerate the deposit scripts from the package
so the 17 drifted duplicates collapse to one source.

**Phase 5 — the paper repo becomes thin.** `notes/analysis` scripts keep their
docstrings, gates and data paths, but call the package instead of re-deriving the
maths. That is what stops the next 41st copy of the hat matrix.

## 6. What not to move

- One-off diagnostics that answered a question once and are documented in
  `notes/`: `cs12_mismatch_*`, `cs12_last_atom`, `fix_009_chelate`, `h_arm_azimuth`,
  `frame_turn_1deg`. Keep them in the paper repo as the record.
- Manuscript plumbing: `wordcount.py`, `audit_refs.py`, `make_overleaf_package.py`,
  `prep_acs_figures.py`, `move_sampler_si.py`.
- Anything that hardcodes a dataset (`tools_*`, `extract_doyle_features.py`) unless it
  is rewritten against a generic loader.
