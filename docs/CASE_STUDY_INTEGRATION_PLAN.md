# Plan: bring the case-study code into DescriPyTor

2026-09-27. This supersedes `docs/CASE_STUDY_CODE_INTEGRATION.md` of 2026-09-16 (in the old
MolFeatures working copy), keeping its structure and updating the state, the counts and the
gates.

## Goal

The **package** holds every method the paper uses. The **paper repository** keeps its data
paths, its asserts and its one-off diagnostics, and calls the package for the maths.

Two rules hold throughout:

1. **The paper stays reproducible.** Tag `paper-v3` is the code behind the submission. Every
   ported method lands with a test that reproduces a number the deposit asserts. Nothing in
   `zenodo_archive/` changes without Eden's explicit approval.
2. **One implementation per method.** A port is not finished until the copies it replaces call
   it or are retired.

## Where things stand

| done since 2026-09-16 | where |
|---|---|
| Closed-form ridge LOO, exhaustive k-search, nested LOO, constraint null, gated on the deposit (0.813; 30/30 folds; nested 0.832; 384 of 1653) | `M3_modeler/ridge_search.py`, `tests/test_ridge_search.py` |
| The two `sterimol_utils.py` copies reconciled; all Sterimol fixes recorded | `docs/STERIMOL_FIXES.md` |
| One bond rule shared by the platform, the standalone copies and the explorer | `docs/FEATURE_FIXES.md` #1–2 |
| The explorer is in git, gated (`gate.js`, `gate_pos.js`, `smoke.js`) | `M2_data_extractor/theta_explorer/` |
| `_tmp_`, `_old` and `_backup` files and the case-study notebooks moved out of the package | `../DescriPyTor_workspace/` |

What the paper still runs outside the package (counted 2026-09-27):

| | files |
|---|---|
| paper scripts: `notes/analysis` 105, `zenodo_archive/scripts` 39, `figure_studio` 14 | 158 |
| case-study notebooks and scripts (workspace) | 33 |
| filenames in both `notes/analysis` and the deposit | 30, of which **23 have drifted apart** |
| scripts that re-implement closed-form LOO | 33 |
| scripts with their own exhaustive combination loop | 33 |
| scripts with a permutation null | 31 |
| scripts with their own ridge solve | 16 |
| scripts with their own Sterimol or bond code | 8 |

## Target: what goes where

| capability | module | functions | replaces | gate (the deposit's own number) |
|---|---|---|---|---|
| Search and nested LOO | `M3_modeler/ridge_search.py` | `loo_predictions`, `adj_q2`, `nested_loo`, `constraint_null` (**done**); add a `stability(selected)` summary | `cs3_search.py`, `cs3_constraint_null.py` and ~30 LOO copies | done: 0.813, 30/30, 0.832, 384/1653 |
| Controls | `M3_modeler/validation.py` | `selection_null` (permute y only; reuse each model's hat matrix, so 1000 permutations cost one einsum), `family_null`, `paired_sign_flip`, `group_holdout(draws=...)`, `baselines` | `selection_null_cs1_dipole.py`, `family_null*.py`, `holdout.py`, the sign-flip blocks | CS3 vs published: closer on 14/30, p 0.62, residual r 0.64; BOX→PyOx transfer against the deposited draws; CS1 random-vector null (`models/null_mazet_cs1_dipole.json`) |
| Frame comparison | `M2_data_extractor/frames.py` | `candidate_axes(mol)` (every backbone axis a structure defines), `features_per_frame(mols, axes, feature)`, `rank_frames(table, y, criterion)` | `cs1_backbone_frames.py`, `cs1_frame_control.py`, `cs2_frame_ablation.py`, `cs2_heldout_frames.py`, `cs3_frames.py` | CS1 13 axes (`cs1_backbone_frames_axes.csv`); CS2 23 frames (`cs2_frames.csv`); CS3 6 anchors (`cp_frame_ranking.csv`) |
| Complex building and xTB | `M1_pre_calculations/complex_builder.py`, `xtb_runner.py` | `complex_from_ligand(xyz, metal, ancillary)`, `xtb_opt`, `xtb_sp`, parsers | `build_general.py`, `rebuild_single_structure.py`, three xTB drivers | the rebuilt CS3 geometries reproduce the deposited matrix (needs xTB: run on the cluster) |
| Explorer from Python | `M2_data_extractor/theta_explorer/render.py` | `session(...)`, `shot(session, png)`, `values(spec)` (the `f1_values.js` path) | `figure_studio/figure_frames.py` helpers, `f1_values.js`, `panels_png.py` plumbing | the 27 Figure 1 numbers; goldens |
| Case-study pipelines | `Getting_started_with_examples/case_studies/` | one notebook per case study, calling the package on the deposited data and asserting the paper's numbers | the deposit notebooks' maths and the workspace notebooks | the table in `zenodo_archive/README.md`, row for row |
| Reproduction switch | `descripytor/compat.py` | `paper_v3()`: a context manager that sets `threshold_distance=1.82`, `STERIMOL_FRAME="lab"` and the pre-September keys | the per-script legacy flags | every gate above passes inside `paper_v3()` on `main` |

## Decide first (blocks phase 2 onward)

0. **Standardisation in the transfer test.** The draft's BOX→PyOx number (0.330) fits ridge on columns
   scaled with `ddof=0`; every other CS3 fit uses `ddof=1`. Pick one before `group_holdout` is ported.

1. **Deposit scripts: thin wrappers or frozen?** Either regenerate the deposit scripts from the
   package (one source, but the deposit changes), or freeze them as they are and point only
   `notes/analysis` at the package. Recommended: freeze the deposit, since its scripts
   already assert the paper's numbers standalone.
2. **The four definitions in `FEATURE_FIXES.md` "Found, not changed"**: the stretch score,
   metal neighbours in CPK typing, signed dihedrals, and grouped dipole input. The pipelines
   of phase 4 should be written after these are settled, not before.
3. **Which CS1/CS2 matrices the examples use**: deposited, or `*_sterimol_regenerated.csv`.
   The paper uses the deposited ones.

## Phases

Each phase ends with its gates passing on `main` and the replaced copies either calling the
package or listed as retired.

| # | phase | output | size |
|---|---|---|---|
| 1 | **Controls** — done 2026-09-27: `validation.selection_null`, `paired_sign_flip` (group_holdout waits for the ddof decision) | `validation.py` + tests; `ridge_search.stability` | about a day; every gate is a deposited number |
| 2 | **Reproduction switch** — done 2026-09-27: `descripytor.compat.paper_v3()` | `compat.paper_v3()`; the whole gate set runs on `main` inside it | half a day; removes the need to check out the tag |
| 3 | **Frames** | `frames.py`; CS1/CS2/CS3 frame tables reproduced | 1–2 days; three different frame conventions to unify |
| 4 | **Case-study notebooks** | CS1, CS2, CS3 notebooks on the package, asserting the README table | 2 days, after decision 2 |
| 5 | **Explorer from Python** | `render.py`; `figure_studio` becomes data plus calls | 1 day |
| 6 | **Builder and xTB** | `complex_builder.py`, `xtb_runner.py`; the CS3 rebuild runs from the package | 1–2 days plus cluster time |
| 7 | **One Sterimol kernel** | the platform, standalone, cube and `metal_complex` kernels behind one function with a frame option; optionally exact B1 (removes the ±0.01 Å grid limit) | 2 days; every Sterimol gate must stay at 0 |
| 8 | **Thin paper repo** | `notes/analysis` calls the package; the 23 drifted duplicates collapse | half a day per 20 scripts |

Phases 1–2 are independent of every decision above and can start now.

## Not moving

- One-off diagnostics documented in `notes/`: `cs12_mismatch_*`, `cs12_last_atom`,
  `cs12_angle_provenance`, `fix_009_chelate`, `h_arm_azimuth`, `frame_turn_1deg`,
  `oa_009_variants`. They stay in the paper repository as the record.
- Manuscript plumbing: `prep_acs_figures.py`, `build_manifest.py`, `wordcount.py`,
  `make_overleaf_package.py`.
- Scripts tied to one dataset that are not part of a method (`tools_*`,
  `extract_doyle_features.py`), unless rewritten against a generic loader.
- `zenodo_archive/scripts/common/sterimol.py` (a fourth, partial Sterimol: no B1, Bondi radii)
  and `cs3/tools_extract_metal_descriptors.py` (reads a sandbox path, cannot run). Flag both
  in the deposit README when the deposit is next revised.
