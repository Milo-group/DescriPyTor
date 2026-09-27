# Sterimol: fixes, conventions and open issues

DescriPyTor has three Sterimol implementations. They share the Verloop definition (CPK radii,
B1 = narrowest width over a rotation scan, B5 = widest, L = length along the axis), but not
their code:

| path | used by | frame of the B1 scan |
|---|---|---|
| `extractor_utils/sterimol_utils.py` via `data_extractor.get_sterimol_df` | the platform, Case Studies 1 and 2 | three atoms (base, axis, a neighbour) |
| `metal_complex.sterimol` via `MetalComplex.geometric_features` | Case Study 3 (`fromM_*`, `sub_*`) | see fix 5 |
| `theta_explorer/sterimol.js` | the explorer and every paper figure | line-for-line port of the platform path |

`sterimol_standalone.py` and `cube_sterimol.py` carry their own copies of the platform scan.

The bug-free version is `main` from commit `0857485` on. Tag **`paper-v3`** is the state that
reproduces the paper's deposited tables.

## Fixed

| # | bug | where | effect | fixed in | pinned by |
|---|---|---|---|---|---|
| 1 | **Mirrored B1 normal.** The scan rotates the plane, t = R p, so the rotated axis in the original frame is Rᵀx = [cos, −sin]. The code stored R x = [cos, sin], its mirror image. | `sterimol_utils.b1s_for_loop_function`, `scan_b1_over_angles` | `B1_B5_angle` (θ) changed when the frame was turned about the axis (29.4 → 35.2° under a 17° turn). B1, B5, L and loc_B5 were unaffected. The deposited CS1/CS2 θ columns hold the mirrored values; neither case study's reported equation uses them. | `5713c51` (entered in `dd938be`, 2026-06-30) | `tests/test_sterimol_b1_normal.py` |
| 2 | **18° B1 scan** in the standalone copies (`scans=90//5`) | `sterimol_standalone.get_b1s_list`, `cube_sterimol.get_b1s_list` | B1 misses the narrowest plane when it falls between grid points | `5713c51`: every degree, as `sterimol_utils` already did | — |
| 3 | **18° B1 scan** in `metal_complex` (coarse 18, 36 … 90°, refined ±18°) | `metal_complex.sterimol` | missed the narrowest plane on 55/708 GOAT arms (metal frame up to 0.38 Å) | `d892a8a`: `STERIMOL_SCAN_STEP = 1` (set 18 to rebuild pre-September tables) | `test_default_scan_finds_the_plane_the_18_degree_grid_misses` |
| 4 | **φ reported for CS3.** The in-plane azimuth φ is undefined when the B5 atom lies on the axis, which happens on every lone-H arm (9 of 30 cp ligands; the tie-break alone moved cp 0.789 ↔ 0.812). | `MetalComplex.geometric_features` | convention-dependent columns | `d892a8a`: reports θ (`STERIMOL_KEYS`); add `"angle"` back to rebuild old tables | `test_theta_of_a_lone_hydrogen_is_zero` |
| 5 | **Lab-frame scan start.** e1 = axis × lab-z, so a rigid rotation of the input file moved θ by up to 25°, φ by up to 83° and B1 by up to 0.004 Å | `metal_complex.sterimol` | CS3 angle columns depended on file orientation | `0857485`: e1 is the fragment atom furthest off the axis. `STERIMOL_FRAME = "lab"` restores the old start, and with it the deposited CS3 tables rebuild exactly. | `tests/test_metal_complex_frame.py` (fails on `paper-v3`) |
| 6 | **`origin` argument ignored.** `preform_coordination_transformation(origin=…)` centred on the basis atoms instead of the atoms passed | `sterimol_utils` | any caller passing `origin` got the wrong centre | `c78d2ca` | — |
| 7 | Circle sampling in the plot and the scan disagreed | `sterimol_utils` (`STERIMOL_CIRCLE_POINTS`) | B1 plot slices misaligned | `5713c51` | — |
| 8 | steriplot endpoint unpacked a DataFrame into two names | `gui_server` | endpoint crashed | `ae36eb0` | — |
| 9 | **Flat 1.82 Å bond cutoff** dropped S–CF₃ (p-OTf, L 7.82 → 10.32 Å) and P–C bonds, and changed the CPK type of the atom left behind | `utils.help_functions.extract_connectivity` and its copies | Sterimol of any group containing a long single bond | after `paper-v3`: covalent radii × 1.15; `threshold_distance=1.82` restores the old rule (see FEATURE_FIXES.md #1) | `tests/test_connectivity.py` |

### Explorer port (`theta_explorer/sterimol.js`)

| bug | fix | gate |
|---|---|---|
| φ folded to ≤ 90° (unsigned B1 normal): 33° where the table says 147° | signed contact direction | `gate.js`, φ within 0.6° of Python |
| no fragment blocking: on a chelate, an N → C* axis walked round the ring and measured the chloride as B5 | `opt.block` (stop the walk) and `opt.drop` (remove atoms), as in `metal_complex` | `smoke.js` |
| sessions dropped `thr` (and before that `labShift`), so figures rebuilt from a saved session used 1.82 Å on phosphines | both are serialised and restored | `smoke.js` goldens |

## Added (no effect on existing columns)

- `loc_B1`, `B1_tangent_atoms`, `B1_B5_azimuth` (directed in-plane angle, 0–180°) and
  `B1_B5_acute` (folded to 0–90°) next to `B1_B5_angle`. sin(`B1_B5_angle`) = |cos(azimuth)| ·
  |p|/|v|, so θ is partly confounded with loc_B5; prefer `B1_B5_acute` unless the near/far gap has
  been checked for the set.
- `metal_complex.sterimol` also returns its frame (`origin`, `axis`, `normal`, `b5_atom`, `atoms`)
  for drawing.
- `BITE_WINDOW` warning when a donor has come off the metal (009_lig).

## Conventions (behaviour, not bugs)

- **Third frame atom follows Python set order** in `direction_atoms_for_sterimol`, not ascending
  order. On 20 of 118 CS1/CS2 axes this picks a different atom, and on 3 of those B1 changes. The
  JS port emulates CPython 3.10 set order to match.
- **Measure the whole molecule, draw the group.** Cutting the fragment out changes the frame's
  third atom and moves B1 by up to 0.09 Å (026_lig). Use `block`/`drop` to restrict the
  measurement, or the explorer's `fadeOff` to restrict only the drawing.
- **Substituent-frame arms whose attachment atom sits on the axis** (the quinoline arms of cp
  037–039, oa 009's phosphine) have a flat B1 minimum: B1 is stable but its plane, and so θ, can
  flip by 21–41° under a frame turn. Report them as a band.
- `metal_complex` and the platform path use different scan frames, so their B1 on the same axis
  can differ by a few thousandths of an Å (041_lig: 1.8987 vs 1.8923).

## Open

- **B1 resolution about ±0.01 Å.** B1 is a minimum that usually sits on a kink where two atoms tie, and the
  scan samples it every 1°, so the scan's starting direction moves B1 by up to 0.01 Å. An exact
  minimisation (rotating calipers on the projected discs) would remove it.
- **Deposited CS1/CS2 Sterimol tables** were computed with each molecule's last atom missing from
  the bond table. This is a defect in how the data was generated, not in current code, and its
  origin is not identified. `*_sterimol_regenerated.csv` holds the corrected values; they are not
  adopted in the paper.

## Reproducing the paper

```bash
git checkout paper-v3
```

Or on `main`: `metal_complex.STERIMOL_FRAME = "lab"` (and `STERIMOL_SCAN_STEP = 18`,
`STERIMOL_KEYS = ("B1", "B5", "L", "angle")` for the pre-September 18° tables).

Checked 2026-09-27: CS1/CS2 regenerated Sterimol blocks, every cell exact. CS3 single-structure
Sterimol/θ block, 16 columns to 7e-15 from the rebuilt geometries. Explorer gates: gate.js 0 of
102, gate_pos.js 0 of 14, smoke.js 0 of 1383.

Background: `notes/RESCAN_AND_THETA_FIX.md` and `notes/STERIMOL_AUDIT_SUMMARY.md` in the paper
repository.
