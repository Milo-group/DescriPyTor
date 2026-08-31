"""
descriptor_extractor.py
=======================

One script to run *many* descriptor-extraction options over a molecule set,
choosing which packages/engines to use and giving the right config to each.

Built on top of the DescriPytor `Molecules` / `Molecules_xyz` API used in
`Practical_Notebook_Features.ipynb` and `Full_Analysis_Notebook.ipynb`.

What it does
------------
1. Loads a *set* of molecules - either a directory of `.feather` files
   (DescriPytor native) or a directory of `.xyz` files.
2. Lets you turn on/off independent extraction ENGINES, each with its own
   config block:

     - descripytor_full   -> Molecules.get_molecules_features_set(answers_dict)
     - descripytor_steric -> Molecules.get_sterimol_dict(pairs)
     - xyz_sterimol       -> Molecules_xyz.get_sterimol_df(pairs, radii)
     - xyz_geometry       -> Molecules_xyz.get_angles_df / get_bond_lengths_df
     - xyz_buried_volume  -> Molecules_xyz.get_buried_volume_df(metal_index,...)
     - morfeus_sterimol   -> morfeus.Sterimol (independent reference engine)

3. Runs every enabled engine, merges all outputs column-wise into ONE
   features table, and writes it to CSV.

Atom selection
--------------
The atom indices each engine needs can be:
  * typed straight into the config (as in the notebooks), or
  * built visually by clicking atoms in the 3Dmol.js picker
    (`atom_picker.html`) and loaded here from the exported JSON.

Usage
-----
    # As a script, driven by a JSON config:
    python descriptor_extractor.py --config my_run.json

    # Print a starter config you can edit:
    python descriptor_extractor.py --print-template > my_run.json

    # From a notebook / REPL:
    from descriptor_extractor import run_from_config, ExtractionConfig
    df = run_from_config("my_run.json")

The config file format is documented in `config_example.json` and in
`build_template_config()` below.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
from pathlib import Path

import pandas as pd

# --------------------------------------------------------------------------- #
# 0.  Make the DescriPytor packages importable.
#     Defaults to the repository root (this file lives in
#     <root>/Getting_started_with_examples/descriptor_extraction_toolkit).
#     Can be overridden by the env var DESCRIPYTOR_ROOT or the config's
#     "root_dir" key.
# --------------------------------------------------------------------------- #
DEFAULT_ROOT_DIR = os.environ.get(
    "DESCRIPYTOR_ROOT",
    str(Path(__file__).resolve().parents[2]),
)


def _add_descripytor_to_path(root_dir: str) -> None:
    """Append the DescriPytor sub-packages to sys.path (idempotent)."""
    root_dir = str(root_dir)
    for sub in ("", "M3_modeler", "M2_data_extractor", "utils", "MolAlign"):
        p = os.path.join(root_dir, sub) if sub else root_dir
        if p not in sys.path and os.path.isdir(p):
            sys.path.append(p)


# --------------------------------------------------------------------------- #
# 1.  The canonical DescriPytor answers_dict keys.
#     These MUST match the strings the GUI / get_molecules_features_set expects,
#     so the picker output drops straight in.  (Taken verbatim from the
#     Practical_Notebook_Features.ipynb cells.)
# --------------------------------------------------------------------------- #
ANSWERS_KEYS = {
    "ring_vibration": "Ring Vibration atoms - by order -> Pick one  atom from a six member ring\n example: 13",
    "stretch_threshold": "Stretch Threshold",
    "stretching_vibration": "Stretching Vibration atoms- enter bonded atom pairs: \n example: 1,2 4,5",
    "bend_threshold": "Bend Threshold",
    "bending_vibration": "Bending Vibration atoms - enter atom pairs that have a common atom: \n example: 4,7",
    "center_atoms_dipole": "Center_Atoms Dipole",
    "dipole": "Dipole atoms - indices for coordination transformation: \n example: 4,5,6 - origin, y-axis, new xy plane",
    "charges": "charges values - Insert atoms to show charge: \n example: 1,2,3,4",
    "charge_diff": "charge_diff - Insert atoms to show charge difference: \n example: 1,2 3,4",
    "sterimol": "Sterimol atoms - Primary axis along: \n example: 7,8",
    "drop_atoms": "drop_atoms - Atoms to drop: \n example: 1,2,3",
    "bond_length": "Bond_length - Atom pairs to calculate difference: \n example: 1,2 4,5",
    "bond_angle": "Bond_angle - Atom triplets to calculate difference: \n example: 1,2,3 4,5,6",
}

# Short keys whose values are scalar thresholds (wrapped in a 1-element list).
_THRESHOLD_KEYS = {"stretch_threshold", "bend_threshold"}


def build_answers_dict(short: dict) -> dict:
    """
    Convert a compact, picker-friendly dict (short keys -> lists) into the
    verbose answers_dict that `Molecules.get_molecules_features_set` expects.

    Example
    -------
    >>> build_answers_dict({
    ...     "stretch_threshold": 1600,
    ...     "sterimol": [[1, 23], [1, 3]],
    ...     "bond_length": [[1, 2], [1, 23]],
    ... })
    {'Sterimol atoms ...': [[1, 23], [1, 3]], ...}
    """
    out = {}
    for short_key, full_key in ANSWERS_KEYS.items():
        val = short.get(short_key, [])
        if short_key in _THRESHOLD_KEYS:
            if val in (None, [], ""):
                val = []
            elif not isinstance(val, list):
                val = [val]
        elif val is None:
            val = []
        out[full_key] = val
    return out


# --------------------------------------------------------------------------- #
# 2.  Engine implementations.  Each takes (loaded_set, engine_config, ctx)
#     and returns a DataFrame indexed by molecule name (or None if skipped).
# --------------------------------------------------------------------------- #
class ExtractionContext:
    """Holds the lazily-loaded molecule sets so engines can share them."""

    def __init__(self, cfg: dict, progress=None):
        self.cfg = cfg
        self.root_dir = cfg.get("root_dir") or DEFAULT_ROOT_DIR
        _add_descripytor_to_path(self.root_dir)
        self.progress = progress
        self._feather_mols = None
        self._xyz_mols = None
        self._xyz_map = None
        self._smiles_map = None
        self._effective_xyz_dir = None

    # -- DescriPytor feather set --------------------------------------------
    @property
    def feather_mols(self):
        if self._feather_mols is None:
            from data_extractor import Molecules  # noqa: late import on purpose

            path = self.cfg["feather_dir"]
            print(f"[load] feather set <- {path}")
            self._feather_mols = Molecules(
                path,
                progress=self.progress,
                file_limit=self.cfg.get("file_limit"),
            )
            names = [m.molecule_name for m in self._feather_mols.molecules]
            print(f"[load] {len(names)} molecules: {names}")
        return self._feather_mols

    # -- xyz set (sterimol_standalone) --------------------------------------
    @property
    def xyz_mols(self):
        if self._xyz_mols is None:
            from sterimol_standalone import Molecules_xyz

            path = self.effective_xyz_dir
            print(f"[load] xyz set <- {path}")
            self._xyz_mols = Molecules_xyz(path)
            print(
                f"[load] {len(self._xyz_mols.molecules)} molecules | "
                f"failed: {self._xyz_mols.failed_molecules}"
            )
        return self._xyz_mols

    # -- name -> xyz path map (used by feature_packages engines) -------------
    @property
    def xyz_map(self):
        if self._xyz_map is None:
            import feature_packages as fp
            self._xyz_map = fp.xyz_map_from_dir(self.effective_xyz_dir)
            print(f"[load] xyz_map: {len(self._xyz_map)} files")
        return self._xyz_map

    @property
    def effective_xyz_dir(self) -> str:
        """
        Directory used by all xyz-based engines.

        Normally this is config['xyz_dir'].  When `derive_xyz_from_feathers` is
        true, it is generated once from Molecule.xyz_df so feather-only projects
        can still run xyz/RDKit/morfeus descriptors.
        """
        if self._effective_xyz_dir is None:
            if self.cfg.get("derive_xyz_from_feathers"):
                self._effective_xyz_dir = self.export_xyz_from_feathers(
                    self.cfg.get("derived_xyz_dir")
                )
            else:
                self._effective_xyz_dir = self.cfg["xyz_dir"]
        return self._effective_xyz_dir

    def export_xyz_from_feathers(self, out_dir: str | None = None) -> str:
        """Write one .xyz file per feather molecule and return the output dir."""
        if not out_dir:
            feather_dir = Path(self.cfg["feather_dir"]).resolve()
            out_dir = str(feather_dir.parent / f"{feather_dir.name}_xyz")
        out = Path(out_dir).resolve()
        out.mkdir(parents=True, exist_ok=True)
        print(f"[prep] exporting xyz from feathers -> {out}")

        count = 0
        for mol in self.feather_mols.molecules:
            df = getattr(mol, "xyz_df", None)
            if df is None or df.empty:
                print(f"[prep] skip {getattr(mol, 'molecule_name', '?')}: no xyz_df")
                continue
            name = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(getattr(mol, "molecule_name", f"mol_{count}")))
            path = out / f"{name}.xyz"
            rows = df[["atom", "x", "y", "z"]]
            with open(path, "w", encoding="utf-8") as f:
                f.write(f"{len(rows)}\n{name}\n")
                for _, r in rows.iterrows():
                    f.write(f"{str(r['atom']):2s} {float(r['x']):14.8f} {float(r['y']):14.8f} {float(r['z']):14.8f}\n")
            count += 1
        print(f"[prep] wrote {count} xyz files")
        return str(out)

    # -- name -> SMILES (perceived once, cached, shared by rdkit/mordred/...) -
    @property
    def smiles_map(self):
        if self._smiles_map is None:
            import feature_packages as fp
            print("[load] perceiving SMILES from xyz (once)...")
            self._smiles_map = fp.extract_smiles(self.xyz_map)
            ok = sum(v is not None for v in self._smiles_map.values())
            print(f"[load] SMILES ok for {ok}/{len(self._smiles_map)} molecules")
        return self._smiles_map


def _natural_sort(df: pd.DataFrame) -> pd.DataFrame:
    """Sort rows so conf_2 comes before conf_10 etc. (matches the notebooks)."""
    def key(idx):
        return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", str(idx))]

    try:
        return df.reindex(sorted(df.index, key=key))
    except Exception:
        return df


def engine_descripytor_full(ctx: ExtractionContext, ecfg: dict) -> pd.DataFrame:
    """Full DescriPytor descriptor set (IR, dipole, charges, sterimol...)."""
    mols = ctx.feather_mols
    answers = build_answers_dict(ecfg.get("atoms", {}))
    out_dir = ecfg.get("out_dir")
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
        os.chdir(out_dir)
    df = mols.get_molecules_features_set(
        entry_widgets=answers,
        answers_list=None,
        save_as=bool(ecfg.get("save_as", False)),
        csv_file_name=ecfg.get("csv_file_name", "descripytor_full_features"),
    )
    return _natural_sort(df)


def engine_descripytor_steric(ctx: ExtractionContext, ecfg: dict) -> pd.DataFrame:
    """Quick DescriPytor Sterimol-only table from the feather set."""
    mols = ctx.feather_mols
    pairs = ecfg["pairs"]
    df = mols.get_sterimol_dict(pairs)
    return _natural_sort(df)


def engine_xyz_sterimol(ctx: ExtractionContext, ecfg: dict) -> pd.DataFrame:
    """Sterimol straight from a directory of .xyz files (sterimol_standalone)."""
    mols = ctx.xyz_mols
    pairs = ecfg["pairs"]
    radii = ecfg.get("radii", "CPK")  # 'CPK' or 'bondi'
    df = mols.get_sterimol_df(pairs, radii=radii)
    return _natural_sort(df)


def engine_xyz_geometry(ctx: ExtractionContext, ecfg: dict) -> pd.DataFrame:
    """Angles, dihedrals and bond lengths from the .xyz set."""
    mols = ctx.xyz_mols
    parts = []
    angles = ecfg.get("angles")            # triplets and/or quadruplets
    if angles:
        parts.append(mols.get_angles_df(angles))
    bonds = ecfg.get("bond_lengths")       # 2-atom pairs
    if bonds:
        parts.append(mols.get_bond_lengths_df(bonds))
    if not parts:
        return pd.DataFrame()
    df = pd.concat(parts, axis=1)
    return _natural_sort(df)


def engine_xyz_buried_volume(ctx: ExtractionContext, ecfg: dict) -> pd.DataFrame:
    """Buried volume around an anchor/metal atom (sterimol_standalone)."""
    mols = ctx.xyz_mols
    df = mols.get_buried_volume_df(
        metal_index=ecfg["metal_index"],
        radius=ecfg.get("radius", 3.5),
    )
    return _natural_sort(df)


def engine_morfeus_sterimol(ctx: ExtractionContext, ecfg: dict) -> pd.DataFrame:
    """
    Independent Sterimol reference engine using the `morfeus` package, run over
    every .xyz file in xyz_dir.  Useful for cross-checking the DescriPytor /
    standalone numbers (see notebook cells that compare 'Extractor 1' vs '2').
    """
    from glob import glob
    from morfeus import read_xyz, Sterimol

    xyz_dir = ctx.effective_xyz_dir
    pairs = ecfg["pairs"]  # 1-indexed [origin, attached]
    rows = []
    for fp in glob(os.path.join(xyz_dir, "*.xyz")):
        conf = os.path.basename(fp).replace(".xyz", "")
        try:
            elements, coords = read_xyz(fp)
        except Exception as e:  # noqa
            print(f"[morfeus] skip {conf}: {e}")
            continue
        row = {}
        for i, j in pairs:
            try:
                s = Sterimol(elements, coords, i, j)
                tag = f"{i}-{j}"
                row[f"morfeus_L_{tag}"] = s.L_value
                row[f"morfeus_B1_{tag}"] = s.B_1_value
                row[f"morfeus_B5_{tag}"] = s.B_5_value
            except Exception:
                continue
        rows.append(pd.Series(row, name=conf))
    if not rows:
        return pd.DataFrame()
    df = pd.DataFrame(rows)
    return _natural_sort(df)


# --------------------------------------------------------------------------- #
# 2b.  External descriptor packages (ported from
#      feature_extraction_and_regression.ipynb -> feature_packages.py).
#      These operate on the xyz set (and SMILES derived from it), or .log files.
# --------------------------------------------------------------------------- #
def engine_rafbl(ctx, ecfg):
    import feature_packages as fp
    return _natural_sort(fp.extract_rafbl(ctx.xyz_map, mult=ecfg.get("mult", 1.2)))


def engine_rdkit(ctx, ecfg):
    import feature_packages as fp
    return _natural_sort(fp.extract_rdkit(ctx.smiles_map, embed=ecfg.get("embed", True)))


def engine_rdkit_fp(ctx, ecfg):
    import feature_packages as fp
    return _natural_sort(fp.extract_rdkit_fingerprints(
        ctx.smiles_map, n_bits=ecfg.get("n_bits", 256),
        include_usrcat=ecfg.get("include_usrcat", True)))


def engine_mordred(ctx, ecfg):
    import feature_packages as fp
    return _natural_sort(fp.extract_mordred(ctx.smiles_map, ignore_3d=ecfg.get("ignore_3d", True)))


def engine_deepchem(ctx, ecfg):
    import feature_packages as fp
    return _natural_sort(fp.extract_deepchem(ctx.smiles_map, n_bits=ecfg.get("n_bits", 256)))


def engine_qm(ctx, ecfg):
    import feature_packages as fp
    log_dir = ecfg.get("log_dir") or ctx.cfg.get("log_dir")
    if not log_dir:
        print("[qm] no log_dir set - skipping")
        return pd.DataFrame()
    return _natural_sort(fp.extract_qm(log_dir))


def engine_aqme_qdescp(ctx, ecfg):
    import feature_packages as fp
    work_dir = ecfg.get("work_dir") or os.path.join(
        os.path.dirname(os.path.abspath(ctx.effective_xyz_dir)), "aqme_qdescp")
    os.makedirs(work_dir, exist_ok=True)
    return _natural_sort(fp.extract_aqme_qdescp(
        ctx.xyz_map, work_dir,
        level=ecfg.get("level", "full"), atoms=ecfg.get("atoms"),
        charge=ecfg.get("charge", 0), mult=ecfg.get("mult", 1),
        solvent=ecfg.get("solvent"), xtb_opt=ecfg.get("xtb_opt", True),
        nprocs=ecfg.get("nprocs", 1), overwrite=ecfg.get("overwrite", False)))


def engine_morfeus_suite(ctx, ecfg):
    import feature_packages as fp
    return _natural_sort(fp.extract_morfeus_suite(
        ctx.xyz_map,
        descriptors=ecfg.get("descriptors"),
        sterimol_pairs=ecfg.get("sterimol_pairs"),
        metal_index=ecfg.get("metal_index"),
        cone_atoms=ecfg.get("cone_atoms"),
        bv_radius=ecfg.get("bv_radius", 3.5),
        pyramid_atoms=ecfg.get("pyramid_atoms")))


# Registry: engine name -> (function, which set it needs)
ENGINES = {
    # DescriPytor native
    "descripytor_full":    (engine_descripytor_full,    "feather"),
    "descripytor_steric":  (engine_descripytor_steric,  "feather"),
    "xyz_sterimol":        (engine_xyz_sterimol,        "xyz"),
    "xyz_geometry":        (engine_xyz_geometry,        "xyz"),
    "xyz_buried_volume":   (engine_xyz_buried_volume,   "xyz"),
    "morfeus_sterimol":    (engine_morfeus_sterimol,    "xyz"),
    # External packages (feature_extraction_and_regression.ipynb)
    "rafbl":               (engine_rafbl,               "xyz"),
    "rdkit":               (engine_rdkit,               "xyz"),
    "rdkit_fp":            (engine_rdkit_fp,            "xyz"),
    "mordred":             (engine_mordred,             "xyz"),
    "deepchem":            (engine_deepchem,            "xyz"),
    "qm":                  (engine_qm,                  "log"),
    "aqme_qdescp":         (engine_aqme_qdescp,         "xyz"),
    "morfeus_suite":       (engine_morfeus_suite,       "xyz"),
}


# --------------------------------------------------------------------------- #
# 3.  Runner
# --------------------------------------------------------------------------- #
def run_from_config(config, stop_check=None, progress=None) -> pd.DataFrame:
    """
    Run all enabled engines described in `config` (a dict or path to a JSON
    file) and return the merged features DataFrame.  Also writes the merged
    CSV if config['output_csv'] is set.

    stop_check : optional callable, checked before each engine starts.  If it
    returns truthy, the run stops early (whatever engines already finished are
    still merged and returned/written) - used by the webapp's Cancel button.
    Cancellation is only checked *between* engines, not mid-engine.

    progress : optional callable(event_dict).  Molecules.__init__ reports
    per-file load events; each engine start is reported as phase "engine".
    """
    if isinstance(config, (str, Path)):
        with open(config, "r", encoding="utf-8") as f:
            cfg = json.load(f)
    else:
        cfg = dict(config)

    # Allow an external atoms file (the picker's export) to feed the full engine.
    atoms_file = cfg.get("atoms_file")
    if atoms_file:
        with open(atoms_file, "r", encoding="utf-8") as f:
            picked = json.load(f)
        # picked may be either the compact short-key form or the verbose form.
        cfg.setdefault("engines", {}).setdefault("descripytor_full", {})
        cfg["engines"]["descripytor_full"].setdefault("atoms", {})
        cfg["engines"]["descripytor_full"]["atoms"].update(
            picked.get("atoms", picked)
        )

    ctx = ExtractionContext(cfg, progress=progress)
    engines_cfg = cfg.get("engines", {})

    results = []
    for name, ecfg in engines_cfg.items():
        if stop_check is not None and stop_check():
            print(f"\n[stop] cancelled before engine '{name}' - "
                 f"merging {len(results)} engine(s) that already finished.")
            break
        if not ecfg.get("enabled", True):
            print(f"[skip] {name} (disabled)")
            continue
        if name not in ENGINES:
            print(f"[warn] unknown engine '{name}' - skipping")
            continue
        fn, _need = ENGINES[name]
        print(f"\n=== engine: {name} ===")
        if progress:
            progress({"event": "progress", "phase": "engine", "name": name, "i": 0, "n": 0})
        try:
            df = fn(ctx, ecfg)
        except Exception as e:  # noqa - report and keep going
            print(f"[error] engine '{name}' failed: {e}")
            continue
        if df is None or len(df) == 0:
            print(f"[info] engine '{name}' produced no columns")
            continue
        # prefix columns so engines never collide
        prefix = ecfg.get("prefix", "")
        if prefix:
            df = df.add_prefix(prefix)
        print(f"[ok] {name}: {df.shape[0]} rows x {df.shape[1]} cols")
        results.append(df)

    if not results:
        print("\n[done] no engine produced output.")
        return pd.DataFrame()

    merged = pd.concat(results, axis=1)
    merged = _natural_sort(merged)

    out_csv = cfg.get("output_csv")
    if out_csv:
        os.makedirs(os.path.dirname(os.path.abspath(out_csv)) or ".", exist_ok=True)
        merged.to_csv(out_csv)
        print(f"\n[done] merged features -> {out_csv}  ({merged.shape})")
    else:
        print(f"\n[done] merged features ready  ({merged.shape}) - no output_csv set")
    return merged


# --------------------------------------------------------------------------- #
# 4.  Starter config
# --------------------------------------------------------------------------- #
def build_template_config() -> dict:
    """A fully-populated example config users can trim down."""
    return {
        "root_dir": DEFAULT_ROOT_DIR,
        "feather_dir": "path/to/your/feathers",
        "xyz_dir": "path/to/your/xyz",
        # If True, xyz-based engines use .xyz files exported from feather xyz_df.
        # Optional derived_xyz_dir controls where those files are written.
        "derive_xyz_from_feathers": False,
        "derived_xyz_dir": None,
        "output_csv": "merged_features.csv",
        # Optional: load atom selections exported by atom_picker.html
        "atoms_file": None,
        "engines": {
            "descripytor_full": {
                "enabled": True,
                "save_as": True,
                "csv_file_name": "descripytor_full_features",
                "prefix": "",
                "atoms": {
                    "ring_vibration": [],
                    "stretch_threshold": 1600,
                    "stretching_vibration": [],
                    "bend_threshold": 3000,
                    "bending_vibration": [],
                    "center_atoms_dipole": [],
                    "dipole": [[23, 1, 3], [2, 1, 23]],
                    "charges": [8, 1],
                    "charge_diff": [[8, 3], [8, 1]],
                    "sterimol": [[1, 23], [1, 3]],
                    "drop_atoms": [],
                    "bond_length": [[1, 2], [1, 23]],
                    "bond_angle": [[23, 1, 3], [23, 1, 2, 3]],
                },
            },
            "descripytor_steric": {
                "enabled": False,
                "pairs": [[2, 19], [2, 4], [6, 5]],
                "prefix": "steric_",
            },
            "xyz_sterimol": {
                "enabled": False,
                "pairs": [[6, 7], [6, 4]],
                "radii": "CPK",
                "prefix": "xyz_",
            },
            "xyz_geometry": {
                "enabled": False,
                "angles": [[2, 1, 5], [5, 1, 8, 6]],
                "bond_lengths": [[1, 5], [5, 7]],
                "prefix": "geom_",
            },
            "xyz_buried_volume": {
                "enabled": False,
                "metal_index": 18,
                "radius": 3.5,
                "prefix": "bv_",
            },
            "morfeus_sterimol": {
                "enabled": False,
                "pairs": [[6, 4]],
                "prefix": "",
            },
            # ---- external descriptor packages (need their pip installs) ----
            "rdkit": {
                "enabled": False,
                "embed": True,
                "prefix": "rdkit_",
            },
            "rdkit_fp": {
                "enabled": False,
                "n_bits": 256,
                "include_usrcat": True,
                "prefix": "",
            },
            "mordred": {
                "enabled": False,
                "ignore_3d": True,
                "prefix": "",
            },
            "deepchem": {
                "enabled": False,
                "n_bits": 256,
                "prefix": "",
            },
            "rafbl": {
                "enabled": False,
                "mult": 1.2,
                "prefix": "rafbl_",
            },
            "qm": {
                "enabled": False,
                "log_dir": "path/to/your/gaussian_logs",
                "prefix": "qm_",
            },
            "aqme_qdescp": {
                "enabled": False,
                "level": "full",
                "atoms": [],
                "charge": 0,
                "mult": 1,
                "solvent": None,
                "xtb_opt": True,
                "nprocs": 1,
                "prefix": "aqme_",
            },
            "morfeus_suite": {
                "enabled": False,
                "descriptors": ["sterimol", "buried_volume", "sasa", "dispersion"],
                "sterimol_pairs": [[6, 4]],
                "metal_index": None,
                "cone_atoms": None,
                "bv_radius": 3.5,
                "pyramid_atoms": None,
                "prefix": "",
            },
        },
    }


# --------------------------------------------------------------------------- #
# 5.  CLI
# --------------------------------------------------------------------------- #
def _main(argv=None):
    ap = argparse.ArgumentParser(description="DescriPytor multi-engine descriptor extractor")
    ap.add_argument("-c", "--config", help="Path to a JSON run config.")
    ap.add_argument("--print-template", action="store_true",
                    help="Print a starter JSON config and exit.")
    ap.add_argument("--list-engines", action="store_true",
                    help="List available engines and exit.")
    # ---- AQME CSEARCH preprocessing (SMILES csv -> xyz ensemble) ----
    ap.add_argument("--csearch", metavar="SMILES_CSV",
                    help="Run AQME CSEARCH on a SMILES CSV to generate an xyz set, then exit.")
    ap.add_argument("--csearch-out", metavar="DIR", help="Output dir for the generated .xyz files.")
    ap.add_argument("--csearch-program", default="rdkit",
                    choices=["rdkit", "summ", "fullmonte", "crest"], help="CSEARCH backend.")
    ap.add_argument("--csearch-sample", type=int, default=10, help="Conformers per molecule.")
    ap.add_argument("--name-col", default="name", help="Name column in the SMILES CSV.")
    ap.add_argument("--smiles-col", default="smiles", help="SMILES column in the SMILES CSV.")
    ap.add_argument("--keep", default="lowest", choices=["lowest", "all"],
                    help="Keep the lowest-energy conformer per molecule, or all of them.")
    ap.add_argument("--list-aqme-descriptors", action="store_true",
                    help="Print the AQME QDESCP descriptor names per level and exit.")
    args = ap.parse_args(argv)

    if args.print_template:
        print(json.dumps(build_template_config(), indent=2))
        return
    if args.list_engines:
        for name, (_, need) in ENGINES.items():
            print(f"{name:22s} (needs {need} set)")
        return
    if args.list_aqme_descriptors:
        _add_descripytor_to_path(DEFAULT_ROOT_DIR)
        import feature_packages as fp
        fp.list_aqme_descriptors()
        return
    if args.csearch:
        if not args.csearch_out:
            ap.error("--csearch also needs --csearch-out DIR")
        import feature_packages as fp
        out = fp.run_aqme_csearch(
            args.csearch, args.csearch_out,
            name_col=args.name_col, smiles_col=args.smiles_col,
            program=args.csearch_program, sample=args.csearch_sample, keep=args.keep)
        if out:
            print(f"\nNext: set \"xyz_dir\": \"{out}\" in your run config and enable the xyz engines.")
        return
    if not args.config:
        ap.error("provide --config FILE (or --print-template to start one)")
    run_from_config(args.config)


if __name__ == "__main__":
    _main()
