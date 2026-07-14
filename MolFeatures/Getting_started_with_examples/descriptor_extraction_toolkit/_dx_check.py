"""
Preflight checks for descriptor_extractor.py runs.

This script is intentionally lightweight: it validates paths, optional packages,
file-name alignment, Gaussian log contents, and atom-index ranges before a run.
It does not extract descriptors.
"""

from __future__ import annotations

import argparse
import contextlib
import glob
import importlib.util
import io
import json
import os
import re
import sys
from pathlib import Path

pd = None
PANDAS_IMPORT_ERROR = None


ENGINE_PACKAGES = {
    "rdkit": ["rdkit", "ase"],
    "rdkit_fp": ["rdkit", "ase"],
    "mordred": ["mordred", "rdkit", "ase"],
    "deepchem": ["deepchem"],
    "rafbl": ["moltop"],
    "qm": ["autoqchem"],
    "aqme_qdescp": ["aqme"],
    "morfeus_sterimol": ["morfeus"],
    "morfeus_suite": ["morfeus"],
}

XYZ_ENGINES = {
    "xyz_sterimol",
    "xyz_geometry",
    "xyz_buried_volume",
    "morfeus_sterimol",
    "morfeus_suite",
    "rdkit",
    "rdkit_fp",
    "mordred",
    "deepchem",
    "rafbl",
    "aqme_qdescp",
}
FEATHER_ENGINES = {"descripytor_full", "descripytor_steric"}
LOG_ENGINES = {"qm"}


class Report:
    def __init__(self):
        self.errors = []
        self.warnings = []
        self.ok = []

    def add(self, level, message):
        getattr(self, level).append(message)

    def print(self):
        for title, items, mark in (
            ("OK", self.ok, "[ok]"),
            ("Warnings", self.warnings, "[warn]"),
            ("Errors", self.errors, "[error]"),
        ):
            if not items:
                continue
            print(f"\n{title}")
            for item in items:
                print(f"  {mark} {item}")
        print(f"\nSummary: {len(self.errors)} errors, {len(self.warnings)} warnings, {len(self.ok)} ok")

    @property
    def exit_code(self):
        return 2 if self.errors else (1 if self.warnings else 0)


def have(module):
    try:
        return importlib.util.find_spec(module) is not None
    except Exception:
        return False


def get_pandas():
    global pd, PANDAS_IMPORT_ERROR
    if pd is not None or PANDAS_IMPORT_ERROR is not None:
        return pd
    try:
        with contextlib.redirect_stderr(io.StringIO()):
            import pandas as _pd
        pd = _pd
    except Exception as e:  # pragma: no cover - environment-specific
        PANDAS_IMPORT_ERROR = e
        pd = None
    return pd


def load_config(path):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def enabled_engines(cfg):
    return {
        name: ecfg
        for name, ecfg in cfg.get("engines", {}).items()
        if ecfg.get("enabled", True)
    }


def list_stems(folder, pattern):
    if not folder or not os.path.isdir(folder):
        return []
    return sorted(Path(p).stem for p in glob.glob(os.path.join(folder, pattern)))


def xyz_atom_counts(xyz_dir):
    counts = {}
    elements = {}
    for path in glob.glob(os.path.join(xyz_dir, "*.xyz")):
        stem = Path(path).stem
        try:
            with open(path, "r", encoding="utf-8", errors="ignore") as f:
                lines = [x.strip() for x in f.readlines() if x.strip()]
            n = int(lines[0])
            els = [line.split()[0] for line in lines[2 : 2 + n]]
            counts[stem] = len(els)
            elements[stem] = els
        except Exception:
            counts[stem] = 0
            elements[stem] = []
    return counts, elements


def feather_summaries(feather_dir):
    counts = {}
    elements = {}
    capabilities = {}
    pandas = get_pandas()
    if pandas is None or not have("pyarrow"):
        return counts, elements, capabilities
    feather_read_disabled = False
    for path in glob.glob(os.path.join(feather_dir, "*.feather")):
        if feather_read_disabled:
            break
        stem = Path(path).stem
        try:
            with contextlib.redirect_stderr(io.StringIO()):
                df = pandas.read_feather(path)
            cols = [str(c).strip() for c in df.columns]
            df.columns = cols
            if {"atom", "x", "y", "z"}.issubset(cols):
                xyz = df[["atom", "x", "y", "z"]].dropna()
                elements[stem] = xyz["atom"].astype(str).tolist()
                counts[stem] = len(elements[stem])
            else:
                xyz = df.iloc[:, 0:4].dropna()
                elements[stem] = xyz.iloc[:, 0].astype(str).tolist()
                counts[stem] = len(elements[stem])
            capabilities[stem] = {
                "dipole": all(c in cols for c in ("dip_x", "dip_y", "dip_z")),
                "freq": "Frequency" in cols,
                "ir": "IR" in cols,
                "nbo": "nbo_charge" in cols,
                "hirshfeld": "hirshfeld_charge" in cols,
                "cm5": "cm5_charge" in cols,
            }
        except Exception:
            feather_read_disabled = True
            counts.clear()
            elements.clear()
            capabilities.clear()
    return counts, elements, capabilities


def flatten_indices(value):
    out = []
    if value is None:
        return out
    if isinstance(value, int):
        return [value]
    if isinstance(value, (list, tuple)):
        for item in value:
            out.extend(flatten_indices(item))
    return out


def max_picked_index(cfg):
    mx = 0
    engines = cfg.get("engines", {})
    full = engines.get("descripytor_full", {}).get("atoms", {})
    for key, val in full.items():
        if key in {"stretch_threshold", "bend_threshold"}:
            continue
        ints = [x for x in flatten_indices(val) if isinstance(x, int)]
        if ints:
            mx = max(mx, max(ints))
    for name, ecfg in engines.items():
        for key in ("pairs", "angles", "bond_lengths", "sterimol_pairs", "cone_atoms", "pyramid_atoms"):
            ints = [x for x in flatten_indices(ecfg.get(key)) if isinstance(x, int)]
            if ints:
                mx = max(mx, max(ints))
        for key in ("metal_index",):
            if isinstance(ecfg.get(key), int):
                mx = max(mx, ecfg[key])
    return mx


def scan_logs(log_dir):
    summary = {}
    for path in glob.glob(os.path.join(log_dir, "*.log")) + glob.glob(os.path.join(log_dir, "*.out")):
        stem = Path(path).stem
        try:
            text = Path(path).read_text(encoding="utf-8", errors="ignore")
            summary[stem] = {
                "normal_termination": "Normal termination" in text,
                "freq": "Frequencies --" in text,
                "dipole": "Dipole moment" in text,
                "standard_orientation": "Standard orientation" in text,
                "nbo": bool(re.search(r"\bNBO\b|Natural Population|NPA", text, re.I)),
            }
        except Exception:
            summary[stem] = {}
    return summary


def compare_name_sets(report, label_a, a, label_b, b):
    if not a or not b:
        return
    set_a, set_b = set(a), set(b)
    missing_b = sorted(set_a - set_b)[:10]
    missing_a = sorted(set_b - set_a)[:10]
    if not missing_a and not missing_b:
        report.add("ok", f"{label_a} names match {label_b} names ({len(set_a & set_b)} shared)")
    else:
        report.add("warnings", f"name mismatch: {label_a} not in {label_b}: {missing_b or 'none'}; {label_b} not in {label_a}: {missing_a or 'none'}")


def run_checks(cfg):
    report = Report()
    engines = enabled_engines(cfg)
    report.add("ok", f"enabled engines: {', '.join(engines) if engines else 'none'}")

    root_dir = cfg.get("root_dir")
    feather_dir = cfg.get("feather_dir")
    xyz_dir = cfg.get("xyz_dir")
    log_dir = cfg.get("log_dir") or engines.get("qm", {}).get("log_dir")
    derive_xyz = bool(cfg.get("derive_xyz_from_feathers"))

    if root_dir and os.path.isdir(root_dir):
        report.add("ok", f"root_dir exists: {root_dir}")
    elif root_dir:
        report.add("warnings", f"root_dir not found: {root_dir}")

    needs_feather = bool(FEATHER_ENGINES & set(engines)) or derive_xyz
    needs_xyz = bool(XYZ_ENGINES & set(engines)) and not derive_xyz
    needs_log = bool(LOG_ENGINES & set(engines))

    feather_counts, feather_elements, feather_caps = {}, {}, {}
    if needs_feather:
        if not feather_dir or not os.path.isdir(feather_dir):
            report.add("errors", f"feather_dir is required but not found: {feather_dir}")
        else:
            stems = list_stems(feather_dir, "*.feather")
            report.add("ok", f"found {len(stems)} feather files")
            feather_counts, feather_elements, feather_caps = feather_summaries(feather_dir)
            if get_pandas() is None:
                report.add("warnings", "pandas is not importable; feather contents were not scanned")
            elif not have("pyarrow"):
                report.add("warnings", "pyarrow is not importable; feather contents were not scanned")
            elif stems and not feather_counts:
                report.add("warnings", "feather metadata scan failed; check pyarrow/numpy compatibility in this environment")

    xyz_counts, xyz_elements = {}, {}
    if needs_xyz:
        if not xyz_dir or not os.path.isdir(xyz_dir):
            report.add("errors", f"xyz_dir is required but not found: {xyz_dir}")
        else:
            stems = list_stems(xyz_dir, "*.xyz")
            report.add("ok", f"found {len(stems)} xyz files")
            xyz_counts, xyz_elements = xyz_atom_counts(xyz_dir)

    if needs_log:
        if not log_dir or not os.path.isdir(log_dir):
            report.add("errors", f"log_dir is required for qm but not found: {log_dir}")
        else:
            logs = scan_logs(log_dir)
            report.add("ok", f"found {len(logs)} Gaussian log/out files")
            bad = [k for k, v in logs.items() if not v.get("normal_termination")]
            if bad:
                report.add("warnings", f"logs without Normal termination: {bad[:10]}")
            if logs and not any(v.get("freq") for v in logs.values()):
                report.add("warnings", "no frequency blocks found in logs; vibration features may be empty")
            if logs and not any(v.get("dipole") for v in logs.values()):
                report.add("warnings", "no dipole blocks found in logs")

    if feather_counts and xyz_counts:
        compare_name_sets(report, "feather", feather_counts, "xyz", xyz_counts)
        shared = sorted(set(feather_elements) & set(xyz_elements))
        mismatched = [s for s in shared if feather_elements[s] != xyz_elements[s]]
        if mismatched:
            report.add("warnings", f"atom element/order mismatch between feather and xyz for: {mismatched[:10]}")
        elif shared:
            report.add("ok", f"atom element/order matches for {len(shared)} shared molecules")

    picked_max = max_picked_index(cfg)
    if picked_max:
        counts = feather_counts or xyz_counts
        too_small = [name for name, n in counts.items() if n and picked_max > n]
        if too_small:
            report.add("errors", f"picked atom index {picked_max} exceeds atom count for: {too_small[:10]}")
        elif counts:
            report.add("ok", f"picked atom indices fit scanned molecules (max index {picked_max})")

    full_atoms = engines.get("descripytor_full", {}).get("atoms", {})
    if full_atoms.get("dipole") and feather_caps:
        missing = [k for k, v in feather_caps.items() if not v.get("dipole")]
        if missing:
            report.add("warnings", f"dipole selected but dipole columns missing in feathers: {missing[:10]}")
    if (full_atoms.get("stretching_vibration") or full_atoms.get("bending_vibration") or full_atoms.get("ring_vibration")) and feather_caps:
        missing = [k for k, v in feather_caps.items() if not v.get("freq")]
        if missing:
            report.add("warnings", f"vibration selections present but Frequency column missing: {missing[:10]}")

    for engine, packages in ENGINE_PACKAGES.items():
        if engine not in engines:
            continue
        missing = [pkg for pkg in packages if not have(pkg)]
        if missing:
            report.add("warnings", f"{engine} enabled but missing import(s): {', '.join(missing)}")
        else:
            report.add("ok", f"{engine} optional imports available")

    return report


def print_install_report(cfg):
    engines = enabled_engines(cfg)
    needed = []
    for name in engines:
        needed.extend(ENGINE_PACKAGES.get(name, []))
    missing = sorted({pkg for pkg in needed if not have(pkg)})
    if not missing:
        print("All optional imports for enabled engines are available.")
        return
    pip_names = {
        "morfeus": "morfeus-ml",
    }
    print("Missing optional packages:")
    for pkg in missing:
        print(f"  - {pkg}")
    print("\nSuggested pip command:")
    print("  pip install " + " ".join(pip_names.get(pkg, pkg) for pkg in missing))


def main(argv=None):
    ap = argparse.ArgumentParser(description="Preflight descriptor extraction config.")
    ap.add_argument("--config", default="run_config.json", help="Path to run_config.json")
    ap.add_argument("--install-report", action="store_true", help="Show missing optional packages for enabled engines")
    args = ap.parse_args(argv)

    if not os.path.exists(args.config):
        print(f"[error] config not found: {args.config}")
        return 2
    cfg = load_config(args.config)

    if args.install_report:
        print_install_report(cfg)
        return 0

    report = run_checks(cfg)
    report.print()
    return report.exit_code


if __name__ == "__main__":
    raise SystemExit(main())
