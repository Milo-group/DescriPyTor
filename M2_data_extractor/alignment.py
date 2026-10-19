"""Atom renumbering with Conformaty (github.com/Milo-group/Conformaty): every molecule of a folder numbered like a
reference, so a descriptor picked on the reference reads the same atom everywhere.

The result per molecule is the 1-based {old: new} dict ``Molecule.renumber_atoms`` takes. Conformaty matches the
largest common substructure, fits the molecule onto the reference and assigns the rest by position, so equivalent
atoms (the two H of a CH2) are told apart by where they sit. The reference must be the smaller molecule (it is the
common part); a molecule larger than it keeps its extra atoms after the reference's.

    from M2_data_extractor.alignment import map_folder, write_renumbered
    results = map_folder(mols, "basic")            # mols: a Molecules
    write_renumbered(mols, results, "renumbered")
"""
from __future__ import annotations

import os
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

INSTALL = "pip install git+https://github.com/Milo-group/Conformaty"   # shortcut: a git URL until Conformaty is on PyPI, then the [align] extra


def _conformaty():
    try:
        import conformaty
    except ImportError as exc:
        raise ImportError(f"Atom renumbering needs Conformaty: {INSTALL}") from exc
    return conformaty


def _worker_threads(n):
    import torch

    torch.set_num_threads(n)          # the cores shared out: each process's torch would otherwise take them all


# Conformaty's losses and Kabsch. rmsd is the default here: on the example set it converged on 9 of 9, as msd and Kabsch;
# mixed (Conformaty's own default) and mad stalled about 6 A away on 5 of 9 from its principal-axis seed, sig on 9.
FITS = ("rmsd", "msd", "mixed", "mad", "sig", "kabsch")


def kabsch(moving, fixed, rows):
    """moving, rotated and moved onto fixed by the closed-form least-squares fit of `rows` (the smallest RMSD)."""
    import numpy as np

    a, b = moving[rows], fixed[rows]
    ca, cb = a.mean(0), b.mean(0)
    u, _s, vt = np.linalg.svd((a - ca).T @ (b - cb))
    return (moving - ca) @ (u @ np.diag([1.0, 1.0, np.sign(np.linalg.det(u @ vt))]) @ vt) + cb


def fit_onto(ref, symbols, coords, fit="rmsd", sigma=30, rows=None):
    """The molecule (symbols, coords in the reference's numbering) fitted onto the reference (symbols, coords) over
    `rows` (0-based; default the heavy atoms holding the reference's element, as Conformaty's own align fits: with the
    hydrogens in, its principal-axis seed can start upside down, and mixed, mad and sig then stall about 6 A away).
    fit: one of FITS, Conformaty's rigid-body optimiser with that loss, or kabsch. Returns the fitted (n, 3) array, or
    the coordinates unmoved if fewer than three rows match."""
    import numpy as np

    rs, rx = ref
    rx, coords = np.asarray(rx, float), np.asarray(coords, float)
    if rows is None:
        rows = [i for i in range(len(rs)) if symbols[i] == rs[i]]
        heavy = [i for i in rows if rs[i] != "H"]
        rows = heavy if len(heavy) >= 3 else rows
    if len(rows) < 3:
        return coords
    if fit == "kabsch":
        return kabsch(coords, rx, rows)
    cf = _conformaty()
    moved = cf.fit_frames(cf.frame_from_atoms(rs, rx), cf.frame_from_atoms(symbols, coords), {i: i for i in rows},
                          loss_type=fit, sigma=sigma)
    return moved.iloc[2:][["X", "Y", "Z"]].to_numpy(float)


def compare_fits(ref, symbols, coords, sigma=30):
    """Every fit (FITS) of one molecule, already in the reference's numbering, onto the reference: the RMSD over the
    positions with the reference's element and over its heavy atoms, how many of them end within 0.25 A, the seconds,
    and the fitted xyz. Kabsch has the smallest RMSD by construction; mixed favours atoms in near-exact overlap."""
    import numpy as np

    rs, rx = ref[0], np.asarray(ref[1], float)
    same = [i for i in range(len(rs)) if symbols[i] == rs[i]]
    heavy = [i for i in same if rs[i] != "H"]
    out = []
    for fit in FITS:
        t0 = time.perf_counter()
        xyz = fit_onto(ref, symbols, coords, fit, sigma)
        d = np.linalg.norm(xyz[same] - rx[same], axis=1)
        out.append({"fit": fit, "rmsd": float(np.sqrt((d ** 2).mean())) if same else None,
                    "rmsd_heavy": float(np.sqrt((np.linalg.norm(xyz[heavy] - rx[heavy], axis=1) ** 2).mean())) if heavy else None,
                    "within": int((d < 0.25).sum()), "n": len(same), "seconds": time.perf_counter() - t0,
                    "xyz": f"{len(symbols)}\nfitted ({fit})\n" + "\n".join(
                        f"{s} {x:.6f} {y:.6f} {z:.6f}" for s, (x, y, z) in zip(symbols, xyz))})
    return out


def map_atoms(ref, query, anchors=None, stages=None, match_charge=False, fit="rmsd", sigma=30, fit_atoms=None):
    """Number one molecule like the reference. ref, query: (symbols, coords) of each.

    anchors: 1-based (reference atom, query atom) pairs that must correspond. stages: Conformaty's middle stages in
    order (complete_rings, extend_mcs, coords_proximity_any; default all three). fit, sigma: the rigid-body fit of
    the overlay (FITS); the numbering keeps Conformaty's own internal fit. fit_atoms: 1-based reference atoms the
    overlay is fitted on (default its heavy atoms with the reference's element, as Conformaty's align).
    Returns plain values (they cross a process boundary): old_to_new {old: new} (1-based), stage (the Conformaty stage
    that matched each atom, in the new order), xyz (the molecule fitted onto the reference, new order), differs (reference
    positions holding another element), rmsd over the positions with the same element (rmsd_kabsch: the same for
    the Kabsch fit, the least-squares optimum, to tell a stalled fit), seconds.
    """
    import numpy as np

    cf = _conformaty()
    (rs, rx), (qs, qx) = ref, query
    if len(qs) < len(rs):
        raise ValueError(f"it has fewer atoms ({len(qs)}) than the reference ({len(rs)}): choose the smaller molecule "
                         "as the reference")
    t0 = time.perf_counter()
    state = cf.map_frames(cf.frame_from_atoms(rs, rx), cf.frame_from_atoms(qs, qx),
                          anchors=[(r - 1, q - 1) for r, q in anchors] if anchors else None,
                          stage_order=stages, match_formal_charge=match_charge)
    order = cf.query_order(state)
    found = {q: entry.name for entry in state.log for q in entry.added}
    syms = [str(qs[q]) for q in order]
    same = [i for i in range(len(rs)) if syms[i] == rs[i]]
    rows = [a - 1 for a in fit_atoms if a - 1 in same] if fit_atoms else None
    rx = np.asarray(rx, float)
    moving = np.asarray(qx, float)[order]

    def dev(xyz):
        return float(np.sqrt(((xyz[same] - rx[same]) ** 2).sum(1).mean())) if same else None

    fitted = fit_onto((rs, rx), syms, moving, fit, sigma, rows)
    rmsd = dev(fitted)
    best = rmsd if fit == "kabsch" else dev(fit_onto((rs, rx), syms, moving, "kabsch", rows=rows))
    return {
        "old_to_new": {int(q) + 1: k + 1 for k, q in enumerate(order)},
        "stage": [found.get(q, "") for q in order],
        "xyz": f"{len(syms)}\nfitted onto the reference\n" + "\n".join(
            f"{s} {x:.6f} {y:.6f} {z:.6f}" for s, (x, y, z) in zip(syms, fitted)),
        "differs": [i + 1 for i in range(len(rs)) if syms[i] != rs[i]],
        "rmsd": rmsd,
        "rmsd_kabsch": best,
        "fit": fit,
        "seconds": time.perf_counter() - t0,
    }


def _atoms(mol):
    x = mol.xyz_df
    return x["atom"].astype(str).tolist(), x[["x", "y", "z"]].to_numpy(float).tolist()


def map_folder(mols, reference, anchors=None, stages=None, match_charge=False, workers=None, progress=None,
               cancelled=lambda: False, fit="rmsd", sigma=30, fit_atoms=None, isolate=False):
    """Number every molecule of a Molecules like `reference` (one of its names). anchors: {name: [(ref atom, atom)]}.

    workers > 1 runs that many molecules at once in separate processes, the cores shared between them; 1 (the default)
    runs them here, one after the other, unless isolate (a server: RDKit's substructure search holds the interpreter
    lock, which would stall its other threads). progress(event) gets one dict per molecule as it finishes; cancelled()
    is asked between molecules. Returns {name: result of map_atoms, or {"error": text}};
    the reference itself maps onto itself unchanged.
    """
    _conformaty()
    by = {m.molecule_name: m for m in mols.molecules}
    if reference not in by:
        raise ValueError(f"no molecule named {reference!r} in this folder")
    ref = _atoms(by[reference])
    n = len(ref[0])
    results = {reference: {"old_to_new": {i: i for i in range(1, n + 1)}, "stage": ["reference"] * n, "differs": [],
                           "rmsd": 0.0, "seconds": 0.0, "xyz": None, "reference": True}}
    todo = [name for name in by if name != reference]
    workers = max(1, min(len(todo), workers or 1))
    if progress:
        progress({"event": "start", "n": len(todo), "names": todo, "reference": reference, "workers": workers})
    if fit not in FITS:
        raise ValueError(f"fit must be one of {', '.join(FITS)}")
    args = {name: (ref, _atoms(by[name]), (anchors or {}).get(name), stages, match_charge, fit, sigma, fit_atoms)
            for name in todo}

    def done(i, name, run):
        try:
            results[name] = run()
        except Exception as exc:              # one molecule that cannot be mapped does not stop the folder
            results[name] = {"error": f"{type(exc).__name__}: {exc}"}
        if progress:
            progress({"event": "molecule", "i": i, "n": len(todo), "name": name, **_summary(results[name]),
                      "result": results[name]})

    if workers == 1 and not isolate:
        for i, name in enumerate(todo, 1):
            if cancelled():
                break
            done(i, name, lambda: map_atoms(*args[name]))
        return results
    import threading

    with ProcessPoolExecutor(workers, initializer=_worker_threads,
                             initargs=(max(1, (os.cpu_count() or 2) // workers),)) as pool:
        jobs = {pool.submit(map_atoms, *args[name]): name for name in todo}

        def watch():                       # one molecule can take minutes: Stop ends its process instead of waiting
            while not all(j.done() for j in jobs):
                if cancelled():
                    for other in jobs:
                        other.cancel()
                    for proc in list(getattr(pool, "_processes", {}).values()):
                        proc.terminate()
                    return
                time.sleep(0.3)

        threading.Thread(target=watch, daemon=True).start()
        for i, job in enumerate(as_completed(jobs), 1):
            if cancelled():
                break
            done(i, jobs[job], job.result)
    return results


def _summary(result):
    if "error" in result:
        return {"error": result["error"]}
    return {"moved": sum(o != n for o, n in result["old_to_new"].items()), "differs": result["differs"],
            "rmsd": result["rmsd"], "rmsd_kabsch": result.get("rmsd_kabsch"), "seconds": result["seconds"]}


def write_renumbered(mols, results, out_dir, reference=""):
    """Write each mapped molecule, renumbered, to out_dir/<name>.feather (the names unchanged, so outcomes.csv still
    joins; it is copied along), and renumbering.json: the {old: new} dict of every molecule. Molecules without one
    (failed) are not written. Returns (written names, skipped names)."""
    import copy
    import json
    import shutil

    from M2_data_extractor.data_extractor import _molecules_folder_files

    os.makedirs(out_dir, exist_ok=True)
    folder, files = _molecules_folder_files(mols.molecules_path)
    xyz_in = {os.path.splitext(f)[0] for f, _s, _z in files if f.lower().endswith(".xyz")}
    written, skipped = [], []
    for mol in mols.molecules:
        result = results.get(mol.molecule_name) or {}
        if "old_to_new" not in result:
            skipped.append(mol.molecule_name)
            continue
        m = copy.deepcopy(mol)
        if mol.molecule_name in xyz_in:          # coordinates stay coordinates; .log, .out and .json become .feather
            m.renumber_atoms(result["old_to_new"], verbose=False, save_feather=False)
            x = m.xyz_df
            with open(os.path.join(out_dir, f"{mol.molecule_name}.xyz"), "w", encoding="utf-8") as fh:
                fh.write(f"{len(x)}\nrenumbered like {reference}\n" + "".join(
                    f"{a} {float(px):.6f} {float(py):.6f} {float(pz):.6f}\n" for a, px, py, pz in x[["atom", "x", "y", "z"]].itertuples(index=False)))
        else:
            m.renumber_atoms(result["old_to_new"], verbose=False, save_feather=True,
                             out_file=os.path.join(out_dir, f"{mol.molecule_name}.feather"))
        written.append(mol.molecule_name)
    outcomes = os.path.join(mols.molecules_path, "outcomes.csv")
    if os.path.isfile(outcomes):
        shutil.copy(outcomes, out_dir)
    keep = ("old_to_new", "stage", "differs", "rmsd", "rmsd_kabsch", "fit", "seconds", "error")
    with open(os.path.join(out_dir, "renumbering.json"), "w", encoding="utf-8") as fh:
        json.dump({"reference": reference, "folder": mols.molecules_path,
                   "molecules": {k: {f: v[f] for f in keep if f in v} for k, v in results.items()}}, fh, indent=1)
    return written, skipped
