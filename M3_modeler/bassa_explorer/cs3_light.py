"""A BOX-type ligand on the case-study-3 complexes by the light route: SMILES -> RDKit conformers with the donors held syn -> build_general (the
metal and its ancillaries on the free ligand) -> one GFN2-xTB optimisation at default settings -> a single point -> the
columns of the cs3 cases. A rough first look, not a reproduction: from the SMILES the measured ligands move by up to
0.37 (cp), 0.99 (cc), 0.09 (oa) and 0.13 (da_f) kcal/mol from their case predictions (the paper's
notes/RESULT_cs3_light_route_v2.md); from the published free ligands the same code reproduces them to 3e-4.

    python cs3_light.py    reads {"smiles"} as JSON on stdin; writes {"rows": {case: {column: value}}, "info", "errors"}
                           (the three builds, CuCl for cp and oa, Ni(H)2 for cc, CuCl2 for da_f, run at once) or {"error"}

Runs under the reader Python (RDKit and the package's readers). xtb: $XTB, else on PATH, else a conda env's xtb.exe.
build_general.py is a verbatim copy of zenodo_archive/scripts/cs3/build_general.py (sha256 bb56b31e...a5493b8c5a2).
"""
import glob, json, os, shutil, subprocess, sys, tempfile, time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np

H = Path(__file__).resolve().parent
sys.path[:0] = [str(H), str(H.parents[1])]
import build_general as BG
from M2_data_extractor import metal_complex as mc, placement_search as ps
from M2_data_extractor.ligand_topology import LigandTopology
from M2_data_extractor.xtb_singlepoint import XtbSinglePoint, read_charges, read_wbo

# build: (metal, ancillary, how many, uhf); the case's reaction decides the build (build_wide.BUILD)
BUILDS = {"CuCl": ("Cu", "Cl", 1, 0), "Ni": ("Ni", "H", 2, 0), "CuCl2": ("Cu", "Cl", 2, 0)}
CASE_BUILD = {"cs3_cp": "CuCl", "cs3_oa": "CuCl", "cs3_cc": "Ni", "cs3_da_f": "CuCl2"}
PLACE = ("CsR", "DCs", "MD", "MCs", "MR", "DR", "RRp"); Q = ("B1", "B5", "L", "loc_B5", "theta", "phi")
XTB = os.environ.get("XTB") or shutil.which("xtb") or next(iter(sorted(glob.glob(str(Path.home() / "anaconda3/envs/*/Library/bin/xtb.exe")))), None)
mc.STERIMOL_KEYS = ("B1", "B5", "L", "theta", "angle")                       # as build_pools.py: the angle is written as phi


def conformers(smiles, n=10, seed=42):
    """free-ligand starts (symbols, coords), lazily: for each graph donor pair in the pipeline's order, ETKDGv3 + MMFF94 with
    the two donors held 2.5-3.0 A apart (pipeline _embed), lowest MMFF energy first (notes/PREREG_cs3_light_route_v2.md)"""
    from rdkit import Chem
    from M1_pre_calculations.pipeline.build import BuildError, _chelate_donor_pairs, _embed
    m = Chem.MolFromSmiles(smiles)
    if m is None: raise ValueError("RDKit cannot read this SMILES")
    for pair in _chelate_donor_pairs(Chem.AddHs(m), {}):
        try: mol, ids = _embed(smiles, n, seed, pair=pair, pair_range=(2.5, 3.0))
        except BuildError: continue                                        # this pair cannot chelate: the next one
        sym = [a.GetSymbol() for a in mol.GetAtoms()]
        for i in ids: yield sym, mol.GetConformer(i).GetPositions()


def xtb(work, args, timeout):
    if not XTB: raise RuntimeError("no xtb found: set XTB to the xtb executable (xtb 6.7.1, conda-forge)")
    env = dict(os.environ, OMP_NUM_THREADS="1", PATH=os.path.dirname(XTB) + os.pathsep + os.environ["PATH"])
    r = subprocess.run([XTB, *args], cwd=work, env=env, capture_output=True, encoding="utf-8", errors="replace", timeout=timeout)
    if r.returncode: raise RuntimeError("xtb %s failed (exit %d): %s" % (args[1], r.returncode, (r.stderr.strip().splitlines() or ["no message"])[-1]))
    return r.stdout


def single_point(work, out):
    """charges, Wiberg orders, dipole (a.u.) and HOMO/LUMO (eV) of an xtb single point, as the pipeline's xtb stage reads them"""
    lines = out.splitlines(); k = next(i for i, l in enumerate(lines) if "molecular dipole" in l)
    dip = next(l.split()[1:4] for l in lines[k:] if l.strip().startswith("full:"))
    homo = [float(l.split()[-2]) for l in lines if "(HOMO)" in l][-1]; lumo = [float(l.split()[-2]) for l in lines if "(LUMO)" in l][0]
    return XtbSinglePoint(charges=read_charges(work / "charges"), wbo=read_wbo(work / "wbo"), dipole=[float(v) for v in dip], homo=homo, lumo=lumo, dipole_unit="au")


def features(path, sp, smiles):
    """every column of the cs3 cases for one relaxed complex (build_pools.py + build_wide.py)"""
    m = mc.MetalComplex.from_xyz(str(path))
    row = {k.replace("_angle_", "_phi_"): v for k, v in m.geometric_features().items()}; row.update(m.electronic_features(sp))
    row.update(LigandTopology(smiles).size_normalized_features())
    F, _ = ps.compute(ps.Dataset.from_xyz([str(path)], ["x"]), ps.chelate_placements())
    for key in PLACE:
        for q in Q:
            a, b = F.get(f"{key}_a_{q}"), F.get(f"{key}_b_{q}")
            if a is None or b is None: continue
            row[f"{key}_{q}_sym"] = float((a.iloc[0] + b.iloc[0]) / 2); row[f"{key}_{q}_asym"] = float(abs(a.iloc[0] - b.iloc[0]))
    return {k: float(v) for k, v in row.items() if v is not None and np.isfinite(v)}


def featurize(build, smiles, work, start=None, timeout=240):
    """(row, info) for one ligand on one build. start: a free-ligand xyz to use instead of the RDKit conformers."""
    metal, anc, n_anc, uhf = BUILDS[build]; work = Path(work); work.mkdir(parents=True, exist_ok=True); t0 = time.time()
    starts = [BG.read_xyz(start)] if start else conformers(smiles); info = "no N, P or S donor pair could be embedded at a chelating distance"
    for k, (sym, X) in enumerate(starts):
        BG.write_xyz(work / "ligand.xyz", sym, X); res, info = BG.build(str(work / "ligand.xyz"), metal, anc, n_anc)
        if res is not None: break
    else: raise ValueError("no donor pair to put the metal on (%s): is this a bidentate N,N or N,P ligand?" % info)
    BG.write_xyz(work / "in.xyz", *res, "built")
    opt = xtb(work, ["in.xyz", "--gfn", "2", "--chrg", "0", "--uhf", str(uhf), "--opt"], timeout)
    if not (work / "xtbopt.xyz").exists(): raise RuntimeError("xtb wrote no optimised structure")
    sp = single_point(work, xtb(work, ["xtbopt.xyz", "--gfn", "2", "--chrg", "0", "--uhf", str(uhf), "--sp"], timeout))
    info = dict(build=build, conformer=None if start else k, donors=list(info["donors"]), bite_built=round(info["bite"], 1), converged="GEOMETRY OPTIMIZATION CONVERGED" in opt,
                energy=float(open(work / "xtbopt.xyz").read().splitlines()[1].split()[1]), seconds=round(time.time() - t0, 1))
    return features(work / "xtbopt.xyz", sp, smiles), info


def one(build, smi, w):
    try:
        row, info = featurize(build, smi, Path(w) / build)
        return (row, info) if info["converged"] else "the xtb optimisation did not converge"
    except subprocess.TimeoutExpired: return "xtb took too long and was stopped"
    except Exception as e: return "%s: %s" % (type(e).__name__, e)


if __name__ == "__main__":
    try:
        smi = str(json.loads(sys.stdin.read()).get("smiles", "")).strip()
        if not smi or len(smi) > 400: raise ValueError("give one SMILES (up to 400 characters)")
        from rdkit import Chem
        if Chem.MolFromSmiles(smi) is None: raise ValueError("RDKit cannot read this SMILES")
        if not XTB: raise RuntimeError("no xtb found: set XTB to the xtb executable (xtb 6.7.1, conda-forge)")
        with tempfile.TemporaryDirectory(prefix="cs3_light_") as w, ThreadPoolExecutor(3) as pool:          # the three builds at once
            got = dict(zip(BUILDS, pool.map(lambda b: one(b, smi, w), BUILDS)))
        print(json.dumps(dict(rows={c: got[b][0] for c, b in CASE_BUILD.items() if isinstance(got[b], tuple)}, info={b: g[1] for b, g in got.items() if isinstance(g, tuple)},
                              errors={c: got[b] for c, b in CASE_BUILD.items() if isinstance(got[b], str)})))
    except Exception as e: print(json.dumps(dict(error="%s: %s" % (type(e).__name__, e))))
