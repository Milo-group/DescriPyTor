"""
feature_packages.py
====================

Open-access / external descriptor *packages*, ported from
`feature_extraction_and_regression.ipynb` (and `xyz_feature_extraction.ipynb`)
so the toolkit can mix DescriPytor's own descriptors with:

    rafbl         - moltop topological graph descriptors (from .xyz)
    rdkit         - full RDKit 2D descriptor list (from SMILES)
    rdkit_fp      - Morgan/ECFP, FCFP, MACCS, AtomPair, TopologicalTorsion, USRCAT
    mordred       - Mordred 2D/3D descriptor families
    deepchem      - DeepChem CircularFingerprint + RDKitDescriptors
    qm            - autoqchem gaussian_log_extractor descriptors (from .log)
    aqme_qdescp   - AQME QDESCP xTB / MORFEUS descriptors (from .xyz, needs xtb)
    morfeus_suite - direct morfeus: Sterimol, BuriedVolume, ConeAngle, SASA,
                    Dispersion, Pyramidalization  (atom-anchored)

Every extractor:
  * returns a pandas DataFrame indexed by molecule name (df.index.name='name'),
  * is guarded by an availability check - if the package isn't installed it
    prints a one-line skip and returns an empty DataFrame, so a run never dies
    because one optional dependency is missing.

Install hints
-------------
    pip install rdkit moltop mordred deepchem morfeus-ml ase networkx
    # qm     -> pip install autoqchem
    # aqme   -> pip install aqme   (and a working xTB)
"""

from __future__ import annotations

import glob
import importlib.util
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def have(pkg: str) -> bool:
    """True if an importable module/spec exists for `pkg`."""
    try:
        return importlib.util.find_spec(pkg) is not None
    except (ImportError, ValueError):
        return False


def _clean(df: pd.DataFrame) -> pd.DataFrame:
    """Drop all-NaN and constant columns (matches the notebook's _clean)."""
    if df is None or len(df) == 0:
        return pd.DataFrame()
    df = df.dropna(axis=1, how="all")
    df = df.loc[:, df.nunique() > 1]
    return df


def xyz_map_from_dir(xyz_dir: str, name_fn=None) -> dict:
    """{molecule_name: path} for every .xyz in xyz_dir."""
    fn = name_fn or (lambda s: s)
    return {fn(Path(p).stem): p for p in sorted(glob.glob(os.path.join(xyz_dir, "*.xyz")))}


# --------------------------------------------------------------------------- #
# SMILES from xyz (ase + RDKit bond perception). Many engines need this.
# --------------------------------------------------------------------------- #
def extract_smiles(xyz_map: dict) -> dict:
    """{name: SMILES or None} - perceives bonds with several charge guesses."""
    if not (have("rdkit") and have("ase")):
        print("[smiles] need rdkit + ase; skipping SMILES perception")
        return {name: None for name in xyz_map}

    import copy
    from ase.io import read as ase_read
    from rdkit import Chem, RDLogger
    from rdkit.Chem import rdDetermineBonds
    RDLogger.DisableLog("rdApp.*")

    out = {}
    for name, path in xyz_map.items():
        mol_out, last_err = None, None
        try:
            atoms = ase_read(path)
            rw = Chem.RWMol()
            for a in atoms:
                rw.AddAtom(Chem.Atom(int(a.number)))
            conf = Chem.Conformer(rw.GetNumAtoms())
            for i, pos in enumerate(atoms.get_positions()):
                conf.SetAtomPosition(i, pos.tolist())
            rw.AddConformer(conf)
            base = rw.GetMol()
            for charge in (0, 1, -1, 2, -2):
                try:
                    mol = copy.deepcopy(base)
                    rdDetermineBonds.DetermineBonds(mol, charge=charge)
                    Chem.SanitizeMol(mol)
                    mol_out = Chem.MolToSmiles(mol)
                    break
                except Exception as e:
                    last_err = e
        except Exception as e:
            last_err = e
        if mol_out is None:
            print(f"[smiles] FAIL {name}: {type(last_err).__name__}")
        out[name] = mol_out
    return out


def _mol_from_smiles(smi, add_hs=False, embed_3d=False):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    if not smi:
        return None
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    if add_hs or embed_3d:
        mol = Chem.AddHs(mol)
    if embed_3d:
        try:
            params = AllChem.ETKDGv3()
            params.randomSeed = 42
            if AllChem.EmbedMolecule(mol, params) == 0:
                AllChem.UFFOptimizeMolecule(mol, maxIters=200)
            else:
                return None
        except Exception:
            return None
    return mol


# --------------------------------------------------------------------------- #
# rafbl - moltop topological graph descriptors (from .xyz)
# --------------------------------------------------------------------------- #
def extract_rafbl(xyz_map: dict, mult: float = 1.2,
                  max_heavy_t: int = 30, max_heavy_chi: int = 40) -> pd.DataFrame:
    if not have("moltop"):
        print("[rafbl] moltop not installed - skipping")
        return pd.DataFrame()
    from moltop import descriptors as md_top
    from moltop.adjacency import g_from_file
    from moltop import utils as mu

    rows = []
    for name, path in xyz_map.items():
        row = {"name": name}
        try:
            G = g_from_file(path, mult=mult)
            g_an = mu.no_hydrogen(G)
            na, nb = G.number_of_nodes(), G.number_of_edges()
            na_an, nb_an = g_an.number_of_nodes(), g_an.number_of_edges()
            row.update({"na": na, "nb": nb, "bpa": nb / na if na else 0,
                        "na_an": na_an, "nb_an": nb_an,
                        "bpa_an": nb_an / na_an if na_an else 0})
            for fn, kw, key in [("estrada", {}, "estrada"), ("wiener", {}, "wiener"),
                                ("global_eff", {}, "global_eff"),
                                ("zagreb", {"n": 1}, "zagreb1"), ("zagreb", {"n": 2}, "zagreb2"),
                                ("molecular_shannon_i", {}, "shannon")]:
                try:
                    row[key] = getattr(md_top, fn)(g_an, **kw)
                except Exception:
                    pass
            for fn in ("balaban_j", "hosoya_z", "redundancy"):
                try:
                    row[fn] = getattr(md_top, fn)(g_an)
                except Exception:
                    pass
            for valence in (False, True):
                sfx = "_v" if valence else ""
                try:
                    row[f"kier_alpha{sfx}"] = md_top.kier_alpha(g_an, valence=valence)
                except Exception:
                    pass
                try:
                    for ki, kv in enumerate(md_top.kier_mkappa(g_an, valence=valence), 1):
                        row[f"kier_mkappa{ki}{sfx}"] = kv
                except Exception:
                    pass
            for fn, base in [("kier_phi", "kier_phi"), ("kier_xi", "kier_xi")]:
                try:
                    for ki, kv in enumerate(getattr(md_top, fn)(g_an), 1):
                        row[f"{base}{ki}"] = kv
                except Exception:
                    pass
            if na_an <= max_heavy_chi:
                try:
                    lim = min(na_an - 1, 8)
                    for order in range(lim + 1):
                        for valence in (False, True):
                            sfx = "_v" if valence else ""
                            row[f"chi{order}{sfx}"] = md_top.hall_kier_mchi(G, order=order, valence=valence)
                except Exception:
                    pass
            if na_an <= max_heavy_t:
                try:
                    row["T"] = md_top.global_top_state(G)
                except Exception:
                    pass
        except Exception as e:
            print(f"[rafbl] FAIL {name}: {e}")
        rows.append(row)
    return _clean(pd.DataFrame(rows).set_index("name"))


# --------------------------------------------------------------------------- #
# rdkit - full 2D descriptor list (from SMILES, embedded)
# --------------------------------------------------------------------------- #
def extract_rdkit(smiles_map: dict, embed: bool = True) -> pd.DataFrame:
    if not have("rdkit"):
        print("[rdkit] rdkit not installed - skipping")
        return pd.DataFrame()
    if not smiles_map:
        print("[rdkit] smiles_map is empty - check xyz_dir or derive_xyz_from_feathers")
        return pd.DataFrame()
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors
    from rdkit.ML.Descriptors import MoleculeDescriptors

    names = [d[0] for d in Descriptors.descList]
    calc = MoleculeDescriptors.MolecularDescriptorCalculator(names)
    rows = []
    for name, smi in smiles_map.items():
        row = {"name": name}
        if smi:
            try:
                mol = Chem.MolFromSmiles(smi)
                if mol:
                    if embed:
                        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
                    row.update(dict(zip(names, calc.CalcDescriptors(mol))))
            except Exception as e:
                print(f"[rdkit] FAIL {name}: {e}")
        rows.append(row)
    return _clean(pd.DataFrame(rows).set_index("name"))


# --------------------------------------------------------------------------- #
# rdkit_fp - fingerprints + USRCAT
# --------------------------------------------------------------------------- #
def extract_rdkit_fingerprints(smiles_map: dict, n_bits: int = 256,
                               include_usrcat: bool = True) -> pd.DataFrame:
    if not have("rdkit"):
        print("[rdkit_fp] rdkit not installed - skipping")
        return pd.DataFrame()
    if not smiles_map:
        print("[rdkit_fp] smiles_map is empty - check xyz_dir or derive_xyz_from_feathers")
        return pd.DataFrame()
    from rdkit import DataStructs
    from rdkit.Chem import rdMolDescriptors, MACCSkeys

    def bv(fp):
        arr = np.zeros((fp.GetNumBits(),), dtype=int)
        DataStructs.ConvertToNumpyArray(fp, arr)
        return arr

    def sc(fp, nb):
        arr = np.zeros((nb,), dtype=float)
        for bit, val in fp.GetNonzeroElements().items():
            arr[int(bit) % nb] += float(val)
        return arr

    rows = []
    for name, smi in smiles_map.items():
        row = {"name": name}
        mol = _mol_from_smiles(smi)
        if mol is not None:
            try:
                row.update({f"ecfp4_{i}": int(v) for i, v in enumerate(
                    bv(rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2, nBits=n_bits)))})
            except Exception as e:
                print(f"[rdkit_fp] ECFP {name}: {e}")
            try:
                row.update({f"fcfp4_{i}": int(v) for i, v in enumerate(
                    bv(rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2, nBits=n_bits, useFeatures=True)))})
            except Exception as e:
                print(f"[rdkit_fp] FCFP {name}: {e}")
            try:
                row.update({f"maccs_{i}": int(v) for i, v in enumerate(bv(MACCSkeys.GenMACCSKeys(mol)))})
            except Exception as e:
                print(f"[rdkit_fp] MACCS {name}: {e}")
            try:
                row.update({f"atom_pair_{i}": v for i, v in enumerate(
                    sc(rdMolDescriptors.GetHashedAtomPairFingerprint(mol, nBits=n_bits), n_bits))})
            except Exception as e:
                print(f"[rdkit_fp] AtomPair {name}: {e}")
            try:
                row.update({f"topo_torsion_{i}": v for i, v in enumerate(
                    sc(rdMolDescriptors.GetHashedTopologicalTorsionFingerprint(mol, nBits=n_bits), n_bits))})
            except Exception as e:
                print(f"[rdkit_fp] Torsion {name}: {e}")
        if include_usrcat:
            mol3d = _mol_from_smiles(smi, embed_3d=True)
            if mol3d is not None:
                try:
                    for i, v in enumerate(rdMolDescriptors.GetUSRCAT(mol3d)):
                        row[f"usrcat_{i}"] = float(v)
                except Exception as e:
                    print(f"[rdkit_fp] USRCAT {name}: {e}")
        rows.append(row)
    return _clean(pd.DataFrame(rows).set_index("name"))


# --------------------------------------------------------------------------- #
# mordred
# --------------------------------------------------------------------------- #
def extract_mordred(smiles_map: dict, ignore_3d: bool = True) -> pd.DataFrame:
    if not (have("mordred") and have("rdkit")):
        print("[mordred] mordred/rdkit not installed - skipping")
        return pd.DataFrame()
    import scipy
    if not hasattr(scipy, "Inf"):
        scipy.Inf = float("inf")  # mordred compat with new scipy
    from mordred import Calculator, descriptors
    calc = Calculator(descriptors, ignore_3D=ignore_3d)
    names, mols = [], []
    for name, smi in smiles_map.items():
        mol = _mol_from_smiles(smi)
        if mol is not None:
            names.append(name)
            mols.append(mol)
    if not mols:
        return pd.DataFrame()
    df = calc.pandas(mols, quiet=True)
    df.index = names
    df.index.name = "name"
    df = df.apply(pd.to_numeric, errors="coerce").add_prefix("mordred_")
    return _clean(df)


# --------------------------------------------------------------------------- #
# deepchem
# --------------------------------------------------------------------------- #
def extract_deepchem(smiles_map: dict, n_bits: int = 256) -> pd.DataFrame:
    if not have("deepchem"):
        print("[deepchem] deepchem not installed - skipping")
        return pd.DataFrame()
    import deepchem as dc
    featurizers = []
    try:
        featurizers.append(("dc_circular", dc.feat.CircularFingerprint(size=n_bits, radius=2)))
    except Exception as e:
        print(f"[deepchem] CircularFingerprint: {e}")
    try:
        featurizers.append(("dc_rdkit2d", dc.feat.RDKitDescriptors()))
    except Exception as e:
        print(f"[deepchem] RDKitDescriptors: {e}")
    pairs = [(n, s) for n, s in smiles_map.items() if s]
    if not featurizers or not pairs:
        return pd.DataFrame()
    names, smiles = zip(*pairs)
    blocks = []
    for prefix, fz in featurizers:
        try:
            arr = np.asarray(fz.featurize(list(smiles)))
            if arr.ndim == 1:
                arr = np.vstack(arr)
            df = pd.DataFrame(arr, index=names)
            df.columns = [f"{prefix}_{i}" for i in range(df.shape[1])]
            blocks.append(df.apply(pd.to_numeric, errors="coerce"))
        except Exception as e:
            print(f"[deepchem] {prefix} FAIL: {e}")
    if not blocks:
        return pd.DataFrame()
    out = pd.concat(blocks, axis=1)
    out.index.name = "name"
    return _clean(out)


# --------------------------------------------------------------------------- #
# qm - autoqchem gaussian_log_extractor (from .log files)
# --------------------------------------------------------------------------- #
_MOL_KEYS = ["dipole", "homo_energy", "lumo_energy", "electronegativity", "hardness",
             "E", "G", "ZPE", "H", "E_scf", "zero_point_correction", "E_zpe",
             "electronic_spatial_extent", "molar_mass", "molar_volume"]
_ATOM_KEYS = ["Mulliken_charge", "APT_charge", "NPA_charge", "NPA_core", "NPA_valence",
              "NPA_Rydberg", "NPA_total", "VBur", "NMR_shift", "NMR_anisotropy"]


def extract_qm(log_dir: str) -> pd.DataFrame:
    if not have("autoqchem"):
        print("[qm] autoqchem not installed - skipping")
        return pd.DataFrame()
    import autoqchem.gaussian_log_extractor as gle
    log_map = {Path(p).stem: p for p in sorted(glob.glob(os.path.join(log_dir, "*.log")))}
    rows = []
    for name, path in log_map.items():
        row = {"name": name}
        try:
            res = gle.gaussian_log_extractor(path).get_descriptors()
            desc = res.get("descriptors", {})
            row.update({k: desc[k] for k in _MOL_KEYS if k in desc and desc[k] is not None})
            for key in _ATOM_KEYS:
                vals = res.get("atom_descriptors", {}).get(key)
                if vals:
                    arr = np.array(vals, dtype=float)
                    row[f"{key}_mean"] = float(np.nanmean(arr))
                    row[f"{key}_std"] = float(np.nanstd(arr))
                    row[f"{key}_max"] = float(np.nanmax(arr))
                    row[f"{key}_min"] = float(np.nanmin(arr))
        except Exception as e:
            print(f"[qm] FAIL {name}: {e}")
        rows.append(row)
    return _clean(pd.DataFrame(rows).set_index("name"))


# --------------------------------------------------------------------------- #
# aqme_qdescp - AQME QDESCP xTB / MORFEUS descriptors (from .xyz)
# --------------------------------------------------------------------------- #
def extract_aqme_qdescp(xyz_map: dict, work_dir: str, level: str = "full",
                        atoms=None, charge: int = 0, mult: int = 1,
                        solvent=None, xtb_opt: bool = True, nprocs: int = 1,
                        overwrite: bool = False) -> pd.DataFrame:
    if not have("aqme"):
        print("[aqme_qdescp] aqme not installed - skipping (needs aqme + xTB)")
        return pd.DataFrame()
    import shutil
    from aqme.qdescp import qdescp

    def _read_table(wd):
        cands = [os.path.join(wd, f"AQME-Descriptors_{level}_AQME_run.csv"),
                 os.path.join(wd, f"AQME-ROBERT_{level}_AQME_run.csv"),
                 os.path.join(wd, "QDESCP", f"AQME-Descriptors_{level}_AQME_run.csv"),
                 os.path.join(wd, "QDESCP", f"AQME-ROBERT_{level}_AQME_run.csv")]
        for path in cands:
            if os.path.exists(path):
                df = pd.read_csv(path)
                ncol = "code_name" if "code_name" in df.columns else df.columns[0]
                df[ncol] = df[ncol].astype(str).str.replace(r"_conf_\d+$", "", regex=True)
                df = df.rename(columns={ncol: "name"}).set_index("name")
                df = df.drop(columns=[c for c in ["SMILES"] if c in df.columns])
                return _clean(df.apply(pd.to_numeric, errors="coerce"))
        raise FileNotFoundError(f"no AQME QDESCP {level} CSV in {wd}")

    input_dir = os.path.join(work_dir, "input_xyz")
    os.makedirs(input_dir, exist_ok=True)
    if not overwrite:
        try:
            cached = _read_table(work_dir)
            print(f"[aqme_qdescp] using cached table {cached.shape}")
            return cached
        except FileNotFoundError:
            pass
    safe = []
    for name, src in xyz_map.items():
        dst = os.path.join(input_dir, f"{name}.xyz")
        shutil.copy2(src, dst)
        safe.append(dst)
    old = os.getcwd()
    os.chdir(work_dir)
    try:
        qdescp(files=safe, destination=os.path.join(work_dir, "QDESCP"),
               program="xtb", charge=charge, mult=mult, nprocs=nprocs,
               qdescp_atoms=atoms or [], qdescp_solvent=solvent,
               xtb_opt=xtb_opt, boltz=True, robert=False)
    finally:
        os.chdir(old)
    return _read_table(work_dir)


# --------------------------------------------------------------------------- #
# aqme CSEARCH - generate a conformer ensemble from SMILES, write .xyz files.
# This is a *preprocessing* step (not a descriptor engine): it produces an
# xyz_dir you then point the run config / picker at.  End-to-end AQME path:
#     SMILES csv --(CSEARCH)--> xyz ensemble --(aqme_qdescp / morfeus / ...)--> features
# --------------------------------------------------------------------------- #
def _rdkit_mol_to_xyz(mol, name: str) -> str:
    conf = mol.GetConformer()
    lines = [str(mol.GetNumAtoms()), name]
    for atom in mol.GetAtoms():
        p = conf.GetAtomPosition(atom.GetIdx())
        lines.append(f"{atom.GetSymbol():2s} {p.x:14.8f} {p.y:14.8f} {p.z:14.8f}")
    return "\n".join(lines) + "\n"


def run_aqme_csearch(smiles_csv: str, out_xyz_dir: str,
                     name_col: str = "name", smiles_col: str = "smiles",
                     program: str = "rdkit", sample: int = 10,
                     charge: int = 0, mult: int = 1, keep: str = "lowest"):
    """
    AQME CSEARCH conformer search from a SMILES CSV, converting the SDF output to
    .xyz files in `out_xyz_dir`.

    program : 'rdkit' | 'summ' | 'fullmonte' | 'crest'   (crest needs CREST installed)
    sample  : conformers to generate per molecule
    keep    : 'lowest' -> one representative xyz per molecule (lowest energy),
              'all'    -> name_conf_1.xyz, name_conf_2.xyz, ...

    Returns out_xyz_dir on success, else None.
    """
    if not have("aqme"):
        print("[csearch] aqme not installed - skipping (pip install aqme)")
        return None
    if not have("rdkit"):
        print("[csearch] rdkit needed to convert conformers - skipping")
        return None
    from aqme.csearch import csearch
    from rdkit import Chem

    os.makedirs(out_xyz_dir, exist_ok=True)
    work = os.path.join(out_xyz_dir, "_csearch_work")
    os.makedirs(work, exist_ok=True)

    # AQME expects a CSV with 'code_name' + 'SMILES' columns.
    df = pd.read_csv(smiles_csv)
    if name_col not in df.columns or smiles_col not in df.columns:
        print(f"[csearch] CSV needs columns '{name_col}' and '{smiles_col}'; got {list(df.columns)}")
        return None
    tmp_csv = os.path.join(work, "csearch_input.csv")
    (df.rename(columns={name_col: "code_name", smiles_col: "SMILES"})[["code_name", "SMILES"]]
       .to_csv(tmp_csv, index=False))

    old = os.getcwd()
    os.chdir(work)
    try:
        csearch(input=tmp_csv, program=program, sample=sample,
                charge=charge, mult=mult, destination=os.path.join(work, "CSEARCH"))
    except Exception as e:
        os.chdir(old)
        print(f"[csearch] AQME csearch failed: {e}")
        return None
    finally:
        os.chdir(old)

    sdfs = glob.glob(os.path.join(work, "CSEARCH", "**", "*.sdf"), recursive=True)
    written = []
    for sdf in sdfs:
        try:
            mols = [m for m in Chem.SDMolSupplier(sdf, removeHs=False) if m is not None]
        except Exception as e:
            print(f"[csearch] read {os.path.basename(sdf)}: {e}")
            continue
        if not mols:
            continue
        stem = os.path.splitext(os.path.basename(sdf))[0]
        name = re.sub(r"_(rdkit|crest|summ|fullmonte).*$", "", stem)
        chosen = mols if keep == "all" else mols[:1]   # AQME orders confs by energy
        for k, m in enumerate(chosen, 1):
            label = name if keep != "all" else f"{name}_conf_{k}"
            fname = f"{label}.xyz"
            with open(os.path.join(out_xyz_dir, fname), "w") as fh:
                fh.write(_rdkit_mol_to_xyz(m, label))
            written.append(fname)
    print(f"[csearch] wrote {len(written)} xyz files to {out_xyz_dir}")
    return out_xyz_dir if written else None


def list_aqme_descriptors():
    """Print the QDESCP descriptor names available at each level (helps pick `level`)."""
    if not have("aqme"):
        print("[aqme] not installed (pip install aqme)")
        return None
    from aqme.qdescp_utils import collect_descp_lists
    d = collect_descp_lists()
    print("QDESCP molecular descriptors:")
    print("  denovo   :", d.get("denovo_mols"))
    print("  interpret:", [x for x in d.get("interpret_mols", []) if x not in d.get("denovo_mols", [])])
    print("QDESCP atom-centered descriptors:")
    print("  denovo   :", d.get("denovo_atoms"))
    print("  interpret:", [x for x in d.get("interpret_atoms", []) if x not in d.get("denovo_atoms", [])])
    return d


# --------------------------------------------------------------------------- #
# morfeus_suite - direct morfeus descriptors, atom-anchored (from .xyz)
# --------------------------------------------------------------------------- #
def extract_morfeus_suite(xyz_map: dict, descriptors=None,
                          sterimol_pairs=None, metal_index=None,
                          cone_atoms=None, bv_radius: float = 3.5,
                          pyramid_atoms=None) -> pd.DataFrame:
    """
    descriptors: subset of
        ['sterimol','buried_volume','cone_angle','sasa','dispersion','pyramidalization']
    Atom indices are 1-indexed (morfeus convention == picker convention).
    """
    if not have("morfeus"):
        print("[morfeus_suite] morfeus-ml not installed - skipping")
        return pd.DataFrame()
    from morfeus import read_xyz
    import morfeus as mf

    want = set(descriptors or ["sterimol", "buried_volume", "sasa", "dispersion"])
    rows = []
    for name, path in xyz_map.items():
        row = {"name": name}
        try:
            elements, coords = read_xyz(path)
        except Exception as e:
            print(f"[morfeus_suite] read FAIL {name}: {e}")
            rows.append(row)
            continue

        if "sterimol" in want and sterimol_pairs:
            for i, j in sterimol_pairs:
                try:
                    s = mf.Sterimol(elements, coords, i, j)
                    tag = f"{i}-{j}"
                    row[f"mf_L_{tag}"] = float(s.L_value)
                    row[f"mf_B1_{tag}"] = float(s.B_1_value)
                    row[f"mf_B5_{tag}"] = float(s.B_5_value)
                except Exception as e:
                    print(f"[morfeus_suite] sterimol {name} {i}-{j}: {e}")

        if "buried_volume" in want and metal_index:
            for mi in (metal_index if isinstance(metal_index, list) else [metal_index]):
                try:
                    bv = mf.BuriedVolume(elements, coords, mi, radius=bv_radius)
                    row[f"mf_vbur_{mi}"] = float(bv.fraction_buried_volume)
                except Exception as e:
                    print(f"[morfeus_suite] buried_volume {name} {mi}: {e}")

        if "cone_angle" in want and cone_atoms:
            for ci in (cone_atoms if isinstance(cone_atoms, list) else [cone_atoms]):
                try:
                    ca = mf.ConeAngle(elements, coords, ci)
                    row[f"mf_cone_angle_{ci}"] = float(ca.cone_angle)
                except Exception as e:
                    print(f"[morfeus_suite] cone_angle {name} {ci}: {e}")

        if "sasa" in want:
            try:
                sa = mf.SASA(elements, coords)
                row["mf_sasa_area"] = float(sa.area)
                row["mf_sasa_volume"] = float(sa.volume)
            except Exception as e:
                print(f"[morfeus_suite] sasa {name}: {e}")

        if "dispersion" in want:
            try:
                dp = mf.Dispersion(elements, coords)
                row["mf_disp_area"] = float(getattr(dp, "area", np.nan))
                row["mf_disp_volume"] = float(getattr(dp, "volume", np.nan))
                if hasattr(dp, "p_int"):
                    row["mf_disp_p_int"] = float(dp.p_int)
            except Exception as e:
                print(f"[morfeus_suite] dispersion {name}: {e}")

        if "pyramidalization" in want and pyramid_atoms:
            for pi in (pyramid_atoms if isinstance(pyramid_atoms, list) else [pyramid_atoms]):
                try:
                    py = mf.Pyramidalization(coords, pi)
                    row[f"mf_pyr_P_{pi}"] = float(getattr(py, "P", np.nan))
                except Exception as e:
                    print(f"[morfeus_suite] pyramidalization {name} {pi}: {e}")

        rows.append(row)
    return _clean(pd.DataFrame(rows).set_index("name"))
