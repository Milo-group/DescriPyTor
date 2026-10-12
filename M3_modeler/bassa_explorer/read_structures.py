"""Structures for an uploaded table, read under the environment that has the DescriPyTor readers (the server may run under another one, the one with bassa_reg).

    <anaconda python> read_structures.py <folder>      compound names as a JSON list on stdin; JSON on stdout: {found: {name: {xyz, dipole, charges}}, missing: [...]}

A compound is matched to the file whose name is exactly its name, or else contains it as a whole word (L1 matches FL_lig_L1_opt, not L12),
or else is the only file carrying its number (L1 matches FL_lig_1_str when no other file has a 1).
A .feather or .log brings the dipole and the charges with the geometry; an .xyz brings the geometry only.
"""
import json, re, sys
from pathlib import Path

from settings import S

sys.path.insert(0, str(Path(S["descripytor"]) / "M2_data_extractor" / "feature_explorer"))
KINDS = (".feather", ".log", ".xyz")                                            # in order of preference: what carries more comes first


def match(name, files):
    word = re.compile(r"(?<![A-Za-z0-9])" + re.escape(name) + r"(?![A-Za-z0-9])")
    for same in (True, False):
        hit = [f for f in files if (f.stem == name if same else word.search(f.stem))]
        if hit: return min(hit, key=lambda f: (KINDS.index(f.suffix.lower()), len(f.stem), f.name))
    num = re.search(r"\d+$", name)                                              # last resort, L1 against FL_lig_1_str: the compound's number as a whole number, and only when one molecule has it
    if num:
        hit = [f for f in files if re.search(r"(?<!\d)" + num.group() + r"(?!\d)", f.stem)]
        if len({f.stem for f in hit}) == 1: return min(hit, key=lambda f: KINDS.index(f.suffix.lower()))
    return None


def main(folder, names):
    files = [f for f in Path(folder).iterdir() if f.suffix.lower() in KINDS]; found, missing = {}, []
    for n in names:
        f = match(n, files)
        if f is None: missing.append(n); continue
        try:
            if f.suffix.lower() == ".xyz": found[n] = dict(xyz=f.read_text(encoding="utf-8", errors="replace").replace("\r", ""), dipole=None, charges=None)
            else:
                import data as TD
                p = TD.payload(str(f), n); found[n] = dict(xyz=p["xyz"], dipole=p.get("dipole"), charges=p.get("charges") or None)
        except Exception as e: missing.append("%s (%s: %s)" % (n, f.name, e))
    return dict(found=found, missing=missing)


if __name__ == "__main__":
    F = [Path(x) for x in ("L1.xyz", "FL_lig_L1_opt.feather", "L12.xyz", "FL_lig_1_str.feather")]
    assert match("L1", F).name == "L1.xyz" and match("L12", F).name == "L12.xyz" and match("1", F).name == "FL_lig_1_str.feather" and match("L3", F) is None
    sys.stdout.reconfigure(encoding="utf-8"); print(json.dumps(main(sys.argv[1], json.loads(sys.stdin.read()))))
