"""SGE array scripts for each stage, rendered from one skeleton.

Every task: pick its molecule from ids.txt, take the previous stage's ``out/<id>.xyz`` (the build
for stage 1), refuse it if the element order differs from the manifest, run in node scratch
(falling back to $TMPDIR), write ``out/<id>.xyz`` and the marker ``status/<id>.<stage>``:
``running``, ``done`` or ``failed: <reason>``. A task whose output is already done exits at once,
so resubmitting a stage only redoes what is missing.

Placeholders are ``@@KEY@@``, so bash's own ``$`` and braces need no escaping.
"""
from __future__ import annotations

import shlex

from .protocol import Protocol, Stage

SKELETON = r"""#!/bin/bash
#$ -N @@JOB@@
#$ -S /bin/bash
#$ -q @@QUEUE@@
#$ -pe shared @@CORES@@
#$ -j y
#$ -o @@ROOT@@/@@STAGE@@/logs
#$ -t 1-@@N@@
@@DIRECTIVES@@
# @@PROTOCOL@@, stage @@STAGE@@ (@@KIND@@). Written by descripytor pipeline.
set -u
ROOT=@@ROOT@@
STAGE=@@STAGE@@
NAME=$(sed -n "${SGE_TASK_ID}p" "$ROOT/ids.txt")
[ -n "$NAME" ] || { echo "no molecule for task $SGE_TASK_ID"; exit 1; }
OUT="$ROOT/$STAGE/out"; WORK="$ROOT/$STAGE/work"; MARK="$ROOT/status/$NAME.$STAGE"
mkdir -p "$OUT" "$WORK" "$ROOT/status"
fail() { echo "failed: $*" > "$MARK"; echo "FAILED $NAME: $*"; exit 1; }
if [ "$(cat "$MARK" 2>/dev/null)" = "done" ]; then echo "skip $NAME: done"; exit 0; fi
IN="@@PREV@@/$NAME.xyz"
[ -f "$IN" ] || fail "missing input $IN"
got=$(awk 'NR>2 && NF>=4 {print $1}' "$IN" | xargs)
want=$(xargs < "$ROOT/elements/$NAME.elements")
[ "$got" = "$want" ] || fail "element order of $IN differs from the manifest"
echo running > "$MARK"
SCR=@@SCRATCH@@
if mkdir -p "$SCR/$USER" 2>/dev/null; then SCR="$SCR/$USER"; else SCR="${TMPDIR:-/tmp}"; fi
W="$SCR/@@JOB@@_${JOB_ID}.${SGE_TASK_ID}"
rm -rf "$W"; mkdir -p "$W"; cd "$W" || fail "no scratch $W"
trap 'cd /; rm -rf "$W"' EXIT
echo "=== $NAME $STAGE host=$(hostname) scratch=$W start=$(date) ==="
cp "$IN" start.xyz
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
@@BODY@@
echo done > "$MARK"
echo "=== $NAME end=$(date) ==="
"""

ORCA_ENV = r"""export PATH="@@ORCA@@:@@MPI@@/bin:$PATH"
export LD_LIBRARY_PATH="@@ORCA@@:@@ORCA@@/lib:@@MPI@@/lib:${LD_LIBRARY_PATH:-}"
export OMPI_MCA_btl=self,vader,tcp XTBPATH="@@ORCA@@"
"""

GOAT_BODY = ORCA_ENV + r"""cat > "$NAME.inp" <<'ORCAINPUT'
%pal nprocs @@CORES@@ end
%maxcore @@MAXCORE@@
@@KEYWORDS@@
* xyzfile @@CHARGE@@ @@MULT@@ start.xyz
ORCAINPUT
"@@ORCA@@/orca" "$W/$NAME.inp" > "$WORK/$NAME.out" 2>&1 || fail "orca exit $?"
grep -q "ORCA TERMINATED NORMALLY" "$WORK/$NAME.out" || fail "orca did not terminate normally"
[ -f "$NAME.globalminimum.xyz" ] || fail "no globalminimum.xyz"
cp -f "$NAME".finalensemble.xyz "$WORK/" 2>/dev/null
awk 'NR<=2 || NF>=4' "$NAME.globalminimum.xyz" > "$OUT/$NAME.xyz"
"""

ORCA_BODY = ORCA_ENV + r"""cat > "$NAME.inp" <<'ORCAINPUT'
%pal nprocs @@CORES@@ end
%maxcore @@MAXCORE@@
@@KEYWORDS@@
* xyzfile @@CHARGE@@ @@MULT@@ start.xyz
ORCAINPUT
"@@ORCA@@/orca" "$W/$NAME.inp" > "$WORK/$NAME.out" 2>&1 || fail "orca exit $?"
grep -q "ORCA TERMINATED NORMALLY" "$WORK/$NAME.out" || fail "orca did not terminate normally"
cp -f "$NAME".property.txt "$NAME".hess "$WORK/" 2>/dev/null
nimag=$(grep -c "\*\*\*imaginary mode\*\*\*" "$WORK/$NAME.out")
[ "$nimag" -eq 0 ] || fail "$nimag imaginary frequencies"
if [ -f "$NAME.xyz" ]; then cp "$NAME.xyz" "$OUT/$NAME.xyz"; else cp start.xyz "$OUT/$NAME.xyz"; fi
"""

XTB_BODY = r"""export PATH="@@ORCA@@:$PATH" LD_LIBRARY_PATH="@@ORCA@@:@@ORCA@@/lib:${LD_LIBRARY_PATH:-}" XTBPATH="@@ORCA@@"
export OMP_NUM_THREADS=@@CORES@@
"@@ORCA@@/otool_xtb" start.xyz --opt @@LEVEL@@ --chrg @@CHARGE@@ --uhf @@UHF@@ > opt.out 2>&1
cp opt.out "$WORK/$NAME.opt.out"
[ -f xtbopt.xyz ] || fail "no xtbopt.xyz"
grep -q "GEOMETRY OPTIMIZATION CONVERGED" opt.out || fail "xtb optimization did not converge"
"@@ORCA@@/otool_xtb" xtbopt.xyz --sp --chrg @@CHARGE@@ --uhf @@UHF@@ > sp.out 2>&1 || fail "xtb single point exit $?"
cp sp.out "$WORK/$NAME.sp.out"; cp charges "$OUT/$NAME.q"; cp wbo "$OUT/$NAME.wbo"
dip=$(awk '/molecular dipole/{g=1} g&&/full:/{print $2, $3, $4; exit}' sp.out)
hl=$(awk '/\(HOMO\)/{h=$(NF-1)} /\(LUMO\)/{if (l=="") l=$(NF-1)} END{print h, l}' sp.out)
echo "$NAME $dip $hl" > "$OUT/$NAME.props"
cp xtbopt.xyz "$OUT/$NAME.xyz"
"""

UMA_BODY = r"""export HF_TOKEN=$(cat @@HF_TOKEN_FILE@@) OMP_NUM_THREADS=@@CORES@@
"@@UMA_ENV@@/bin/python" "$ROOT/$STAGE/uma_one.py" start.xyz "$OUT/$NAME.xyz" @@CHARGE@@ @@MULT@@ @@FMAX@@ @@STEPS@@ @@MODEL@@ @@TASK@@ \
    > "$WORK/$NAME.uma.log" 2>&1 || fail "uma exit $?"
cat "$WORK/$NAME.uma.log"
"""

UMA_ONE = r'''"""Reoptimize one structure with a FAIRChem UMA model. Written by descripytor pipeline.
usage: uma_one.py in.xyz out.xyz charge multiplicity fmax steps model task"""
import sys

from ase.io import read
from ase.optimize import LBFGS
from fairchem.core import FAIRChemCalculator, pretrained_mlip

src, dst, charge, mult, fmax, steps, model, task = sys.argv[1:9]
atoms = read(src)
atoms.info.update(charge=int(charge), spin=int(mult))
atoms.calc = FAIRChemCalculator(pretrained_mlip.get_predict_unit(model, device="cpu"), task_name=task)
opt = LBFGS(atoms, logfile="-")
converged = opt.run(fmax=float(fmax), steps=int(steps))
e = atoms.get_potential_energy()
f = float(((atoms.get_forces() ** 2).sum(axis=1) ** 0.5).max())
with open(dst, "w", newline="\n") as fh:
    fh.write(f"{len(atoms)}\nUMA {model} E={e:.8f} eV fmax={f:.4f} steps={opt.get_number_of_steps()} converged={bool(converged)}\n")
    for s, (x, y, z) in zip(atoms.get_chemical_symbols(), atoms.get_positions()):
        fh.write(f"{s:<2s}{x:20.10f}{y:20.10f}{z:20.10f}\n")
print("converged" if converged else "NOT converged", f, flush=True)
'''

GAUSSIAN_BODY = r"""export g16root=@@G16ROOT@@ GAUSS_SCRDIR="$W"
source "$g16root/g16/bsd/g16.profile"
{
  printf '%%nprocshared=%s\n%%mem=%sGB\n%%chk=%s.chk\n' @@CORES@@ @@MEM@@ "$NAME"
  cat <<'ROUTE'
@@ROUTE@@
ROUTE
  printf '\n%s\n\n%s %s\n' "$NAME" @@CHARGE@@ @@MULT@@
  awk 'NR>2 && NF>=4 {printf "%-2s %16.8f %16.8f %16.8f\n", $1, $2, $3, $4}' start.xyz
  printf '\n'
  cat <<'TAIL'
@@TAIL@@
TAIL
  printf '\n'
} > "$NAME.com"
cp "$NAME.com" "$WORK/"
"$g16root/g16/g16" < "$NAME.com" > "$OUT/$NAME.log" 2>&1
tail -n 3 "$OUT/$NAME.log" | grep -q "Normal termination" || fail "gaussian did not terminate normally"
nimag=$(grep -c "imaginary frequencies (negative Signs)" "$OUT/$NAME.log")
[ "$nimag" -eq 0 ] || fail "imaginary frequencies"
"""

STATUS = r"""#!/bin/bash
# Pipeline status from the markers on disk. Written by descripytor pipeline. No python.
ROOT=@@ROOT@@
STAGES="@@STAGES@@"
printf 'id'; for s in $STAGES; do printf '\t%s' "$s"; done; printf '\n'
while IFS= read -r id || [ -n "$id" ]; do
  [ -z "$id" ] && continue
  printf '%s' "$id"
  for s in $STAGES; do
    m="$ROOT/status/$id.$s"
    if [ -f "$m" ]; then v=$(head -n 1 "$m" | cut -c1-70); else v=waiting; fi
    printf '\t%s' "$v"
  done
  printf '\n'
done < "$ROOT/ids.txt"
echo "== queue"
qstat -u "$USER" 2>/dev/null | awk -v p="@@PREFIX@@" 'NR>2 && index($3, p) == 1 {n[$3" "$5]++} END {for (k in n) print k, n[k]}'
"""

DEFAULTS = {
    "goat": dict(keywords="! GOAT XTB", maxcore=2000),
    "orca": dict(keywords="! r2SCAN-3c Opt Freq", maxcore=3000),
    "xtb": dict(level="--gfn 2"),
    "uma": dict(fmax=0.05, steps=300, model="uma-s-1p1", task="omol"),
    "gaussian": dict(route="", mem_gb=32, tail=""),          # route is required: the level is a choice
}
BODIES = {"goat": GOAT_BODY, "orca": ORCA_BODY, "xtb": XTB_BODY, "uma": UMA_BODY, "gaussian": GAUSSIAN_BODY}


def job_prefix(p: Protocol) -> str:
    return (p.name[:6] or "pipe").replace("-", "_")


def _fill(text: str, values: dict) -> str:
    for k, v in values.items():
        text = text.replace(f"@@{k}@@", str(v))
    if "@@" in text:
        raise ValueError(f"unfilled placeholder in stage script: {text[text.index('@@'):][:40]}")
    return text


def render_stage(p: Protocol, index: int, n: int) -> dict[str, str]:
    """{filename: text} for stage ``index`` (1-based) of ``p`` over ``n`` molecules."""
    stage: Stage = p.stages[index - 1]
    sdir = stage.dirname(index)
    opt = dict(DEFAULTS[stage.kind], **stage.options)
    if stage.kind == "gaussian" and not opt["route"].strip():
        raise ValueError("a gaussian stage needs its route line, e.g. \"route\": \"#p M062X/def2TZVP opt freq pop=nbo\"")
    c = p.cluster
    root = p.remote_root
    prev = f"{root}/build" if index == 1 else f"{root}/{p.stages[index - 2].dirname(index - 1)}/out"
    directives = "\n".join(f"#$ -l {r}" for r in opt.get("resources", []))
    values = dict(
        JOB=f"{job_prefix(p)}{index}{stage.kind[:3]}", QUEUE=c["queue"], CORES=stage.cores, ROOT=root, STAGE=sdir,
        N=n, DIRECTIVES=directives, PROTOCOL=p.name, KIND=stage.kind, PREV=prev, SCRATCH=c["scratch"],
        ORCA=c["orca"], MPI=c["mpi"], CHARGE=p.charge, MULT=p.multiplicity, UHF=p.multiplicity - 1,
        MAXCORE=opt.get("maxcore", 2000), KEYWORDS=opt.get("keywords", ""), LEVEL=opt.get("level", ""),
        FMAX=opt.get("fmax", 0.05), STEPS=opt.get("steps", 300), MODEL=shlex.quote(str(opt.get("model", ""))),
        TASK=shlex.quote(str(opt.get("task", ""))), UMA_ENV=c.get("uma_env", ""), HF_TOKEN_FILE=c.get("hf_token_file", ""),
        G16ROOT=c.get("g16root", ""), MEM=opt.get("mem_gb", 32), ROUTE=opt.get("route", ""), TAIL=opt.get("tail", ""),
    )
    body = _fill(BODIES[stage.kind], values)
    files = {f"{sdir}/run.sh": _fill(SKELETON.replace("@@BODY@@", body.rstrip("\n")), values)}
    if stage.kind == "uma":
        files[f"{sdir}/uma_one.py"] = UMA_ONE
    return files


def render_all(p: Protocol, n: int) -> dict[str, str]:
    files = {}
    for i in range(1, len(p.stages) + 1):
        files.update(render_stage(p, i, n))
    files["status.sh"] = _fill(STATUS, dict(ROOT=p.remote_root, STAGES=" ".join(p.stage_dirs()), PREFIX=job_prefix(p)))
    return files
