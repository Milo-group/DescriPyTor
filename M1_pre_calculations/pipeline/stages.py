"""SGE job scripts for each stage, rendered from one skeleton.

Every job: pick its molecule from ids.txt (line PIPE_TASK), take the previous stage's ``out/<id>.xyz`` (the build
for stage 1), refuse it if the element order differs from the manifest, run in node scratch
(falling back to $TMPDIR), write ``out/<id>.xyz`` and the marker ``status/<id>.<stage>``:
``running``, ``done`` or ``failed: <reason>``. A task whose output is already done exits at once,
so resubmitting a stage only redoes what is missing.

Placeholders are ``@@KEY@@``, so bash's own ``$`` and braces need no escaping.
"""
from __future__ import annotations

import shlex

from .checks import check_awk
from .protocol import Protocol, Stage

SKELETON = r"""#!/bin/bash
#$ -N @@JOB@@
#$ -S /bin/bash
#$ -q @@QUEUE@@
#$ -pe shared @@CORES@@
#$ -j y
#$ -o @@ROOT@@/@@STAGE@@/logs
@@DIRECTIVES@@
# @@PROTOCOL@@, stage @@STAGE@@ (@@KIND@@). Written by descripytor pipeline.
set -u
ROOT=@@ROOT@@
STAGE=@@STAGE@@
# the molecule: its line in ids.txt, from qsub -v PIPE_TASK=n (N1GE 6.0 rejects qsub -t),
# or the array task id if this script is ever submitted as an array
TASK=${PIPE_TASK:-${SGE_TASK_ID:-}}
NAME=$(sed -n "${TASK}p" "$ROOT/ids.txt")
[ -n "$NAME" ] || { echo "no molecule for task $TASK"; exit 1; }
OUT="$ROOT/$STAGE/out"; WORK="$ROOT/$STAGE/work"; MARK="$ROOT/status/$NAME.$STAGE"
CHECKS="$ROOT/checks/$NAME.$STAGE"
mkdir -p "$OUT" "$WORK" "$ROOT/status" "$ROOT/checks"
fail() { echo "failed: $*" > "$MARK"; echo "FAILED $NAME: $*"; exit 1; }
reject() { mkdir -p "$OUT/rejected"; mv -f "$OUT/$NAME".* "$OUT/rejected/" 2>/dev/null; fail "$*"; }
if [ "$(cat "$MARK" 2>/dev/null)" = "done" ]; then echo "skip $NAME: done"; exit 0; fi
PREVMARK="@@PREVMARK@@"
if [ -n "$PREVMARK" ] && [ "$(cat "$ROOT/status/$NAME.$PREVMARK" 2>/dev/null)" != "done" ]; then
  fail "previous stage $PREVMARK is not done ($(head -c 80 "$ROOT/status/$NAME.$PREVMARK" 2>/dev/null || echo waiting))"
fi
rm -f "$CHECKS"
IN="@@PREV@@/$NAME.xyz"
[ -f "$IN" ] || fail "missing input $IN"
got=$(awk 'NR>2 && NF>=4 {print $1}' "$IN" | xargs)
want=$(xargs < "$ROOT/elements/$NAME.elements")
[ "$got" = "$want" ] || fail "element order of $IN differs from the manifest"
echo running > "$MARK"
SCR=@@SCRATCH@@
if mkdir -p "$SCR/$USER" 2>/dev/null; then SCR="$SCR/$USER"; else SCR="${TMPDIR:-/tmp}"; fi
W="$SCR/@@JOB@@_${JOB_ID}.${TASK}"
rm -rf "$W"; mkdir -p "$W"; cd "$W" || fail "no scratch $W"
trap 'cd /; rm -rf "$W"' EXIT
echo "=== $NAME $STAGE host=$(hostname) scratch=$W start=$(date) ==="
cp "$IN" start.xyz
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
@@BODY@@
@@STRUCTCHECK@@
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
grep -Eiq "aborting the run|error termination" "$WORK/$NAME.out" && fail "orca: $(grep -Ei 'aborting the run|error termination' "$WORK/$NAME.out" | head -n 1)"
grep -q "ORCA TERMINATED NORMALLY" "$WORK/$NAME.out" || fail "orca did not terminate normally"
[ -f "$NAME.globalminimum.xyz" ] || fail "no globalminimum.xyz"
grep -qi "did not converge" "$WORK/$NAME.out" && echo "WARN some GOAT optimizations did not converge" >> "$CHECKS"
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
ERR="optimization did not converge|scf not converged|scf is not converged|aborting the run|error termination"
grep -Eiq "$ERR" "$WORK/$NAME.out" && fail "orca: $(grep -Ei "$ERR" "$WORK/$NAME.out" | head -n 1)"
grep -q "ORCA TERMINATED NORMALLY" "$WORK/$NAME.out" || fail "orca did not terminate normally"
cp -f "$NAME".property.txt "$NAME".hess "$WORK/" 2>/dev/null
if grep -qiE '^!.*\b(num)?freq\b' "$NAME.inp" && ! grep -q "VIBRATIONAL FREQUENCIES" "$WORK/$NAME.out"; then
  fail "frequencies were requested but none were printed"
fi
imag=$(awk -v t=@@IMAG_TOL@@ '/imaginary mode/ {f = $2 + 0; if (-f > t) printf "%s ", $2}' "$WORK/$NAME.out")
[ -z "$imag" ] || fail "imaginary frequencies (cm-1): $imag"
awk -v t=@@IMAG_TOL@@ '/imaginary mode/ {f = $2 + 0; if (-f <= t) printf "WARN small imaginary frequency %s cm-1 (below the %s tolerance)\n", $2, t}' "$WORK/$NAME.out" >> "$CHECKS"
if [ -f "$NAME.xyz" ]; then cp "$NAME.xyz" "$OUT/$NAME.xyz"; else cp start.xyz "$OUT/$NAME.xyz"; fi
"""

XTB_BODY = r"""export PATH="@@ORCA@@:$PATH" LD_LIBRARY_PATH="@@ORCA@@:@@ORCA@@/lib:${LD_LIBRARY_PATH:-}" XTBPATH="@@ORCA@@"
export OMP_NUM_THREADS=@@CORES@@
ERR="abnormal termination|failed to converge|#ERROR|convergence criteria cannot be satisfied"
xtb_ok() {  # $1 output file, $2 label
  grep -Eiq "$ERR" "$1" && fail "xtb $2: $(grep -Ei "$ERR" "$1" | head -n 1)"
  grep -q "normal termination of xtb" "$1" || fail "xtb $2 did not terminate normally"
}
cp start.xyz cur.xyz
tries=0
while :; do
  rm -f xtbopt.xyz xtbhess.xyz vibspectrum
  "@@ORCA@@/otool_xtb" cur.xyz @@XTBMODE@@ @@LEVEL@@ --chrg @@CHARGE@@ --uhf @@UHF@@ > opt.out 2>&1
  cp opt.out "$WORK/$NAME.opt.out"
  xtb_ok opt.out optimization
  [ -f xtbopt.xyz ] || fail "no xtbopt.xyz"
  grep -q "GEOMETRY OPTIMIZATION CONVERGED" opt.out || fail "xtb optimization did not converge"
  imag=""
  if [ "@@XTBMODE@@" = "--ohess" ]; then
    [ -f vibspectrum ] || fail "no vibspectrum from --ohess"
    cp vibspectrum "$WORK/$NAME.vibspectrum"
    # rows are "mode [symmetry] wavenumber intensity selection": the first decimal after the mode
    imag=$(awk -v t=@@IMAG_TOL@@ '!/^[$#]/ && NF >= 4 {
        f = ($2 ~ /^-?[0-9]+\.[0-9]+$/) ? $2 : $3
        if (f !~ /^-?[0-9]+\.[0-9]+$/) { bad++; next }
        n++; if (f + 0 < -t) printf "%s ", f
      } END { if (n == 0 || bad) printf "UNREADABLE" }' vibspectrum)
    case "$imag" in *UNREADABLE*) fail "could not read the frequencies in vibspectrum";; esac
  fi
  [ -z "$imag" ] && break
  if [ "$tries" -lt @@IMAG_RETRY@@ ] && [ -f xtbhess.xyz ]; then
    tries=$((tries + 1)); echo "imaginary $imag cm-1: restart $tries from xtbhess.xyz (displaced along the mode)"
    cp xtbhess.xyz cur.xyz; continue
  fi
  fail "imaginary frequencies (cm-1): $imag after $tries restart(s)"
done
[ "$tries" -gt 0 ] && echo "WARN imaginary frequency removed by $tries restart(s) along the mode" >> "$CHECKS"
"@@ORCA@@/otool_xtb" xtbopt.xyz --sp --chrg @@CHARGE@@ --uhf @@UHF@@ > sp.out 2>&1 || fail "xtb single point exit $?"
cp sp.out "$WORK/$NAME.sp.out"
xtb_ok sp.out "single point"
cp charges "$OUT/$NAME.q"; cp wbo "$OUT/$NAME.wbo"
dip=$(awk '/molecular dipole/{g=1} g&&/full:/{print $2, $3, $4; exit}' sp.out)
hl=$(awk '/\(HOMO\)/{h=$(NF-1)} /\(LUMO\)/{if (l=="") l=$(NF-1)} END{print h, l}' sp.out)
echo "$NAME $dip $hl" > "$OUT/$NAME.props"
cp xtbopt.xyz "$OUT/$NAME.xyz"
"""

UMA_BODY = r"""export HF_TOKEN=$(cat @@HF_TOKEN_FILE@@) OMP_NUM_THREADS=@@CORES@@
"@@UMA_ENV@@/bin/python" "$ROOT/$STAGE/uma_one.py" start.xyz "$OUT/$NAME.xyz" @@CHARGE@@ @@MULT@@ @@FMAX@@ @@STEPS@@ @@MODEL@@ @@TASK@@ \
    > "$WORK/$NAME.uma.log" 2>&1 || fail "uma exit $? ($(tail -n 1 "$WORK/$NAME.uma.log"))"
cat "$WORK/$NAME.uma.log"
if sed -n 2p "$OUT/$NAME.xyz" | grep -q "converged=False"; then
  if [ "@@ALLOW_UNCONV@@" = "yes" ]; then echo "WARN UMA stopped at the step limit: $(sed -n 2p "$OUT/$NAME.xyz")" >> "$CHECKS"
  else reject "UMA did not converge within the step limit: $(sed -n 2p "$OUT/$NAME.xyz" | cut -c1-90)"; fi
fi
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
grep -q "Error termination" "$OUT/$NAME.log" && reject "gaussian: $(grep -B 3 'Error termination' "$OUT/$NAME.log" | head -n 1 | cut -c1-80)"
tail -n 3 "$OUT/$NAME.log" | grep -q "Normal termination" || reject "gaussian did not terminate normally"
if grep -qi "freq" "$NAME.com" && ! grep -q "Frequencies --" "$OUT/$NAME.log"; then
  reject "frequencies were requested but none were printed"
fi
imag=$(awk -v t=@@IMAG_TOL@@ '/Frequencies --/ {for (i = 3; i <= NF; i++) if ($i + 0 < -t) printf "%s ", $i}' "$OUT/$NAME.log")
[ -z "$imag" ] || reject "imaginary frequencies (cm-1): $imag"
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
    [ "$v" = "done" ] && [ -s "$ROOT/checks/$id.$s" ] && v="done*"
    printf '\t%s' "$v"
  done
  printf '\n'
done < "$ROOT/ids.txt"
echo "== checks"
for f in "$ROOT"/checks/*; do [ -s "$f" ] && printf '%s: %s\n' "$(basename "$f")" "$(paste -sd ';' - < "$f" | cut -c1-200)"; done
echo "== queue"
qstat -u "$USER" 2>/dev/null | awk -v p="@@PREFIX@@" 'NR>2 && index($3, p) == 1 {n[$3" "$5]++} END {for (k in n) print k, n[k]}'
"""

CHAIN = r"""#!/bin/bash
# submit_chain.sh PLAN - written by descripytor pipeline.
# Each PLAN line is "task stage stage ...": that molecule's stages are submitted as single-task
# jobs, each held on the one before (-hold_jid), so every molecule moves on as soon as its own
# previous step ends. N1GE 6.0 has no -terse, no -hold_jid_ad and rejects qsub -t, so the
# molecule goes in as -v PIPE_TASK=n.
ROOT=@@ROOT@@
while read -r t rest; do
  [ -z "$t" ] && continue
  prev=""
  for s in $rest; do
    out=$(qsub -v PIPE_TASK="$t" ${prev:+-hold_jid "$prev"} "$ROOT/$s/run.sh" < /dev/null 2>&1)
    id=$(printf '%s\n' "$out" | sed -n 's/^Your job[-a-z]* \([0-9][0-9]*\).*/\1/p' | head -n 1)
    if [ -z "$id" ]; then echo "QSUB_FAIL $t $s $(printf '%s' "$out" | tr '\n' ' ')"; break; fi
    echo "JOB $t $s $id"
    prev=$id
  done
done < "$1"
"""

DEFAULTS = {
    "goat": dict(keywords="! GOAT XTB", maxcore=2000),
    "orca": dict(keywords="! r2SCAN-3c Opt Freq", maxcore=3000, imag_tol=20),
    "xtb": dict(level="--gfn 2", hess=True, imag_tol=20, imag_retry=1),
    "uma": dict(fmax=0.05, steps=300, model="uma-s-1p1", task="omol", allow_unconverged=False),
    "gaussian": dict(route="", mem_gb=32, tail="", imag_tol=20),          # route is required: the level is a choice
}
BODIES = {"goat": GOAT_BODY, "orca": ORCA_BODY, "xtb": XTB_BODY, "uma": UMA_BODY, "gaussian": GAUSSIAN_BODY}


STRUCTCHECK = '# structure check on the output: bonds against the build, clashes, donors and bite angle\nCHK=$(awk -v ref="$ROOT/refs/$NAME.bonds" -v donors="$(cat "$ROOT/refs/$NAME.donors" 2>/dev/null)" \\\n    -v anc="$(cat "$ROOT/refs/$NAME.anc" 2>/dev/null)" -f "$ROOT/check_structure.awk" "$OUT/$NAME.xyz")\nprintf \'%s\\n\' "$CHK" | grep \'^WARN\' >> "$CHECKS"\nif printf \'%s\\n\' "$CHK" | grep -q \'^FAIL\'; then\n  printf \'%s\\n\' "$CHK" | grep \'^FAIL\' >> "$CHECKS"\n  reject "structure: $(printf \'%s\\n\' "$CHK" | grep \'^FAIL\' | head -n 3 | cut -c6- | paste -sd \';\' -)"\nfi\n[ -s "$CHECKS" ] || rm -f "$CHECKS"'


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
        IMAG_TOL=float(opt.get("imag_tol", 20)), IMAG_RETRY=int(opt.get("imag_retry", 1)),
        XTBMODE="--ohess" if opt.get("hess", True) else "--opt", ALLOW_UNCONV="yes" if opt.get("allow_unconverged") else "no",
        PREVMARK="" if index == 1 else p.stages[index - 2].dirname(index - 1),
    )
    body = _fill(BODIES[stage.kind], values)
    check = "" if stage.kind == "gaussian" else STRUCTCHECK          # gaussian writes a log, not an xyz
    files = {f"{sdir}/run.sh": _fill(SKELETON.replace("@@BODY@@", body.rstrip("\n")).replace("@@STRUCTCHECK@@", check), values)}
    if stage.kind == "uma":
        files[f"{sdir}/uma_one.py"] = UMA_ONE
    return files


def render_all(p: Protocol, n: int) -> dict[str, str]:
    files = {}
    for i in range(1, len(p.stages) + 1):
        files.update(render_stage(p, i, n))
    files["submit_chain.sh"] = _fill(CHAIN, dict(ROOT=p.remote_root))
    files["check_structure.awk"] = check_awk()
    files["status.sh"] = _fill(STATUS, dict(ROOT=p.remote_root, STAGES=" ".join(p.stage_dirs()), PREFIX=job_prefix(p)))
    return files
