"""Bring finished outside calculations into a pipeline run as a completed stage.

A structure computed before the pipeline existed (say a GOAT minimum from a hand-written job)
becomes ``<stage>/out/<id>.xyz`` with the stage marked done, so the pipeline continues from the
next stage. The outside run started from its own input file, possibly in another atom order. The
permutation is found by matching the pipeline's build to that input atom by atom on coordinates,
which only works when both were built identically; anything else is refused rather than guessed.

    sources.csv: id, remote_source, local_start
        id            pipeline id (m001 ...)
        remote_source finished structure on the cluster, in the outside run's atom order
        local_start   the outside run's input xyz (same atom order), used for the permutation
"""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from .cluster import ClusterError, remote_for, stage_files
from .protocol import Protocol

ADOPT = r"""#!/bin/bash
# adopt.sh STAGE TSV - written by descripytor pipeline. Each TSV line: id <TAB> remote source.
# The source is reordered with adopt/<id>.perm (source atom numbers in pipeline order), checked
# against elements/<id>.elements, written to STAGE/out/<id>.xyz and marked done.
ROOT=@@ROOT@@
STAGE=$1
mkdir -p "$ROOT/$STAGE/out" "$ROOT/status"
while IFS=$'\t' read -r id src; do
  [ -z "$id" ] && continue
  mark="$ROOT/status/$id.$STAGE"
  if [ "$(cat "$mark" 2>/dev/null)" = "done" ]; then echo "ALREADY $id"; continue; fi
  if [ ! -s "$src" ]; then echo "MISSING $id $src"; continue; fi
  tmp="$ROOT/$STAGE/out/.$id.adopt"
  awk -v permfile="$ROOT/adopt/$id.perm" -v src="$src" '
    BEGIN { while ((getline l < permfile) > 0) perm[++np] = l }
    NR == 1 { n = $1; next }
    NR == 2 { comment = $0; next }
    NR <= n + 2 { line[NR - 2] = sprintf("%-2s %16.8f %16.8f %16.8f", $1, $2, $3, $4) }
    END {
      if (np != n) { print "PERM_MISMATCH" > "/dev/stderr"; exit 1 }
      print n; print "adopted from " src " | " comment
      for (i = 1; i <= np; i++) print line[perm[i]]
    }' "$src" > "$tmp" || { echo "FAILED $id perm"; echo "failed: adopt permutation" > "$mark"; rm -f "$tmp"; continue; }
  got=$(awk 'NR>2 && NF>=4 {print $1}' "$tmp" | xargs)
  want=$(xargs < "$ROOT/elements/$id.elements")
  if [ "$got" != "$want" ]; then echo "FAILED $id elements"; echo "failed: adopted structure has the wrong element order" > "$mark"; rm -f "$tmp"; continue; fi
  mv "$tmp" "$ROOT/$STAGE/out/$id.xyz"
  echo done > "$mark"
  echo "ADOPTED $id"
done < "$2"
"""


def _read_xyz(path):
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    n = int(lines[0].split()[0])
    rows = [l.split() for l in lines[2:2 + n]]
    return [r[0] for r in rows], np.array([[float(v) for v in r[1:4]] for r in rows])


def permutation(build_xyz, start_xyz, tol=1e-4) -> list[int]:
    """1-based atom numbers of ``start_xyz`` in the order of ``build_xyz`` (matched on coordinates)."""
    sb, xb = _read_xyz(build_xyz)
    ss, xs = _read_xyz(start_xyz)
    if sorted(sb) != sorted(ss):
        raise ClusterError(f"{Path(start_xyz).name}: different atoms from the pipeline build")
    perm = []
    for s, x in zip(sb, xb):
        d = np.linalg.norm(xs - x, axis=1)
        j = int(d.argmin())
        if d[j] > tol or ss[j] != s:
            raise ClusterError(f"{Path(start_xyz).name}: no atom at the build's position (off by {d[j]:.3g} A); "
                               "the outside run did not start from this build")
        perm.append(j + 1)
    if len(set(perm)) != len(perm):
        raise ClusterError(f"{Path(start_xyz).name}: two build atoms matched one start atom")
    return perm


def read_sources(path) -> list[dict]:
    """Rows of the sources CSV; a relative local_start is read from the CSV's folder."""
    base = Path(path).resolve().parent
    with open(path, newline="", encoding="utf-8-sig") as fh:
        rows = [dict(r) for r in csv.DictReader(fh)]
    for r in rows:
        if not Path(r["local_start"]).is_absolute():
            r["local_start"] = str((base / r["local_start"]).resolve())
    return rows


def adopt(p: Protocol, stage: int, sources: list[dict]) -> dict:
    """Adopt finished structures as stage ``stage`` (1-based). Uploads the run files (no jobs are
    submitted), writes the permutations, runs adopt.sh. Returns {adopted, already, missing, failed};
    adopted ids are also appended to <workdir>/adopted.txt."""
    sdir = p.stage_dirs()[stage - 1]
    tree = p.workdir / "remote"
    for name, text in stage_files(p).items():
        f = tree / name
        f.parent.mkdir(parents=True, exist_ok=True)
        f.write_text(text, encoding="utf-8", newline="\n")
    for sub in ("build", "elements"):
        (tree / sub).mkdir(exist_ok=True)
        for f in (p.workdir / sub).glob("*"):
            (tree / sub / f.name).write_bytes(f.read_bytes())
    for name in ("ids.txt", "manifest.json"):
        (tree / name).write_bytes((p.workdir / name).read_bytes())
    (tree / "protocol.json").write_bytes(p.path.read_bytes())
    (tree / "adopt").mkdir(exist_ok=True)
    ids = set((p.workdir / "ids.txt").read_text(encoding="utf-8").split())
    rows = []
    for s in sources:
        if s["id"] not in ids:
            raise ClusterError(f"adopt: {s['id']} is not a molecule of this run")
        perm = permutation(p.workdir / "build" / f"{s['id']}.xyz", s["local_start"])
        (tree / "adopt" / f"{s['id']}.perm").write_text("\n".join(map(str, perm)) + "\n", encoding="utf-8", newline="\n")
        rows.append(f"{s['id']}\t{s['remote_source']}\n")
    (tree / "adopt" / f"{sdir}.tsv").write_text("".join(rows), encoding="utf-8", newline="\n")
    (tree / "adopt.sh").write_text(ADOPT.replace("@@ROOT@@", p.remote_root), encoding="utf-8", newline="\n")

    r = remote_for(p)
    root = p.remote_root
    r.upload_tree(tree, root)
    r.run(f"cd {root} && chmod +x *.sh */run.sh && mkdir -p status " + " ".join(f"{d}/logs {d}/out {d}/work" for d in p.stage_dirs()))
    out = r.run(f"bash {root}/adopt.sh {sdir} {root}/adopt/{sdir}.tsv", tries=1)
    res = {k: [] for k in ("adopted", "already", "missing", "failed")}
    key = {"ADOPTED": "adopted", "ALREADY": "already", "MISSING": "missing", "FAILED": "failed"}
    for line in out.splitlines():
        f = line.split()
        if f and f[0] in key:
            res[key[f[0]]].append(f[1])
    done_file = p.workdir / "adopted.txt"
    old = set(done_file.read_text(encoding="utf-8").split()) if done_file.exists() else set()
    done_file.write_text("\n".join(sorted(old | set(res["adopted"]) | set(res["already"]))) + "\n", encoding="utf-8")
    return res
