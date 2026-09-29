"""The cluster pipeline without a cluster: numbering, metal placement, stage scripts, the job chain."""
import json
import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

from M1_pre_calculations import pipeline as pl

EXAMPLE = Path(pl.__file__).parent / "examples" / "protocol_bgu.json"


def _protocol(tmp_path, csv_rows, **over):
    raw = json.loads(EXAMPLE.read_text(encoding="utf-8"))
    raw.update(over)
    (tmp_path / "mols.csv").write_text("ligand,smiles\n" + "\n".join(csv_rows) + "\n", encoding="utf-8")
    raw["molecules"] = "mols.csv"
    f = tmp_path / "protocol.json"
    f.write_text(json.dumps(raw), encoding="utf-8")
    return pl.load(f)


def _xyz(path):
    lines = path.read_text().splitlines()
    n = int(lines[0])
    sym = [l.split()[0] for l in lines[2:2 + n]]
    X = np.array([[float(v) for v in l.split()[1:4]] for l in lines[2:2 + n]])
    return sym, X


def test_core_smarts_gives_every_molecule_the_same_numbers(tmp_path):
    p = _protocol(tmp_path, ["pOMe,COc1ccc(cc1)C(=O)O", "mCl,Clc1cccc(c1)C(=O)O"],
                  build={"type": "organic", "core_smarts": "OC(=O)c1ccccc1"})
    pl.build_all(p)
    for mid in ("m001", "m002"):
        sym, _ = _xyz(p.workdir / "build" / f"{mid}.xyz")
        assert sym[:9] == ["O", "C", "O", "C", "C", "C", "C", "C", "C"]      # the core, in SMARTS order
    m = json.loads((p.workdir / "manifest.json").read_text())["molecules"]
    assert m[0]["core_1based"] == list(range(1, 10))


def test_monodentate_metal_sits_on_the_lone_pair(tmp_path):
    p = _protocol(tmp_path, ["PPh3,P(c1ccccc1)(c1ccccc1)c1ccccc1"])
    rows = pl.build_all(p)
    sym, X = _xyz(p.workdir / "build" / "m001.xyz")
    assert sym[:2] == ["Ni", "P"] and rows[0]["donor_1based"] == [2]
    assert abs(np.linalg.norm(X[0] - X[1]) - 2.2) < 1e-6
    assert (p.workdir / "elements" / "m001.elements").read_text().split() == sym


def test_chelate_order_is_metal_ancillaries_donors(tmp_path):
    p = _protocol(tmp_path, ["bipy,c1ccc(nc1)-c1ccccn1"],
                  build={"type": "metal_chelate", "metal": "Ni", "ancillary": "H", "n_ancillary": 2})
    rows = pl.build_all(p)
    sym, X = _xyz(p.workdir / "build" / "m001.xyz")
    assert sym[:5] == ["Ni", "H", "H", "N", "N"]
    assert 65 < rows[0]["bite"] < 105
    assert all(1.9 < np.linalg.norm(X[0] - X[k]) < 2.3 for k in (3, 4))


def test_every_stage_script_renders_and_targets_the_queue(tmp_path):
    stages = [{"kind": "goat"}, {"kind": "uma"}, {"kind": "xtb"}, {"kind": "orca"},
              {"kind": "gaussian", "route": "#p M062X/def2TZVP opt freq pop=(nbo,hirshfeld)"}]
    p = _protocol(tmp_path, ["PMe3,CP(C)C"], stages=stages)
    files = pl.render_all(p, 7)
    assert set(files) >= {"s1_goat/run.sh", "s2_uma/run.sh", "s2_uma/uma_one.py", "s5_gaussian/run.sh", "status.sh"}
    for name, text in files.items():
        assert "@@" not in text
        if name.endswith("run.sh"):
            assert "#$ -q fairshare.q" in text and "#$ -t" not in text and "PIPE_TASK" in text
    assert "/s1_goat/out/$NAME.xyz" in files["s2_uma/run.sh"]            # stage 2 reads stage 1
    assert "/build/$NAME.xyz" in files["s1_goat/run.sh"]
    if shutil.which("bash"):
        for name, text in files.items():
            if name.endswith(".sh"):
                r = subprocess.run(["bash", "-n"], input=text.encode(), capture_output=True)   # bytes: no CRLF on Windows
                assert r.returncode == 0, (name, r.stderr.decode())


def test_submit_dry_run_plans_every_molecule_through_every_stage(tmp_path):
    p = _protocol(tmp_path, ["PMe3,CP(C)C", "PPh3,P(c1ccccc1)(c1ccccc1)c1ccccc1"])
    pl.build_all(p)
    out = pl.submit(p, dry_run=True)
    assert out["n"] == 2 and out["plan"] == ["1 s1_goat s2_uma s3_xtb", "2 s1_goat s2_uma s3_xtb"]
    chain = (Path(out["local"]) / "submit_chain.sh").read_text()
    qsub_line = next(l for l in chain.splitlines() if "out=$(qsub" in l)
    assert '-hold_jid "$prev"' in qsub_line and '-v PIPE_TASK="$t"' in qsub_line
    assert "-hold_jid_ad" not in qsub_line and "-terse" not in qsub_line and " -t " not in qsub_line
    assert pl.plan_lines(p, [(2, 2)]) == "2 s2_uma s3_xtb\n"


def test_submit_chain_script_reads_sge6_output(tmp_path):
    """submit_chain.sh against a fake qsub that answers the way SGE 6 does."""
    if not shutil.which("bash"):
        pytest.skip("no bash")
    p = _protocol(tmp_path, ["PMe3,CP(C)C"])
    chain = pl.render_all(p, 1)["submit_chain.sh"]
    fake = ('qsub() { n=$(( $(cat counter 2>/dev/null || echo 100) + 1 )); echo $n > counter; '
            'echo "$*" >> args; echo "Your job-array $n.1-1:1 (\\"x\\") has been submitted"; }\n')
    script = fake + chain.replace('done < "$1"', "done <<'PLAN'\n1 s1_goat s2_uma s3_xtb\nPLAN")
    r = subprocess.run(["bash"], input=script.encode(), capture_output=True, cwd=tmp_path)
    out = r.stdout.decode()
    assert r.returncode == 0 and out.split().count("JOB") == 3, out + r.stderr.decode()
    args = (tmp_path / "args").read_text().splitlines()
    assert "-hold_jid" not in args[0] and "-hold_jid 101" in args[1] and "-hold_jid 102" in args[2]
    assert all("PIPE_TASK=1" in a for a in args)
    assert "JOB 1 s1_goat 101" in out and "JOB 1 s3_xtb 103" in out


def test_protocol_rejects_bad_input(tmp_path):
    with pytest.raises(pl.ProtocolError):
        _protocol(tmp_path, ["PMe3,CP(C)C"], stages=[{"kind": "gaussian", "route": "#p hf"}, {"kind": "xtb"}])
    with pytest.raises(pl.ProtocolError):
        _protocol(tmp_path, ["PMe3,CP(C)C"], stages=[{"kind": "crest"}])
    p = _protocol(tmp_path, ["PMe3,CP(C)C"], stages=[{"kind": "gaussian"}])
    with pytest.raises(ValueError, match="route"):
        pl.render_all(p, 1)


def test_adopt_permutation_maps_an_outside_run_onto_the_build(tmp_path):
    """An outside run kept RDKit order (Ni first, P where RDKit put it); the build puts P second."""
    p = _protocol(tmp_path, ["PMe3,CP(C)C"])
    pl.build_all(p)
    sym, X = _xyz(p.workdir / "build" / "m001.xyz")
    old_order = [0, 2, 3, 1] + list(range(4, len(sym)))           # P moved to fourth place
    start = tmp_path / "outside.xyz"
    start.write_text(f"{len(sym)}\nx\n" + "".join(f"{sym[k]} {X[k][0]:.8f} {X[k][1]:.8f} {X[k][2]:.8f}\n" for k in old_order))
    perm = pl.permutation(p.workdir / "build" / "m001.xyz", start)
    assert [old_order[j - 1] for j in perm] == list(range(len(sym)))
    X2 = X.copy(); X2[5] += 0.5                                    # a different start structure
    bad = tmp_path / "other.xyz"
    bad.write_text(f"{len(sym)}\nx\n" + "".join(f"{s} {x:.8f} {y:.8f} {z:.8f}\n" for s, (x, y, z) in zip(sym, X2)))
    with pytest.raises(pl.ClusterError):
        pl.permutation(p.workdir / "build" / "m001.xyz", bad)


def test_adopt_script_reorders_and_checks_elements(tmp_path):
    if not shutil.which("bash"):
        pytest.skip("no bash")
    from M1_pre_calculations.pipeline.adopt import ADOPT
    root = tmp_path / "run"
    for d in ("adopt", "elements", "s1_goat/out", "status"):
        (root / d).mkdir(parents=True)
    files = {"elements/m001.elements": "Ni P C\n", "adopt/m001.perm": "1\n3\n2\n",
             "src.xyz": "3\ngoat\nNi 0 0 0\nC 1 0 0\nP 0 1 0\n", "adopt/s1.tsv": "m001\tsrc.xyz\nm002\tnone.xyz\n"}
    for name, text in files.items():
        (root / name).write_bytes(text.encode())                   # LF only, as the pipeline writes them
    script = ADOPT.replace("@@ROOT@@", ".")
    r = subprocess.run(["bash", "-s", "s1_goat", "adopt/s1.tsv"], input=script.encode(), capture_output=True, cwd=root)
    out = r.stdout.decode()
    assert "ADOPTED m001" in out and "MISSING m002" in out, out + r.stderr.decode()
    got = (root / "s1_goat" / "out" / "m001.xyz").read_text().splitlines()[2:5]
    assert [l.split()[0] for l in got] == ["Ni", "P", "C"]
    assert (root / "status" / "m001.s1_goat").read_text().strip() == "done"


# ---------------------------------------------------------------- structure and frequency checks
def _bash_path(path):
    """The folder as the first bash on PATH sees it (Git Bash: /c/..., WSL: /mnt/c/...)."""
    return subprocess.run(["bash"], input=b"pwd" + bytes([10]), capture_output=True, cwd=path).stdout.decode().strip()


def _run_check(tmp_path, p, mid, sym, X):
    from M1_pre_calculations.pipeline.checks import check_awk
    f = tmp_path / "probe.xyz"
    f.write_bytes((f"{len(sym)}\nprobe\n" + "".join(f"{s} {x:.6f} {y:.6f} {z:.6f}\n" for s, (x, y, z) in zip(sym, X))).encode())
    (tmp_path / "check.awk").write_bytes(check_awk().encode())
    for ext in ("bonds", "donors", "anc"):
        (tmp_path / f"ref.{ext}").write_bytes((p.workdir / "refs" / f"{mid}.{ext}").read_bytes())
    cmd = b'awk -v ref=ref.bonds -v donors="$(cat ref.donors)" -v anc="$(cat ref.anc)" -f check.awk probe.xyz' + bytes([10])
    r = subprocess.run(["bash"], input=cmd, capture_output=True, cwd=tmp_path)     # stdin: no Windows quoting
    assert not r.stderr, r.stderr.decode()
    return r.stdout.decode().splitlines()


def test_structure_check_passes_the_build_and_catches_each_fault(tmp_path):
    if not shutil.which("bash"):
        pytest.skip("no bash")
    p = _protocol(tmp_path, ["PMe3,CP(C)C"])
    pl.build_all(p)
    sym, X = _xyz(p.workdir / "build" / "m001.xyz")
    assert _run_check(tmp_path, p, "m001", sym, X) == []                      # awk agrees with the python rule
    c = next(i for i, s in enumerate(sym) if s == "C")
    far = X.copy(); far[1] = X[0] + (X[1] - X[0]) * 1.8                         # P pulled 4 A off the Ni
    assert any("donor P2" in l and l.startswith("FAIL") for l in _run_check(tmp_path, p, "m001", sym, far))
    broken = X.copy(); broken[c] = X[c] + (X[c] - X[1]) * 1.5                   # a P-C bond stretched
    assert any(l.startswith("FAIL broken bond") for l in _run_check(tmp_path, p, "m001", sym, broken))
    clash = X.copy(); clash[c] = X[c + 1] + 0.1                                 # two atoms on top of each other
    assert any("apart" in l and l.startswith("FAIL") for l in _run_check(tmp_path, p, "m001", sym, clash))


def test_structure_check_flags_a_bad_bite_and_warns_on_a_contact(tmp_path):
    if not shutil.which("bash"):
        pytest.skip("no bash")
    p = _protocol(tmp_path, ["bipy,c1ccc(nc1)-c1ccccn1"],
                  build={"type": "metal_chelate", "metal": "Ni", "ancillary": "H", "n_ancillary": 2})
    pl.build_all(p)
    sym, X = _xyz(p.workdir / "build" / "m001.xyz")
    assert _run_check(tmp_path, p, "m001", sym, X) == []
    lost_h = X.copy(); lost_h[1] = X[0] + (X[1] - X[0]) * 3                     # a hydride left the metal
    assert any("ancillary H2" in l for l in _run_check(tmp_path, p, "m001", sym, lost_h))
    extra = np.vstack([X, X[0] + np.array([0.0, 0.0, 2.0])])                     # a carbon parked on the metal
    ref = p.workdir / "refs"
    out = _run_check(tmp_path, p, "m001", sym + ["C"], extra)
    assert any(l.startswith("WARN metal contact C") for l in out)


def test_xtb_stage_restarts_once_along_an_imaginary_mode(tmp_path):
    """The whole xtb run.sh against a fake otool_xtb: the first Hessian has a -85 cm-1 mode and
    writes xtbhess.xyz; the restart is clean. The stage must end done, with a warning."""
    if not shutil.which("bash"):
        pytest.skip("no bash")
    fake = tmp_path / "orca"
    fake.mkdir()
    (fake / "otool_xtb").write_bytes(b"""#!/bin/bash
in=$1; mode=$2
if [ "$mode" = "--sp" ]; then
  printf 'molecular dipole:\n full:  0.1 0.2 0.3 0.4\n  -8.1 (HOMO)\n -6.2 (LUMO)\n normal termination of xtb\n'
  echo "q" > charges; echo "w" > wbo; exit 0
fi
cp "$in" xtbopt.xyz
printf ' *** GEOMETRY OPTIMIZATION CONVERGED AFTER 5 ITERATIONS ***\n normal termination of xtb\n'
n=$(cat "$CALLS" 2>/dev/null || echo 0); echo $((n+1)) > "$CALLS"
if [ "$n" -eq 0 ]; then f=-85.0; cp "$in" xtbhess.xyz; else f=45.0; fi
printf '$vibrational spectrum\n#  mode     symmetry     wave number   IR intensity    selection rules\n#                         cm**(-1)      (km*mol-1)        IR\n     1                      -0.00         0.00000          -\n     7        a            %s         0.63990         YES\n     8        a             28.06         0.10190         YES\n$end\n' "$f" > vibspectrum
""")
    stages = [{"kind": "xtb", "cores": 1}]
    p = _protocol(tmp_path, ["PMe3,CP(C)C"], stages=stages)
    bt = _bash_path(tmp_path)
    raw = json.loads(p.path.read_text())
    raw["cluster"].update(root=f"{bt}/remote_root", orca=f"{bt}/orca", scratch=f"{bt}/scr")
    p.path.write_text(json.dumps(raw)); p = pl.load(p.path)
    pl.build_all(p)
    run = tmp_path / "remote_root" / p.name
    files = pl.render_all(p, 1)
    for sub in ("build", "elements", "refs"):
        (run / sub).mkdir(parents=True, exist_ok=True)
        for f in (p.workdir / sub).glob("*"):
            (run / sub / f.name).write_bytes(f.read_bytes())
    (run / "ids.txt").write_bytes(b"m001\n")
    for name, text in files.items():
        (run / name).parent.mkdir(parents=True, exist_ok=True)
        (run / name).write_bytes(text.encode())
    script = f"export PIPE_TASK=1 JOB_ID=7 USER=t CALLS; chmod +x {bt}/orca/otool_xtb; cd {bt}/remote_root/{p.name}; CALLS=$PWD/.calls bash s1_xtb/run.sh"
    r = subprocess.run(["bash"], input=script.encode(), capture_output=True, cwd=tmp_path)
    log = r.stdout.decode() + r.stderr.decode()
    assert (run / "status" / "m001.s1_xtb").read_text().strip() == "done", log
    assert "restart 1 from xtbhess.xyz" in log
    assert "imaginary frequency removed" in (run / "checks" / "m001.s1_xtb").read_text()
    assert (run / "s1_xtb" / "out" / "m001.props").read_text().split()[-2:] == ["-8.1", "-6.2"]


def test_orca_and_gaussian_imaginary_checks_use_the_tolerance(tmp_path):
    """The imag=$(awk ...) lines of the rendered scripts, run on sample output: -12 cm-1 is under
    the 20 cm-1 tolerance, -85 and -312 are not."""
    if not shutil.which("bash"):
        pytest.skip("no bash")
    stages = [{"kind": "orca"}, {"kind": "gaussian", "route": "#p hf/sto-3g freq"}]
    files = pl.render_all(_protocol(tmp_path, ["PMe3,CP(C)C"], stages=stages), 1)
    samples = {
        "s1_orca/run.sh": ("VIBRATIONAL FREQUENCIES\n   6:     -12.10 cm**-1 ***imaginary mode***\n"
                           "   7:     -85.40 cm**-1 ***imaginary mode***\n   8:      55.00 cm**-1\n"),
        "s2_gaussian/run.sh": (" Frequencies --   -312.1480    -12.0000     45.2000\n"
                               " Frequencies --     88.1000    120.0000    150.0000\n"),
    }
    for name, text in samples.items():
        line = next(l for l in files[name].splitlines() if l.startswith("imag=$(awk"))
        (tmp_path / "sample.out").write_bytes(text.encode())
        cmd = line.replace('"$WORK/$NAME.out"', "sample.out").replace('"$OUT/$NAME.log"', "sample.out") + '; echo "[$imag]"'
        r = subprocess.run(["bash"], input=cmd.encode(), capture_output=True, cwd=tmp_path)
        r.stdout, r.stderr = r.stdout.decode(), r.stderr.decode()
        assert r.stdout.strip() == ("[-85.40 ]" if "orca" in name else "[-312.1480 ]"), (name, r.stdout, r.stderr)


def test_xtb_frequency_parser_reads_the_real_vibspectrum(tmp_path):
    """Rows 1-6 carry no symmetry label, later rows do; one selection column. A file with no
    readable frequency must fail the stage, not pass it."""
    if not shutil.which("bash"):
        pytest.skip("no bash")
    run = pl.render_all(_protocol(tmp_path, ["PMe3,CP(C)C"], stages=[{"kind": "xtb"}]), 1)["s1_xtb/run.sh"]
    lines = run.splitlines()
    a = next(i for i, l in enumerate(lines) if l.strip().startswith("imag=$(awk"))
    b = next(i for i in range(a, len(lines)) if "UNREADABLE" in lines[i] and "printf" in lines[i])
    snippet = "\n".join(lines[a:b + 1]) + '\necho "[$imag]"\n'
    real = ("$vibrational spectrum\n#  mode     symmetry     wave number   IR intensity    selection rules\n"
            "#                         cm**(-1)      (km*mol-1)        IR\n"
            "     1                      -0.00         0.00000          -\n"
            "     7        a            -34.07         0.63990         YES\n"
            "     8        a            -12.88         0.59398         YES\n"
            "     9        a             28.06         0.10190         YES\n$end\n")
    for text, want in ((real, "[-34.07 ]"), ("$vibrational spectrum\n     7   a   x   y   YES\n$end\n", "[UNREADABLE]")):
        (tmp_path / "vibspectrum").write_bytes(text.encode())
        r = subprocess.run(["bash"], input=snippet.encode(), capture_output=True, cwd=tmp_path)
        assert r.stdout.decode().strip() == want, (r.stdout.decode(), r.stderr.decode())


def test_xtb_opt_level_reaches_the_command_line(tmp_path):
    run = pl.render_all(_protocol(tmp_path, ["PMe3,CP(C)C"], stages=[{"kind": "xtb", "opt_level": "vtight", "imag_retry": 0}]), 1)["s1_xtb/run.sh"]
    assert "cur.xyz --ohess vtight --gfn 2" in run and '[ "$tries" -lt 0 ]' in run
    with pytest.raises(ValueError):
        pl.render_all(_protocol(tmp_path, ["PMe3,CP(C)C"], stages=[{"kind": "xtb", "opt_level": "supertight"}]), 1)
