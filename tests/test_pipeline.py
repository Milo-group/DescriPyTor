"""The cluster pipeline without a cluster: numbering, metal placement, stage scripts, the job chain."""
import json
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
