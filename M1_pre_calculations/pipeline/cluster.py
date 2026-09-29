"""Talk to an SGE cluster over ssh: upload a run, chain its stages, report, retry, fetch.

Everything runs as bash on the login node (no python there). Each molecule's stages are plain jobs
chained with ``qsub -hold_jid`` and told their molecule with ``-v PIPE_TASK=n`` (the cluster's
N1GE 6.0 has no ``-terse``, no per-task ``-hold_jid_ad``, and rejects ``qsub -t``), so nothing has
to poll: a molecule moves on the moment its previous step ends. A failed step writes ``failed: <reason>`` and the steps after it
fail at once with "missing input"; ``retry`` resubmits exactly those.

File transfer is base64 over ssh, because the login banner corrupts scp. The login shell (tcsh)
swallows ssh's stdin and rejects any command-line word over about 8000 characters, so an upload
goes as 4000-character parts, each written to its own numbered file (a retry cannot duplicate
data), joined remotely and checked by md5. Downloads are a tar stream between markers.
"""
from __future__ import annotations

import base64
import hashlib
import io
import json
import os
import re
import subprocess
import tarfile
import time
from pathlib import Path

from .protocol import Protocol
from .stages import job_prefix, render_all

BEGIN, END = "___DPT_BEGIN___", "___DPT_END___"
PART = 4000            # base64 characters per ssh call; tcsh rejects words over ~8000


class ClusterError(RuntimeError):
    pass


class Remote:
    def __init__(self, host: str, key: str, timeout: int = 180):
        self.host, self.key, self.timeout = host, os.path.expanduser(key), timeout

    def _ssh(self, command: str) -> str:
        cmd = ["ssh", "-n", "-T", "-o", "BatchMode=yes", "-o", "ConnectTimeout=30", "-o", "IdentitiesOnly=yes",
               "-i", self.key, self.host, command]
        r = subprocess.run(cmd, capture_output=True, timeout=self.timeout)
        return (r.stdout or b"").decode("utf-8", "replace")

    def run(self, script: str, tries: int = 3) -> str:
        """Run a bash script remotely; returns its stdout. Retries when the reply comes back empty
        or cut short (the end marker is missing). Scripts must be safe to run twice."""
        payload = base64.b64encode(f"echo {BEGIN}\n{script}\necho {END}\n".encode()).decode()
        if len(payload) > 7000:
            raise ClusterError("remote script too long for the tcsh login shell; upload it as a file instead")
        last = ""
        for _ in range(tries):
            last = self._ssh(f"echo {payload} | base64 -d | bash")
            if BEGIN in last and END in last:
                return last.split(BEGIN, 1)[1].rsplit(END, 1)[0].strip("\n")
            time.sleep(2)
        raise ClusterError(f"no complete reply from {self.host}: {last[-300:]!r}")

    def upload_tree(self, local: Path, remote: str) -> None:
        buf = io.BytesIO()
        with tarfile.open(fileobj=buf, mode="w") as tar:
            for f in sorted(local.rglob("*")):
                if f.is_file():
                    tar.add(f, arcname=f.relative_to(local).as_posix())
        data = buf.getvalue()
        digest = hashlib.md5(data).hexdigest()
        b64 = base64.b64encode(data).decode()
        self.run(f"mkdir -p {remote} && rm -f {remote}/.upload.part.*")
        for k, i in enumerate(range(0, len(b64), PART)):
            for _ in range(3):
                if "PART_OK" in self._ssh(f"printf %s {b64[i:i + PART]} > {remote}/.upload.part.{k:06d} && echo PART_OK"):
                    break
            else:
                raise ClusterError(f"upload part {k} failed three times")
        out = self.run(f"cd {remote} && cat .upload.part.* > .upload.b64 && rm -f .upload.part.* "
                       f"&& base64 -d .upload.b64 > .upload.tar && md5sum .upload.tar")
        if digest not in out:
            raise ClusterError(f"upload corrupted: local md5 {digest}, remote {out!r}")
        self.run(f"cd {remote} && tar xf .upload.tar && rm -f .upload.tar .upload.b64")

    def download_tree(self, remote: str, paths: list[str], local: Path) -> None:
        out = self.run(f"cd {remote} && tar cf - {' '.join(paths)} 2>/dev/null | base64 -w0; echo")
        data = base64.b64decode("".join(out.split()))
        local.mkdir(parents=True, exist_ok=True)
        with tarfile.open(fileobj=io.BytesIO(data)) as tar:
            try:
                tar.extractall(local, filter="data")
            except TypeError:                    # Python < 3.12 has no filter argument
                tar.extractall(local)


def remote_for(p: Protocol) -> Remote:
    return Remote(p.cluster["host"], p.cluster["key"])


def _ids(p: Protocol) -> list[str]:
    f = p.workdir / "ids.txt"
    if not f.exists():
        raise ClusterError("nothing built yet: run `descripytor pipeline build` first")
    return [x for x in f.read_text(encoding="utf-8").split() if x]


def stage_files(p: Protocol) -> dict[str, str]:
    return render_all(p, len(_ids(p)))


def plan_lines(p: Protocol, tasks) -> str:
    """Plan text for submit_chain.sh: one line per (task, first stage)."""
    stages = p.stage_dirs()
    return "".join(f"{t} {' '.join(stages[first - 1:])}\n" for t, first in tasks)


def _put_file(r: Remote, remote_dir: str, name: str, text: str) -> None:
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        Path(tmp, name).write_text(text, encoding="utf-8", newline="\n")
        r.upload_tree(Path(tmp), remote_dir)


def _run_chain(r: Remote, p: Protocol, plan: str, plan_name: str) -> dict:
    """Upload a plan and submit it. Returns {task: {stage: job id}}. Raises on any qsub failure,
    reporting what was submitted, and never resubmits blindly."""
    _put_file(r, p.remote_root, plan_name, plan)
    out = r.run(f"bash {p.remote_root}/submit_chain.sh {p.remote_root}/{plan_name}", tries=1)
    jobs, failed = {}, []
    for line in out.splitlines():
        f = line.split()
        if f[:1] == ["JOB"] and len(f) == 4:
            jobs.setdefault(f[1], {})[f[2]] = f[3]
        elif f[:1] == ["QSUB_FAIL"]:
            failed.append(line)
    if failed:
        raise ClusterError(f"qsub failed ({len(failed)}); submitted so far: {jobs}; first failure: {failed[0]}")
    return jobs


def submit(p: Protocol, dry_run: bool = False, force: bool = False, from_stage: int = 1,
           only: list[str] | None = None) -> dict:
    """Upload the run and submit each molecule's chain from ``from_stage`` (1-based). ``only``
    limits it to those ids (e.g. the ones adopted at the stage before). Job ids are merged into
    jobs.json; a molecule whose chain from that stage is already on record is refused unless
    ``force``."""
    ids = _ids(p)
    files = stage_files(p)
    stage_dir = p.workdir / "remote"
    for name, text in files.items():
        f = stage_dir / name
        f.parent.mkdir(parents=True, exist_ok=True)
        f.write_text(text, encoding="utf-8", newline="\n")
    for sub in ("build", "elements"):
        for f in (p.workdir / sub).glob("*"):
            (stage_dir / sub).mkdir(exist_ok=True)
            (stage_dir / sub / f.name).write_bytes(f.read_bytes())
    for name in ("ids.txt", "manifest.json"):
        (stage_dir / name).write_bytes((p.workdir / name).read_bytes())
    (stage_dir / "protocol.json").write_bytes(p.path.read_bytes())
    if only is not None:
        unknown = set(only) - set(ids)
        if unknown:
            raise ClusterError(f"not molecules of this run: {sorted(unknown)}")
    tasks = [(t, from_stage) for t, mid in enumerate(ids, 1) if only is None or mid in only]
    first = p.stage_dirs()[from_stage - 1]
    rec_file = p.workdir / "jobs.json"
    record = json.loads(rec_file.read_text(encoding="utf-8")) if rec_file.exists() else dict(root=p.remote_root, n=len(ids), jobs={})
    clash = [ids[t - 1] for t, _ in tasks if first in record["jobs"].get(str(t), {})]
    if clash and not force:
        raise ClusterError(f"already submitted at {first}: {clash}; use status / retry, or --force")
    plan = plan_lines(p, tasks)
    if dry_run:
        return dict(dry_run=True, n=len(tasks), files=sorted(files), plan=plan.splitlines(), local=str(stage_dir))
    if not tasks:
        return dict(record, note="nothing to submit")

    r = remote_for(p)
    root = p.remote_root
    r.upload_tree(stage_dir, root)
    r.run(f"cd {root} && chmod +x *.sh */run.sh && mkdir -p status " + " ".join(f"{d}/logs {d}/out {d}/work" for d in p.stage_dirs()))
    jobs = _run_chain(r, p, plan, "plan.txt")
    for t, stages in jobs.items():
        record["jobs"].setdefault(t, {}).update(stages)
    record.setdefault("submissions", []).append(dict(at=time.strftime("%Y-%m-%d %H:%M:%S"), from_stage=first,
                                                     molecules=[ids[int(t) - 1] for t in jobs]))
    rec_file.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    _put_file(r, root, "jobs.json", json.dumps(record, indent=2) + "\n")
    return record


def status(p: Protocol) -> tuple[list[dict], str]:
    """Rows {id, <stage>: marker} and the queue summary."""
    out = remote_for(p).run(f"bash {p.remote_root}/status.sh")
    table, _, queue = out.partition("== queue")
    lines = [l for l in table.strip().splitlines() if l.strip()]
    head = lines[0].split("\t")
    rows = [dict(zip(head, l.split("\t"))) for l in lines[1:]]
    return rows, queue.strip()


def summarize(p: Protocol, rows: list[dict]) -> dict:
    out = {}
    for s in p.stage_dirs():
        c = {}
        for r in rows:
            v = r.get(s, "waiting")
            v = "failed" if v.startswith("failed") else v
            c[v] = c.get(v, 0) + 1
        out[s] = c
    return out


def retry(p: Protocol, dry_run: bool = False) -> list[dict]:
    """Resubmit every molecule from its first failed stage onwards. A stage left 'running' with
    none of this pipeline's jobs in the queue (a killed task) counts as failed."""
    rows, queue = status(p)
    idle = not queue.strip()
    ids = _ids(p)
    stages = p.stage_dirs()
    plan = []
    for r in rows:
        for k, s in enumerate(stages, 1):
            v = r.get(s, "waiting")
            if v.startswith("failed") or (idle and v == "running"):
                plan.append(dict(id=r["id"], task=ids.index(r["id"]) + 1, from_stage=k, reason=v))
                break
    if dry_run or not plan:
        return plan
    rem = remote_for(p)
    for item in plan:
        clear = " ".join(f"{p.remote_root}/status/{item['id']}.{st}" for st in stages[item["from_stage"] - 1:])
        rem.run(f"rm -f {clear}")
    jobs = _run_chain(rem, p, plan_lines(p, [(i["task"], i["from_stage"]) for i in plan]), "retry_plan.txt")
    for item in plan:
        item["jobs"] = jobs.get(str(item["task"]), {})
    return plan


def watch(p: Protocol, interval: int = 600, echo=print) -> dict:
    """Report every ``interval`` seconds until every molecule is done or failed at some stage."""
    last = p.stage_dirs()[-1]
    while True:
        rows, queue = status(p)
        summary = summarize(p, rows)
        echo(time.strftime("%H:%M ") + "  ".join(f"{s} {c}" for s, c in summary.items()))
        finished = all(r.get(last) == "done" or any(r.get(s, "").startswith("failed") for s in p.stage_dirs()) for r in rows)
        if finished and not queue.strip():
            return summary
        time.sleep(interval)


def fetch(p: Protocol, all_stages: bool = False) -> Path:
    """Copy the last stage's outputs (or every stage's, with ``all_stages``) and the manifest into
    <workdir>/fetched/."""
    dirs = p.stage_dirs() if all_stages else p.stage_dirs()[-1:]
    dest = p.workdir / "fetched"
    remote_for(p).download_tree(p.remote_root, ["manifest.json", "ids.txt", "status"] + [f"{d}/out" for d in dirs], dest)
    return dest

