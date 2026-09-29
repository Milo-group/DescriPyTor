"""Talk to an SGE cluster over ssh: upload a run, chain its stages, report, retry, fetch.

Everything runs as bash on the login node (no python there). Each molecule's stage k waits only
for its own stage k-1 (``qsub -hold_jid_ad``), so nothing has to poll: a molecule moves on the
moment its previous step ends. A failed step writes ``failed: <reason>`` and the steps after it
fail at once with "missing input"; ``retry`` resubmits exactly those.

File transfer is base64 over ssh, because the login banner corrupts scp: uploads go in chunks
appended remotely and are checked by md5; downloads are a tar stream between markers.
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
CHUNK = 20000          # base64 characters per ssh call; stays under the Windows command-line limit


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
        self.run(f"mkdir -p {remote} && rm -f {remote}/.upload.b64")
        for i in range(0, len(b64), CHUNK):
            out = self._ssh(f"printf %s {b64[i:i + CHUNK]} >> {remote}/.upload.b64 && echo CHUNK_OK")
            if "CHUNK_OK" not in out:
                raise ClusterError(f"upload chunk {i // CHUNK} failed")
        out = self.run(f"cd {remote} && base64 -d .upload.b64 > .upload.tar && md5sum .upload.tar")
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


def qsub_commands(p: Protocol, first_stage: int = 1, task: int | None = None) -> list[tuple[str, str]]:
    """(stage dir, qsub command) in submission order. ``{prev}`` is replaced by the job id of the
    stage before, so each task waits for the same task of that stage (array dependency)."""
    out = []
    for i in range(first_stage, len(p.stages) + 1):
        sdir = p.stages[i - 1].dirname(i)
        t = f" -t {task}-{task}" if task else ""
        hold = " -hold_jid_ad {prev}" if i > first_stage else ""
        out.append((sdir, f"qsub -terse{t}{hold} {p.remote_root}/{sdir}/run.sh"))
    return out


def submit(p: Protocol, dry_run: bool = False, force: bool = False) -> dict:
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
    cmds = qsub_commands(p)
    if dry_run:
        return dict(dry_run=True, n=len(ids), files=sorted(files), commands=[c for _, c in cmds], local=str(stage_dir))

    r = remote_for(p)
    root = p.remote_root
    if not force and "EXISTS" in r.run(f"[ -f {root}/jobs.json ] && echo EXISTS || true"):
        raise ClusterError(f"{root} already has a submitted run; use status / retry, or --force to overwrite")
    r.upload_tree(stage_dir, root)
    r.run(f"cd {root} && chmod +x status.sh */run.sh && mkdir -p status " + " ".join(f"{d}/logs {d}/out {d}/work" for d in p.stage_dirs()))
    jobs, prev = {}, None
    for sdir, cmd in cmds:
        out = r.run(cmd.replace("{prev}", str(prev)) if prev else cmd)
        m = re.search(r"^(\d+)", out.strip().splitlines()[-1] if out.strip() else "")
        if not m:
            raise ClusterError(f"qsub for {sdir} gave no job id: {out!r}")
        jobs[sdir] = prev = m.group(1)
    record = dict(submitted=time.strftime("%Y-%m-%d %H:%M:%S"), root=root, n=len(ids), jobs=jobs)
    (p.workdir / "jobs.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    r.run(f"cat > {root}/jobs.json <<'EOF'\n{json.dumps(record, indent=2)}\nEOF")
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
        clear = " ".join(f"{p.remote_root}/status/{item['id']}.{s}" for s in stages[item["from_stage"] - 1:])
        rem.run(f"rm -f {clear}")
        prev = None
        for sdir, cmd in qsub_commands(p, item["from_stage"], item["task"]):
            out = rem.run(cmd.replace("{prev}", str(prev)) if prev else cmd)
            m = re.search(r"^(\d+)", out.strip().splitlines()[-1] if out.strip() else "")
            if not m:
                raise ClusterError(f"retry qsub for {item['id']} {sdir} gave no job id: {out!r}")
            prev = m.group(1)
        item["job"] = prev
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

