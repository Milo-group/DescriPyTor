"""The pipeline protocol: one JSON file per dataset.

Everything a run needs is in this file, so the file is also the record of how the structures
were made. Relative paths are read from the protocol file's folder. See docs/PIPELINE.md.
"""
from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from pathlib import Path

STAGE_KINDS = ("goat", "uma", "xtb", "orca", "gaussian")
DEFAULT_CORES = {"goat": 8, "uma": 8, "xtb": 4, "orca": 16, "gaussian": 16}
BUILD_TYPES = ("organic", "metal_mono", "metal_chelate")
CLUSTER_KEYS = ("host", "key", "root", "queue", "orca", "mpi", "scratch")


class ProtocolError(ValueError):
    pass


@dataclass
class Stage:
    kind: str
    cores: int
    options: dict = field(default_factory=dict)

    def dirname(self, index: int) -> str:
        return f"s{index}_{self.kind}"


@dataclass
class Protocol:
    path: Path
    name: str
    molecules: Path
    smiles_column: str
    name_column: str | None
    id_column: str | None
    charge: int
    multiplicity: int
    build: dict
    stages: list[Stage]
    cluster: dict

    @property
    def folder(self) -> Path:
        return self.path.parent

    @property
    def workdir(self) -> Path:
        """Local folder that mirrors the remote run root."""
        return self.folder / f"{self.name}_pipeline"

    @property
    def remote_root(self) -> str:
        return f"{self.cluster['root'].rstrip('/')}/{self.name}"

    def stage_dirs(self) -> list[str]:
        return [s.dirname(i) for i, s in enumerate(self.stages, 1)]


def load(path) -> Protocol:
    path = Path(path).resolve()
    raw = json.loads(path.read_text(encoding="utf-8"))
    name = raw.get("name", "")
    if not re.fullmatch(r"[A-Za-z0-9_-]{1,40}", name):
        raise ProtocolError("'name' must be 1-40 characters of letters, digits, _ or -")
    if "molecules" not in raw:
        raise ProtocolError("'molecules' (a CSV with a SMILES column) is required")
    build = dict(raw.get("build", {"type": "organic"}))
    if build.get("type") not in BUILD_TYPES:
        raise ProtocolError(f"build.type must be one of {BUILD_TYPES}")
    if build["type"] != "organic" and not build.get("metal"):
        raise ProtocolError("build.metal is required for a metal build")

    stages = []
    for s in raw.get("stages", []):
        s = dict(s)
        kind = s.pop("kind", None)
        if kind not in STAGE_KINDS:
            raise ProtocolError(f"stage kind {kind!r} is not one of {STAGE_KINDS}")
        stages.append(Stage(kind, int(s.pop("cores", DEFAULT_CORES[kind])), s))
    if not stages:
        raise ProtocolError("at least one stage is required")
    for s in stages[:-1]:
        if s.kind == "gaussian":
            raise ProtocolError("a gaussian stage must be the last stage (it writes a log, not an xyz)")

    cluster = dict(raw.get("cluster", {}))
    missing = [k for k in CLUSTER_KEYS if k not in cluster]
    if missing:
        raise ProtocolError(f"cluster section is missing {missing}; see the example protocol")
    kinds = {s.kind for s in stages}
    need = {"uma": ("uma_env", "hf_token_file"), "gaussian": ("g16root",)}
    for kind, keys in need.items():
        if kind in kinds:
            for k in keys:
                if k not in cluster:
                    raise ProtocolError(f"a {kind} stage needs cluster.{k}")

    return Protocol(
        path=path, name=name, molecules=(path.parent / raw["molecules"]).resolve(),
        smiles_column=raw.get("smiles_column", "smiles"), name_column=raw.get("name_column"),
        id_column=raw.get("id_column"), charge=int(raw.get("charge", 0)),
        multiplicity=int(raw.get("multiplicity", 1)), build=build, stages=stages, cluster=cluster,
    )
