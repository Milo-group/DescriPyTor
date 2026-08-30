"""CLI surface a first-time user actually types."""
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _run(*args: str) -> subprocess.CompletedProcess[str]:
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join([str(ROOT), env.get("PYTHONPATH", "")])
    return subprocess.run(
        [sys.executable, "-m", "descripytor", *args],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
        env=env,
    )


def test_descripytor_help_lists_visual():
    result = _run("--help")
    assert result.returncode == 0, result.stderr
    assert "visual" in result.stdout


def test_visual_help():
    result = _run("visual", "--help")
    assert result.returncode == 0, result.stderr
    assert "--no-browser" in result.stdout
