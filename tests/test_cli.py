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


def test_help_functions_does_not_import_ipywidgets():
    """Extractor imports help_functions; ipywidgets is notebook-only."""
    import ast

    tree = ast.parse((ROOT / "utils" / "help_functions.py").read_text(encoding="utf-8"))
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and (node.module or "").split(".", 1)[0] == "ipywidgets":
            raise AssertionError("utils.help_functions must not import ipywidgets at module level")
        if isinstance(node, ast.Import):
            for alias in node.names:
                if alias.name.split(".", 1)[0] == "ipywidgets":
                    raise AssertionError("utils.help_functions must not import ipywidgets at module level")
