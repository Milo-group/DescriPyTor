"""Bundled example feathers and JSON ship with the package."""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from descripytor.examples import (
    config_example_json,
    feather_example_dir,
    input_example_json,
)


def test_twenty_six_feathers():
    directory = feather_example_dir()
    feathers = sorted(directory.glob("*.feather"))
    assert len(feathers) == 26
    assert (directory / "basic.feather").is_file()


def test_extractor_jsons_exist():
    assert input_example_json().is_file()
    assert config_example_json().is_file()
    assert "Sterimol" in input_example_json().read_text(encoding="utf-8")
