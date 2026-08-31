"""Bundled example feathers and JSON ship with the package."""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from descripytor.examples import (
    baptiste_example_dir,
    config_example_json,
    feather_example_dir,
    input_example_json,
)


def test_baptiste_example_set():
    directory = baptiste_example_dir()
    feathers = sorted(directory.glob("*.feather"))
    assert len(feathers) == 18
    assert (directory / "unsub.feather").is_file()
    outcomes = directory / "outcomes.csv"
    assert outcomes.is_file()
    text = outcomes.read_text(encoding="utf-8")
    assert "unsub" in text
    assert "output" in text.splitlines()[0]
    presets = directory / "presets.json"
    assert presets.is_file()
    assert "sterimol" in presets.read_text(encoding="utf-8")
    names = {p.stem for p in feathers}
    outcome_names = {
        line.split(",", 1)[0].strip()
        for line in outcomes.read_text(encoding="utf-8").splitlines()[1:]
        if line.strip()
    }
    assert names == outcome_names


def test_benzene_example_set():
    directory = feather_example_dir()
    feathers = sorted(directory.glob("*.feather"))
    assert len(feathers) == 10
    assert (directory / "basic.feather").is_file()
    presets = directory / "presets.json"
    assert presets.is_file()
    text = presets.read_text(encoding="utf-8")
    assert "sterimol" in text
    assert "basic.feather" in text
    assert directory != baptiste_example_dir()


def test_extractor_jsons_exist():
    assert input_example_json().is_file()
    assert config_example_json().is_file()
    assert "Sterimol" in input_example_json().read_text(encoding="utf-8")
