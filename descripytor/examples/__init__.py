"""Bundled example molecules and extractor JSON.

After ``pip install descripytor``::

    from descripytor.examples import feather_example_dir, input_example_json
    molset = Molecules(str(feather_example_dir()), threshold=1.82)
"""

from __future__ import annotations

from pathlib import Path

_HERE = Path(__file__).resolve().parent
_REPO_FEATHER = (
    Path(__file__).resolve().parents[2]
    / "Getting_started_with_examples"
    / "feather_example"
)


def feather_example_dir() -> Path:
    """Directory of the 26 substituted-benzene ``.feather`` files."""
    bundled = _HERE / "feather_example"
    if (bundled / "basic.feather").is_file():
        return bundled
    if (_REPO_FEATHER / "basic.feather").is_file():
        return _REPO_FEATHER
    raise FileNotFoundError(
        "Example feathers not found. Expected basic.feather under "
        f"{bundled} or {_REPO_FEATHER}."
    )


def input_example_json() -> Path:
    """Extractor input JSON that matches the benzene example set."""
    for candidate in (
        _HERE / "input_example.json",
        feather_example_dir() / "input_example.json",
        _REPO_FEATHER.parent / "input_example.json",
    ):
        if candidate.is_file():
            return candidate
    raise FileNotFoundError("input_example.json not found")


def config_example_json() -> Path:
    """Toolkit ``run_config`` JSON example."""
    bundled = _HERE / "config_example.json"
    if bundled.is_file():
        return bundled
    clone = (
        Path(__file__).resolve().parents[2]
        / "Getting_started_with_examples"
        / "descriptor_extraction_toolkit"
        / "config_example.json"
    )
    if clone.is_file():
        return clone
    raise FileNotFoundError("config_example.json not found")
