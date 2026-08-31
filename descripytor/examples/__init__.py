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
_REPO_BAPTISTE = (
    Path(__file__).resolve().parents[2]
    / "Getting_started_with_examples"
    / "baptiste_products"
)


def _first_existing(*candidates: Path, marker: str) -> Path:
    for path in candidates:
        if (path / marker).is_file():
            return path
    raise FileNotFoundError(
        f"{marker} not found under: " + ", ".join(str(p) for p in candidates)
    )


def feather_example_dir() -> Path:
    """Directory of the 26 substituted-benzene ``.feather`` files."""
    return _first_existing(
        _HERE / "feather_example",
        _REPO_FEATHER,
        marker="basic.feather",
    )


def baptiste_example_dir() -> Path:
    """Directory of the Baptiste product ``.feather`` files (paper case, not the GUI default)."""
    return _first_existing(
        _HERE / "baptiste_products",
        _REPO_BAPTISTE,
        marker="unsub.feather",
    )


def input_example_json() -> Path:
    """Extractor input JSON that ships with the package."""
    for candidate in (
        _HERE / "input_example.json",
        Path(__file__).resolve().parents[2]
        / "Getting_started_with_examples"
        / "input_example.json",
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
