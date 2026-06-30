"""
Feature extraction layer for MolFeatures.

The public API centers on:

- ``Molecule`` for one parsed molecule.
- ``Molecules`` for batch descriptor extraction.
- ``logs_to_feather`` for Gaussian log conversion.
- ``Molecules_xyz`` for standalone XYZ Sterimol extraction.

Exports are loaded lazily so importing lightweight utility modules does not
eagerly import pandas, PyArrow, RDKit, or other heavy scientific dependencies.
"""

__all__ = [
    "Molecule",
    "Molecules",
    "Molecules_xyz",
    "logs_to_feather",
    "show_highly_correlated_pairs",
]


def __getattr__(name):
    if name in {"Molecule", "Molecules", "show_highly_correlated_pairs"}:
        from .data_extractor import Molecule, Molecules, show_highly_correlated_pairs

        exports = {
            "Molecule": Molecule,
            "Molecules": Molecules,
            "show_highly_correlated_pairs": show_highly_correlated_pairs,
        }
        return exports[name]

    if name == "logs_to_feather":
        from .feather_extractor import logs_to_feather

        return logs_to_feather

    if name == "Molecules_xyz":
        from .sterimol_standalone import Molecules_xyz

        return Molecules_xyz

    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
