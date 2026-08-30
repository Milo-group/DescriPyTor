"""
Feature extraction layer for DescriPyTor.

The public API centers on:

- ``Molecule`` for one parsed molecule.
- ``Molecules`` for batch descriptor extraction.
- ``logs_to_feather`` for Gaussian log conversion.
- ``Molecules_xyz`` for standalone XYZ Sterimol extraction.
- ``MetalComplex`` / ``MetalComplexEnsemble`` for GOAT/xTB metal-complex features.
- ``XtbSinglePoint`` for xTB charge / WBO / dipole dumps.
- ``LigandTopology`` for SMILES graph indices.
- ``XYZEnsemble`` / ``GoatEnsemble`` for multi-structure XYZ files.

Exports are loaded lazily so importing lightweight utility modules does not
eagerly import pandas, PyArrow, RDKit, or other heavy scientific dependencies.
"""

from importlib import import_module

__all__ = [
    "Molecule",
    "Molecules",
    "Molecules_xyz",
    "logs_to_feather",
    "show_highly_correlated_pairs",
    "XYZEnsemble",
    "GoatEnsemble",
    "XtbSinglePoint",
    "MetalComplex",
    "MetalComplexEnsemble",
    "MetalComplexSet",
    "LigandTopology",
    "parse_props_file",
]


_LAZY = {
    "Molecule": (".data_extractor", "Molecule"),
    "Molecules": (".data_extractor", "Molecules"),
    "show_highly_correlated_pairs": (".data_extractor", "show_highly_correlated_pairs"),
    "logs_to_feather": (".feather_extractor", "logs_to_feather"),
    "Molecules_xyz": (".sterimol_standalone", "Molecules_xyz"),
    "XYZEnsemble": (".xyz_io", "XYZEnsemble"),
    "GoatEnsemble": (".xyz_io", "GoatEnsemble"),
    "XtbSinglePoint": (".xtb_singlepoint", "XtbSinglePoint"),
    "parse_props_file": (".xtb_singlepoint", "parse_props_file"),
    "MetalComplex": (".metal_complex", "MetalComplex"),
    "MetalComplexEnsemble": (".metal_complex", "MetalComplexEnsemble"),
    "MetalComplexSet": (".metal_complex", "MetalComplexSet"),
    "LigandTopology": (".ligand_topology", "LigandTopology"),
}


def __getattr__(name):
    if name not in _LAZY:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    module_name, attr = _LAZY[name]
    module = import_module(module_name, __name__)
    value = getattr(module, attr)
    globals()[name] = value
    return value
