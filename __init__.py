"""
DescriPyTor package.

DescriPyTor provides tools for preparing molecular calculations, extracting
quantum-chemistry descriptors.

Most users enter through the command-line interface (``descripytor`` /
``python -m descripytor``) or through the high-level classes:

- ``M2_data_extractor.data_extractor.Molecule``
- ``M2_data_extractor.data_extractor.Molecules``
- ``M2_data_extractor.metal_complex.MetalComplex``
"""

from descripytor import __version__

__all__ = ["__version__"]
