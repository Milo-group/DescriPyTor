"""
DescriPyTor package.

DescriPyTor provides tools for preparing molecular calculations, extracting
quantum-chemistry descriptors, and running molecular modeling workflows.

Most users enter through the command-line interface (``descripytor`` /
``python -m descripytor``) or through the high-level classes:

- ``M2_data_extractor.data_extractor.Molecule``
- ``M2_data_extractor.data_extractor.Molecules``
- ``M2_data_extractor.metal_complex.MetalComplex``
- ``M3_modeler.modeling.LinearRegressionModel``
- ``M3_modeler.modeling.ClassificationModel``
"""

from descripytor import __version__

__all__ = ["__version__"]
