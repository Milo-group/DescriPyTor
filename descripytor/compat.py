"""Reproduce the numbers of an earlier DescriPyTor release without checking it out.

    from descripytor.compat import paper_v3

    with paper_v3():
        mols = Molecules("feathers/")           # flat 1.82 A bonds
        mc = MetalComplex.from_xyz("cx.xyz")    # B1 scan started from the lab frame

Inside the block the package computes what tag ``paper-v3`` (the DescriPyTor paper) computes;
outside it, the current definitions apply. The settings are module-level, so the block covers
everything that runs inside it, and they are restored on exit even if an error is raised.

What ``paper_v3`` changes (see docs/STERIMOL_FIXES.md and docs/FEATURE_FIXES.md):

- bonds: the flat 1.82 A cutoff for non-metal pairs instead of covalent radii x 1.15
  (fix #9 / FEATURE_FIXES #1: S-CF3 and P-C bonds are dropped again);
- ``metal_complex.STERIMOL_FRAME = "lab"``: the B1 scan starts from the lab frame
  (STERIMOL_FIXES #5);
- ``metal_complex.STERIMOL_THETA_RULE = "scan"``: theta at the one B1 direction the scan lands
  on, not the soft average over tied directions (STERIMOL_FIXES #11);
- ``metal_complex.SUB_FRAGMENT_BOUND = "donor"``: the C*->R fragment is blocked only at the
  donor and the metal, so a fused substituent can run into the backbone (STERIMOL_FIXES #10);
- ``metal_complex.STEREOCENTRE_RULE = "first"``: among equally sized substituents the first in
  atom order is the arm's (C*, R), and a P-aryl arm is measured along ipso -> ortho;
- ``metal_complex.DIPOLE_ABS = False``: signed mu_desym / mu_outofplane (FEATURE_FIXES #12).

It does not undo the feature fixes that only removed wrong results (ring positions, ring angles,
missing stretch windows, tuple inputs); the paper does not use those features.
"""
from __future__ import annotations

import contextlib

__all__ = ["paper_v3"]


def _settings():
    """(module, attribute, paper-v3 value) for every module that can be imported here."""
    from utils import help_functions
    from M2_data_extractor import metal_complex
    out = [(help_functions, "DEFAULT_BOND_THRESHOLD", 1.82), (metal_complex, "STERIMOL_FRAME", "lab"),
           (metal_complex, "STERIMOL_THETA_RULE", "scan"), (metal_complex, "SUB_FRAGMENT_BOUND", "donor"),
           (metal_complex, "DIPOLE_ABS", False), (metal_complex, "STEREOCENTRE_RULE", "first")]
    try:                                    # needs morfeus; skip it where morfeus is absent
        from M2_data_extractor import sterimol_standalone
        out.append((sterimol_standalone, "DEFAULT_BOND_THRESHOLD", 1.82))
    except ImportError:
        pass
    return out


@contextlib.contextmanager
def paper_v3():
    """Run the enclosed code with the definitions of tag paper-v3."""
    settings = _settings()
    saved = [(module, name, getattr(module, name)) for module, name, _ in settings]
    try:
        for module, name, value in settings:
            setattr(module, name, value)
        yield
    finally:
        for module, name, value in saved:
            setattr(module, name, value)
