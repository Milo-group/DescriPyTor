"""Multi-structure XYZ ensembles and Boltzmann weights.

This module does not import morfeus. It is the geometry-file layer used by
:class:`~M2_data_extractor.metal_complex.MetalComplex` and
:class:`~M2_data_extractor.metal_complex.MetalComplexEnsemble`.

Energy conventions
------------------
CREST writes the GFN2 energy as the **first** float on the comment line.
ORCA GOAT often writes an RMSD first and the energy last; the Case Study 3
extractors therefore take the **last** float. Using the first float on a
GOAT file that starts with an RMSD silently Boltzmann-weights the wrong
number.

The GOAT files used for the Corminboeuf tables happen to start with the
energy (``-72.91 converged=true``), so last-float still recovers it because
``true`` is not a float. Last-float is the safe default for GOAT.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np

HARTREE_TO_KCAL = 627.5095          # matches the CS3 extract_cb.py constant
GAS_CONSTANT_KCAL = 0.0019872041    # kcal / (mol K)


def energy_from_comment(comment: str, convention: str = "last") -> float:
    """Parse a Hartree energy from an XYZ comment line.

    Parameters
    ----------
    comment : str
        Line 2 of an XYZ block.
    convention : {'last', 'first', 'auto'}
        ``last`` is GOAT-safe. ``first`` matches CREST. ``auto`` takes the last
        float whose absolute value is at least 1 Eh when one exists, otherwise
        the last float (CREST relative zeros stay ~0).
    """
    if comment is None:
        return float("nan")
    tokens = comment.replace("=", " ").split()
    floats: list[float] = []
    for tok in tokens:
        try:
            floats.append(float(tok))
        except ValueError:
            continue
    if not floats:
        return float("nan")
    conv = (convention or "last").lower()
    if conv == "first":
        return floats[0]
    if conv == "last":
        return floats[-1]
    if conv == "auto":
        abs_eh = [v for v in floats if abs(v) >= 1.0]
        return abs_eh[-1] if abs_eh else floats[-1]
    raise ValueError(f"Unknown energy convention {convention!r}")


def boltzmann_weights(energies_hartree, temperature: float = 298.15) -> np.ndarray:
    """Normalized Boltzmann weights from absolute energies in Hartree.

    Uses the CS3 constants (``HARTREE_TO_KCAL = 627.5095``) so ensemble averages
    match the published Case Study 3 tables. NaN energies get a NaN weight and
    are left out of the normalization.
    """
    energies = np.asarray(energies_hartree, dtype=float)
    weights = np.full(energies.shape, np.nan)
    valid = ~np.isnan(energies)
    if not valid.any():
        return weights
    rel_kcal = (energies[valid] - energies[valid].min()) * HARTREE_TO_KCAL
    exponent = -rel_kcal / (GAS_CONSTANT_KCAL * temperature)
    exponent = exponent - exponent.max()
    w = np.exp(exponent)
    w = w / w.sum()
    weights[valid] = w
    return weights


def boltzmann_average(
    rows: Sequence[dict],
    energies_hartree,
    temperature: float = 298.15,
) -> dict:
    """Boltzmann-average numeric dicts. Columns with any NaN are dropped.

    A single row is returned unchanged (no energy required), matching
    ``confsets.average``.
    """
    if not rows:
        return {}
    if len(rows) == 1:
        return dict(rows[0])
    cols = sorted(set().union(*(r.keys() for r in rows)))
    values = np.array([[r.get(c, np.nan) for c in cols] for r in rows], dtype=float)
    weights = boltzmann_weights(energies_hartree, temperature=temperature)
    valid = ~np.isnan(weights)
    if not valid.any():
        return {}
    w = weights[valid]
    w = w / w.sum()
    keep = ~np.isnan(values[valid]).any(axis=0)
    out = {}
    for j, col in enumerate(cols):
        if not keep[j]:
            continue
        out[col] = float((values[valid, j] * w).sum())
    return out


@dataclass
class XYZFrame:
    """One structure from a (possibly multi-block) XYZ file."""

    symbols: list
    coords: np.ndarray
    energy: float = float("nan")
    comment: str = ""
    index: int = 0

    def __post_init__(self):
        self.coords = np.asarray(self.coords, dtype=float)
        if self.coords.ndim != 2 or self.coords.shape[1] != 3:
            raise ValueError(f"coords must be (n_atoms, 3), got {self.coords.shape}")
        if len(self.symbols) != len(self.coords):
            raise ValueError("symbols and coords length mismatch")


def parse_xyz_ensemble(
    filepath,
    energy_convention: str = "last",
    relative_zero_threshold: float | None = None,
) -> list[XYZFrame]:
    """Parse a single- or multi-structure XYZ file into :class:`XYZFrame` objects.

    Parameters
    ----------
    filepath : str
        Path to the XYZ file.
    energy_convention : {'last', 'first', 'auto'}
        See :func:`energy_from_comment`.
    relative_zero_threshold : float or None
        If set, and every parsed energy is missing or has absolute value below
        this cutoff (CREST's ``0.00000000`` relative dump), every frame energy
        is set to NaN so the ensemble is not pooled with absolute GOAT energies.
    """
    with open(filepath, encoding="utf-8") as handle:
        lines = handle.read().splitlines()
    while lines and not lines[-1].strip():
        lines.pop()

    frames: list[XYZFrame] = []
    i = 0
    n_lines = len(lines)
    while i < n_lines:
        header = lines[i].strip()
        if not header:
            i += 1
            continue
        try:
            natoms = int(header.split()[0])
        except (ValueError, IndexError) as exc:
            raise ValueError(
                f"Malformed XYZ block header at line {i + 1} in {filepath}: {lines[i]!r}"
            ) from exc
        comment = lines[i + 1] if i + 1 < n_lines else ""
        energy = energy_from_comment(comment, convention=energy_convention)
        atom_lines = lines[i + 2:i + 2 + natoms]
        symbols = []
        rows = []
        for line in atom_lines:
            parts = line.split()
            if len(parts) < 4:
                continue
            symbols.append(parts[0])
            rows.append([float(parts[1]), float(parts[2]), float(parts[3])])
        if len(rows) != natoms:
            raise ValueError(
                f"Block starting at line {i + 1} in {filepath} declares {natoms} atoms "
                f"but {len(rows)} coordinate line(s) were parsed."
            )
        frames.append(XYZFrame(
            symbols=symbols,
            coords=np.asarray(rows, dtype=float),
            energy=float(energy),
            comment=comment.strip(),
            index=len(frames),
        ))
        i += 2 + natoms

    if not frames:
        raise ValueError(f"No structures parsed from {filepath}")

    if relative_zero_threshold is not None:
        energies = [f.energy for f in frames]
        if all(np.isnan(e) or abs(e) < relative_zero_threshold for e in energies):
            for frame in frames:
                frame.energy = float("nan")
    return frames


class XYZEnsemble:
    """Parsed multi-XYZ file with Boltzmann weights.

    Parameters
    ----------
    filepath : str
        Ensemble or single-structure XYZ.
    energy_convention : {'last', 'first', 'auto'}
        Default ``last`` (GOAT-safe). Pass ``first`` for CREST.
    temperature : float
        Kelvin, default 298.15.
    relative_zero_threshold : float or None
        See :func:`parse_xyz_ensemble`. Use ``1.0`` to match the CS3 CREST
        pooling rule.
    molecule_name : str, optional
        Defaults to the file stem with ensemble suffixes stripped.
    """

    ENSEMBLE_SUFFIXES = (
        ".finalensemble.xyz",
        ".ens.xyz",
        "_conformers.xyz",
        "_rotamers.xyz",
        "_best.xyz",
        ".xyz",
    )

    def __init__(
        self,
        filepath,
        energy_convention: str = "last",
        temperature: float = 298.15,
        relative_zero_threshold: float | None = None,
        molecule_name: str | None = None,
    ):
        self.filepath = filepath
        self.energy_convention = energy_convention
        self.temperature = temperature
        self.frames = parse_xyz_ensemble(
            filepath,
            energy_convention=energy_convention,
            relative_zero_threshold=relative_zero_threshold,
        )
        self.molecule_name = molecule_name or self._name_from_path(filepath)
        self.energies_hartree = np.array([f.energy for f in self.frames], dtype=float)
        self.weights = boltzmann_weights(self.energies_hartree, temperature=temperature)

    @staticmethod
    def _name_from_path(filepath) -> str:
        import os

        name = os.path.basename(filepath)
        lower = name.lower()
        for suffix in XYZEnsemble.ENSEMBLE_SUFFIXES:
            if lower.endswith(suffix):
                return name[: -len(suffix)] or name
        return os.path.splitext(name)[0]

    def __len__(self) -> int:
        return len(self.frames)

    @property
    def n_conformers(self) -> int:
        return len(self.frames)

    @property
    def lowest_index(self) -> int:
        energies = self.energies_hartree.copy()
        energies[np.isnan(energies)] = np.inf
        return int(np.argmin(energies))

    def lowest_frame(self) -> XYZFrame:
        return self.frames[self.lowest_index]

    def average(self, rows: Sequence[dict]) -> dict:
        """Boltzmann-average descriptor dicts aligned with ``self.frames``."""
        if len(rows) != len(self.frames):
            raise ValueError("rows must line up 1:1 with ensemble frames")
        keep_rows, keep_e = [], []
        for row, energy in zip(rows, self.energies_hartree):
            if not row:
                continue
            keep_rows.append(row)
            keep_e.append(energy)
        return boltzmann_average(keep_rows, keep_e, temperature=self.temperature)

    @classmethod
    def from_file(cls, filepath, **kwargs) -> "XYZEnsemble":
        return cls(filepath, **kwargs)


class GoatEnsemble(XYZEnsemble):
    """GOAT multi-XYZ with last-float comment energies.

    For metal-referenced CS3 descriptors (``fromM_*``, ``q_anc_sum``, …) use
    :class:`~M2_data_extractor.metal_complex.MetalComplexEnsemble`. This class
    is the file/energy layer only.
    """

    def __init__(self, filepath, energy_convention: str = "last", **kwargs):
        super().__init__(filepath, energy_convention=energy_convention, **kwargs)
