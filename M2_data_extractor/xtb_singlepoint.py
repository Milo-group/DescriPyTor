"""Parse xTB single-point dumps (charges, Wiberg pairs, dipole, HOMO/LUMO).

Three on-disk layouts are supported, all of which appear in the Case Study 3
scratchpad:

1. Native xtb files: ``name.q`` (Mulliken charges), ``name.wbo`` (1-based
   ``i j value`` triples), optional ``name.dip``.
2. A props table (``props.txt``): ``name dx dy dz homo lumo`` with the dipole
   already in Debye.
3. The compact per-frame dump used for cheap/ensemble electronics::

       F <source> <frame_index>
       Q  q0 q1 ...
       W  j1 w1 j2 w2 ...     # metal-centric, 1-based partner index
       D  dx dy dz            # atomic units
       E  homo lumo
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass, field
import numpy as np

AU_TO_DEBYE = 2.541746


def read_charges(path) -> np.ndarray:
    """Read an xtb ``.q`` file (whitespace-separated Mulliken charges)."""
    return np.array([float(x) for x in open(path, encoding="utf-8").read().split()], dtype=float)


def read_wbo(path) -> dict:
    """Read an xtb ``.wbo`` file into ``{(i, j): value}`` with 0-based ``i < j``."""
    out = {}
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                i, j, value = int(parts[0]) - 1, int(parts[1]) - 1, float(parts[2])
            except ValueError:
                continue
            out[(min(i, j), max(i, j))] = value
    return out


def read_dip_file(path) -> np.ndarray:
    """Read an xtb-style ``.dip`` file (``full: dx dy dz |mu|``)."""
    text = open(path, encoding="utf-8").read()
    numbers = re.findall(r"[-+]?\d+\.\d+", text)
    if len(numbers) < 3:
        raise ValueError(f"Could not parse dipole components from {path}")
    return np.array([float(x) for x in numbers[:3]], dtype=float)


def parse_props_file(path) -> dict:
    """Parse ``name dx dy dz homo lumo`` rows.

    Returns
    -------
    dict
        ``{name: {'dipole': ndarray(3), 'homo': float, 'lumo': float}}``
        Dipole is Debye, matching the CS3 ``props*.txt`` files.
    """
    out = {}
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            parts = line.split()
            if len(parts) < 6:
                continue
            name = parts[0]
            out[name] = {
                "dipole": np.array([float(x) for x in parts[1:4]], dtype=float),
                "homo": float(parts[4]),
                "lumo": float(parts[5]),
            }
    return out


def _metal_wbo_from_compact(values: list[float]) -> dict:
    """Compact ``W j1 w1 j2 w2 ...`` lines are metal-centric (atom 0)."""
    out = {}
    for i in range(0, len(values) - 1, 2):
        partner = int(values[i]) - 1
        out[(min(0, partner), max(0, partner))] = float(values[i + 1])
    return out


@dataclass
class XtbSinglePoint:
    """One xTB single point: charges, Wiberg pairs, dipole, frontier orbitals.

    Parameters
    ----------
    charges : array-like
        Mulliken charges, one per atom, same order as the XYZ.
    wbo : dict
        ``{(i, j): order}`` with 0-based ``i < j``. Missing metal-donor pairs
        are treated as 0.0 by :meth:`MetalComplex.electronic_features` — that
        is the CS3 convention when a donor has swung off.
    dipole : array-like
        Cartesian dipole. See ``dipole_unit``.
    homo, lumo : float
        Orbital energies as printed by xtb (eV in the CS3 dumps).
    dipole_unit : {'debye', 'au'}
        Compact dumps store atomic units; ``props.txt`` and ``.dip`` store Debye.
    """

    charges: np.ndarray
    wbo: dict = field(default_factory=dict)
    dipole: np.ndarray | None = None
    homo: float | None = None
    lumo: float | None = None
    dipole_unit: str = "debye"
    source: str | None = None
    frame_index: int | None = None

    def __post_init__(self):
        self.charges = np.asarray(self.charges, dtype=float)
        if self.dipole is not None:
            self.dipole = np.asarray(self.dipole, dtype=float).reshape(3)

    @property
    def dipole_debye(self) -> np.ndarray:
        if self.dipole is None:
            raise ValueError("This single point has no dipole")
        if self.dipole_unit == "au":
            return self.dipole * AU_TO_DEBYE
        if self.dipole_unit == "debye":
            return self.dipole
        raise ValueError(f"Unknown dipole_unit {self.dipole_unit!r}")

    @property
    def gap(self) -> float:
        if self.homo is None or self.lumo is None:
            raise ValueError("This single point has no HOMO/LUMO")
        return self.lumo - self.homo

    @classmethod
    def from_files(
        cls,
        charges_path,
        wbo_path=None,
        dip_path=None,
        homo=None,
        lumo=None,
        dipole=None,
        dipole_unit: str = "debye",
    ) -> "XtbSinglePoint":
        """Load native xtb ``.q`` / ``.wbo`` / optional ``.dip``."""
        charges = read_charges(charges_path)
        wbo = read_wbo(wbo_path) if wbo_path and os.path.exists(wbo_path) else {}
        if dipole is None and dip_path and os.path.exists(dip_path):
            dipole = read_dip_file(dip_path)
            dipole_unit = "debye"
        return cls(
            charges=charges,
            wbo=wbo,
            dipole=dipole,
            homo=homo,
            lumo=lumo,
            dipole_unit=dipole_unit,
        )

    @classmethod
    def from_xtb_dir(cls, directory, stem, props=None) -> "XtbSinglePoint":
        """Load ``{stem}.q`` / ``.wbo`` / ``.dip`` from a charges directory.

        ``props`` is an optional dict with ``dipole``, ``homo``, ``lumo`` as
        produced by :func:`parse_props_file`.
        """
        props = props or {}
        q_path = os.path.join(directory, f"{stem}.q")
        wbo_path = os.path.join(directory, f"{stem}.wbo")
        dip_path = os.path.join(directory, f"{stem}.dip")
        return cls.from_files(
            charges_path=q_path,
            wbo_path=wbo_path if os.path.exists(wbo_path) else None,
            dip_path=dip_path if os.path.exists(dip_path) else None,
            homo=props.get("homo"),
            lumo=props.get("lumo"),
            dipole=props.get("dipole"),
            dipole_unit="debye",
        )

    @classmethod
    def parse_dump(cls, path) -> dict:
        """Parse a compact per-frame dump into ``{source: {frame: XtbSinglePoint}}``.

        Dipoles in this format are atomic units and are stored with
        ``dipole_unit='au'``.
        """
        out: dict[str, dict[int, XtbSinglePoint]] = {}
        source = None
        idx = None
        cur: dict = {}

        def _flush():
            if source is None or not cur:
                return
            if "q" not in cur:
                return
            sp = cls(
                charges=cur["q"],
                wbo=cur.get("w", {}),
                dipole=cur.get("d"),
                homo=None if "e" not in cur else cur["e"][0],
                lumo=None if "e" not in cur else cur["e"][1],
                dipole_unit="au",
                source=source,
                frame_index=idx,
            )
            out.setdefault(source, {})[idx] = sp

        with open(path, encoding="utf-8") as handle:
            for line in handle:
                parts = line.split()
                if not parts:
                    continue
                tag = parts[0]
                if tag == "F":
                    _flush()
                    source, idx, cur = parts[1], int(parts[2]), {}
                elif tag == "Q":
                    cur["q"] = np.array([float(x) for x in parts[1:]], dtype=float)
                elif tag == "W":
                    cur["w"] = _metal_wbo_from_compact([float(x) for x in parts[1:]])
                elif tag == "D":
                    cur["d"] = np.array([float(x) for x in parts[1:4]], dtype=float)
                elif tag == "E":
                    cur["e"] = (float(parts[1]), float(parts[2]))
        _flush()
        return out
