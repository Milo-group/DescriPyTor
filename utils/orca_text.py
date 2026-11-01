"""Read the last geometry, frequencies, normal modes and IR intensities of an ORCA .out file."""
from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from typing import Optional

import numpy as np

IMAG_CUT = -1.0
BOHR = 0.529177210903


@dataclass
class ViewMol:
    name: str
    source: str
    symbols: list
    coordinates: np.ndarray
    frequencies: np.ndarray = field(default_factory=lambda: np.zeros(0))
    modes: Optional[np.ndarray] = None  # (3N, n_modes)
    intensities: Optional[np.ndarray] = None
    trajectory_xyz: str = ""
    note: str = ""


def _last_header(lines: list[str], title: str) -> int:
    hit = -1
    t = title.lower()
    for i, line in enumerate(lines):
        if t in line.lower():
            hit = i
    return hit


def _orca_last_xyz(lines: list[str]) -> tuple[list[str], np.ndarray]:
    idx = _last_header(lines, "CARTESIAN COORDINATES (ANGSTROEM)")
    if idx < 0:
        idx = _last_header(lines, "CARTESIAN COORDINATES (A.U.)")
        au = True
    else:
        au = "A.U." in lines[idx].upper()
    if idx < 0:
        raise ValueError("No CARTESIAN COORDINATES block in ORCA output")
    symbols, coords = [], []
    for line in lines[idx + 1:]:
        s = line.strip()
        if not s or s.startswith("-"):
            if symbols:
                break
            continue
        if s.lower().startswith("cartesian") or s.lower().startswith("internal"):
            break
        parts = s.split()
        if len(parts) < 4:
            if symbols:
                break
            continue
        try:
            xyz = [float(parts[1]), float(parts[2]), float(parts[3])]
        except ValueError:
            if symbols:
                break
            continue
        symbols.append(parts[0])
        coords.append(xyz)
    if not symbols:
        raise ValueError("Empty coordinate block in ORCA output")
    arr = np.asarray(coords, dtype=float)
    if au:
        arr = arr * BOHR
    return symbols, arr


def _orca_frequencies(lines: list[str]) -> tuple[np.ndarray, list[bool]]:
    idx = _last_header(lines, "VIBRATIONAL FREQUENCIES")
    if idx < 0:
        return np.zeros(0), []
    freqs, flags = [], []
    pat = re.compile(r"^\s*(\d+):\s+(-?\d+\.\d+)\s+cm")
    for line in lines[idx + 1:]:
        s = line.strip()
        if s.startswith("NORMAL") or s.startswith("IR SPECTRUM") or s.startswith("RAMAN"):
            break
        m = pat.match(line)
        if not m:
            continue
        v = float(m.group(2))
        imag = "imaginary" in line.lower()
        if imag and v > 0:
            v = -v
        freqs.append(v)
        flags.append(imag or v < IMAG_CUT)
    return np.asarray(freqs, dtype=float), flags


def _orca_normal_modes(lines: list[str], n_cart: int, n_modes: int) -> Optional[np.ndarray]:
    idx = _last_header(lines, "NORMAL MODES")
    if idx < 0 or n_modes <= 0:
        return None
    matrix = np.zeros((n_cart, n_modes), dtype=float)
    columns: list[int] = []
    started = False
    for line in lines[idx + 1:]:
        s = line.strip()
        if not s:
            continue
        if s.startswith("IR SPECTRUM") or s.startswith("RAMAN") or s.startswith("VIBRATIONAL"):
            break
        if s.startswith("-") or s.lower().startswith("these values"):
            continue
        parts = s.split()
        if not parts:
            continue
        header = all("." not in p and p.lstrip("+-").isdigit() for p in parts)
        if header:
            columns = [int(p) for p in parts]
            started = True
            continue
        if not started or not columns:
            continue
        if not parts[0].lstrip("+-").isdigit():
            break
        row = int(parts[0])
        vals = [float(p) for p in parts[1:]]
        if row >= n_cart or len(vals) != len(columns):
            continue
        for col, val in zip(columns, vals):
            if 0 <= col < n_modes:
                matrix[row, col] = val
    if not started:
        return None
    return matrix


def _orca_ir(lines: list[str], n_modes: int) -> np.ndarray:
    idx = _last_header(lines, "IR SPECTRUM")
    out = np.full(n_modes, math.nan)
    if idx < 0 or n_modes <= 0:
        return out
    pat = re.compile(r"^\s*(\d+):\s+")
    for line in lines[idx + 1:]:
        s = line.strip()
        if not s:
            if np.isfinite(out).any():
                break
            continue
        if s.startswith("RAMAN") or s.startswith("NORMAL") or s.startswith("THE"):
            break
        m = pat.match(line)
        if not m:
            continue
        parts = s.replace(":", " ").split()
        try:
            k = int(parts[0])
            # ORCA: Mode freq eps Int ...  Int is km/mol, typically 4th number
            nums = [float(p) for p in parts[1:] if re.match(r"^-?\d", p)]
            if k < n_modes and len(nums) >= 3:
                out[k] = nums[2]
            elif k < n_modes and len(nums) >= 1:
                out[k] = nums[-1]
        except (ValueError, IndexError):
            continue
    return out


def parse_orca_out_text(text: str, name: str = "orca") -> ViewMol:
    lines = text.splitlines()
    symbols, coords = _orca_last_xyz(lines)
    freqs, _flags = _orca_frequencies(lines)
    n_cart = 3 * len(symbols)
    n_modes = int(len(freqs)) if len(freqs) else n_cart
    modes = _orca_normal_modes(lines, n_cart, n_modes)
    intensities = _orca_ir(lines, n_modes) if n_modes else None
    note = "ORCA output"
    if freqs.size:
        n_imag = int((freqs < IMAG_CUT).sum())
        note += f" · {n_modes} frequencies · {n_imag} imaginary"
        if modes is None:
            note += " · no NORMAL MODES block (cannot animate)"
    else:
        note += " · no VIBRATIONAL FREQUENCIES (OptTS did not print a Hessian)"
    return ViewMol(
        name=name, source=name, symbols=symbols, coordinates=coords,
        frequencies=freqs, modes=modes, intensities=intensities, note=note,
    )
