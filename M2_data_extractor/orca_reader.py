"""ORCA output (.out) -> the structured molecule JSON the package reads (utils.help_functions.json_payload_handler).

What it carries, from the last block of each kind in the file: the final geometry, the energy (FINAL SINGLE POINT
ENERGY), the dipole (Total Dipole Moment, converted to Debye as the Gaussian feathers store it), every charge
population ORCA printed (Mulliken, Loewdin, Hirshfeld), the isotropic polarizability, and the vibrations of a
frequency job (normal modes, wavenumbers, IR intensities; the six zero modes left out, as Gaussian prints none).
NBO and CM5 charges are Gaussian's: an ORCA file has none, so those descriptors are not offered for it.

    from M2_data_extractor.orca_reader import read_orca, is_orca_output
    payload = read_orca("job.out")      # {"molecule": {"atoms": [...]}, "properties": {...}, "charges": {...}, ...}
"""
from __future__ import annotations

import re

import numpy as np

AU_TO_DEBYE = 2.541746
_POPULATIONS = (("MULLIKEN ATOMIC CHARGES", "mulliken"), ("LOEWDIN ATOMIC CHARGES", "loewdin"))


def is_orca_output(path, head=None):
    """An ORCA output, told by its banner (an .out can be anything)."""
    if head is None:
        try:
            with open(path, "r", encoding="utf-8", errors="replace") as fh:
                head = fh.read(20000)
        except OSError:
            return False
    return "O   R   C   A" in head


def _last(lines, title):
    hits = [i for i, ln in enumerate(lines) if ln.strip().startswith(title)]
    return hits[-1] if hits else None


def _population(lines, title, n):
    """Per-atom charges of the last `title` block ("   0 C :   -0.198087" lines), or None."""
    i = _last(lines, title)
    if i is None:
        return None
    out = []
    for ln in lines[i + 2:i + 2 + n]:
        m = re.match(r"\s*\d+\s+[A-Za-z]+\s*:\s*(-?\d+\.\d+)", ln)
        if not m:
            return None
        out.append(float(m.group(1)))
    return out if len(out) == n else None


def _hirshfeld(lines, n):
    i = _last(lines, "HIRSHFELD ANALYSIS")
    if i is None:
        return None
    for j in range(i, min(i + 40, len(lines))):
        if lines[j].split()[:2] == ["ATOM", "CHARGE"]:
            vals = [float(ln.split()[2]) for ln in lines[j + 1:j + 1 + n] if len(ln.split()) >= 3]
            return vals if len(vals) == n else None
    return None


def _number(lines, title, field=-1):
    i = _last(lines, title)
    if i is None:
        return None
    try:
        return float(lines[i].split()[field])
    except (ValueError, IndexError):
        return None


def read_orca(path):
    """The payload json_payload_handler reads; raises ValueError when the file has no geometry."""
    from utils.mode_viewer import parse_orca_out_text

    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        text = fh.read()
    if not is_orca_output(path, text[:20000]):
        raise ValueError(f"{path} is not an ORCA output")
    view = parse_orca_out_text(text, str(path))
    if not view.symbols:
        raise ValueError(f"{path}: no CARTESIAN COORDINATES (ANGSTROEM) block")
    lines = text.splitlines()
    n = len(view.symbols)
    atoms = [{"atom": s, "x": float(x), "y": float(y), "z": float(z)} for s, (x, y, z) in zip(view.symbols, view.coordinates)]

    props = {}
    energy = _number(lines, "FINAL SINGLE POINT ENERGY")
    if energy is not None:
        props["energies"] = {"final_single_point_energy_hartree_orca": energy}
    i = _last(lines, "Total Dipole Moment")
    if i is not None:
        d = [float(v) * AU_TO_DEBYE for v in lines[i].split(":")[1].split()[:3]]
        props["dipole"] = {"x": d[0], "y": d[1], "z": d[2], "total": float(np.linalg.norm(d))}
    iso = _number(lines, "Isotropic polarizability")
    if iso is not None:
        props["polarizability"] = {"iso": iso, "aniso": None}

    charges = {key: q for title, key in _POPULATIONS if (q := _population(lines, title, n)) is not None}
    if (h := _hirshfeld(lines, n)) is not None:
        charges["hirshfeld"] = h

    ir = {}                                       # IR SPECTRUM rows are "mode: freq eps Int ...", by ORCA's mode number
    i = _last(lines, "IR SPECTRUM")
    if i is not None:
        for ln in lines[i + 1:i + 6 + 3 * n]:
            m = re.match(r"\s*(\d+):\s+-?\d+\.\d+\s+\S+\s+(\d+\.\d+)", ln)
            if m:
                ir[int(m.group(1))] = float(m.group(2))
    modes = []
    if view.frequencies.size and view.modes is not None:
        for k, f in enumerate(view.frequencies):
            if abs(f) < 1e-3:                     # translations and rotations
                continue
            disp = view.modes[:, k].reshape(n, 3)
            modes.append({"frequency_cm-1": float(f), "ir_intensity_km_per_mol": ir.get(k, 0.0),
                          "displacements": [{"dx": float(a), "dy": float(b), "dz": float(c)} for a, b, c in disp]})
    return {"molecule": {"atoms": atoms}, "properties": props, "charges": charges,
            "vibrations": {"modes": modes}, "source": {"program": "ORCA", "file": str(path)}}
