"""descripytor.compat.paper_v3: main reproduces the paper inside the block, and only there."""
import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from descripytor.compat import paper_v3
from M2_data_extractor import metal_complex as mc
from utils.help_functions import extract_connectivity

ROOT = Path(__file__).resolve().parents[1]
FIX = ROOT / "M2_data_extractor" / "theta_explorer" / "fixtures"
PAPER = Path(os.environ["DESCRIPYTOR_PAPER"]) if os.environ.get("DESCRIPYTOR_PAPER") else None


def _xyz(path):
    lines = Path(path).read_text().splitlines()
    n = int(lines[0].split()[0])
    df = pd.DataFrame([l.split()[:4] for l in lines[2:2 + n]], columns=["atom", "x", "y", "z"])
    df[["x", "y", "z"]] = df[["x", "y", "z"]].astype(float)
    return df


def _has_s_cf3():
    return (34, 37) in {tuple(sorted(r)) for r in extract_connectivity(_xyz(FIX / "p-Otf.xyz")).values.tolist()}


def test_switches_on_and_back_off():
    assert _has_s_cf3() and mc.STERIMOL_FRAME == "fragment"
    with paper_v3():
        assert not _has_s_cf3() and mc.STERIMOL_FRAME == "lab"
    assert _has_s_cf3() and mc.STERIMOL_FRAME == "fragment"


def test_restores_after_an_error():
    with pytest.raises(RuntimeError):
        with paper_v3():
            raise RuntimeError
    assert _has_s_cf3() and mc.STERIMOL_FRAME == "fragment"


@pytest.mark.skipif(PAPER is None, reason="set DESCRIPYTOR_PAPER to the paper repository")
def test_cs3_sterimol_block_is_the_deposited_one():
    """single_structure/cp.csv, 16 Sterimol/theta columns, from the geometries they came from."""
    ref = pd.read_csv(PAPER / "zenodo_archive" / "feature_matrices" / "single_structure" / "cp.csv").set_index("name")
    cols = [f"{fr}_{k}_{p}" for fr in ("fromM", "sub") for k in ("B1", "B5", "L", "theta") for p in ("sym", "asym")]
    with paper_v3():
        got = pd.DataFrame({n: mc.MetalComplex.from_xyz(str(PAPER / "notes" / "analysis" / "single_structure_rebuilt" / "cu1" / f"{n}.xyz")).geometric_features()
                            for n in ref.index}).T
    assert float((got[cols].astype(float) - ref[cols]).abs().max().max()) < 1e-9
