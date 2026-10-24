"""The feather writer stores polarizability as [iso, aniso] at positions 8 and 9 (feather_extractor's pol_df). Both positional
readers named position 8 'aniso' until 2026-10-18, so every table read through them had iso and aniso swapped (CS2's deposited
'aniso' column is the isotropic polarizability). These tests pin the order on an old-style feather with unnamed columns."""
import numpy as np
import pandas as pd

from M2_data_extractor import feather_extractor as FE
from M2_data_extractor import gaussian_handler as GH

ISO, ANISO = 180.521, 118.05                            # L15 H-bonded: the log prints iso 0.180521D+03, aniso 0.118050D+03


def old_style_frame(n=3):
    """unnamed columns as the pre-2026 writer left them: xyz 0-3, dipole 4-7, polarizability 8-9 (iso first), charges 10-12, freq/IR 13-14, one vibration 15-17"""
    rows = max(n, 2); d = pd.DataFrame(np.nan, index=range(rows), columns=range(18), dtype=object)
    for i in range(n): d.iloc[i, 0:4] = ["C", float(i), 0.0, 0.0]; d.iloc[i, 10:13] = [-0.1, -0.05, -0.08]; d.iloc[i, 15:18] = [0.1, 0.0, 0.0]
    d.iloc[0, 4:8] = [1.0, 0.0, 0.0, 1.0]; d.iloc[0, 8:10] = [ISO, ANISO]; d.iloc[0, 13:15] = [100.0, 5.0]
    return d


def test_df_file_handler_reads_iso_first():
    df_dict = FE.df_file_handler(old_style_frame())[0]
    pol = next(v for v in df_dict.values() if isinstance(v, pd.DataFrame) and "iso" in v.columns)
    assert np.isclose(float(pol["iso"].iloc[0]), ISO) and np.isclose(float(pol["aniso"].iloc[0]), ANISO)


def test_feather_file_handler_positional_branch_reads_iso_first(tmp_path):
    d = old_style_frame(); d.columns = [str(c) for c in d.columns]
    for c in d.columns: d[c] = d[c].astype(str).where(d[c].notna(), None)
    f = tmp_path / "old.feather"; d.to_feather(f)
    df_dict = GH.feather_file_handler(str(f))[0]
    pol = next(v for v in df_dict.values() if isinstance(v, pd.DataFrame) and "iso" in v.columns)
    assert np.isclose(float(pol["iso"].iloc[0]), ISO) and np.isclose(float(pol["aniso"].iloc[0]), ANISO)
