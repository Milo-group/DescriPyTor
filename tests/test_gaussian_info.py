"""process_gaussian_info: the IR column is 'IR Inten', not 'Frc consts' (Gaussian prints both per block)."""
from M2_data_extractor.feather_extractor import process_gaussian_info, process_gaussian_vibs_string

LOG = """
 Harmonic frequencies (cm**-1), IR intensities (KM/Mole), Raman scattering
                      1                      2                      3
                      A                      A                      A
 Frequencies --    -31.1985                46.4283               100.5681
 Red. masses --      1.9973                 6.8825                 3.8991
 Frc consts  --      0.0011                 0.0087                 0.0232
 IR Inten    --      0.5408                 0.1564                 2.6504
  Atom  AN      X      Y      Z        X      Y      Z        X      Y      Z
     1   6     0.00   0.00   0.10     0.00   0.10   0.00     0.10   0.00   0.00
                      4
                      A
 Frequencies --   3193.2336
 Red. masses --      1.0892
 Frc consts  --      6.5432
 IR Inten    --     12.3456
  Atom  AN      X      Y      Z
     1   6     0.00   0.00   0.10
 -------------------
 - Thermochemistry -
"""


def test_ir_column_is_intensity():
    df = process_gaussian_info(process_gaussian_vibs_string(LOG))
    assert df["Frequency"].tolist() == [-31.1985, 46.4283, 100.5681, 3193.2336]
    assert df["IR"].tolist() == [0.5408, 0.1564, 2.6504, 12.3456]
