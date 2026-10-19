"""Molecule files of every kind the package reads: the folder scan keeps one file per name and skips anything else."""
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_xyz_and_feather_folders_load(tmp_path):
    from M2_data_extractor.data_extractor import Molecules, molecule_file_kind

    shutil.copy(ROOT / "tests" / "fixtures" / "072_lig.finalensemble.xyz", tmp_path / "a.xyz")
    shutil.copy(ROOT / "descripytor" / "examples" / "feather_example" / "basic.feather", tmp_path / "b.feather")
    shutil.copy(ROOT / "descripytor" / "examples" / "feather_example" / "basic.feather", tmp_path / "a.feather")  # beats a.xyz
    (tmp_path / "run.job.log").write_text("SGE job log\n")
    (tmp_path / "run_config.json").write_text("{}")
    assert molecule_file_kind(str(tmp_path / "run.job.log")) is None
    ms = Molecules(str(tmp_path))
    assert sorted(m.molecule_name for m in ms.molecules) == ["a", "b"] and not ms.failed_molecules
    assert all(len(m.xyz_df) == 29 for m in ms.molecules)            # a.feather, not a.xyz (27 atoms)
    alone = tmp_path / "x"
    alone.mkdir()
    shutil.copy(ROOT / "tests" / "fixtures" / "072_lig.finalensemble.xyz", alone)
    m = Molecules(str(alone)).molecules[0]
    assert len(m.xyz_df) == 27 and all(df.empty for df in m.charge_dict.values())


ORCA = """                                 * O   R   C   A *
---------------------------------
CARTESIAN COORDINATES (ANGSTROEM)
---------------------------------
  O      0.000000    0.000000    0.117300
  H      0.000000    0.757200   -0.469200
  H      0.000000   -0.757200   -0.469200

-----------------------
MULLIKEN ATOMIC CHARGES
-----------------------
   0 O :   -0.640000
   1 H :    0.320000
   2 H :    0.320000
Sum of atomic charges:   -0.0000000

----------------------
LOEWDIN ATOMIC CHARGES
----------------------
   0 O :   -0.420000
   1 H :    0.210000
   2 H :    0.210000

FINAL SINGLE POINT ENERGY       -76.400000000000
                        -----------------------------------------
Total Dipole Moment    :      0.000000000       0.000000000       0.780000000
                        -----------------------------------------
Magnitude (a.u.)       :      0.780000000
Magnitude (Debye)      :      1.982561880
"""


def test_an_orca_output_reads_like_a_feather(tmp_path):
    from M2_data_extractor.data_extractor import Molecule, folder_inputs

    (tmp_path / "water.out").write_text(ORCA)
    m = Molecule(str(tmp_path / "water.out"))
    assert m.xyz_df["atom"].tolist() == ["O", "H", "H"]
    assert m.charge_dict["mulliken"]["charge"].tolist() == [-0.64, 0.32, 0.32]
    assert abs(m.charge_dict["loewdin"]["charge"].iloc[0] + 0.42) < 1e-9 and m.charge_dict["nbo"].empty
    d = folder_inputs(str(tmp_path))
    assert d["kinds"] == {"orca_out": 1}
    assert d["shared"] == {"dipole": True, "polarizability": False, "vibrations": False, "charges": ["loewdin", "mulliken"]}


def test_renumbered_xyz_stays_xyz(tmp_path):
    from M2_data_extractor.alignment import write_renumbered
    from M2_data_extractor.data_extractor import Molecules

    shutil.copy(ROOT / "tests" / "fixtures" / "072_lig.finalensemble.xyz", tmp_path / "a.xyz")
    mols = Molecules(str(tmp_path))
    n = len(mols.molecules[0].xyz_df)
    old_to_new = {i: i for i in range(1, n + 1)}
    old_to_new[1], old_to_new[2] = 2, 1
    write_renumbered(mols, {"a": {"old_to_new": old_to_new}}, str(tmp_path / "out"), "a")
    lines = (tmp_path / "out" / "a.xyz").read_text().splitlines()
    first = mols.molecules[0].xyz_df.iloc[0]
    assert lines[0] == str(n) and lines[3].split()[0] == first["atom"] and abs(float(lines[3].split()[1]) - float(first["x"])) < 1e-6
