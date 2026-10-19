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
