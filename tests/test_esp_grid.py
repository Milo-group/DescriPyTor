"""The ESP sampling grid must never place a point inside an atom.

esp.f accumulates the nuclear term as z(j)/(1.0d-14 + r), so a point sitting on
a nucleus returns ~1e14 Hartree rather than an error. One such point makes the
normalised angular moments meaningless for that ligand, silently. The offset
support-function construction rules it out by geometry, and this pins that.

Unlike test_sterimol_loc_b1, these import normally: esp_grid_utils keeps the
GeneralConstants import inside get_radii precisely so the geometry does not drag
in rdkit and plotly.
"""
import numpy as np
import pandas as pd
import pytest

from M2_data_extractor.extractor_utils.esp_grid_utils import (
    ANGSTROM_TO_BOHR,
    angular_moments,
    build_cylindrical_grid,
    cross_section_support,
    read_xtb_esp,
    write_esp_coord,
)

BONDI = {'P': 1.80, 'C': 1.70, 'H': 1.10, 'Cl': 1.75}


@pytest.fixture
def ligand():
    """Narrow near the origin, one bulky arm far out along the axis, so no
    single cylinder radius can stay outside the arm and outside the neck."""
    return pd.DataFrame({
        'atom': ['P', 'C', 'C', 'C', 'C', 'C', 'H', 'H', 'H', 'Cl'],
        'x': [0.0, 0.0, 1.4, 2.1, 1.4, 0.0, -0.9, 3.2, 1.9, -2.6],
        'y': [0.0, 1.5, 2.2, 3.4, 4.6, 4.0, 1.9, 3.6, 5.5, 4.2],
        'z': [0.0, 0.2, 0.6, 0.9, 0.4, -0.3, 0.1, 1.2, 0.7, -1.1],
    })


def radii_of(df):
    return df['atom'].map(BONDI).to_numpy(float)


def min_clearance(points, df):
    """Smallest (distance to nucleus - vdW radius). Negative means inside."""
    centres = df[['x', 'y', 'z']].to_numpy(float)
    d = np.linalg.norm(points[:, None, :] - centres[None, :, :], axis=2)
    return (d - radii_of(df)[None, :]).min()


def test_no_grid_point_lands_inside_an_atom(ligand):
    grid = build_cylindrical_grid(ligand, n_phi=36, n_axial=20,
                                  probe_offset=0.5, radii=radii_of(ligand))
    pts = grid[['x', 'y', 'z']].to_numpy(float)
    assert min_clearance(pts, ligand) > 0


def test_fixed_radius_shell_does_land_inside_atoms(ligand):
    """The alternative sampling definition, kept as a test so the reason the
    grid follows the surface stays visible rather than becoming folklore."""
    levels = np.linspace(-1.0, 5.5, 20)
    phi = np.linspace(0, 2 * np.pi, 36, endpoint=False)
    yy, pp = np.meshgrid(levels, phi, indexing='ij')
    for rho0 in (2.0, 3.0, 4.0):
        shell = np.stack([rho0 * np.cos(pp).ravel(), yy.ravel(),
                          rho0 * np.sin(pp).ravel()], axis=1)
        assert min_clearance(shell, ligand) < 0


def test_probe_offset_must_be_positive(ligand):
    with pytest.raises(ValueError):
        build_cylindrical_grid(ligand, probe_offset=0.0, radii=radii_of(ligand))


def test_support_function_reaches_b5_in_the_plane_of_the_b5_atom(ligand):
    """max over phi of h is Sterimol B5, so the grid radius and the Sterimol
    widths come from one construction rather than two.

    The cross-section at height y sees atom i with the sliced radius
    sqrt(r_i**2 - (y - y_i)**2), which equals r_i only in the plane through the
    atom's own centre. So B5 is recovered by evaluating at the atom centres, not
    on an arbitrary level grid -- a grid that misses y_i underestimates it.
    """
    phi = np.linspace(0, 2 * np.pi, 720, endpoint=False)
    r = radii_of(ligand)
    h = np.nanmax([cross_section_support(ligand, y, phi, radii=r)
                   for y in ligand['y'].to_numpy(float)], axis=0)
    assert 0 < h.min() < h.max()

    xz = np.linalg.norm(ligand[['x', 'z']].to_numpy(float), axis=1)
    assert h.max() == pytest.approx((xz + r).max(), abs=1e-3)


def test_esp_coord_is_written_in_bohr(ligand, tmp_path):
    grid = build_cylindrical_grid(ligand, n_phi=12, n_axial=6,
                                  radii=radii_of(ligand))
    path = tmp_path / 'esp_coord'
    n = write_esp_coord(grid, str(path))
    assert n == len(grid)
    written = np.loadtxt(path)
    assert np.allclose(written,
                       grid[['x', 'y', 'z']].to_numpy(float) * ANGSTROM_TO_BOHR)


def test_truncated_esp_file_is_rejected(ligand, tmp_path):
    """A stale esp_coord makes xtb sample a different grid than the one asked
    for. That has to raise, not return numbers for the wrong molecule."""
    grid = build_cylindrical_grid(ligand, n_phi=12, n_axial=6,
                                  radii=radii_of(ligand))
    bohr = grid[['x', 'y', 'z']].to_numpy(float) * ANGSTROM_TO_BOHR
    dat = tmp_path / 'xtb_esp.dat'
    np.savetxt(dat, np.column_stack([bohr, np.zeros(len(bohr))])[:-3])
    with pytest.raises(ValueError, match='stale esp_coord'):
        read_xtb_esp(grid, str(dat))


def test_shifted_coordinates_are_rejected(ligand, tmp_path):
    """Right number of rows, wrong points -- the case a row count alone misses."""
    grid = build_cylindrical_grid(ligand, n_phi=12, n_axial=6,
                                  radii=radii_of(ligand))
    bohr = grid[['x', 'y', 'z']].to_numpy(float) * ANGSTROM_TO_BOHR
    dat = tmp_path / 'xtb_esp.dat'
    np.savetxt(dat, np.column_stack([bohr + 0.1, np.zeros(len(bohr))]))
    with pytest.raises(ValueError, match='do not match'):
        read_xtb_esp(grid, str(dat))


def test_isotropic_field_gives_vanishing_moment(ligand, tmp_path):
    """The property the moments exist for: no angular structure, no direction.
    An argmin would still return some angle here; |Z| reports that there is
    nothing to report."""
    grid = build_cylindrical_grid(ligand, n_phi=72, n_axial=6,
                                  radii=radii_of(ligand))
    # Constant potential, and a circular cross-section, so both fields are flat.
    grid = grid.assign(rho=1.0)
    bohr = grid[['x', 'y', 'z']].to_numpy(float) * ANGSTROM_TO_BOHR
    dat = tmp_path / 'xtb_esp.dat'
    np.savetxt(dat, np.column_stack([bohr, np.zeros(len(bohr))]))

    sampled = read_xtb_esp(grid, str(dat))
    moments = angular_moments(sampled)
    assert moments['Z_S_1'].abs().max() < 1e-8
    assert (moments['A_1'].abs() < 1e-8).all()
