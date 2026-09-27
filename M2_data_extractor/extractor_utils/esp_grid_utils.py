"""
Cylindrical ESP sampling grids for xtb --esp.

The grid follows the molecular surface rather than sitting at a fixed radius.
That is not a stylistic choice: esp.f accumulates the nuclear term as

    espe(i) = espe(i) + z(j)/(1.0d-14 + r)

so a grid point that lands inside an atom returns ~1e14 Hartree instead of an
error, and a single such point destroys any normalised angular moment computed
from the field. A fixed-radius shell punches through atoms for any ligand whose
bulk varies along the axis; the surface-following construction below cannot.

Frame convention matches sterimol_utils: y is the Sterimol axis (L), and x-z is
the transverse plane, so azimuth is measured as atan2(z, x). Feed this the same
dataframe preform_coordination_transformation produces.
"""

import numpy as np
import pandas as pd

try:
    from numpy.typing import ArrayLike, NDArray
except ImportError:  # numpy < 1.20, as on the cluster's psi4conda build
    ArrayLike = NDArray = np.ndarray

ANGSTROM_TO_BOHR = 1.8897261246257702

# Clearance between the vdW surface and the sampling shell. 0.5 A keeps the
# probe in the region a substrate contact would occupy without pushing it so
# far out that the angular contrast washes away.
DEFAULT_PROBE_OFFSET = 0.5


def get_radii(coordinates_df: pd.DataFrame, radii: str = 'bondi') -> NDArray:
    """
    Van der Waals radii per atom, keyed on element symbol.

    Bondi is the default here rather than CPK because CPK needs the perceived
    atom types from a bonds dataframe, and the sampling shell does not need
    that resolution -- only a surface to stand off from.

    GeneralConstants is imported here rather than at module scope so the
    geometry above stays importable without pulling in the sterimol_utils
    import chain, which reaches plotly, dash and rdkit.
    """
    from .sterimol_utils import GeneralConstants

    radii_map = (GeneralConstants.BONDI_RADII.value if radii == 'bondi'
                 else GeneralConstants.CPK_RADII.value)
    missing = set(coordinates_df['atom']) - set(radii_map)
    if missing:
        raise KeyError(f'no {radii} radius for: {sorted(missing)}')
    return coordinates_df['atom'].map(radii_map).to_numpy(dtype=float)


def _resolve_radii(coordinates_df: pd.DataFrame, radii) -> NDArray:
    """Accept either a named radii table or an explicit per-atom array."""
    if isinstance(radii, str):
        return get_radii(coordinates_df, radii)
    return np.asarray(radii, dtype=float)


def cross_section_support(coordinates_df: pd.DataFrame, y_level: float,
                          phi: ArrayLike, radii='bondi') -> NDArray:
    """
    Support function h(y, phi) of the molecular cross-section at one axial level.

    Atom i cuts the plane at height y in a disc of radius
    sqrt(r_i**2 - (y - y_i)**2), so the support function of the union of those
    discs is

        h(y, phi) = max_i [ (x_i, z_i) . u_phi + r_i(y) ]

    taken over the atoms the plane actually intersects. Evaluating the Sterimol
    support function this way rather than over the whole molecule is what makes
    the grid track the surface along the axis.

    `radii` is either the name of a table or an explicit per-atom array.
    Returns NaN where no atom reaches the level, which the caller skips.
    """
    xz = coordinates_df[['x', 'z']].to_numpy(dtype=float)
    y = coordinates_df['y'].to_numpy(dtype=float)
    r = _resolve_radii(coordinates_df, radii)

    # Discs of the atoms this plane cuts through.
    dy = y_level - y
    cuts = np.abs(dy) < r
    if not cuts.any():
        return np.full(np.shape(phi), np.nan)

    r_slice = np.sqrt(r[cuts] ** 2 - dy[cuts] ** 2)
    u = np.stack([np.cos(phi), np.sin(phi)], axis=-1)

    # (n_atoms, n_phi): projection of each disc centre plus its own radius.
    reach = xz[cuts] @ u.T + r_slice[:, None]
    return reach.max(axis=0)


def build_cylindrical_grid(coordinates_df: pd.DataFrame, n_phi: int = 36,
                           n_axial: int = 20, probe_offset: float = DEFAULT_PROBE_OFFSET,
                           axial_range: tuple = None, radii='bondi') -> pd.DataFrame:
    """
    Surface-following cylindrical sampling grid.

    Each point sits at radius h(y, phi) + probe_offset along direction phi, at
    axial level y. Because h is a support function, every point of the
    cross-section satisfies p.u_phi <= h, so a point at h + offset lies strictly
    outside the supporting half-plane and therefore outside every atom in that
    plane. No point can be singular for offset > 0.

    The construction is convex-hull-following, so it stands off the outer
    envelope and does not descend into concave grooves between substituents.
    That is the conservative choice for ESP sampling -- a groove-following grid
    would have to re-admit the singularity risk this exists to remove.

    Returns a frame of x, y, z in Angstrom, with the (axial, phi) indices kept
    so the sampled potential can be folded back into V(y, phi).
    """
    if probe_offset <= 0:
        raise ValueError('probe_offset must be positive, or points may fall inside atoms')

    y = coordinates_df['y'].to_numpy(dtype=float)
    r = _resolve_radii(coordinates_df, radii)
    if axial_range is None:
        axial_range = (float((y - r).min()), float((y + r).max()))

    # Half-step inset: the extreme levels are tangent planes where the
    # cross-section degenerates to a point and the support function is noise.
    edges = np.linspace(*axial_range, n_axial + 1)
    levels = 0.5 * (edges[:-1] + edges[1:])
    phi = np.linspace(0.0, 2.0 * np.pi, n_phi, endpoint=False)

    rows = []
    for i_axial, y_level in enumerate(levels):
        h = cross_section_support(coordinates_df, y_level, phi, radii)
        if np.isnan(h).all():
            continue
        rho = h + probe_offset
        rows.append(pd.DataFrame({
            'i_axial': i_axial,
            'i_phi': np.arange(n_phi),
            'y_level': y_level,
            'phi': phi,
            'rho': rho,
            'x': rho * np.cos(phi),
            'y': y_level,
            'z': rho * np.sin(phi),
        }))

    if not rows:
        raise ValueError('no axial level intersects the molecule')
    return pd.concat(rows, ignore_index=True)


def write_esp_coord(grid_df: pd.DataFrame, path: str = 'esp_coord') -> int:
    """
    Write the grid in the format esp.f expects: xyz triples in BOHR.

    The reader tries five columns (coordinates plus two surfac weights) and
    falls back to three, so a bare triple file is accepted. Returns the point
    count, which the caller should check against xtb_esp.dat.
    """
    xyz_bohr = grid_df[['x', 'y', 'z']].to_numpy(dtype=float) * ANGSTROM_TO_BOHR
    np.savetxt(path, xyz_bohr, fmt='%20.12E')
    return len(xyz_bohr)


def read_xtb_esp(grid_df: pd.DataFrame, path: str = 'xtb_esp.dat') -> pd.DataFrame:
    """
    Attach the sampled potential to the grid it was requested on.

    xtb_esp.dat is written as '(3E18.10,F14.8)': coordinates in Bohr, potential
    in Hartree, in the order the grid file supplied them. The row count is
    checked rather than assumed -- a mismatch means xtb read a different grid
    than the one written, which is the failure mode a stale esp_coord produces.
    """
    data = np.loadtxt(path)
    if len(data) != len(grid_df):
        raise ValueError(
            f'{path} has {len(data)} points but the grid has {len(grid_df)}; '
            'xtb likely read a stale esp_coord from a previous molecule'
        )

    out = grid_df.copy()
    out['esp'] = data[:, 3]

    # Coordinates round-trip through the file, so disagreement beyond write
    # precision means the rows are not the points that were asked for.
    written = out[['x', 'y', 'z']].to_numpy(dtype=float) * ANGSTROM_TO_BOHR
    if not np.allclose(data[:, :3], written, atol=1e-6):
        raise ValueError(f'{path} coordinates do not match the written grid')
    return out


def angular_moments(esp_df: pd.DataFrame, orders=(1, 2), beta: float = 1.0) -> pd.DataFrame:
    """
    Circular moments of the steric and electronic angular fields, per axial level.

        Z^(m) = sum_phi f(phi) exp(i m phi) / sum_phi f(phi)

    Steric openness is weighted exp(-beta * h), so open directions dominate;
    electronic favourability uses the potential shifted to be non-negative,
    since a signed weight would not give a normalisable circular distribution.

    |Z| is the point of using moments at all: a near-isotropic field returns
    |Z| ~ 0 and contributes nothing, instead of an argmin picking an arbitrary
    direction out of shallow competing minima.
    """
    records = []
    for i_axial, level in esp_df.groupby('i_axial'):
        phi = level['phi'].to_numpy(dtype=float)
        openness = np.exp(-beta * (level['rho'].to_numpy(dtype=float)
                                   - level['rho'].min()))
        favour = level['esp'].max() - level['esp'].to_numpy(dtype=float)

        row = {'i_axial': i_axial, 'y_level': level['y_level'].iloc[0]}
        for m in orders:
            phase = np.exp(1j * m * phi)
            z_s = (openness * phase).sum() / openness.sum()
            z_e = (favour * phase).sum() / favour.sum() if favour.sum() > 0 else 0j
            row[f'Z_S_{m}'] = z_s
            row[f'Z_E_{m}'] = z_e
            row[f'A_{m}'] = float(np.real(z_e * np.conj(z_s)))
            row[f'H_{m}'] = float(np.imag(z_e * np.conj(z_s)))
        records.append(row)

    return pd.DataFrame(records)
