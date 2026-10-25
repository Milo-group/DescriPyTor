import numpy as np
import pandas as pd 
import numpy.typing as npt
import matplotlib.pyplot as plt
plt.ion() 

from utils.help_functions import *
try:
    from extractor_utils.bond_graph import (
        diameter_vertex_path,
        graph_from_bonds_df,
    )
except ImportError:
    from bond_graph import (
        diameter_vertex_path,
        graph_from_bonds_df,
    )


class GeneralConstants(Enum):
    """
    Holds constants for calculations and conversions
    1. covalent radii from Alvarez (2008) DOI: 10.1039/b801115j
    2. atomic numbers
    2. atomic weights
    """
    
    PYYKKO_RADII= {
            'H': 0.31, 'He': 0.28, 'Li': 1.28,
            'Be': 0.96, 'B': 0.84, 'C': 0.76, 
            'N': 0.71, 'O': 0.66, 'F': 0.57, 'Ne': 0.58,
            'Na': 1.66, 'Mg': 1.41, 'Al': 1.21, 'Si': 1.11, 
            'P': 1.07, 'S': 1.05, 'Cl': 1.02, 'Ar': 1.06,
            'K': 2.03, 'Ca': 1.76, 'Sc': 1.70, 'Ti': 1.60, 
            'V': 1.53, 'Cr': 1.39, 'Mn': 1.61, 'Fe': 1.52, 
            'Co': 1.50, 'Ni': 1.24, 'Cu': 1.32, 'Zn': 1.22, 
            'Ga': 1.22, 'Ge': 1.20, 'As': 1.19, 'Se': 1.20, 
            'Br': 1.20, 'Kr': 1.16, 'Rb': 2.20, 'Sr': 1.95,
            'Y': 1.90, 'Zr': 1.75, 'Nb': 1.64, 'Mo': 1.54,
            'Tc': 1.47, 'Ru': 1.46, 'Rh': 1.42, 'Pd': 1.39,
            'Ag': 1.45, 'Cd': 1.44, 'In': 1.42, 'Sn': 1.39,
            'Sb': 1.39, 'Te': 1.38, 'I': 1.39, 'Xe': 1.40,
            'Cs': 2.44, 'Ba': 2.15, 'La': 2.07, 'Ce': 2.04,
            'Pr': 2.03, 'Nd': 2.01, 'Pm': 1.99, 'Sm': 1.98,
            'Eu': 1.98, 'Gd': 1.96, 'Tb': 1.94, 'Dy': 1.92,
            'Ho': 1.92, 'Er': 1.89, 'Tm': 1.90, 'Yb': 1.87,
            'Lu': 1.87, 'Hf': 1.75, 'Ta': 1.70, 'W': 1.62,
            'Re': 1.51, 'Os': 1.44, 'Ir': 1.41, 'Pt': 1.36,
            'Au': 1.36, 'Hg': 1.32, 'Tl': 1.45, 'Pb': 1.46,
            'Bi': 1.48, 'Po': 1.40, 'At': 1.50, 'Rn': 1.50, 
            'Fr': 2.60, 'Ra': 2.21, 'Ac': 2.15, 'Th': 2.06,
            'Pa': 2.00, 'U': 1.96, 'Np': 1.90, 'Pu': 1.87,
            'Am': 1.80, 'Cm': 1.69
    }
    
    BONDI_RADII={
        'H': 1.10, 'C': 1.70, 'F': 1.47,
        'S': 1.80, 'B': 1.92, 'I': 1.98,
        'N': 1.55, 'O': 1.52,
        'Br': 1.83, 'Si': 2.10,
        'P': 1.80, 'Cl': 1.75,
        # alkali / alkaline-earth (Alvarez 2013, DOI: 10.1039/c3dt50599e)
        'Li': 1.82, 'Be': 1.53, 'Na': 2.27, 'Mg': 1.73,
        'K': 2.75, 'Ca': 2.31, 'Rb': 3.03, 'Sr': 2.49,
        'Cs': 3.43, 'Ba': 2.68,
        # p-block metals / heavier metalloids
        'Al': 1.84, 'Ga': 1.87, 'Ge': 2.11, 'As': 1.85,
        'In': 1.93, 'Sn': 2.17, 'Sb': 2.06, 'Te': 2.06,
        'Tl': 1.96, 'Pb': 2.02, 'Bi': 2.07,
        # first-row transition metals
        'Sc': 2.18, 'Ti': 2.11, 'V': 2.07, 'Cr': 2.06,
        'Mn': 2.05, 'Fe': 2.04, 'Co': 2.00, 'Ni': 1.97,
        'Cu': 1.96, 'Zn': 2.01,
        # second-row transition metals
        'Y': 2.32, 'Zr': 2.23, 'Nb': 2.18, 'Mo': 2.17,
        'Tc': 2.16, 'Ru': 2.13, 'Rh': 2.10, 'Pd': 2.10,
        'Ag': 2.11, 'Cd': 2.18,
        # third-row transition metals
        'Hf': 2.23, 'Ta': 2.22, 'W': 2.18, 'Re': 2.16,
        'Os': 2.16, 'Ir': 2.13, 'Pt': 2.09, 'Au': 2.14,
        'Hg': 2.23,
        # lanthanides
        'La': 2.43, 'Ce': 2.42, 'Pr': 2.40, 'Nd': 2.39,
        'Pm': 2.38, 'Sm': 2.36, 'Eu': 2.35, 'Gd': 2.34,
        'Tb': 2.33, 'Dy': 2.31, 'Ho': 2.30, 'Er': 2.29,
        'Tm': 2.27, 'Yb': 2.26, 'Lu': 2.24,
        # actinides
        'Ac': 2.47, 'Th': 2.45, 'Pa': 2.43, 'U': 2.41,
        'Np': 2.39, 'Pu': 2.37, 'Am': 2.35, 'Cm': 2.35,
    }
    CPK_RADII = {
    'C': 1.50,
    'C3': 1.60,
    'C6/N6': 1.70,
    'H': 1.00,
    'N': 1.50,
    'N4': 1.45,
    'O': 1.35,
    'O2': 1.35,
    'P': 1.40,
    'S': 1.70,
    'S1': 1.00,
    'F': 1.35,
    'Cl': 1.80,
    'S4': 1.40,
    'Br': 1.95,
    'I': 2.15,
    'X': 1.92,
    # metals — element symbols used directly (Alvarez 2013 VdW values)
    'Li': 1.82, 'Be': 1.53, 'Na': 2.27, 'Mg': 1.73,
    'K': 2.75, 'Ca': 2.31, 'Rb': 3.03, 'Sr': 2.49,
    'Cs': 3.43, 'Ba': 2.68,
    'Al': 1.84, 'Ga': 1.87, 'Ge': 2.11, 'In': 1.93,
    'Sn': 2.17, 'Tl': 1.96, 'Pb': 2.02, 'Bi': 2.07,
    'Sc': 2.18, 'Ti': 2.11, 'V': 2.07, 'Cr': 2.06,
    'Mn': 2.05, 'Fe': 2.04, 'Co': 2.00, 'Ni': 1.97,
    'Cu': 1.96, 'Zn': 2.01,
    'Y': 2.32, 'Zr': 2.23, 'Nb': 2.18, 'Mo': 2.17,
    'Tc': 2.16, 'Ru': 2.13, 'Rh': 2.10, 'Pd': 2.10,
    'Ag': 2.11, 'Cd': 2.18,
    'Hf': 2.23, 'Ta': 2.22, 'W': 2.18, 'Re': 2.16,
    'Os': 2.16, 'Ir': 2.13, 'Pt': 2.09, 'Au': 2.14,
    'Hg': 2.23,
    'La': 2.43, 'Ce': 2.42, 'Pr': 2.40, 'Nd': 2.39,
    'Pm': 2.38, 'Sm': 2.36, 'Eu': 2.35, 'Gd': 2.34,
    'Tb': 2.33, 'Dy': 2.31, 'Ho': 2.30, 'Er': 2.29,
    'Tm': 2.27, 'Yb': 2.26, 'Lu': 2.24,
    'Ac': 2.47, 'Th': 2.45, 'Pa': 2.43, 'U': 2.41,
    'Np': 2.39, 'Pu': 2.37, 'Am': 2.35, 'Cm': 2.35,
}
    

# Points sampled per atom circle when building the xz plane. plot_b1_visualization
# slices the plane back into circles with this number, so the two must agree.
STERIMOL_CIRCLE_POINTS = 100


def generate_circle(center_x, center_y, radius, n_points=20):
    """
    Generate circle coordinates given a center and radius.
    Returns a DataFrame with columns 'x' and 'y'.
    """
    theta = np.linspace(0, 2 * np.pi, n_points)
    x = center_x + radius * np.cos(theta)
    y = center_y + radius * np.sin(theta)
    return np.column_stack((x, y))



    
def calc_new_base_atoms(coordinates_array: npt.ArrayLike, atom_indices: npt.ArrayLike):  #help function for calc_coordinates_transformation
    """
    a function that calculates the new base atoms for the transformation of the coordinates.
    optional: if the atom_indices[0] is list, compute the new origin as the middle of the first atoms.
    """
   
    if isinstance(atom_indices[0], list):
        new_origin=np.mean(coordinates_array[atom_indices[0]], axis=0)
    else:
        new_origin=coordinates_array[atom_indices[0]]
    new_y=(coordinates_array[atom_indices[1]]-new_origin)/np.linalg.norm((coordinates_array[atom_indices[1]]-new_origin))
    coplane=((coordinates_array[atom_indices[2]]-new_origin)/np.linalg.norm((coordinates_array[atom_indices[2]]-new_origin)+0.00000001))
    
    return (new_origin,new_y,coplane)

def np_cross_and_vstack(plane_1, plane_2):
    cross_plane=np.cross(plane_1, plane_2)
    united_results=np.vstack([plane_1, plane_2, cross_plane])
    return united_results

from numpy.typing import ArrayLike

def calc_basis_vector(origin, y: ArrayLike, coplane: ArrayLike):
    """
    Calculate the new basis vector.
    
    Parameters
    ----------
    origin : array-like
        The origin of the new basis.
    y : array-like
        The new basis's y direction.
    coplane : array-like
        A vector coplanar with the new y direction.
    
    Returns
    -------
    new_basis : np.array
        The computed new basis matrix.
    """
    coef_mat = np_cross_and_vstack(coplane, y)
    angle_new_y_coplane = calc_angle(coplane, y)
    cop_ang_x = angle_new_y_coplane - (np.pi/2)
    result_vector = [np.cos(cop_ang_x), 0, 0]
    new_x, _, _, _ = np.linalg.lstsq(coef_mat, result_vector, rcond=None)
    new_basis = np_cross_and_vstack(new_x, y)
    return new_basis

def transform_row(row_array, new_basis, new_origin, round_digits):
    """
    Transform a single row of coordinates.
    
    Parameters
    ----------
    row_array : array-like
        A row of coordinates.
    new_basis : array-like
        New basis matrix.
    new_origin : array-like
        The new origin.
    round_digits : int
        Number of decimal places to round to.
    
    Returns
    -------
    np.array
        The transformed coordinates.
    """
    row_array = np.squeeze(row_array)
    translocated_row = row_array - new_origin
    ## make sure or change to shape (3,) if needed
    if translocated_row.shape != (3,):
        try:
            translocated_row = np.reshape(translocated_row, (3,))
        except Exception as e:
            print("transform_row: Error reshaping translocated_row:", e)
    
    result = np.dot(new_basis, translocated_row).round(round_digits)
    return result

def calc_coordinates_transformation(coordinates_array: ArrayLike, base_atoms_indices: ArrayLike, round_digits: int = 4, origin: ArrayLike = None) -> ArrayLike:
    """
    Transform a coordinates array to a new basis defined by the atoms in base_atoms_indices.
    
    Parameters
    ----------
    coordinates_array : np.array
        The original xyz molecule coordinates.
    base_atoms_indices : list or array-like
        Indices of atoms to define the new basis.
    round_digits : int, optional
        Number of digits to round the result.
    origin : array-like, optional
        A new origin for the transformation. If None, the first base atom is used.
    
    Returns
    -------
    transformed_coordinates : np.array
        The transformed coordinates.
    """
    # indices = adjust_indices(base_atoms_indices)
    # Calculate new basis using helper function calc_new_base_atoms (assumed to return origin, y, and coplane)
    new_basis = calc_basis_vector(*calc_new_base_atoms(coordinates_array, base_atoms_indices))
    if origin is None:
        new_origin = coordinates_array[base_atoms_indices[0]]
    else:
        new_origin = origin
   
    transformed_coordinates = np.apply_along_axis(lambda x: transform_row(x, new_basis, new_origin, round_digits), 1, coordinates_array)
  
    return transformed_coordinates


def get_center_of_mass(xyz_df):
    coordinates=np.array(xyz_df[['x','y','z']].values,dtype=float)
    atoms_symbol=np.array(xyz_df['atom'].values)
    masses=np.array([GeneralConstants.ATOMIC_WEIGHTS.value[symbol] for symbol in atoms_symbol])
    center_of_mass=np.sum(coordinates*masses[:,None],axis=0)/np.sum(masses)
    return center_of_mass

def get_closest_atom_to_center(xyz_df,center_of_mass):
    distances = np.sqrt((xyz_df['x'] - center_of_mass[0]) ** 2 + (xyz_df['y'] - center_of_mass[1]) ** 2 + (xyz_df['z'] - center_of_mass[2]) ** 2)
    idx_closest = np.argmin(distances)
    center_atom = xyz_df.loc[idx_closest]
    return center_atom

def get_sterimol_base_atoms(center_atom, bonds_df):
    
    center_atom_id=int(center_atom.name)+1
    base_atoms = [center_atom_id]
    if (any(bonds_df[0] == center_atom_id)):
        base_atoms.append(int(bonds_df[(bonds_df[0]==center_atom_id)][1].iloc[0]))
    else:
        base_atoms.append(int(bonds_df[(bonds_df[1]==center_atom_id)][0].iloc[0]))
    
    return base_atoms

def center_substructure(coordinates_array,atom_indices):
    atom_indices=adjust_indices(atom_indices)
    substructure=coordinates_array[atom_indices]
    center_substructure=np.mean(substructure,axis=0)
    return center_substructure



def get_sterimol_indices(coordinates,bonds_df):
    center=get_center_of_mass(coordinates)
    center_atom=get_closest_atom_to_center(coordinates,center)
    base_atoms=get_sterimol_base_atoms(center_atom,bonds_df)
    return base_atoms

def filter_atoms_for_sterimol(bonded_atoms_df,coordinates_df):
    """
   
    """

    allowed_bonds_indices= pd.concat([bonded_atoms_df['index_1'],bonded_atoms_df['index_2']],axis=1).reset_index(drop=True)
    atom_filter=adjust_indices(np.unique([atom for sublist in allowed_bonds_indices.values.tolist() for atom in sublist]))
    # save current index as a list
    index_list=coordinates_df.loc[atom_filter].index.tolist()
    edited_coordinates_df=coordinates_df.loc[atom_filter].reset_index(drop=True)

    return edited_coordinates_df, index_list


def get_extended_df_for_sterimol(coordinates_df, bonds_df, radii='CPK'):
    radii_map = GeneralConstants.CPK_RADII.value if radii == 'CPK' else GeneralConstants.BONDI_RADII.value
    df = coordinates_df.copy()

    if radii == 'bondi':
        df['atype'] = df['atom']
    else:
        df['atype'] = nob_atype(coordinates_df, bonds_df)

    df['magnitude'] = np.linalg.norm(df[['x', 'z']].astype(float), axis=1)
    df['radius'] = df['atype'].map(radii_map)
    df['B5'] = df['radius'] + df['magnitude']
    df['L'] = df['y'] + df['radius']
    df['loc_B5'] = df['y']

    return df

def get_transfomed_plane_for_sterimol(plane, degree):
    """
    Rotate a 2D plane by a given angle in degrees.
    plane: np.ndarray of shape (n_points, 2)
    """
    theta = np.deg2rad(degree)
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)

    rot_matrix = np.array([
        [cos_theta, -sin_theta],
        [sin_theta,  cos_theta]
    ])

    return plane @ rot_matrix.T


def b1s_for_loop_function(degree_list, plane):
    """
    Scan rotated planes over a list of angles and compute B1-related metrics.

    Every angle at once: the planes are rotated in one batched product (bit for bit the
    per-angle ``plane @ R.T``), and each angle's extents, B5 point and angles are read off
    along the angle axis.

    Parameters
    ----------
    degree_list : iterable of float
        Angles in degrees to scan.
    plane : np.ndarray, shape (n_points, 2)
        Original 2D steric point cloud.

    Returns
    -------
    pd.DataFrame
        Columns: degree, B1, B1_B5_angle, b1_coords, b5_value, plane, b1_normal.
    """
    degrees = list(degree_list)
    if not degrees:
        return pd.DataFrame(columns=["degree", "B1", "B1_B5_angle", "b1_coords", "b5_value", "plane", "b1_normal"])
    theta = np.deg2rad(np.asarray(degrees, dtype=float))
    c, s = np.cos(theta), np.sin(theta)
    rot_t = np.stack([np.stack([c, s], axis=1), np.stack([-s, c], axis=1)], axis=1)     # R.T for each angle
    planes = np.matmul(plane[None, :, :], rot_t)                                        # (angles, points, 2)

    max_x, min_x = planes[:, :, 0].max(axis=1), planes[:, :, 0].min(axis=1)
    max_y, min_y = planes[:, :, 1].max(axis=1), planes[:, :, 1].min(axis=1)
    extents = np.abs(np.stack([max_x, min_x, max_y, min_y], axis=1))
    min_index = np.argmin(extents, axis=1)
    b1_value = extents[np.arange(len(degrees)), min_index]
    on_x = min_index < 2
    b1_x = np.where(min_index == 0, max_x, np.where(min_index == 1, min_x, 0.0))
    b1_y = np.where(min_index == 2, max_y, np.where(min_index == 3, min_y, 0.0))

    norms_sq = np.sum(planes ** 2, axis=2)
    b5_point = planes[np.arange(len(degrees)), np.argmax(norms_sq, axis=1)]
    b5_value = np.linalg.norm(b5_point, axis=1)

    angle_b1 = np.arctan2(b1_y, b1_x) % (2 * np.pi)
    angle_b5 = np.arctan2(b5_point[:, 1], b5_point[:, 0]) % (2 * np.pi)
    angle_diff = np.abs(angle_b5 - angle_b1)
    angle_diff = np.where(angle_diff > np.pi, 2 * np.pi - angle_diff, angle_diff)

    # B1 plane normal in the original (pre-rotation) xz frame.
    # The scan rotates the plane, t = R p, so a rotated axis is the original
    # direction R^T applied to it: rotated x is [cos θ, -sin θ], rotated z is
    # [sin θ, cos θ]. (R x̂ = [cos θ, sin θ] is the mirror image across x, and
    # makes B1_B5_angle depend on how the frame is turned about the axis.)
    # B1 minimum falls along the rotated x-axis (min_index 0/1) or z-axis (2/3).
    normals = np.where(on_x[:, None], np.stack([c, -s], axis=1), np.stack([s, c], axis=1))

    return pd.DataFrame({
        "degree": degrees,
        "B1": b1_value,
        "B1_B5_angle": np.degrees(angle_diff),
        "b1_coords": list(zip(b1_x.tolist(), b1_y.tolist())),
        "b5_value": b5_value,
        "plane": list(planes),
        "b1_normal": list(normals),
    })

def preform_coordination_transformation(xyz_df, indices=None, origin=None):
    """
    Perform a coordination transformation on the xyz DataFrame.

    Parameters
    ----------
    xyz_df : pd.DataFrame
        DataFrame containing columns 'x', 'y', 'z'.
    indices : array-like, optional
        Atom indices defining the new basis. If None, [1,2,3] is used.
    origin : array-like, optional
        Atom indices whose centroid becomes the new origin. If None, the first
        entry of `indices` is used, which is calc_coordinates_transformation's
        default.

        Note: this argument used to take the centroid of `indices` rather than
        of `origin`, so it acted as a flag ("centre on the base atoms") instead
        of naming the atoms to centre on. Callers that relied on that are
        responsible for passing the atoms they actually want.

    Returns
    -------
    xyz_copy : pd.DataFrame
        DataFrame with transformed coordinates.
    """
    xyz_copy = xyz_df.copy()
    coordinates = np.array(xyz_copy[['x', 'y', 'z']].values)

    if origin is not None and len(np.atleast_1d(origin)) > 0:
        origin_idx = [int(i) for i in np.atleast_1d(origin)]
        geometric_origin = np.mean(xyz_df.iloc[origin_idx][['x', 'y', 'z']].values, axis=0)
    else:
        geometric_origin = None

    basis_indices = [1, 2, 3] if indices is None else indices
    transformed = calc_coordinates_transformation(coordinates, basis_indices,
                                                  origin=geometric_origin)

    xyz_copy[['x', 'y', 'z']] = transformed
    return xyz_copy



def direction_atoms_for_sterimol(bonds_df, base_atoms) -> list:
    """
    Return [origin, direction, third] for Sterimol coordinate setup.
    - bonds_df: DataFrame with two columns (0,1), undirected bonds between integer atom indices.
    - base_atoms: [origin, direction] or [origin, direction, third].
       If third == origin, choose a neighbor of 'direction' different from origin.
    Raises ValueError with a clear message if it can't determine 'third'.
    """
    import pandas as pd

    if not isinstance(bonds_df, pd.DataFrame) or not set(bonds_df.columns[:2]) >= {0, 1}:
        raise ValueError("bonds_df must be a DataFrame with columns 0 and 1 listing bonds.")
    
    # if base atoms list of list do base_atoms=base_atoms[0]
    if isinstance(base_atoms[0], list):
        base_atoms=base_atoms[0]

    if len(base_atoms) < 2:
        raise ValueError(f"base_atoms must contain at least two indices [origin, direction]; got: {base_atoms}")

    # Normalize to ints
    origin, direction = int(base_atoms[0]), int(base_atoms[1])

    # Build neighbor sets
    # neighbors(n): all atoms bonded to n (undirected)
    col0, col1 = bonds_df[0].astype(int), bonds_df[1].astype(int)
    nbr_dir = set(col1[col0 == direction].astype(int)).union(set(col0[col1 == direction].astype(int)))
    nbr_org = set(col1[col0 == origin].astype(int)).union(set(col0[col1 == origin].astype(int)))

    # Helper to pick a 'third' atom:
    def pick_third():
        # Prefer neighbor of direction that is not origin
        cand = [a for a in nbr_dir if a != origin]
        if cand:
            return int(cand[0])
        # If direction is terminal (only bonded to origin), try neighbor of origin that is not direction
        cand2 = [a for a in nbr_org if a != direction]
        if cand2:
            return int(cand2[0])
        # Nothing usable
        return None

    # If user supplied a 3rd atom explicitly
    if len(base_atoms) >= 3:
        third = int(base_atoms[2])
        if third == origin:
            third = pick_third()
            if third is None:
                raise ValueError(
                    f"Cannot determine third atom: 'direction'={direction} has no neighbor != origin={origin}, "
                    f"and 'origin' has no neighbor != direction."
                )
            return [origin, direction, third]
        else:
            # If the provided third is bonded to direction, accept it; otherwise try to fix gracefully
            if (third in nbr_dir) and (third != origin):
                return [origin, direction, third]
            # Try automatic pick if given third is not valid
            auto = pick_third()
            if auto is not None:
                return [origin, direction, auto]
            raise ValueError(
                f"Provided third atom {third} is not bonded to 'direction'={direction}; "
                f"also failed to auto-pick a neighbor."
            )

    # Only two atoms given → pick a third automatically
    third = pick_third()
    if third is None:
        raise ValueError(
            f"Cannot determine third atom from bonds: 'direction'={direction} neighbors={sorted(nbr_dir)}, "
            f"'origin'={origin} neighbors={sorted(nbr_org)}."
        )

    return [origin, direction, third]



    

def get_molecule_connections(bonds_df, source, direction, mode='all'):
    """The atoms on the direction side of the source-direction bond, with the source (sorted).

    These are exactly the atoms on simple paths that leave ``source`` through ``direction``:
    everything ``direction`` reaches without passing back through ``source``. One breadth-first
    search finds them, where enumerating the paths grows exponentially with the rings.
    ``mode='shortest'``: the shortest such path, the bond itself.
    """
    source, direction = int(source), int(direction)
    if mode == 'shortest':
        return [source, direction]
    graph = graph_from_bonds_df(bonds_df, directed=False)
    graph.add_edge(source, direction)
    seen, queue = {source, direction}, [direction]
    while queue:
        atom = queue.pop()
        for nxt in graph.adj[atom]:
            if nxt not in seen:
                seen.add(nxt); queue.append(nxt)
    return np.array(sorted(seen), dtype=int)


def find_longest_simple_path(bonds_df):
    """
    bonds_df : pandas.DataFrame with exactly two columns [u, v] (integer atom indices)

    Returns
    -------
    List[int]
        Vertex sequence along a graph-diameter shortest path.
    """
    return diameter_vertex_path(graph_from_bonds_df(bonds_df, directed=False))

def get_specific_bonded_atoms_df(bonds_df, longest_path, coordinates_df):
    """
    Returns a DataFrame of atoms bonded along the longest path.
    Includes debug prints to trace internal state.

    Parameters
    ----------
    bonds_df : pd.DataFrame
        DataFrame of bonds, with atom indices (1-based)
    longest_path : list or None
        List of atom indices to include (1-based)
    coordinates_df : pd.DataFrame
        DataFrame of atomic coordinates, must include 'atom' column

    Returns
    -------
    pd.DataFrame
        A DataFrame with bonded atom names and their original bond indices
    """

    # 1) Filter bonds based on longest_path
    if longest_path is not None:
   
        mask = bonds_df.isin(longest_path)
  
        edited_bonds_df = bonds_df[mask.all(axis=1)].dropna().reset_index(drop=True)
       
    else:
        edited_bonds_df = bonds_df.copy()

    if edited_bonds_df.empty:
        return pd.DataFrame(columns=XYZConstants.BONDED_COLUMNS.value)


    # 2) Build bonds_array (zero-based indices)
    bonds_array = (edited_bonds_df.values - 1).astype(int)
    # 3) Atom names of both ends of every bond, shape (n_bonds, 2)
    atom_bonds = coordinates_df['atom'].to_numpy()[bonds_array]

    # 5) Concatenate with edited_bonds_df
    bonded_atoms_df = pd.concat([pd.DataFrame(atom_bonds), edited_bonds_df], axis=1)
    bonded_atoms_df.columns = XYZConstants.BONDED_COLUMNS.value
    return bonded_atoms_df

def scan_b1_over_angles(plane, degree_list):
    """b1s_for_loop_function with the arguments the other way round (kept for its callers)."""
    return b1s_for_loop_function(degree_list, plane)


def get_b1s_list(extended_df, scans=1, plot_result=False):
    """
    Calculate B1 values by scanning over rotation angles using circle-expanded atoms.

    Parameters
    ----------
    extended_df : pd.DataFrame
        DataFrame with at least the columns: 'x', 'z', 'radius', 'L', 'B5', 'loc_B5'.
    scans : int, optional
        Coarse step size in degrees.
    plot_result : bool, optional
        Kept for compatibility.

    Returns
    -------
    list
        [b1s, b1_b5_angle, planes, b1_normals]

        b1_b5_angle is the DIRECTED azimuthal separation of the B1 and B5
        directions, in [0, 180]. calc_sterimol folds it to [0, 90] before
        reporting it as B1_B5_angle, because the B1 side choice is unstable
        when the two supporting planes are near-equidistant; the directed value
        is passed through as B1_B5_azimuth.
    """
    # every atom's circle at once, point for point what generate_circle draws, atom after atom
    theta = np.linspace(0, 2 * np.pi, STERIMOL_CIRCLE_POINTS)
    cx, cz, r = (extended_df[k].to_numpy(float)[:, None] for k in ("x", "z", "radius"))
    original_plane_xz = np.column_stack(((cx + r * np.cos(theta)).ravel(), (cz + r * np.sin(theta)).ravel()))

    # Coarse scan on original geometry
    degree_list = list(range(18, 108, scans))
    sterimol_df = b1s_for_loop_function(degree_list, original_plane_xz)

    best_idx = int(sterimol_df["B1"].idxmin())
    best_angle = int(sterimol_df.loc[best_idx, "degree"])

    # Fine scan, still on the original geometry
    back_ang = best_angle - scans
    front_ang = best_angle + scans
    new_degree_list = list(range(back_ang, front_ang + 1))

    sterimol_df = b1s_for_loop_function(new_degree_list, original_plane_xz)

    b1s = sterimol_df["B1"].to_numpy()
    b1_b5_angle = sterimol_df["B1_B5_angle"].to_numpy()
    plane_xz = sterimol_df["plane"].to_list()
    b1_normals = np.array(sterimol_df["b1_normal"].to_list())  # shape (n, 2): [x, z]

    return [b1s, b1_b5_angle, plane_xz, b1_normals]


B1_TANGENT_TOL = 0.05


def calc_loc_b1(best_plane, edited_coordinates_df, tol=B1_TANGENT_TOL):
    """Position along the substituent axis of the contact that sets B1.

    loc_B5 is unambiguous because B5 is realised by a single atom, so it has a
    y of its own. B1 is not an atom: it is a plane containing the y-axis, and a
    plane parallel to that axis has no position along it. What does have one is
    the contact - the atoms whose vdW circles touch the plane - so loc_B1 is
    reported as the y of that contact set.

    The set almost never has one member. A minimum-width supporting plane is
    minimal precisely because it lies flush against two atoms; if it touched
    only one, rotating would narrow it further. Measured on the Ackermann set
    the contact is a single atom in 0% of cases (median 2, up to 5), and on the
    Corminboeuf ligands in 16%. Picking "the" tangent atom is therefore an
    arbitrary tiebreak, and the y is averaged over the contact set instead,
    weighted by how much of each atom lies on the plane. B1_tangent_atoms is
    returned alongside so callers can see how sharp the contact was.

    Parameters
    ----------
    best_plane : np.ndarray, shape (n_atoms * STERIMOL_CIRCLE_POINTS, 2)
        The sampled vdW circles rotated to the angle that minimises B1 - the
        same array get_b1s_list returns and calc_sterimol selects from.
    edited_coordinates_df : pd.DataFrame
        The substructure atoms, in the row order the circles were built in,
        carrying 'y' (the position along the Sterimol axis).

    Returns
    -------
    (loc_B1, n_tangent_atoms)
    """
    n_atoms = len(edited_coordinates_df)
    if n_atoms == 0 or len(best_plane) != n_atoms * STERIMOL_CIRCLE_POINTS:
        return float('nan'), 0

    owner = np.repeat(np.arange(n_atoms), STERIMOL_CIRCLE_POINTS)
    y_vals = edited_coordinates_df['y'].to_numpy(dtype=float)

    # Which of the four extents is B1, i.e. which edge the plane sits on.
    extents = [best_plane[:, 0].max(), best_plane[:, 0].min(),
               best_plane[:, 1].max(), best_plane[:, 1].min()]
    min_index = int(np.argmin(np.abs(extents)))
    column = 0 if min_index in (0, 1) else 1
    take_max = min_index in (0, 2)

    values = best_plane[:, column]
    edge = values.max() if take_max else values.min()
    touching = np.abs(values - edge) < tol
    if not touching.any():                      # tolerance too tight
        touching = values == edge

    atoms, counts = np.unique(owner[touching], return_counts=True)
    loc_b1 = float(np.average(y_vals[atoms], weights=counts))
    return loc_b1, int(len(atoms))


def calc_sterimol(bonded_atoms_df, extended_df, visualize_bool=False):
    edited_coordinates_df, index_list = filter_atoms_for_sterimol(bonded_atoms_df, extended_df)

    b1s, b1_b5_angle, planes, b1_normals = get_b1s_list(edited_coordinates_df)

    best_idx = int(np.argmin(b1s))
    best_b1 = float(b1s[best_idx])

    # B5, L, and loc_B5 should come from the original Sterimol dataframe
    best_b5 = float(np.max(edited_coordinates_df["B5"]))
    best_b5_index = edited_coordinates_df["B5"].idxmax()
    best_loc_b5 = float(edited_coordinates_df.loc[best_b5_index, "loc_B5"])
    best_L = float(np.max(edited_coordinates_df["L"]))

    best_loc_b1, n_tangent = calc_loc_b1(planes[best_idx], edited_coordinates_df)

    # Three related angles, reported separately because they are not the same
    # quantity and the literature has used more than one of them.
    #
    #   B1_B5_angle    the documented descriptor: the acute angle between the
    #                  B5 atom's origin-to-atom vector and the B1 plane, in
    #                  [0, 90]. This is what parts/desc_sterimol.tex defines
    #                  and what case studies 1 and 2 contain, so its meaning is
    #                  left unchanged.
    #   B1_B5_azimuth  the DIRECTED separation of the B1 and B5 directions
    #                  measured in the projection plane, in [0, 180].
    #   B1_B5_acute    the same separation folded to [0, 90].
    #
    # They are related but not interchangeable. n lies in the projection plane,
    # so v.n = p.n for the perpendicular part p, giving
    #     sin(B1_B5_angle) = |cos(azimuth)| * (|p| / |v|)
    # -- the in-plane angle scaled by how far along the axis the B5 atom sits,
    # which loc_B5 already reports. The scale factor runs down to 0.46 in
    # practice, so B1_B5_angle is partly confounded with loc_B5 and the two
    # in-plane forms are not.
    #
    # Which of the in-plane forms is safe depends on the ligands. B1 is the
    # closer of two supporting planes, so its normal points to whichever side
    # is nearer; where the two sides are near-equidistant that choice flips and
    # the directed value jumps to its supplement. Measuring the gap
    # |d_far - d_near|: below 0.5 A in 23 of 44 arms on the Ackermann set, but
    # in only 23 of 1178 arm-conformers on the Corminboeuf bidentates. Prefer
    # B1_B5_acute unless the gap has been checked for the set in hand.
    best_azimuth = float(b1_b5_angle[best_idx])
    best_acute = float(min(best_azimuth, 180.0 - best_azimuth))

    b5_row = edited_coordinates_df.loc[best_b5_index]
    v_b5 = np.array([float(b5_row["x"]), float(b5_row["y"]), float(b5_row["z"])])
    n_xz = b1_normals[best_idx]
    n_3d = np.array([n_xz[0], 0.0, n_xz[1]])
    v_norm = np.linalg.norm(v_b5)
    if v_norm > 1e-10:
        sin_val = np.clip(np.abs(np.dot(v_b5, n_3d)) / v_norm, 0.0, 1.0)
        best_angle = float(np.degrees(np.arcsin(sin_val)))
    else:
        best_angle = 0.0
    result_df = pd.DataFrame([{
        "B1": best_b1,
        "B5": best_b5,
        "L": best_L,
        "loc_B1": best_loc_b1,
        "loc_B5": best_loc_b5,
        "B1_tangent_atoms": n_tangent,
        "B1_B5_angle": best_angle,
        "B1_B5_azimuth": best_azimuth,
        "B1_B5_acute": best_acute
    }])

    if visualize_bool:
        from utils.visualize import plot_b1_visualization
        edited_coordinates_df.index = [index_list[i] + 1 for i in range(len(index_list))]
        plot_b1_visualization(planes[best_idx], edited_coordinates_df=edited_coordinates_df,
                              sterimol_df=result_df)

    return result_df

def prepare_sterimol_inputs(coordinates_df, bonds_df, base_atoms, connected_from_direction,
                            radii='CPK', sub_structure=True, drop_atoms=None, mode='all'):
    """
    Align the molecule to the Sterimol frame and build the two tables the
    Sterimol calculation runs on.

    Returns
    -------
    (bonded_atoms_df, extended_df, base_atoms)
        base_atoms is returned because dropping atoms renumbers it.
    """
    # Drop atoms and update indices
    if drop_atoms is not None:

        drop_atoms_adjusted = adjust_indices(drop_atoms)
        coordinates_df = coordinates_df.drop(drop_atoms_adjusted)
        coordinates_df = coordinates_df.reset_index(drop=True)

        # Create a mapping from old to new indices
        old_to_new = {old_idx: new_idx for new_idx, old_idx in enumerate(sorted(coordinates_df.index))}
        # Map bonds_df to new indices, drop any bonds involving dropped atoms
        mask = ~bonds_df[0].isin(drop_atoms) & ~bonds_df[1].isin(drop_atoms)

        bonds_df = bonds_df[mask].copy()
        bonds_df[0] = bonds_df[0].map(old_to_new)
        bonds_df[1] = bonds_df[1].map(old_to_new)
        bonds_df = bonds_df.dropna().reset_index(drop=True)  # Drop any rows with NaN values after mapping
        # Update base_atoms as well
        base_atoms = [old_to_new[a] for a in base_atoms if a in old_to_new]
        if connected_from_direction is not None:
            connected_from_direction = [old_to_new[a] for a in connected_from_direction if a in old_to_new]

    bonds_direction = adjust_indices(direction_atoms_for_sterimol(bonds_df, base_atoms))
    if isinstance(base_atoms[0], list):
        base_atoms=base_atoms[0]
        # print(f'Debug! base_atoms adjusted to: {base_atoms}')
    # verify bonds direction is 3 atoms
    if len(bonds_direction) != 3:
        raise ValueError("Bonds direction must contain exactly 3 atoms, check molecule connectivity.")
    new_coordinates_df = preform_coordination_transformation(coordinates_df, bonds_direction)

    if sub_structure:
        if connected_from_direction is None:
            connected_from_direction = get_molecule_connections(bonds_df, base_atoms[0], base_atoms[1], mode=mode)
    else:
        connected_from_direction = None

    bonded_atoms_df = get_specific_bonded_atoms_df(bonds_df, connected_from_direction, new_coordinates_df)
    extended_df = get_extended_df_for_sterimol(new_coordinates_df, bonds_df, radii)

    return bonded_atoms_df, extended_df, base_atoms


def get_sterimol_df(coordinates_df, bonds_df, base_atoms, connected_from_direction,
                    radii='CPK', sub_structure=True, drop_atoms=None, visualize_bool=False, mode='all'):
    bonded_atoms_df, extended_df, base_atoms = prepare_sterimol_inputs(
        coordinates_df, bonds_df, base_atoms, connected_from_direction,
        radii=radii, sub_structure=sub_structure, drop_atoms=drop_atoms, mode=mode)

    sterimol_df = calc_sterimol(bonded_atoms_df, extended_df, visualize_bool)
    sterimol_df = sterimol_df.rename(index={0: str(base_atoms[0]) + '-' + str(base_atoms[1])})
    sterimol_df = sterimol_df.round(4)

    return sterimol_df


def get_sterimol_plot_data(coordinates_df, bonds_df, base_atoms, connected_from_direction=None,
                           radii='CPK', sub_structure=True, drop_atoms=None, mode='all'):
    """
    The two inputs the steriplots need, from the same pipeline get_sterimol_df runs.

    Returns
    -------
    (edited_coordinates_df, best_plane)
        edited_coordinates_df : substructure atoms in the Sterimol frame, with
            'atom', 'x', 'y', 'z', 'radius' and 'B5' - indexed 1-based to match
            the atom numbering shown to users.
        best_plane : (n_atoms * 100, 2) array of sampled vdW circles in the xz
            plane, rotated to the angle that minimises B1.

    Pass these to plot_b1_visualization(best_plane, edited_coordinates_df) and
    plot_L_B5_plane(edited_coordinates_df, ...).
    """
    bonded_atoms_df, extended_df, _ = prepare_sterimol_inputs(
        coordinates_df, bonds_df, base_atoms, connected_from_direction,
        radii=radii, sub_structure=sub_structure, drop_atoms=drop_atoms, mode=mode)

    edited_coordinates_df, index_list = filter_atoms_for_sterimol(bonded_atoms_df, extended_df)
    b1s, _, planes, _ = get_b1s_list(edited_coordinates_df)
    best_idx = int(np.argmin(b1s))

    # match calc_sterimol's visualize path: label rows with 1-based atom numbers
    edited_coordinates_df = edited_coordinates_df.copy()
    edited_coordinates_df.index = [i + 1 for i in index_list]

    return edited_coordinates_df, planes[best_idx]