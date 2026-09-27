import pandas as pd
import numpy as np
import os
import sys
import re
import glob
import math
from enum import Enum
import networkx as nx

from typing import *

import warnings
from scipy.spatial.distance import pdist, squareform

from sklearn.preprocessing import MinMaxScaler
from morfeus import Sterimol, BuriedVolume, read_xyz
warnings.filterwarnings("ignore", category=RuntimeWarning)


def flatten_list(nested_list_arg: List[list]) -> List:
    """
    Flatten a nested list.
    turn [[1,2],[3,4]] to [1,2,3,4]
    """
    flat_list=[item for sublist in nested_list_arg for item in sublist]
    return flat_list

def split_strings(strings_list):
    split_list = []
    for string in strings_list:
        split_list.extend(string.split())
    return split_list

def get_df_from_file(filename,columns=['atom','x','y','z'],index=None):
    """
    Parameters
    ----------
    filename : str
        full file name to read.
    columns : str , optional
        list of column names for DataFrame. The default is None.
    splitter : str, optional
        input for [.split().] , for csv-',' for txt leave empty. The default is None.
    dtype : type, optional
        type of variables for dataframe. The default is None.

    Returns
    -------
    df : TYPE
        DESCRIPTION.

    """
    with open(filename, 'r') as f:
        lines=f.readlines()[2:]
    splitted_lines=split_strings(lines)
    df=pd.DataFrame(np.array(splitted_lines).reshape(-1,4),columns=columns,index=index)
    df[['x','y','z']]=df[['x','y','z']].astype(float)
    return df


class XYZConstants(Enum):
    """
    Constants related to XYZ file processing
    """
    DF_COLUMNS=['atom','x','y','z']
    STERIMOL_INDEX = ['B1', 'B5', 'L', 'loc_B1','loc_B5']
    DIPOLE_COLUMNS = ['dip_x', 'dip_y', 'dip_z', 'total_dipole']
    RING_VIBRATION_COLUMNS = ['cross', 'cross_angle', 'para', 'para_angle']
    RING_VIBRATION_INDEX=['Product','Frequency','Sin_angle']
    VIBRATION_INDEX = ['Frequency', 'Amplitude']
    BONDED_COLUMNS = ['atom_1', 'atom_2', 'index_1', 'index_2']
    NOF_ATOMS = ['N', 'O', 'F']
    STERIC_PARAMETERS = ['B1', 'B5', 'L', 'loc_B1', 'loc_B5','RMSD']
    ELECTROSTATIC_PARAMETERS = ['dip_x', 'dip_y', 'dip_z', 'total_dipole','energy']

class GeneralConstants(Enum):
    """
    Holds constants for calculations and conversions
    1. covalent radii from Alvarez (2008) DOI: 10.1039/b801115j
    2. atomic numbers
    2. atomic weights
    """
    COVALENT_RADII= {
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
    # CPK_RADII={
    #     'C':1.50,   'H':1.00,   'S.O':1.70,  'Si':2.10,
    #     'C2':1.60,  'N':1.50,   'S1':1.00,   'Co':2.00,
    #     'C3':1.60,  'C66':1.70, 'F':1.35,    'Ni':2.00,
    #     'C4':1.50,  'N4':1.45,  'Cl':1.75,
    #     'C5/N5':1.70, 'O':1.35, 'S4':1.40,
    #     'C6/N6':1.70, 'O2':1.35, 'Br':1.95,
    #     'C7':1.70,    'P':1.40,  'I':2.15,
    #     'C8':1.50,    'S':1.70,  'B':1.92,
    
    # }

    REGULAR_BOND_TYPE = {

        'O.2': 'O', 'N.2': 'N', 'S.3': 'S',
        'O.3': 'O', 'N.1': 'N', 'S.O2': 'S',
        'O.co2': 'O', 'N.3': 'N', 'P.3': 'P',
        'C.1': 'C', 'N.ar': 'N',
        'C.2': 'C', 'N.am': 'N',
        "C.cat": 'C', 'N.pl3': 'N',
        'C.3': 'C', 'N.4': 'N',
        'C.ar': 'C', 'S.2': 'S',
    }

    BOND_TYPE={
        
        'O.2':'O2', 'N.2':'C6/N6','S.3':'S4',
        'O.3':'O', 'N.1':'N', 'S.O2':'S',
        'O.co2':'O', 'N.3':'C6/N6','P.3':'P',
        'C.1':'C', 'N.ar':'C6/N6',
        'C.2':'C3', 'N.am':'C6/N6',
        "C.cat":'C3', 'N.pl3':'C6/N6',
        'C.3':'C', 'N.4':'N4',
        'C.ar':'C6/N6', 'S.2':'S','H':'H' 
        }
    
    ATOMIC_NUMBERS ={
    '1':'H', '5':'B', '6':'C', '7':'N', '8':'O', '9':'F', '14':'Si',
             '15':'P', '16':'S', '17':'Cl', '35':'Br', '53':'I', '27':'Co', '28':'Ni'}
        

    ATOMIC_WEIGHTS = {
            'H' : 1.008,'He' : 4.003, 'Li' : 6.941, 'Be' : 9.012,
            'B' : 10.811, 'C' : 12.011, 'N' : 14.007, 'O' : 15.999,
            'F' : 18.998, 'Ne' : 20.180, 'Na' : 22.990, 'Mg' : 24.305,
            'Al' : 26.982, 'Si' : 28.086, 'P' : 30.974, 'S' : 32.066,
            'Cl' : 35.453, 'Ar' : 39.948, 'K' : 39.098, 'Ca' : 40.078,
            'Sc' : 44.956, 'Ti' : 47.867, 'V' : 50.942, 'Cr' : 51.996,
            'Mn' : 54.938, 'Fe' : 55.845, 'Co' : 58.933, 'Ni' : 58.693,
            'Cu' : 63.546, 'Zn' : 65.38, 'Ga' : 69.723, 'Ge' : 72.631,
            'As' : 74.922, 'Se' : 78.971, 'Br' : 79.904, 'Kr' : 84.798,
            'Rb' : 84.468, 'Sr' : 87.62, 'Y' : 88.906, 'Zr' : 91.224,
            'Nb' : 92.906, 'Mo' : 95.95, 'Tc' : 98.907, 'Ru' : 101.07,
            'Rh' : 102.906, 'Pd' : 106.42, 'Ag' : 107.868, 'Cd' : 112.414,
            'In' : 114.818, 'Sn' : 118.711, 'Sb' : 121.760, 'Te' : 126.7,
            'I' : 126.904, 'Xe' : 131.294, 'Cs' : 132.905, 'Ba' : 137.328,
            'La' : 138.905, 'Ce' : 140.116, 'Pr' : 140.908, 'Nd' : 144.243,
            'Pm' : 144.913, 'Sm' : 150.36, 'Eu' : 151.964, 'Gd' : 157.25,
            'Tb' : 158.925, 'Dy': 162.500, 'Ho' : 164.930, 'Er' : 167.259,
            'Tm' : 168.934, 'Yb' : 173.055, 'Lu' : 174.967, 'Hf' : 178.49,
            'Ta' : 180.948, 'W' : 183.84, 'Re' : 186.207, 'Os' : 190.23,
            'Ir' : 192.217, 'Pt' : 195.085, 'Au' : 196.967, 'Hg' : 200.592,
            'Tl' : 204.383, 'Pb' : 207.2, 'Bi' : 208.980, 'Po' : 208.982,
            'At' : 209.987, 'Rn' : 222.081, 'Fr' : 223.020, 'Ra' : 226.025,
            'Ac' : 227.028, 'Th' : 232.038, 'Pa' : 231.036, 'U' : 238.029,
            'Np' : 237, 'Pu' : 244, 'Am' : 243, 'Cm' : 247
    }

_METAL_ELEMENTS = {
    'Li', 'Be', 'Na', 'Mg', 'K', 'Ca', 'Rb', 'Sr', 'Cs', 'Ba', 'Fr', 'Ra',
    'Al', 'Ga', 'Ge', 'In', 'Sn', 'Sb', 'Tl', 'Pb', 'Bi',
    'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
    'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
    'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
    'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy',
    'Ho', 'Er', 'Tm', 'Yb', 'Lu',
    'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm',
}

import numpy.typing as npt


def adjust_indices(indices: npt.ArrayLike, adjustment_num: int=1) -> npt.ArrayLike:
    """
    adjust indices by adjustment_num
    """
    return np.array(indices)-adjustment_num

def adjust_indices_xyz(indices: npt.ArrayLike) -> npt.ArrayLike:
    """
    adjust indices by adjustment_num
    """
    return adjust_indices(indices, adjustment_num=1)

def calc_angle(p1: npt.ArrayLike, p2: npt.ArrayLike, degrees: bool=False) -> float: ###works, name in R: 'angle' , radians
    dot_product=np.dot(p1, p2)
    norm_p1=np.linalg.norm(p1)
    norm_p2=np.linalg.norm(p2)
    thetha=np.arccos(dot_product/(norm_p1*norm_p2))
    if degrees:
        thetha=np.degrees(thetha)   
    return thetha
    
def calc_new_base_atoms(coordinates_array: npt.ArrayLike, atom_indices: npt.ArrayLike):  #help function for calc_coordinates_transformation
    """
    a function that calculates the new base atoms for the transformation of the coordinates.
    optional: if the atom_indices is 4, the origin will be the middle of the first two atoms.
    """
    new_origin=coordinates_array[atom_indices[0]]
    if (len(atom_indices)==4):
        new_origin=(new_origin+coordinates_array[atom_indices[1]])/2
    new_y=(coordinates_array[atom_indices[-2]]-new_origin)/np.linalg.norm((coordinates_array[atom_indices[-2]]-new_origin))
    coplane=((coordinates_array[atom_indices[-1]]-new_origin)/np.linalg.norm((coordinates_array[atom_indices[-1]]-new_origin)+0.00000001))
    return (new_origin,new_y,coplane)

def np_cross_and_vstack(plane_1, plane_2):
    cross_plane=np.cross(plane_1, plane_2)
    united_results=np.vstack([plane_1, plane_2, cross_plane])
    return united_results

def calc_basis_vector(origin, y: npt.ArrayLike, coplane: npt.ArrayLike):#help function for calc_coordinates_transformation
    """
    origin: origin of the new basis
    y: y direction of the new basis
    coplane: a vector that is coplanar with the new y direction
    """
    # cross_y_plane=np.cross(coplane,y)
    # coef_mat=np.vstack([y, coplane, cross_y_plane])
    coef_mat=np_cross_and_vstack(coplane, y)
    angle_new_y_coplane=calc_angle(coplane,y)
    cop_ang_x=angle_new_y_coplane-(np.pi/2)
    # result_vector=[0,np.cos(cop_ang_x),0]
    result_vector=[np.cos(cop_ang_x), 0, 0]
    new_x,_,_,_=np.linalg.lstsq(coef_mat,result_vector,rcond=None)
    new_basis=np_cross_and_vstack(new_x, y)
    # new_z=np.cross(new_x,y)
    # new_basis=np.vstack([new_x, y, new_z])
    return new_basis

def transform_row(row_array, new_basis, new_origin, round_digits):
    translocated_row = row_array - new_origin
    return np.dot(new_basis, translocated_row).round(round_digits)



def calc_coordinates_transformation(coordinates_array: npt.ArrayLike, base_atoms_indices: npt.ArrayLike, round_digits:int=4 ,origin:npt.ArrayLike=None) -> npt.ArrayLike:#origin_atom, y_direction_atom, xy_plane_atom
    """
    a function that recives coordinates_array and new base_atoms_indices to transform the coordinates by
    and returns a dataframe with the shifted coordinates
    parameters:
    ----------
    coordinates_array: np.array
        xyz molecule array
    base_atoms_indices: list of nums
        indices of new atoms to shift coordinates by.
    origin: in case you want to change the origin of the new basis, middle of the ring for example.
    returns:
        transformed xyz molecule dataframe
    -------
        
    example:
    -------
    calc_coordinates_transformation(coordinates_array,[2,3,4])
    
    Output:
        atom       x       y       z
      0    H  0.3477 -0.5049 -1.3214
      1    B     0.0     0.0     0.0
      2    B    -0.0  1.5257     0.0
    """
    indices=adjust_indices_xyz(base_atoms_indices)
    new_basis=calc_basis_vector(*calc_new_base_atoms(coordinates_array,indices))    
    if origin is None:
        new_origin=coordinates_array[indices[0]]
    else:
        new_origin=origin

    transformed_coordinates = np.apply_along_axis(lambda x: transform_row(x, new_basis, new_origin, round_digits), 1,
                                                  coordinates_array)
    # transformed_coordinates=np.array([np.dot(new_basis,(row-new_origin)) for row in coordinates_array]).round(round_digits)
    return transformed_coordinates

def preform_coordination_transformation(xyz_df, indices=None):
    xyz_copy=xyz_df.copy()
    coordinates=np.array(xyz_copy[['x','y','z']].values)
    if indices is None:
        xyz_copy[['x','y','z']]=calc_coordinates_transformation(coordinates, [1,2,3])
    else:
      
        xyz_copy[['x','y','z']]=calc_coordinates_transformation(coordinates, indices)
 
    return xyz_copy

def calc_dipole_gaussian(coordinates_array, gauss_dipole_array, base_atoms_indices ,geometric_transformation_indices=None):
    """
    a function that recives coordinates and gaussian dipole, transform the coordinates
    by the new base atoms and calculates the dipole in each axis
    """
    if geometric_transformation_indices:
        # Calculate the geometric center of specified indices
        geometric_center = np.mean(coordinates_array[geometric_transformation_indices], axis=0)
        # Translate all coordinates
        coordinates_array -= geometric_center

    indices=adjust_indices(base_atoms_indices)
    basis_vector=calc_basis_vector(*calc_new_base_atoms(coordinates_array, indices))
    gauss_dipole_array[0,0:3]=np.matmul(basis_vector,gauss_dipole_array[0,0:3])
    dipole_df=pd.DataFrame(gauss_dipole_array,columns=['dipole_x','dipole_y','dipole_z','total'])
    # print(geometric_transformation_indices, dipole_df)
    return dipole_df

def check_imaginary_frequency(info_df):##return True if no complex frequency, called ground.state in R
        bool_imaginary=not any([isinstance(frequency, complex) for frequency in info_df['Frequency']])
        return bool_imaginary



def indices_to_coordinates_vector(coordinates_array,indices):
    """
    a function that recives coordinates_array and indices of two atoms
    and returns the bond vector between them
    """

    if  isinstance(indices[0], tuple):
        bond_vector=[(coordinates_array[index[0]]-coordinates_array[index[1]]) for index in indices]
    else:
        bond_vector= coordinates_array[indices[0]]-coordinates_array[indices[1]]

    return bond_vector

def get_bonds_vector_for_calc_angle(coordinates_array,atoms_indices): ##for calc_angle_between_atoms

    indices=adjust_indices(atoms_indices)#three atoms-angle four atoms-dihedral
    augmented_indices=[indices[0],indices[1],indices[1],indices[2]]
    if len(indices)==4:
        augmented_indices.extend([indices[2],indices[3]])
    indices_pairs=list(zip(augmented_indices[::2],augmented_indices[1::2]))
  
    bond_vector=indices_to_coordinates_vector(coordinates_array,indices_pairs)
    return bond_vector
  

def calc_angle_between_atoms(coordinates_array,atoms_indices): #gets a list of atom indices
    """
    a function that gets 3/4 atom indices, and returns the angle between thos atoms.
    Parameters
    ----------
    coordinates_array: np.array
        contains x y z atom coordinates
        
    atoms_indices- list of ints
        a list of atom indices to calculate the angle between- [2,3,4]
   
    Returns
    -------
    angle: float
        the bond angle between the atoms
    """
    bonds_list=get_bonds_vector_for_calc_angle(coordinates_array,atoms_indices)
    if len(atoms_indices)==3:
      
        angle=calc_angle(bonds_list[0], bonds_list[1]*(-1), degrees=True)
    else:
        first_cross=np.cross(bonds_list[0],bonds_list[1]*(-1))
        second_cross=np.cross(bonds_list[2]*(-1),bonds_list[1]*(-1)) 
        angle=calc_angle(first_cross, second_cross, degrees=True)
    return angle

def get_angle_df(coordinates_array, atom_indices):
    """
    a function that gets a list of atom indices, and returns a dataframe of angles between thos atoms.
    Parameters
    ----------
    coordinates_array: np.array
        contains x y z atom coordinates
    
    atom_indices- list of lists of ints
        a list of atom indices to calculate the angle between- [[2,3,4],[2,3,4,5]]
    """
 
    if isinstance(atom_indices, list) and all(isinstance(elem, list) for elem in atom_indices):
        indices_list=['angle_{}'.format(index) if len(index)==3 else 'dihedral_{}'.format(index) for index in atom_indices]
        angle_list=[calc_angle_between_atoms(coordinates_array,index) for index in atom_indices]
        return pd.DataFrame(angle_list,index=indices_list)
    else:
        indices_list=['angle_{}'.format(atom_indices) if len(atom_indices)==3 else 'dihedral_{}'.format(atom_indices)]
        angle=[calc_angle_between_atoms(coordinates_array,atom_indices)]
        return pd.DataFrame(angle,index=indices_list)


def calc_single_bond_length(coordinates_array,atom_indices):
    """
    a function that gets 2 atom indices, and returns the distance between thos atoms.
    Parameters
    ----------
    coordinates_array: np.array
        contains x y z atom coordinates
        
    atom_indices- list of ints
        a list of atom indices to calculate the distance between- [2,3]
   
    Returns
    -------
    distance: float
        the bond distance between the atoms
    """
    indices=adjust_indices(atom_indices)
    distance=np.linalg.norm(coordinates_array[indices[0]]-coordinates_array[indices[1]])
    return distance

def calc_bonds_length(coordinates_array,atom_pairs): 
    """
    a function that calculates the distance between each pair of atoms.
    help function for molecule class
    
    Parameters
    ----------
    coordinates_array: np.array
        xyz coordinates 
        
    atom_pairs : iterable
        list containing atom pairs-(([2,3],[4,5]))
        
    Returns
    -------
    pairs_df : dataframe
        distance between each pair 
        
        Output:
                               0
    bond length[2, 3]            1.525692
    bond length[4, 5]            2.881145
    
    """
    # Check the order of bond list and pairs
    bond_list=[calc_single_bond_length(coordinates_array,pair) for pair in atom_pairs]
    pairs=adjust_indices(atom_pairs)
    index=[('bond_length')+str(pair) for pair in pairs]
    pairs_df=pd.DataFrame(bond_list,index=index)
    return pairs_df



def direction_atoms_for_sterimol(bonds_df,base_atoms)->list: #help function for sterinol
    """
    a function that return the base atom indices for coordination transformation according to the bonded atoms.
    you can insert two atom indicess-[1,2] output [1,2,8] or the second bonded atom
    if the first one repeats-[1,2,1] output [1,2,3]
    """
    
    base_atoms_copy=base_atoms[0:2]
    origin,direction=base_atoms[0],base_atoms[1]
    bonds_df = bonds_df[~((bonds_df[0] == origin) & (bonds_df[1] == direction)) & 
                              ~((bonds_df[0] == direction) & (bonds_df[1] == origin))]
    
    try :
        base_atoms[2]==origin
        if(any(bonds_df[0]==direction)):
            # take the second atom in the bond where the first equeal to the direction, second option
            base_atoms_copy[2]=int(bonds_df[(bonds_df[0]==direction)][1].iloc[1])
        else:
            # take the first atom in the bond where the first equeal to the direction, second option
            base_atoms_copy[2]=int(bonds_df[(bonds_df[1]==direction)][0].iloc[1])
    except: 
       
        for _, row in bonds_df.iterrows():
            if row[0] == direction:
                base_atoms_copy.append(row[1])
                break
            elif row[1] == direction:
                base_atoms_copy.append(row[0])
                break
    return base_atoms_copy

def get_molecule_connections(bonds_df, source, direction):
    source, direction = int(source), int(direction)
    edges = pd.DataFrame({
        0: pd.to_numeric(bonds_df.iloc[:, 0]).astype(int),
        1: pd.to_numeric(bonds_df.iloc[:, 1]).astype(int),
    })
    graph = nx.from_pandas_edgelist(edges, source=0, target=1, create_using=nx.Graph)
    if source not in graph:
        return np.array([], dtype=int)
    paths = []
    for target in list(graph.nodes()):
        if int(target) == source:
            continue
        try:
            paths.extend(nx.all_simple_paths(graph, source, int(target)))
        except (nx.NetworkXError, nx.NodeNotFound):
            continue
    with_direction = [path for path in paths if direction in path]
    if not with_direction:
        return np.array([], dtype=int)
    return np.unique(flatten_list(with_direction))




def get_specific_bonded_atoms_df(bonds_df,longest_path,coordinates_df):
    """
    a function that returns a dataframe of the atoms that are bonded in the longest path.
    bonded_atoms_df: dataframe
        the atom type and the index of each bond.

       atom_1 atom_2 index_1 index_2
0       C      N       1       4
1       C      N       1       5
2       C      N       2       3
    """
    if longest_path is not None:
        edited_bonds_df=bonds_df[(bonds_df.isin(longest_path))].dropna().reset_index(drop=True)
    else:
        edited_bonds_df=bonds_df
    bonds_array=(np.array(edited_bonds_df)-1).astype(int) # adjust indices? 
    atom_bonds=np.vstack([(coordinates_df.iloc[bond]['atom'].values) for bond in bonds_array]).reshape(-1,2)
    bonded_atoms_df=(pd.concat([pd.DataFrame(atom_bonds),edited_bonds_df],axis=1))
    bonded_atoms_df.columns=[XYZConstants.BONDED_COLUMNS.value]
    return bonded_atoms_df


def remove_atom_bonds(bonded_atoms_df,atom_remove='H'):
    atom_bonds_array=np.array(bonded_atoms_df)
    delete_rows_left=np.where(atom_bonds_array[:,0]==atom_remove)[0] #itterrow [0] is index [1] are the values
    delete_rows_right=np.where(atom_bonds_array[:,1]==atom_remove)[0]
    atoms_to_delete=np.concatenate((delete_rows_left,delete_rows_right))
    new_bonded_atoms_df=bonded_atoms_df.drop((atoms_to_delete),axis=0)
    return new_bonded_atoms_df



# A non-metal pair is bonded below BOND_SCALE x (sum of Pyykko covalent radii). Over 462
# structures (22,210 atoms: CS1-3, their ligands, the explorer fixtures) real bonds reach 1.05x
# and the closest non-bonded pair sits at 1.26x (a P...N contact); every scale from 1.10 to 1.25
# gives the same bonds.
BOND_SCALE = 1.15

# The flat cutoff (A) used when a caller passes no threshold_distance. None means the covalent
# rule; descripytor.compat.paper_v3() sets 1.82 for the duration of a block.
DEFAULT_BOND_THRESHOLD = None


def extract_connectivity(xyz_df, threshold_distance=None, metals=None,
                         metal_threshold=2.8, max_coordination=6, scale=BOND_SCALE):
    """
    Build a connectivity table from XYZ coordinates.

    Parameters
    ----------
    threshold_distance : float or None
        None (the default): a non-metal pair is bonded when it is closer than
        ``scale`` x the sum of the two covalent radii, so long single bonds (S-CF3
        1.84 A, P-C 1.82-1.85 A, C-Br, C-I, Si-C, S-S) are kept. A number restores
        the old flat cutoff for non-metal pairs; 1.82 reproduces every table built
        before this change (tag paper-v3), including its halogen window up to 2.6 A.
    scale : float
        Covalent-radius multiple for the default rule.
    metals : None | str | list[str]
        Metal element symbols to treat with relaxed distance rules.
        None → auto-detect from the full periodic-table metal set.
        Pass an empty list [] to disable metal handling entirely.
    metal_threshold : float
        Max bond distance (Å) for metal–ligand pairs.
    max_coordination : int
        Maximum number of bonds kept per metal centre.
    """
    if metals is None:
        active_metals = _METAL_ELEMENTS
    elif isinstance(metals, str):
        active_metals = {metals}
    else:
        active_metals = set(metals)

    coordinates = np.array(xyz_df[['x', 'y', 'z']].values)
    atoms_symbol = np.array(xyz_df['atom'].values)
    distances = pdist(coordinates)
    dist_matrix = squareform(distances)

    dist_df = pd.DataFrame(dist_matrix).stack().reset_index()
    dist_df.columns = ['a1', 'a2', 'value']
    dist_df['first_atom'] = [atoms_symbol[i] for i in dist_df['a1']]
    dist_df['second_atom'] = [atoms_symbol[i] for i in dist_df['a2']]

    remove_list = []
    dist_array = np.array(dist_df)
    special_atoms = {'Cl', 'Br', 'F', 'I'}
    if threshold_distance is None:
        threshold_distance = DEFAULT_BOND_THRESHOLD
    flat = threshold_distance is not None
    radii = GeneralConstants.COVALENT_RADII.value

    for idx, row in enumerate(dist_array):
        i, j, dist, atom1, atom2 = row
        remove_flag = False

        if i == j:
            remove_flag = True

        if ((atom1 == 'H') and (atom2 not in XYZConstants.NOF_ATOMS.value)) or \
           ((atom1 == 'H') and (atom2 == 'H')) or \
           (flat and (atom1 == 'H' or atom2 == 'H') and float(dist) >= 1.5):
            remove_flag = True

        involves_metal = atom1 in active_metals or atom2 in active_metals

        if not involves_metal:
            limit = threshold_distance if flat else \
                scale * (radii.get(atom1, 0.77) + radii.get(atom2, 0.77))
            if float(dist) >= limit or float(dist) == 0:
                remove_flag = True
        else:
            if float(dist) > metal_threshold or float(dist) == 0:
                remove_flag = True
            # A metal-H pair has to be a hydride, not an agostic contact (Cu...H 2.5-2.7 A);
            # the flat rule got that from its 1.5 A limit on every H.
            if not flat and 'H' in (atom1, atom2) and                     float(dist) >= scale * (radii.get(atom1, 0.77) + radii.get(atom2, 0.77)):
                remove_flag = True

        # Flat rule only: halogens bonded between threshold and 2.6 Å are allowed
        # (C–Br, C–I); the covalent rule already reaches them.
        if flat and (atom1 in special_atoms or atom2 in special_atoms) and \
           (threshold_distance <= float(dist) < 2.6) and not involves_metal:
            remove_flag = False

        if remove_flag:
            remove_list.append(idx)

    dist_df = dist_df.drop(remove_list)
    dist_df[['min_col', 'max_col']] = pd.DataFrame(
        np.sort(dist_df[['a1', 'a2']], axis=1), index=dist_df.index
    )
    dist_df = dist_df.drop(columns=['a1', 'a2']).rename(columns={'min_col': 0, 'max_col': 1})
    dist_df = dist_df.drop_duplicates(subset=[0, 1])

    # Keep only the closest bond per halogen (terminal halogens bond to one atom)
    special_atoms_idxs = {}
    for idx, row in dist_df.iterrows():
        if row['first_atom'] in special_atoms:
            special_atoms_idxs.setdefault(row[0], []).append((idx, row['value']))
        if row['second_atom'] in special_atoms:
            special_atoms_idxs.setdefault(row[1], []).append((idx, row['value']))

    special_atoms_to_remove = []
    for atom_idx, bonds in special_atoms_idxs.items():
        if len(bonds) > 1:
            bonds.sort(key=lambda x: x[1])
            special_atoms_to_remove.extend([i for i, _ in bonds[1:]])
    dist_df = dist_df.drop(special_atoms_to_remove)

    # For each metal centre keep only the max_coordination shortest bonds
    metal_mask = dist_df['first_atom'].isin(active_metals) | dist_df['second_atom'].isin(active_metals)
    metal_bonds = dist_df[metal_mask].copy()
    non_metal = dist_df[~metal_mask].copy()

    kept_metal_idx = []
    if not metal_bonds.empty:
        metal_bonds['_metal_idx'] = metal_bonds.apply(
            lambda r: r[0] if r['first_atom'] in active_metals else r[1], axis=1
        )
        for _, group in metal_bonds.groupby('_metal_idx'):
            kept_metal_idx.extend(group.nsmallest(max_coordination, 'value').index)

    metal_kept = metal_bonds.loc[kept_metal_idx, [0, 1]] if kept_metal_idx else pd.DataFrame(columns=[0, 1])

    final = pd.concat([non_metal[[0, 1]], metal_kept], ignore_index=True)
    return pd.DataFrame(final[[0, 1]].apply(pd.to_numeric).astype(int) + 1)


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
    # print(center_atom)
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

def nob_atype(xyz_df, bonds_df):
    
    symbols = xyz_df['atom'].values
    
    list_results=[]
    for index,symbol in enumerate(symbols):
        index+=1
        nob = bonds_df[(bonds_df[0] == index) | (bonds_df[1] == index)].shape[0]
        if symbol == 'H':
            result = 'H'
        elif symbol == 'F':
            result = 'F'
        elif symbol == 'P':
            result = 'P'
        elif symbol == 'Cl':
            result = 'Cl'
        elif symbol == 'Br':
            result = 'Br'
        elif symbol == 'I':
            result = 'I'
        elif symbol == 'O':
            if nob < 1.5:
                result = 'O2'
            elif nob > 1.5:
                result = 'O'
        elif symbol == 'S':
            if nob < 2.5:
                result = 'S'
            elif 2.5 < nob < 5.5:
                result = 'S4'
            elif nob > 5.5:
                result = 'S1'
        elif symbol == 'N':
            if nob < 2.5:
                result = 'C6/N6'
            elif nob > 2.5:
                result = 'N'
        elif symbol == 'C':
            if nob < 2.5:
                result = 'C3'
            elif 2.5 < nob < 3.5:
                result = 'C6/N6'
            elif nob > 3.5:
                result = 'C'
        elif symbol in _METAL_ELEMENTS:
            result = symbol
        else:
            result = 'X'

        list_results.append(result)

    return list_results

def get_sterimol_indices(coordinates,bonds_df):
    center=get_center_of_mass(coordinates)
    center_atom=get_closest_atom_to_center(coordinates,center)
    base_atoms=get_sterimol_base_atoms(center_atom,bonds_df)
    return base_atoms

def filter_atoms_for_sterimol(bonded_atoms_df,coordinates_df):
    """
    a function that filter out NOF bonds and H bonds and returns
     a dataframe of the molecule coordinates without them.
    """
    allowed_bonds_indices= pd.concat([bonded_atoms_df['index_1'],bonded_atoms_df['index_2']],axis=1).reset_index(drop=True)
    atom_filter=adjust_indices(np.unique([atom for sublist in allowed_bonds_indices.values.tolist() for atom in sublist]))
    edited_coordinates_df=coordinates_df.loc[atom_filter].reset_index(drop=True)
   
    return edited_coordinates_df



def get_extended_df_for_sterimol(coordinates_df, bonds_df, radii='CPK'):
    """
    A function that adds information to the regular coordinates_df

    Parameters
    ----------
    coordinates_df : dataframe
    bond_type : str
        The bond type of the molecule
    radii : str, optional
        The type of radii to use ('bondi' or 'CPK'), by default 'bondi'

    Returns
    -------
    dataframe
        The extended dataframe with additional columns

    """
    
    bond_type_map_regular = GeneralConstants.REGULAR_BOND_TYPE.value
    bond_type_map=GeneralConstants.BOND_TYPE.value
    ## if radius is cpk mapping should be done on atype, else on atom
    radii_map = GeneralConstants.CPK_RADII.value if radii == 'CPK' else GeneralConstants.BONDI_RADII.value
    
    df = coordinates_df.copy()  # make a copy of the dataframe to avoid modifying the original
    
    
    
    if radii == 'bondi':
        df['atype']=df['atom']
    else:
        df['atype']=nob_atype(coordinates_df, bonds_df)

    df['magnitude'] = calc_magnitude_from_coordinates_array(df[['x', 'z']].astype(float))
    
    df['radius'] = df['atype'].map(radii_map)
    df['B5'] = df['radius'] + df['magnitude']
    df['L'] = df['y'] + df['radius']
    return df


def get_transfomed_plane_for_sterimol(plane,degree):
    """
    a function that gets a plane and rotates it by a given degree
    in the case of sterimol the plane is the x,z plane.
    Parameters:
    ----------
    plane : np.array
        [x,z] plane of the molecule coordinates.
        example:
            [-0.6868 -0.4964]
    degree : float
    """
  
    cos_deg=np.cos(degree*(np.pi/180))
    sin_deg=np.sin(degree*(np.pi/180))
    rot_matrix=np.array([[cos_deg,-1*sin_deg],[sin_deg,cos_deg]])
    transformed_plane=np.vstack([np.matmul(rot_matrix,row) for row in plane]).round(4)

    
    ## return the inversed rotation matrix to transform the plane back

    return transformed_plane


def calc_B1(transformed_plane,avs,edited_coordinates_df,column_index):
    """
    Parameters
    ----------
    transformed_plane : np.array
        [x,z] plane of the molecule coordinates.
        example:
            [-0.6868 -0.4964]
            [-0.7384 -0.5135]
            [-0.3759 -0.271 ]
            [-1.1046 -0.8966]
            [ 0.6763  0.5885]
    avs : list
        the max & min of the [x,z] columns from the transformed_plane.
        example:[0.6763, -1.1046, 0.5885, -0.8966
                 ]
    edited_coordinates_df : TYPE
        DESCRIPTION.
    column_index : int
        0 or 1 depending- being used for transformed plane.
    """
    
    ## get the index of the min value in the column compared to the avs.min
    idx=np.where(np.isclose(np.abs(transformed_plane[:,column_index]),(avs.min())))[0][0]  ## .round(4)
    # Compute number of points per substituent (assumes evenly divided)
    # print('inside calc_B1')
    # n_total = transformed_plane.shape[0]
    # n_subs = edited_coordinates_df.shape[0]
    # n_points = n_total // n_subs if n_subs > 0 else 1

    # # Print debug info
    # print("calc_B1: transformed_plane.shape =", transformed_plane.shape)
    # print("calc_B1: len(extended_df) =", n_subs)
    # print("calc_B1: n_points per substituent =", n_points)
    
    # # Find first index where the absolute value is close to the minimum of avs
    # idx = np.where(np.isclose(np.abs(transformed_plane[:, column_index]), avs.min()))[0][0]
    # # Map the plane index back to the corresponding DataFrame row
    # index_df = idx // n_points
    # print("calc_B1: raw idx =", idx, "mapped index_df =", index_df)
    # idx=index_df
    if transformed_plane[idx,column_index]<0:
       
        new_idx=np.where(np.isclose(transformed_plane[:,column_index],transformed_plane[:,column_index].min()))[0][0]
        bool_list=np.logical_and(transformed_plane[:,column_index]>=transformed_plane[new_idx,column_index],
                                 transformed_plane[:,column_index]<=transformed_plane[new_idx,column_index]+1)
        
        transformed_plane[:,column_index]=-transformed_plane[:,column_index]
    else:
     
        bool_list=np.logical_and(transformed_plane[:,column_index]>=transformed_plane[idx,column_index]-1,
                                 transformed_plane[:,column_index]<=transformed_plane[idx,column_index])
        
    against,against_loc=[],[]
    B1,B1_loc=[],[]

    ### return the part of the transformed plane that b1 is calculated from
    ### convert it back with the inverse matrix to get the original coordinates of b1 location
    
    for i in range(0,transformed_plane.shape[0]): 
        if bool_list[i]:
            against.append(np.array(transformed_plane[i,column_index]+edited_coordinates_df['radius'].iloc[i]))
            against_loc.append(edited_coordinates_df['L'].iloc[i])
        

        if len(against)>0:
        
            B1.append(max(against))
            B1_loc.append(against_loc[against.index(max(against))])
          
        else:
       
            B1.append(np.abs(transformed_plane[idx,column_index]+edited_coordinates_df['radius'].iloc[idx]))
            B1_loc.append(edited_coordinates_df['radius'].iloc[idx])
            
            
      
    return [B1,B1_loc] 

def generate_circle(center_x, center_y, radius, n_points=20):
    """
    Generate circle coordinates given a center and radius.
    Returns a DataFrame with columns 'x' and 'y'.
    """
    theta = np.linspace(0, 2 * np.pi, n_points)
    x = center_x + radius * np.cos(theta)
    y = center_y + radius * np.sin(theta)
    return np.column_stack((x, y))


def b1s_for_loop_function(extended_df, b1s, b1s_loc, degree_list, plane, b1_planes):
    """
    For each degree in degree_list, rotate the plane and compute B1 values and the B1-B5 angle.
    Instead of returning after the first iteration, this version accumulates all results
    into a DataFrame.
    
    Parameters
    ----------
    extended_df : pd.DataFrame
        DataFrame with at least columns 'x', 'z', 'radius', 'L'.
    b1s : list
        (Unused here; kept for compatibility)
    b1s_loc : list
        (Unused here; kept for compatibility)
    degree_list : list
        List of rotation angles (in degrees) to scan.
    plane : np.array
        Array of shape (n_points_total, 2) that contains the combined circle points.
    b1_planes : list
        List to store the rotated plane from each iteration.
        
    Returns
    -------
    pd.DataFrame
        DataFrame with one row per degree, containing:
         - 'degree': the rotation angle,
         - 'B1': the minimum extreme value from the rotated plane,
         - 'B1_B5_angle': the angle between the B1 and B5 arrows in degrees,
         - 'b1_coords': the coordinate (as a tuple) chosen for B1,
         - 'b5_value': the distance of the farthest point (B5) from the origin.
    """
    results = []  # List to accumulate results
    
    for degree in degree_list:
        transformed_plane = get_transfomed_plane_for_sterimol(plane, degree)  # Rotate the plane
     
        max_x = np.max(transformed_plane[:, 0])
        min_x = np.min(transformed_plane[:, 0])
        max_y = np.max(transformed_plane[:, 1])
        min_y = np.min(transformed_plane[:, 1])
        avs = np.abs([max_x, min_x, max_y, min_y])
       
        min_val = np.min(avs)
        min_index = np.argmin(avs)
        
        # Mimic R's switch to pick the B1 coordinate:
        if min_index == 0:
            b1_coords = (max_x, 0)
        elif min_index == 1:
            b1_coords = (min_x, 0)
        elif min_index == 2:
            b1_coords = (0, max_y)
        else:
            b1_coords = (0, min_y)
        
        # Determine B5 as the farthest point from the origin.
        norms_sq = np.sum(transformed_plane**2, axis=1)
        b5_index = np.argmax(norms_sq)
        b5_point = transformed_plane[b5_index]
        b5_value = np.linalg.norm(b5_point)
        
        # Calculate angles for the arrows.
        angle_b1 = np.arctan2(b1_coords[1], b1_coords[0]) % (2*np.pi)
        angle_b5 = np.arctan2(b5_point[1], b5_point[0]) % (2*np.pi)
        angle_diff = abs(angle_b5 - angle_b1)
        if angle_diff > np.pi:
            angle_diff = 2*np.pi - angle_diff
        # Convert the angle difference to degrees.
        B1_B5 = np.degrees(angle_diff)
        B1 = min_val
        
        # Save the transformed plane for this iteration.
        b1_planes.append(transformed_plane)
        
        # Accumulate the result for this degree.
        results.append({
            'degree': degree,
            'B1': B1,
            'B1_B5_angle': B1_B5,
            'b1_coords': b1_coords,
            'b5_value': b5_value,
            'plane': transformed_plane
        })
    
    # Create a DataFrame from the accumulated results.
    sterimol_df = pd.DataFrame(results)
  
    return sterimol_df

   
    
def get_b1s_list(extended_df, scans=1,plot_result=False):
    """
    Calculate B1 values by scanning over a range of rotation angles.
    Instead of using only the center points, this version generates circle points
    for each substituent from extended_df and then stacks them into one plane.
    
    Parameters
    ----------
    extended_df : pd.DataFrame
        DataFrame with at least the columns: 'x', 'z', 'radius', 'L'.
    scans : int, optional
        Degree step for the initial scan.
        
    Returns
    -------
    tuple
        (np.array of B1 values, np.array of B1 location values, list of rotated planes)
    """
    b1s, b1s_loc, b1_planes = [], [], []
    degree_list = list(range(18, 108, scans))
    
    # Generate circles for each substituent and combine them.
    circles = []
    for idx, row in extended_df.iterrows():
        circle_points = generate_circle(row['x'], row['z'], row['radius'], n_points=100)
        circles.append(circle_points)
    plane = np.vstack(circles)  # All circle points combined.
    
    sterimol_df=b1s_for_loop_function(extended_df, b1s, b1s_loc, degree_list, plane, b1_planes)

    b1s=sterimol_df['B1']
    
    try:
        back_ang=degree_list[np.where(b1s==min(b1s))[0][0]]-scans   
        front_ang=degree_list[np.where(b1s==min(b1s))[0][0]]+scans
        new_degree_list=range(back_ang,front_ang+1)
    except:
        
        back_ang=degree_list[np.where(np.isclose(b1s, min(b1s), atol=1e-8))[0][0]]-scans
        front_ang=degree_list[np.where(np.isclose(b1s, min(b1s), atol=1e-8))[0][0]]+scans
        new_degree_list=range(back_ang,front_ang+1)
    # plane=sterimol_df['plane']
    
    b1s, b1s_loc, b1_planes = [], [], []
    sterimol_df=b1s_for_loop_function(extended_df, b1s, b1s_loc, list(new_degree_list), plane, b1_planes)
  
    b1s=sterimol_df['B1']

    b1_b5_angle=sterimol_df['B1_B5_angle']
    plane=sterimol_df['plane']
    

    return [b1s, b1_b5_angle, plane]

import matplotlib.pyplot as plt

def plot_b1_visualization(rotated_plane, extended_df, n_points=100, title="Rotated Plane Visualization"):
    """
    Visualize the rotated plane by plotting:
      - Complete circles (each generated from a substituent),
      - Dashed lines at extreme x and y values,
      - Arrows for the four extreme directions (with the B1 arrow highlighted),
      - A B5 arrow (the farthest point from the origin),
      - An arc indicating the angle between the B1 and B5 arrows.
    
    Parameters
    ----------
    rotated_plane : np.array
        Rotated plane points (stacked complete circles; shape: [n_total_points, 2]).
    extended_df : pd.DataFrame
        DataFrame with columns 'radius' and 'L' (used for annotations).
    n_points : int, optional
        Number of points per circle (default is 20).
    title : str, optional
        Title for the plot.
    """
    # Compute extreme values from all points
    max_x = np.max(rotated_plane[:, 0])
    min_x = np.min(rotated_plane[:, 0])
    max_y = np.max(rotated_plane[:, 1])
    min_y = np.min(rotated_plane[:, 1])
    avs = np.abs([max_x, min_x, max_y, min_y])
    min_val = np.min(avs)
    min_index = np.argmin(avs)
    
    # Determine B1 arrow coordinates based on the minimum extreme
    if min_index == 0:
        b1_coords = np.array([max_x, 0])
    elif min_index == 1:
        b1_coords = np.array([min_x, 0])
    elif min_index == 2:
        b1_coords = np.array([0, max_y])
    else:
        b1_coords = np.array([0, min_y])
    
    # Determine B5 as the farthest point from the origin
    norms_sq = np.sum(rotated_plane**2, axis=1)
    b5_index = np.argmax(norms_sq)
    b5_point = rotated_plane[b5_index]
    b5_value = np.linalg.norm(b5_point)
  
    # Calculate angles for the arrows
    angle_b1 = np.arctan2(b1_coords[1], b1_coords[0]) % (2 * np.pi)
    angle_b5 = np.arctan2(b5_point[1], b5_point[0]) % (2 * np.pi)
    angle_diff = abs(angle_b5 - angle_b1)
    if angle_diff > np.pi:
        angle_diff = 2 * np.pi - angle_diff
    angle_diff_deg = np.degrees(angle_diff)
    
    plt.figure(figsize=(8, 8))
    
    # Plot complete circles.
    n_total = rotated_plane.shape[0]
    n_circles = n_total // n_points
    for i in range(n_circles):
        circle_points = rotated_plane[i * n_points:(i + 1) * n_points, :]
        # Close the circle by appending the first point to the end
        circle_points = np.vstack([circle_points, circle_points[0]])
        plt.plot(circle_points[:, 0], circle_points[:, 1], color='cadetblue', linewidth=1.5)
    
    # Plot dashed extreme lines
    plt.axvline(x=max_x, color='darkred', linestyle='dashed')
    plt.axvline(x=min_x, color='darkred', linestyle='dashed')
    plt.axhline(y=max_y, color='darkgreen', linestyle='dashed')
    plt.axhline(y=min_y, color='darkgreen', linestyle='dashed')
    
    # Draw arrows for each extreme (all black except the B1 arrow highlighted)
    arrow_colors = ['black'] * 4
    arrow_colors[min_index] = '#8FBC8F'
    plt.arrow(0, 0, max_x, 0, head_width=0.1, length_includes_head=True, color=arrow_colors[0])
    plt.arrow(0, 0, min_x, 0, head_width=0.1, length_includes_head=True, color=arrow_colors[1])
    plt.arrow(0, 0, 0, max_y, head_width=0.1, length_includes_head=True, color=arrow_colors[2])
    plt.arrow(0, 0, 0, min_y, head_width=0.1, length_includes_head=True, color=arrow_colors[3])
    
    # Draw the B5 arrow in red
    plt.arrow(0, 0, b5_point[0], b5_point[1], head_width=0.1, length_includes_head=True, color="#CD3333")
    
    # Annotate B1 and B5 values
    plt.text(b1_coords[0] * 0.5, b1_coords[1] * 0.5, f"B1\n{min_val:.2f}", 
             fontsize=12, ha='center', va='bottom', fontweight='bold')
    plt.text(b5_point[0] * 0.66, b5_point[1] * 0.66, f"B5\n{b5_value:.2f}", 
             fontsize=12, ha='center', va='bottom', fontweight='bold')
    
    # Draw an arc between the B1 and B5 arrows to represent the angle difference.
    arc_theta = np.linspace(min(angle_b1, angle_b5), max(angle_b1, angle_b5), 100)
    arc_x = 0.5 * np.cos(arc_theta)
    arc_y = 0.5 * np.sin(arc_theta)
    plt.plot(arc_x, arc_y, color='gray', linewidth=1.5)
    
    # Annotate the angle in degrees at the midpoint of the arc.
    mid_angle = (min(angle_b1, angle_b5) + max(angle_b1, angle_b5)) / 2
    plt.text(0.8 * np.cos(mid_angle), 0.8 * np.sin(mid_angle), f"{angle_diff_deg:.1f}°",
             fontsize=12, ha='center', va='center', fontweight='bold')
    
    plt.title(title)
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.axis('equal')
    plt.show()



# def get_b1s_list(extended_df, scans=90//5):
    
#     b1s,b1s_loc,b1_planes=[],[],[]
#     scans=scans
#     degree_list=list(range(18,108,scans))
#     plane=np.array(extended_df[['x','z']].astype(float))
#     b1s_for_loop_function(extended_df, b1s, b1s_loc, degree_list, plane,b1_planes)

#     if b1s:
#         try:
#             back_ang=degree_list[np.where(b1s==min(b1s))[0][0]]-scans   
#             front_ang=degree_list[np.where(b1s==min(b1s))[0][0]]+scans
#             degree_list=range(back_ang,front_ang+1)
#         except:
            
#             back_ang=degree_list[np.where(np.isclose(b1s, min(b1s), atol=1e-8))[0][0]]-scans
#             front_ang=degree_list[np.where(np.isclose(b1s, min(b1s), atol=1e-8))[0][0]]+scans
#             degree_list=range(back_ang,front_ang+1)
#     else:
     
#         return [np.array(b1s),np.array(b1s_loc)]
    
#     b1s,b1s_loc,b1_planes=[],[] ,[]
#     b1s_for_loop_function(extended_df, b1s, b1s_loc, degree_list, plane,b1_planes)
   
#     return [np.array(b1s),np.array(b1s_loc),b1_planes]

def calc_sterimol(bonded_atoms_df,extended_df,visualize=False):
    edited_coordinates_df=filter_atoms_for_sterimol(bonded_atoms_df,extended_df)
  
    b1s,b1_b5_angle,plane=get_b1s_list(edited_coordinates_df)
   
    valid_indices = np.where(b1s >= 0)[0]
    best_idx = valid_indices[np.argmin(b1s[valid_indices])]
    best_b1_plane = plane[best_idx]

    max_x = np.max(best_b1_plane[:, 0])
    min_x = np.min(best_b1_plane[:, 0])
    max_y = np.max(best_b1_plane[:, 1])
    min_y = np.min(best_b1_plane[:, 1])
    avs = np.abs([max_x, min_x, max_y, min_y])
    min_val = np.min(avs)
    B1= min_val
    b1_index=np.where(b1s==B1)[0][0]
    angle=b1_b5_angle[b1_index]
    norms_sq = np.sum(best_b1_plane**2, axis=1)
    b5_index = np.argmax(norms_sq)
    b5_point = best_b1_plane[b5_index]
    b5_value = np.linalg.norm(b5_point)
    # loc_B1=max(b1s_loc[np.where(b1s[b1s>=0]==min(b1s[b1s>=0]))])
    # get the idx of the row with the  biggest b5 value from edited
    max_row=edited_coordinates_df['B5'].idxmax()
    max_row=int(max_row)
    L=max(edited_coordinates_df['L'].values)
    loc_B5 = edited_coordinates_df['y'].iloc[np.where(edited_coordinates_df['B5']==max(edited_coordinates_df['B5']))[0][0]]
    
    sterimol_df = pd.DataFrame([B1, b5_value, L ,loc_B5,angle], index=['B1', 'B5', 'L','loc_B5','B1_B5_angle'])
    if visualize:
        plot_b1_visualization(best_b1_plane, edited_coordinates_df)

    return sterimol_df.T 


def get_sterimol_df(coordinates_df, bonds_df, base_atoms,connected_from_direction, radii='bondi', sub_structure=True, drop_atoms=None,visualize=False):

    if drop_atoms is not None:
        drop_atoms=adjust_indices(drop_atoms)
        for atom in drop_atoms:
            bonds_df = bonds_df[~((bonds_df[0] == atom) | (bonds_df[1] == atom))]
            ## drop the rows from coordinates_df
            coordinates_df = coordinates_df.drop(atom)

    bonds_direction = direction_atoms_for_sterimol(bonds_df, base_atoms)
    
    new_coordinates_df = preform_coordination_transformation(coordinates_df, bonds_direction)


    if sub_structure:
        if connected_from_direction is None:
            connected_from_direction = get_molecule_connections(bonds_df, base_atoms[0], base_atoms[1])
        else:
            connected_from_direction = connected_from_direction
    else:
        connected_from_direction = None
    
    
    bonded_atoms_df = get_specific_bonded_atoms_df(bonds_df, connected_from_direction, new_coordinates_df)
    
    extended_df = get_extended_df_for_sterimol(new_coordinates_df, bonds_df, radii)
    
    ###calculations
    
    sterimol_df = calc_sterimol(bonded_atoms_df, extended_df,visualize)
    sterimol_df= sterimol_df.rename(index={0: str(base_atoms[0]) + '-' + str(base_atoms[1])})
   
    sterimol_df = sterimol_df.round(4)
    
    return sterimol_df



def calc_magnitude_from_coordinates_array(coordinates_array: npt.ArrayLike) -> List[float]:
    """
    Calculates the magnitudes of each row in the given coordinates array.

    Parameters
    ----------
    coordinates_array: np.ndarray
        A nx3 array representing the x, y, z coordinates of n atoms.

    Returns
    -------
    magnitudes: List[float]
        A list of magnitudes corresponding to the rows of the input array.

    """
    magnitude = np.linalg.norm(coordinates_array, axis=1)
    return magnitude



class Molecule():

    def __init__(self, molecule_xyz_filename):

   
        self.molecule_name = molecule_xyz_filename.split('.')[0]
        self.molecule_path = os.path.dirname(os.path.abspath(molecule_xyz_filename))
        os.chdir(self.molecule_path)
        
        
        self.xyz_df = get_df_from_file(molecule_xyz_filename)

        with open(molecule_xyz_filename, 'r') as _f:
            _comment = _f.readlines()[1].strip()
        try:
            self.energy = float(_comment)
        except (ValueError, IndexError):
            self.energy = np.nan

        self.coordinates_array = np.array(self.xyz_df[['x', 'y', 'z']].astype(float))
        self.bonds_df = extract_connectivity(self.xyz_df)
        
        
            



    def process_sterimol_atom_group(self, atoms, radii, sub_structure=True, drop_atoms=None,visualize=False) -> pd.DataFrame:

        connected = get_molecule_connections(self.bonds_df, atoms[0], atoms[1])
        
        return get_sterimol_df(self.xyz_df, self.bonds_df, atoms, connected, radii, sub_structure=sub_structure, drop_atoms=drop_atoms, visualize=visualize)

    def get_sterimol(self, base_atoms: Union[None, Tuple[int, int]] = None, radii: str = 'bondi',sub_structure=True, drop_atoms=None,visualize=False) -> pd.DataFrame:
        """
        Returns a DataFrame with the Sterimol parameters calculated based on the specified base atoms and radii.

        Args:
            base_atoms (Union[None, Tuple[int, int]], optional): The indices of the base atoms to use for the Sterimol calculation. Defaults to None.
            radii (str, optional): The radii to use for the Sterimol calculation. Defaults to 'bondi'.

        Returns:
            pd.DataFrame: A DataFrame with the Sterimol parameters.
            
            to add
            - only_sub- sterimol of only one part.
            - drop some atoms.
        """
        if base_atoms is None:
            base_atoms = get_sterimol_indices(self.xyz_df, self.bonds_df)
        
        if isinstance(base_atoms[0], list):
            # If base_atoms is a list of lists, process each group individually and concatenate the results
            sterimol_list = [self.process_sterimol_atom_group(atoms, radii, sub_structure=sub_structure, drop_atoms=drop_atoms,visualize=visualize) for atoms in base_atoms]
            sterimol_df = pd.concat(sterimol_list, axis=0)

        else:
            # If base_atoms is a single group, just process that group
            sterimol_df = self.process_sterimol_atom_group(base_atoms, radii,sub_structure=sub_structure, drop_atoms=drop_atoms,visualize=visualize)
        return sterimol_df


    def swap_atom_pair(self, pair_indices: Tuple[int, int]) -> pd.DataFrame:
        """
        Swaps the positions of two atoms in the molecule and returns a new DataFrame with the updated coordinates.

        Args:
            pair_indices (Tuple[int, int]): The indices of the atoms to swap.

        Returns:
            pd.DataFrame: A new DataFrame with the updated coordinates.
        """
        pairs = adjust_indices(pair_indices)
        xyz_df = self.xyz_df
        temp = xyz_df.iloc[pairs[0]].copy()
        xyz_df.iloc[pairs[0]] = self.coordinates_array[pairs[1]]
        xyz_df.iloc[pairs[1]] = temp
        return xyz_df
  
 
    def get_coordination_transformation_df(self, base_atoms_indices: List[int]) -> pd.DataFrame:
        """
        Returns a new DataFrame with the coordinates transformed based on the specified base atoms.

        Args:
            base_atoms_indices (List[int]): The indices of the base atoms to use for the transformation.

        Returns:
            pd.DataFrame: A new DataFrame with the transformed coordinates.
        """
        new_coordinates_df = preform_coordination_transformation(self.xyz_df, base_atoms_indices)
        return new_coordinates_df

    ## not working after renumbering for some reason
    
    

    def get_bond_angle(self, atom_indices: List[int]) -> pd.DataFrame:
        """
        Returns a DataFrame with the bond angles calculated based on the specified atom indices.

        Args:
            atom_indices (List[int]): The indices of the atoms to use for the bond angle calculation.

        Returns:
            pd.DataFrame: A DataFrame with the bond angles.
        """
        return get_angle_df(self.coordinates_array, atom_indices)
    
    def get_bond_length_single(self, atom_pair):
        bond_length = calc_single_bond_length(self.coordinates_array, atom_pair)
        bond_length_df = pd.DataFrame([bond_length], index=[f'bond_length_{atom_pair[0]}-{atom_pair[1]}'])
        return bond_length_df

    def get_bond_length(self, atom_pairs):
        """
            Returns a DataFrame with the bond lengths calculated based on the specified atom pairs.

            Args:
                atom_pairs (Union[List[Tuple[int, int]], Tuple[int, int]]): The pairs of atoms to use for the bond length calculation.

            Returns:
                pd.DataFrame: A DataFrame with the bond lengths.
            """
        if isinstance(atom_pairs[0], list):
            # If atom_pairs is a list of lists, process each pair individually and concatenate the results
            bond_length_list = [self.get_bond_length_single(pair) for pair in atom_pairs]
            bond_df = pd.concat(bond_length_list, axis=0)
        else:
            # If atom_pairs is a single pair, just process that pair
            bond_df = self.get_bond_length_single(atom_pairs)
        return bond_df

    def get_buried_volume(
        self,
        metal_index: int,
        radii: str = 'bondi',
        radius: float = 3.5,
        radii_scale: float = 1.17,
        include_hs: bool = False,
        z_axis_atoms=None,
        xz_plane_atoms=None,
    ) -> pd.DataFrame:
        """
        Calculate the buried volume around a given metal/anchor atom using
        morfeus BuriedVolume.

        Args:
            metal_index (int): 1-based index of the atom at the sphere centre
                               (typically the metal or key anchor atom).
            radii (str): VdW radii set to use ('bondi' or 'alvarez').
            radius (float): Sphere radius in Å (default 3.5).
            radii_scale (float): Scaling factor for VdW radii (default 1.17).
            include_hs (bool): Whether to include H atoms (default False).
            z_axis_atoms: Atom index (or list of indices) that define the
                          z-axis direction (optional).
            xz_plane_atoms: Atom index (or list of indices) that define the
                            xz-plane (optional).

        Returns:
            pd.DataFrame: Single-row DataFrame with columns:
                percent_buried_volume, fraction_buried_volume,
                buried_volume, free_volume.
        """
        elements = self.xyz_df['atom'].tolist()
        bv = BuriedVolume(
            elements=elements,
            coordinates=self.coordinates_array,
            metal_index=metal_index,
            radii_type=radii,
            radius=radius,
            radii_scale=radii_scale,
            include_hs=include_hs,
            z_axis_atoms=z_axis_atoms,
            xz_plane_atoms=xz_plane_atoms,
        )
        return pd.DataFrame([{
            'percent_buried_volume':   bv.percent_buried_volume,
            'fraction_buried_volume':  bv.fraction_buried_volume,
            'buried_volume':           bv.buried_volume,
            'free_volume':             bv.free_volume,
        }])


def dict_to_horizontal_df(data_dict):
    # Initialize an empty DataFrame to store the transformed data
    df_transformed = pd.DataFrame()
    # Loop through each key-value pair in the original dictionary
    for mol, df in data_dict.items():
        transformed_data = {}
        print(df.head())    
        # Loop through each row and column in the DataFrame
        for index, row in df.iterrows():
            
            index_words = set(index.split('_'))
            for col in df.columns:
                # Create a new key using the format: col_index
                try:
                    col_words = set(col.split('_'))
                except:
                    col_words = []
                 # Check if the index and the column have the same words and remove one
                common_words = index_words.intersection(col_words)
                if col != 0 and '0':
                    if common_words:
                        unique_col_words = col_words - common_words
                        unique_index_words = index_words - common_words
                        new_key_parts = ['_'.join(common_words)] if common_words else []
                        new_key_parts.extend([part for part in ['_'.join(unique_col_words), '_'.join(unique_index_words)] if part])
                        new_key = '_'.join(new_key_parts)
                    else:
                        new_key = f"{col}_{index}"
                else:
                    new_key = f"{index}"
                # Store the corresponding value in the transformed_data dictionary
                transformed_data[new_key] = row[col]
        # Convert the dictionary into a DataFrame row with the molecule name as the index
        df_row = pd.DataFrame([transformed_data], index=[mol])
        # Append the row to df_transformed
        df_transformed = pd.concat([df_transformed, df_row], ignore_index=False)
    return df_transformed




class Molecules_xyz():
    
    def __init__(self,molecules_dir_name, renumber=False):
        self.molecules_path=os.path.abspath(molecules_dir_name)
        os.chdir(self.molecules_path) 
        self.molecules=[]
        self.failed_molecules=[]
        for file in os.listdir(): 
            if file.endswith('.xyz'):
                try:
                    self.molecules.append(Molecule(file))
                except:
                    self.failed_molecules.append(file)
                    print(f'Error: {file} could not be processed')
                    failed_file = file.rsplit('.feather', 1)[0] + '.feather_fail'
                # self.molecules.append(Molecule(log_file))
       
    


    def filter_molecules(self, indices):
        self.molecules = [self.molecules[i] for i in indices]
        self.molecules_names = [self.molecules_names[i] for i in indices]

    def get_energy_df(self) -> pd.DataFrame:
        """
        Returns a DataFrame with the energy for each molecule.

        Energy is read from the comment line (line 2) of each .xyz file.
        Molecules with an empty or non-numeric comment line get NaN.

        Returns:
            pd.DataFrame: index = molecule names, column = 'energy'.
        """
        return pd.DataFrame(
            {'energy': {mol.molecule_name: mol.energy for mol in self.molecules}}
        )

    def get_sterimol_dict(self,atom_indices,radii='CPK'):
        sterimol_dict={}
        for molecule in self.molecules:
            try:
                sterimol_dict[molecule.molecule_name]=molecule.get_sterimol(atom_indices,radii)
                print(f'calculated sterimol for {molecule.molecule_name}')
            except ValueError as e:
                print(f'Error: {e}')
                print(f'failed to calculate for {molecule.molecule_name}')
                sterimol_dict[molecule.molecule_name]=np.nan

        return sterimol_dict
   
    def get_sterimol_df(self,atom_indices,radii='CPK'):
        sterimol_dict=self.get_sterimol_dict(atom_indices,radii)
        sterimol_df=dict_to_horizontal_df(sterimol_dict)
        return sterimol_df

    # ── Angles & dihedrals ─────────────────────────────────────────────────

    def get_angles_df(self, atom_indices_list) -> pd.DataFrame:
        """
        Calculate bond angles and/or dihedral angles for all molecules.

        Pass 3-atom groups for bond angles, 4-atom groups for dihedrals,
        or a mix of both.  All indices are 1-based.

        Args:
            atom_indices_list: a single group (list of ints) or a list of
                groups, e.g. [[1,2,3], [1,2,3,4]].

        Returns:
            pd.DataFrame: molecules × angle/dihedral columns (degrees).
                Column names follow the existing convention:
                ``angle_[i,j,k]`` or ``dihedral_[i,j,k,l]``.
        """
        rows = {}
        for molecule in self.molecules:
            try:
                df = molecule.get_bond_angle(atom_indices_list)
                rows[molecule.molecule_name] = df[0].to_dict()
                print(f'calculated angles for {molecule.molecule_name}')
            except Exception as e:
                print(f'Error ({molecule.molecule_name}): {e}')
                rows[molecule.molecule_name] = np.nan
        return pd.DataFrame(rows).T

    # ── Bond lengths ───────────────────────────────────────────────────────

    def get_bond_lengths_df(self, atom_pairs) -> pd.DataFrame:
        """
        Calculate bond lengths for all molecules.

        Args:
            atom_pairs: a single pair (list of 2 ints) or a list of pairs,
                e.g. [[1,2], [3,4]].  Indices are 1-based.

        Returns:
            pd.DataFrame: molecules × bond-length columns (Å).
                Column names: ``bond_length_i-j``.
        """
        rows = {}
        for molecule in self.molecules:
            try:
                df = molecule.get_bond_length(atom_pairs)
                rows[molecule.molecule_name] = df[0].to_dict()
                print(f'calculated bond lengths for {molecule.molecule_name}')
            except Exception as e:
                print(f'Error ({molecule.molecule_name}): {e}')
                rows[molecule.molecule_name] = np.nan
        return pd.DataFrame(rows).T

    # ── Buried volume ──────────────────────────────────────────────────────

    def get_buried_volume_df(
        self,
        metal_index: int,
        radii: str = 'bondi',
        radius: float = 3.5,
        radii_scale: float = 1.17,
        include_hs: bool = False,
        z_axis_atoms=None,
        xz_plane_atoms=None,
    ) -> pd.DataFrame:
        """
        Calculate buried volume around a given atom for all molecules.

        Args:
            metal_index (int): 1-based index of the sphere-centre atom
                               (metal, P, N, …).
            radii (str): VdW radii set ('bondi' or 'alvarez').
            radius (float): Sphere radius in Å (default 3.5).
            radii_scale (float): VdW scaling factor (default 1.17).
            include_hs (bool): Include H atoms in the calculation
                               (default False).
            z_axis_atoms: Atom index (or list) defining the z-axis
                          (optional, for quadrant/octant analysis).
            xz_plane_atoms: Atom index (or list) defining the xz-plane
                            (optional).

        Returns:
            pd.DataFrame: molecules × {percent_buried_volume,
                fraction_buried_volume, buried_volume, free_volume}.
        """
        rows = {}
        for molecule in self.molecules:
            try:
                df = molecule.get_buried_volume(
                    metal_index=metal_index,
                    radii=radii,
                    radius=radius,
                    radii_scale=radii_scale,
                    include_hs=include_hs,
                    z_axis_atoms=z_axis_atoms,
                    xz_plane_atoms=xz_plane_atoms,
                )
                rows[molecule.molecule_name] = df.iloc[0].to_dict()
                print(f'calculated buried volume for {molecule.molecule_name}')
            except Exception as e:
                print(f'Error ({molecule.molecule_name}): {e}')
                rows[molecule.molecule_name] = np.nan
        return pd.DataFrame(rows).T

# =============================================================================
# CREST conformer ensembles -> Boltzmann-averaged Sterimol
# =============================================================================
#
# Everything above this line works on ONE structure per .xyz file. The classes
# below extend that to CREST's native output: a directory per ligand containing
# a multi-structure conformer ensemble (energy in Hartree on the comment line
# of every block), which is parsed, Sterimol'd conformer-by-conformer, and
# combined into a single Boltzmann-weighted value per ligand.
#
# Layout expected per ligand (standard `crest` output naming):
#
#   L1/
#     L1_conformers.xyz   <- de-duplicated unique conformers (used by default)
#     L1_rotamers.xyz     <- fallback: includes symmetry-equivalent duplicates
#     L1_best.xyz         <- fallback: single lowest-energy structure only
#     raw_crest_run/crest_conformers.xyz  <- fallback: raw CREST working dir
#
# Usage (a directory of ligand subfolders):
#
#   from sterimol_standalone import CrestLigandSet
#   dataset = CrestLigandSet(r"C:\...\doyle_conformers")
#   summary_df, detail_df = dataset.get_sterimol_df(
#       base_atoms=[1, 2],      # 1-indexed anchor atoms; omit for auto-detection
#       radii='CPK',
#       save_dir=r"C:\...\doyle_conformers\sterimol_results",
#   )
#
# Usage (a single ligand's ensemble file):
#
#   from sterimol_standalone import ConformerEnsemble
#   ens = ConformerEnsemble(ligand_dir=r"C:\...\doyle_conformers\L1")
#   per_conformer_df, boltzmann_row = ens.get_sterimol(base_atoms=[1, 2])
#
# =============================================================================

HARTREE_TO_KCAL = 627.5094740631  # CODATA, matches CREST's own conversion
_KB_KCAL = 1.9872041e-3           # gas constant, kcal / (mol K)


def boltzmann_weights(energies_hartree, temperature=298.15):
    """
    Convert a list of absolute energies (Hartree, e.g. straight off a CREST/xtb
    ensemble comment line) into normalized Boltzmann weights at `temperature` (K).

    Uses a log-sum-exp shift for numerical stability, and propagates NaN energies
    as NaN weights (rather than silently treating them as 0 or crashing) so the
    caller can decide how to handle/renormalize around missing values.

    Parameters
    ----------
    energies_hartree : array-like
        Absolute electronic (or free) energies, Hartree.
    temperature : float, default 298.15
        Temperature in Kelvin.

    Returns
    -------
    np.ndarray
        Weights summing to 1 over the valid (non-NaN) entries; NaN energies
        get a NaN weight.
    """
    energies = np.asarray(energies_hartree, dtype=float)
    weights = np.full(energies.shape, np.nan)
    valid = ~np.isnan(energies)
    if not valid.any():
        return weights
    rel_kcal = (energies[valid] - energies[valid].min()) * HARTREE_TO_KCAL
    exponent = -rel_kcal / (_KB_KCAL * temperature)
    exponent = exponent - exponent.max()
    w = np.exp(exponent)
    w = w / w.sum()
    weights[valid] = w
    return weights


def _energy_from_xyz_comment(comment, convention='first'):
    """Parse a Hartree energy from an XYZ comment line.

    ``first`` is the CREST convention (energy is the first token). ``last`` is
    GOAT-safe: ORCA GOAT may write an RMSD before the energy, and taking the
    first float then Boltzmann-weights the RMSD. ``auto`` takes the last float
    whose absolute value is at least 1 Eh when one exists.
    """
    if not comment:
        return np.nan
    floats = []
    for tok in comment.replace('=', ' ').split():
        try:
            floats.append(float(tok))
        except ValueError:
            continue
    if not floats:
        return np.nan
    conv = (convention or 'first').lower()
    if conv == 'first':
        return floats[0]
    if conv == 'last':
        return floats[-1]
    if conv == 'auto':
        abs_eh = [v for v in floats if abs(v) >= 1.0]
        return abs_eh[-1] if abs_eh else floats[-1]
    raise ValueError(f"Unknown energy convention {convention!r}")


def parse_xyz_ensemble(filepath, energy_convention='first'):
    """
    Parse a multi-structure XYZ file (CREST/xtb `*_conformers.xyz`,
    `*_rotamers.xyz`, GOAT `*.ens.xyz` / `*.finalensemble.xyz`, or a
    single-structure `*_best.xyz`) into a list of conformer blocks.

    Each block's comment line (line 2) is interpreted as the absolute energy
    in Hartree. Default ``energy_convention='first'`` is the CREST convention.
    Pass ``'last'`` for GOAT files that may write an RMSD before the energy.

    Returns
    -------
    list[dict]
        One dict per structure: {'energy': float, 'xyz_df': DataFrame[atom,x,y,z]}

    Raises
    ------
    ValueError
        If a block header isn't a valid atom count, or a block is truncated
        (fewer coordinate lines than declared) -- both indicate a corrupted
        ensemble file that should not be silently averaged over.
    """
    with open(filepath, 'r') as f:
        lines = f.read().splitlines()
    while lines and not lines[-1].strip():
        lines.pop()

    blocks = []
    i = 0
    n_lines = len(lines)
    while i < n_lines:
        header = lines[i].strip()
        if not header:
            i += 1
            continue
        try:
            natoms = int(header.split()[0])
        except (ValueError, IndexError):
            raise ValueError(f"Malformed XYZ block header at line {i + 1} in {filepath}: '{lines[i]}'")

        comment = lines[i + 1].strip() if i + 1 < n_lines else ''
        energy = _energy_from_xyz_comment(comment, convention=energy_convention)

        atom_lines = lines[i + 2:i + 2 + natoms]
        rows = [ln.split() for ln in atom_lines if len(ln.split()) >= 4]
        if len(rows) != natoms:
            raise ValueError(
                f"Block starting at line {i + 1} in {filepath} declares {natoms} atoms "
                f"but {len(rows)} coordinate line(s) were parsed -- truncated/corrupted ensemble."
            )

        xyz_df = pd.DataFrame([r[:4] for r in rows], columns=['atom', 'x', 'y', 'z'])
        xyz_df[['x', 'y', 'z']] = xyz_df[['x', 'y', 'z']].astype(float)
        blocks.append({'energy': energy, 'xyz_df': xyz_df})
        i += 2 + natoms

    if not blocks:
        raise ValueError(f"No structures parsed from {filepath}")
    return blocks


class _ConformerFrame:
    """
    Lightweight in-memory analogue of `Molecule`, built directly from an
    already-parsed XYZ block (atom/x/y/z DataFrame + energy) instead of
    reading a file from disk. Mirrors the attributes/methods `Molecule`
    exposes so the existing `get_sterimol_df` machinery works unchanged
    -- deliberately kept separate from `Molecule` so looping over hundreds
    of conformers never has to `os.chdir()` per-structure.
    """

    def __init__(self, xyz_df, energy, molecule_name, threshold=None):
        self.molecule_name = molecule_name
        self.energy = energy
        self.xyz_df = xyz_df.reset_index(drop=True)
        self.coordinates_array = np.array(self.xyz_df[['x', 'y', 'z']].astype(float))
        self.bonds_df = extract_connectivity(self.xyz_df, threshold_distance=threshold)

    def process_sterimol_atom_group(self, atoms, radii, sub_structure=True, drop_atoms=None, visualize=False):
        # Only pay for `get_molecule_connections` (an igraph all-simple-paths
        # enumeration that can blow up on larger/cyclic ligands) when the
        # substructure trim is actually requested. `Molecule` upstream always
        # computes this eagerly; for whole free-ligand ensembles (this class's
        # use case) that made a 28-atom, 138-conformer ensemble hang for
        # minutes on paths nobody asked for.
        connected = get_molecule_connections(self.bonds_df, atoms[0], atoms[1]) if sub_structure else None
        return get_sterimol_df(self.xyz_df, self.bonds_df, atoms, connected, radii,
                                sub_structure=sub_structure, drop_atoms=drop_atoms, visualize=visualize)

    def _auto_base_atom_candidates(self):
        """
        Candidate [origin, direction] pairs for auto-detected Sterimol, most
        preferred first. `get_sterimol_indices` always returns the SAME single
        pair (center-of-mass atom + whichever neighbor happens to sort first
        in `bonds_df`, which is often a hydrogen) -- fine most of the time,
        but on some geometries that particular direction is exactly collinear
        with the coplane vector `calc_basis_vector` builds against, which
        makes its `np.linalg.lstsq` singular ("SVD did not converge"). Seen
        on ligand L25, where the default pick was a stray C-H bond.

        This tries the standard pick first (so normal behavior/results don't
        change), then every other neighbor of that same center atom -- heavy
        atoms before hydrogens, since a substituent bond is chemically the
        more sensible axis anyway -- so a single degenerate geometry doesn't
        have to throw the whole conformer away.
        """
        default_pair = get_sterimol_indices(self.xyz_df, self.bonds_df)
        center_id = default_pair[0]
        neighbors = (
            self.bonds_df[self.bonds_df[0] == center_id][1].tolist()
            + self.bonds_df[self.bonds_df[1] == center_id][0].tolist()
        )
        atoms = self.xyz_df['atom'].values
        seen = {default_pair[1]}
        alt_heavy, alt_h = [], []
        for n in neighbors:
            n = int(n)
            if n in seen:
                continue
            seen.add(n)
            (alt_h if atoms[n - 1] == 'H' else alt_heavy).append(n)
        return [default_pair] + [[center_id, n] for n in alt_heavy + alt_h]

    def get_sterimol(self, base_atoms=None, radii='CPK', sub_structure=True, drop_atoms=None, visualize=False):
        # NOTE: deliberately a single attempt at the default auto-pick here,
        # not a per-conformer retry loop -- that turned out to multiply the
        # cost of every failing conformer by the number of candidates (fine
        # for a 6-conformer ensemble like L25, but on a 138-conformer one
        # like L_10 with ~27 failures it pushed a ~34s ligand well past a
        # minute). `ConformerEnsemble.get_sterimol` retries with an alternate
        # atom pair at the ensemble level instead, and only when the DEFAULT
        # pick fails for literally every conformer -- see there.
        if base_atoms is None:
            base_atoms = get_sterimol_indices(self.xyz_df, self.bonds_df)
        if isinstance(base_atoms[0], list):
            sterimol_list = [
                self.process_sterimol_atom_group(atoms, radii, sub_structure, drop_atoms, visualize)
                for atoms in base_atoms
            ]
            return pd.concat(sterimol_list, axis=0)
        return self.process_sterimol_atom_group(base_atoms, radii, sub_structure, drop_atoms, visualize)

    # ── Angles, dihedrals, bond lengths, buried volume ──────────────────────
    # Mirrors of the equivalent `Molecule` methods, ported over so a single
    # conformer frame (built in-memory, no file I/O) supports the same
    # descriptor set `ConformerEnsemble` Boltzmann-averages across conformers.

    def get_bond_angle(self, atom_indices: List[int]) -> pd.DataFrame:
        """
        Bond angle (3 indices, 1-based) or dihedral (4 indices, 1-based)
        between the given atoms -- `get_angle_df` dispatches on length.
        """
        return get_angle_df(self.coordinates_array, atom_indices)

    def get_bond_length_single(self, atom_pair):
        bond_length = calc_single_bond_length(self.coordinates_array, atom_pair)
        bond_length_df = pd.DataFrame([bond_length], index=[f'bond_length_{atom_pair[0]}-{atom_pair[1]}'])
        return bond_length_df

    def get_bond_length(self, atom_pairs):
        """Bond length(s) between the given 1-based atom pair(s)."""
        if isinstance(atom_pairs[0], list):
            bond_length_list = [self.get_bond_length_single(pair) for pair in atom_pairs]
            bond_df = pd.concat(bond_length_list, axis=0)
        else:
            bond_df = self.get_bond_length_single(atom_pairs)
        return bond_df

    def get_buried_volume(
        self,
        metal_index: int,
        radii: str = 'bondi',
        radius: float = 3.5,
        radii_scale: float = 1.17,
        include_hs: bool = False,
        z_axis_atoms=None,
        xz_plane_atoms=None,
    ) -> pd.DataFrame:
        """Buried volume around `metal_index` (1-based) via morfeus BuriedVolume."""
        elements = self.xyz_df['atom'].tolist()
        bv = BuriedVolume(
            elements=elements,
            coordinates=self.coordinates_array,
            metal_index=metal_index,
            radii_type=radii,
            radius=radius,
            radii_scale=radii_scale,
            include_hs=include_hs,
            z_axis_atoms=z_axis_atoms,
            xz_plane_atoms=xz_plane_atoms,
        )
        return pd.DataFrame([{
            'percent_buried_volume':   bv.percent_buried_volume,
            'fraction_buried_volume':  bv.fraction_buried_volume,
            'buried_volume':           bv.buried_volume,
            'free_volume':             bv.free_volume,
        }])


class ConformerEnsemble:
    """
    One ligand's CREST conformer ensemble: parses the multi-structure XYZ file,
    computes Boltzmann weights from the Hartree energies on each block's
    comment line, and computes Sterimol per-conformer plus the weighted average.

    Auto-discovers which file to parse inside `ligand_dir` (in priority order:
    de-duplicated conformers -> raw CREST working-dir conformers -> rotamers ->
    best-only), or pass `ensemble_file` directly to skip discovery.
    """

    CONFORMER_FILENAME_PRIORITY = (
        '{name}_conformers.xyz',
        os.path.join('raw_crest_run', 'crest_conformers.xyz'),
        '{name}_rotamers.xyz',
        os.path.join('raw_crest_run', 'crest_rotamers.xyz'),
        '{name}_best.xyz',
        os.path.join('raw_crest_run', 'crest_best.xyz'),
    )

    def __init__(self, ligand_dir=None, ensemble_file=None, molecule_name=None,
                 temperature=298.15, threshold=None, energy_convention='first'):
        self.temperature = temperature
        self.threshold = threshold
        self.energy_convention = energy_convention
        self.warnings = []

        if ensemble_file is None:
            if ligand_dir is None:
                raise ValueError("Provide either `ligand_dir` or `ensemble_file`.")
            ligand_dir = os.path.abspath(ligand_dir)
            name = molecule_name or os.path.basename(ligand_dir.rstrip(os.sep).rstrip('/'))
            ensemble_file = self._discover_ensemble_file(ligand_dir, name)
            self.molecule_name = name
        else:
            ensemble_file = os.path.abspath(ensemble_file)
            self.molecule_name = molecule_name or os.path.splitext(os.path.basename(ensemble_file))[0]

        self.ensemble_file = ensemble_file
        self.used_fallback_single_structure = os.path.basename(ensemble_file) in ('L_best.xyz',) or \
            os.path.basename(ensemble_file).endswith(('_best.xyz',)) or \
            os.path.basename(ensemble_file) == 'crest_best.xyz'

        blocks = parse_xyz_ensemble(ensemble_file, energy_convention=energy_convention)

        natoms_set = {len(b['xyz_df']) for b in blocks}
        if len(natoms_set) > 1:
            raise ValueError(
                f"[{self.molecule_name}] Inconsistent atom counts across conformers in "
                f"{ensemble_file}: {sorted(natoms_set)}. Refusing to average a corrupted ensemble."
            )

        self.conformers = [
            _ConformerFrame(b['xyz_df'], b['energy'], f"{self.molecule_name}_conf{i + 1}", threshold=threshold)
            for i, b in enumerate(blocks)
        ]
        self.energies_hartree = np.array([c.energy for c in self.conformers], dtype=float)

        if np.isnan(self.energies_hartree).any():
            self.warnings.append(
                f"[{self.molecule_name}] {int(np.isnan(self.energies_hartree).sum())} of "
                f"{len(self.conformers)} conformer(s) had an unparsable energy comment line and "
                f"will be excluded from the Boltzmann average."
            )
        self.weights = boltzmann_weights(self.energies_hartree, temperature=temperature)

        if self.used_fallback_single_structure:
            self.warnings.append(
                f"[{self.molecule_name}] No conformer/rotamer ensemble found -- fell back to a single "
                f"lowest-energy structure ({os.path.basename(ensemble_file)}). Reported values are "
                f"NOT Boltzmann-averaged (n=1)."
            )

    @classmethod
    def _discover_ensemble_file(cls, ligand_dir, name):
        for pattern in cls.CONFORMER_FILENAME_PRIORITY:
            candidate = os.path.join(ligand_dir, pattern.format(name=name))
            if os.path.isfile(candidate):
                return candidate
        loose_xyz = sorted(glob.glob(os.path.join(ligand_dir, '*.xyz')))
        if loose_xyz:
            return loose_xyz[0]
        raise FileNotFoundError(
            f"No CREST conformer/rotamer/best .xyz file found for ligand '{name}' in {ligand_dir}."
        )

    @classmethod
    def from_multiple(cls, sources, molecule_name, temperature=298.15, threshold=None,
                      energy_convention='first'):
        """
        Pool conformers from several ligand directories (or explicit ensemble
        files) into ONE ensemble and Boltzmann-weight across the union.

        For merging duplicate/rerun CREST jobs for what's chemically the same
        ligand (e.g. `L12` and a redone `L12_T`, or `L20`/`L20b`) into a single
        result -- this re-derives one Boltzmann weighting over every conformer
        from every source, which is the statistically correct way to combine
        them (as opposed to averaging two already-Boltzmann-averaged numbers,
        which double-counts/under-counts depending on how many conformers
        backed each one).

        Parameters
        ----------
        sources : list[str]
            Ligand directories (auto-discovers the conformer file in each,
            same priority order as the normal constructor) or explicit
            ensemble .xyz file paths -- may mix both.
        molecule_name : str
            Name for the merged result (e.g. 'L12').

        Raises
        ------
        ValueError
            If the sources don't all have the same atom count (they aren't
            actually the same ligand/molecule).
        """
        self = cls.__new__(cls)
        self.temperature = temperature
        self.threshold = threshold
        self.energy_convention = energy_convention
        self.warnings = []
        self.molecule_name = molecule_name
        self.used_fallback_single_structure = False
        self.source_files = []

        all_blocks = []
        for src in sources:
            if os.path.isdir(src):
                name_guess = os.path.basename(os.path.normpath(src))
                f = cls._discover_ensemble_file(src, name_guess)
            else:
                f = src
            self.source_files.append(f)
            all_blocks.extend(parse_xyz_ensemble(f, energy_convention=energy_convention))

        natoms_set = {len(b['xyz_df']) for b in all_blocks}
        if len(natoms_set) > 1:
            raise ValueError(
                f"[{molecule_name}] Sources don't agree on atom count -- refusing to merge "
                f"{self.source_files}: {sorted(natoms_set)}."
            )

        self.ensemble_file = self.source_files
        self.conformers = [
            _ConformerFrame(b['xyz_df'], b['energy'], f"{molecule_name}_conf{i + 1}", threshold=threshold)
            for i, b in enumerate(all_blocks)
        ]
        self.energies_hartree = np.array([c.energy for c in self.conformers], dtype=float)
        if np.isnan(self.energies_hartree).any():
            self.warnings.append(
                f"[{molecule_name}] {int(np.isnan(self.energies_hartree).sum())} of "
                f"{len(self.conformers)} merged conformer(s) had an unparsable energy comment "
                f"line and will be excluded from the Boltzmann average."
            )
        self.weights = boltzmann_weights(self.energies_hartree, temperature=temperature)
        return self

    def get_sterimol(self, base_atoms=None, radii='CPK', sub_structure=False, drop_atoms=None,
                      energy_cutoff_kcal=None, min_weight=None):
        """
        Compute Sterimol for every conformer and Boltzmann-average the result.

        `sub_structure` defaults to False here (unlike `Molecule`/`Molecules_xyz`):
        these are whole, standalone CREST-optimized ligands, not a substituent
        hanging off a shared scaffold, so there's usually nothing to trim down
        to and no reason to pay for the underlying path-enumeration search
        (which can be slow-to-hanging on larger/cyclic ligands). Pass
        `sub_structure=True` explicitly if you do want the substituent-only trim.

        Parameters
        ----------
        base_atoms : list[int] | list[list[int]] | None
            1-indexed [origin, direction] atom pair (or list of such pairs),
            same convention as `Molecule.get_sterimol`. If None, the anchor
            atoms are auto-detected per-conformer from the center of mass
            (`get_sterimol_indices`) -- convenient for a first pass, but for
            a publication table you almost certainly want to pass the actual
            attachment-vector atoms explicitly.
        radii : str, default 'CPK'
        sub_structure, drop_atoms : passthrough to the underlying Sterimol calc.
        energy_cutoff_kcal : float | None
            If given, conformers above this relative energy (kcal/mol vs. the
            ensemble minimum) are excluded from the average (weight forced to 0,
            remaining weights renormalized).
        min_weight : float | None
            If given, conformers with raw Boltzmann weight below this threshold
            are excluded the same way.

        Returns
        -------
        (per_conformer_df, boltzmann_row) : (pd.DataFrame, pd.DataFrame)
            `per_conformer_df` has one row per conformer (energy, rel. energy,
            raw weight, weight actually used, and every Sterimol column) --
            keep this around to sanity-check/audit what went into the average.
            `boltzmann_row` is a single-row DataFrame (index = ligand name)
            with the weighted-average Sterimol columns plus `n_conformers`
            and `n_conformers_used`.
        """
        per_conformer_df, boltzmann_row = self._get_sterimol_once(
            base_atoms, radii, sub_structure, drop_atoms, energy_cutoff_kcal, min_weight
        )

        # If auto-detection picked an atom pair that's degenerate for EVERY
        # conformer (seen on ligand L25 -- a stray C-H bond exactly collinear
        # with the coplane vector, singular `np.linalg.lstsq` every time),
        # retry with the next candidate neighbor of that same center atom
        # instead of reporting a flat NaN row. Retrying at the ensemble level
        # (one atom pair applied to every conformer) rather than per-conformer
        # keeps this cheap: it only fires when the WHOLE ensemble came back
        # empty, not on the far more common case of a handful of individual
        # conformers failing (e.g. L_10, 27/138 conformers) where the default
        # pick already gives a perfectly good weighted average.
        if base_atoms is None and boltzmann_row['n_conformers_used'].iloc[0] == 0 and self.conformers:
            candidates = self.conformers[0]._auto_base_atom_candidates()[1:]
            for candidate in candidates:
                retry_per_conf, retry_boltz = self._get_sterimol_once(
                    candidate, radii, sub_structure, drop_atoms, energy_cutoff_kcal, min_weight
                )
                if retry_boltz['n_conformers_used'].iloc[0] > 0:
                    self.warnings.append(
                        f"[{self.molecule_name}] default auto-detected atom pair failed for every "
                        f"conformer; retried successfully with base_atoms={candidate}."
                    )
                    return retry_per_conf, retry_boltz
            self.warnings.append(
                f"[{self.molecule_name}] tried {1 + len(candidates)} auto-detected atom pairs, "
                f"all failed for every conformer -- pass `base_atoms` explicitly for this ligand."
            )

        return per_conformer_df, boltzmann_row

    def _get_sterimol_once(self, base_atoms, radii, sub_structure, drop_atoms,
                            energy_cutoff_kcal, min_weight):
        """One pass of the per-conformer Sterimol + Boltzmann-average computation
        for a single fixed `base_atoms` choice -- see `get_sterimol` for the
        public entry point (which adds the ensemble-level retry-on-total-failure
        behavior on top of this)."""
        rel_kcal = (self.energies_hartree - np.nanmin(self.energies_hartree)) * HARTREE_TO_KCAL
        rows = []

        # Fix the expected Sterimol column schema *before* running anything,
        # purely from `base_atoms` (known regardless of success/failure). This
        # matters for the edge case where EVERY conformer fails: without a
        # pre-fixed schema, `per_conformer_df` would end up with zero value
        # columns, the NaN-based failure check below would have nothing to
        # check against, and a fully-failed ligand would silently get
        # `n_conformers_used == n_conformers` with empty descriptor columns
        # instead of being correctly flagged as a total failure.
        _base_cols = ['B1', 'B5', 'L', 'loc_B5', 'B1_B5_angle']
        if base_atoms is not None and isinstance(base_atoms[0], list):
            expected_cols = [f"{g[0]}-{g[1]}_{c}" for g in base_atoms for c in _base_cols]
        else:
            expected_cols = list(_base_cols)

        for conf, w, rel in zip(self.conformers, self.weights, rel_kcal):
            row = {'energy_hartree': conf.energy, 'rel_energy_kcal': rel, 'weight': w}
            row.update({c: np.nan for c in expected_cols})
            try:
                sdf = conf.get_sterimol(base_atoms, radii=radii, sub_structure=sub_structure, drop_atoms=drop_atoms)
                single_group = len(sdf) == 1
                for idx, srow in sdf.iterrows():
                    for col, val in srow.items():
                        key = col if single_group else f"{idx}_{col}"
                        row[key] = val
            except Exception as e:
                self.warnings.append(f"[{self.molecule_name}] conformer {conf.molecule_name} failed Sterimol: {e}")
            rows.append(row)

        per_conformer_df = pd.DataFrame(rows, index=[c.molecule_name for c in self.conformers])

        # Use the pre-fixed schema (not `per_conformer_df.columns`) so a
        # totally-failed ensemble is still correctly detected below.
        value_cols = expected_cols
        eff_weight = per_conformer_df['weight'].copy()
        eff_weight[per_conformer_df[value_cols].isna().any(axis=1)] = 0.0
        if energy_cutoff_kcal is not None:
            eff_weight[per_conformer_df['rel_energy_kcal'] > energy_cutoff_kcal] = 0.0
        if min_weight is not None:
            eff_weight[per_conformer_df['weight'].fillna(0) < min_weight] = 0.0

        total = eff_weight.sum()
        if not total or np.isnan(total) or total <= 0:
            self.warnings.append(
                f"[{self.molecule_name}] All conformers failed Sterimol or were excluded by the "
                f"energy/weight cutoff -- no Boltzmann average computed."
            )
            eff_weight[:] = 0.0
            boltzmann_row = pd.DataFrame([{c: np.nan for c in value_cols}], index=[self.molecule_name])
        else:
            eff_weight = eff_weight / total
            # NaN * 0.0 is still NaN in IEEE754, so a failed conformer's NaN
            # descriptor values would silently poison the whole ligand's
            # average even though its (correctly zeroed) weight shouldn't
            # contribute anything. Zero those cells out first -- their
            # weight is already 0, so this doesn't change the math for any
            # conformer that actually succeeded.
            avg = per_conformer_df[value_cols].fillna(0.0).mul(eff_weight, axis=0).sum(skipna=False)
            boltzmann_row = pd.DataFrame([avg], index=[self.molecule_name])

        per_conformer_df.insert(3, 'weight_used', eff_weight)
        boltzmann_row['n_conformers'] = len(self.conformers)
        boltzmann_row['n_conformers_used'] = int((eff_weight > 0).sum())
        return per_conformer_df, boltzmann_row

    # ── Generic Boltzmann-averaging engine ──────────────────────────────────
    # Same weighting/failure-handling logic as `_get_sterimol_once`, factored
    # out so angle/dihedral/bond-length/buried-volume can reuse it instead of
    # re-deriving the NaN-poisoning-fix, cutoff handling, and total-failure
    # bookkeeping each time. `_get_sterimol_once` is left untouched (it
    # predates this and is already verified against the real dataset).

    def _boltzmann_average_descriptor(self, descriptor_name, expected_cols, per_conformer_calc,
                                       energy_cutoff_kcal=None, min_weight=None):
        """
        Parameters
        ----------
        descriptor_name : str
            Only used in warning messages (e.g. "angle[2, 1, 5]").
        expected_cols : list[str]
            Fixed column schema the average is computed over -- fixed up
            front (not derived from results) so a total failure across every
            conformer is still detected correctly, same reasoning as in
            `_get_sterimol_once`.
        per_conformer_calc : callable(conf) -> dict
            Called once per `_ConformerFrame`; must return a dict mapping
            (a subset of) `expected_cols` to scalar values. Raising inside
            this callable marks that conformer as failed (NaN, weight forced
            to 0 for the average) without aborting the rest of the ensemble.

        Returns
        -------
        (per_conformer_df, boltzmann_row) -- same shape/convention as
        `get_sterimol`'s return value.
        """
        rel_kcal = (self.energies_hartree - np.nanmin(self.energies_hartree)) * HARTREE_TO_KCAL
        rows = []

        for conf, w, rel in zip(self.conformers, self.weights, rel_kcal):
            row = {'energy_hartree': conf.energy, 'rel_energy_kcal': rel, 'weight': w}
            row.update({c: np.nan for c in expected_cols})
            try:
                values = per_conformer_calc(conf)
                for c in expected_cols:
                    if c in values:
                        row[c] = values[c]
            except Exception as e:
                self.warnings.append(
                    f"[{self.molecule_name}] conformer {conf.molecule_name} failed {descriptor_name}: {e}"
                )
            rows.append(row)

        per_conformer_df = pd.DataFrame(rows, index=[c.molecule_name for c in self.conformers])

        value_cols = expected_cols
        eff_weight = per_conformer_df['weight'].copy()
        eff_weight[per_conformer_df[value_cols].isna().any(axis=1)] = 0.0
        if energy_cutoff_kcal is not None:
            eff_weight[per_conformer_df['rel_energy_kcal'] > energy_cutoff_kcal] = 0.0
        if min_weight is not None:
            eff_weight[per_conformer_df['weight'].fillna(0) < min_weight] = 0.0

        total = eff_weight.sum()
        if not total or np.isnan(total) or total <= 0:
            self.warnings.append(
                f"[{self.molecule_name}] All conformers failed {descriptor_name} or were excluded "
                f"by the energy/weight cutoff -- no Boltzmann average computed."
            )
            eff_weight[:] = 0.0
            boltzmann_row = pd.DataFrame([{c: np.nan for c in value_cols}], index=[self.molecule_name])
        else:
            eff_weight = eff_weight / total
            # See `_get_sterimol_once` -- NaN * 0.0 is still NaN in IEEE754,
            # so failed-but-zero-weighted conformers must be zeroed first.
            avg = per_conformer_df[value_cols].fillna(0.0).mul(eff_weight, axis=0).sum(skipna=False)
            boltzmann_row = pd.DataFrame([avg], index=[self.molecule_name])

        per_conformer_df.insert(3, 'weight_used', eff_weight)
        boltzmann_row['n_conformers'] = len(self.conformers)
        boltzmann_row['n_conformers_used'] = int((eff_weight > 0).sum())
        return per_conformer_df, boltzmann_row

    # ── Angles ───────────────────────────────────────────────────────────────

    def get_angle(self, atoms, energy_cutoff_kcal=None, min_weight=None):
        """
        Boltzmann-average a bond angle across every conformer.

        Parameters
        ----------
        atoms : list[int]
            Exactly 3 atom indices, 1-based, e.g. [2, 1, 5].

        Returns
        -------
        (per_conformer_df, boltzmann_row) -- `boltzmann_row` has a single
        value column 'angle' (degrees) plus n_conformers/n_conformers_used.
        """
        if len(atoms) != 3:
            raise ValueError(f"get_angle expects exactly 3 atom indices, got {atoms}")

        def _calc(conf):
            df = conf.get_bond_angle(atoms)
            return {'angle': float(df.iloc[0, 0])}

        return self._boltzmann_average_descriptor(
            f"angle{atoms}", ['angle'], _calc,
            energy_cutoff_kcal=energy_cutoff_kcal, min_weight=min_weight,
        )

    # ── Dihedrals ────────────────────────────────────────────────────────────

    def get_dihedral(self, atoms, energy_cutoff_kcal=None, min_weight=None):
        """
        Boltzmann-average a dihedral angle across every conformer.

        Parameters
        ----------
        atoms : list[int]
            Exactly 4 atom indices, 1-based, e.g. [5, 1, 8, 6].

        Returns
        -------
        (per_conformer_df, boltzmann_row) -- `boltzmann_row` has a single
        value column 'dihedral' (degrees) plus n_conformers/n_conformers_used.
        """
        if len(atoms) != 4:
            raise ValueError(f"get_dihedral expects exactly 4 atom indices, got {atoms}")

        def _calc(conf):
            # `get_angle_df`/`get_bond_angle` dispatch on index-list length,
            # so a 4-atom group is automatically treated as a dihedral.
            df = conf.get_bond_angle(atoms)
            return {'dihedral': float(df.iloc[0, 0])}

        return self._boltzmann_average_descriptor(
            f"dihedral{atoms}", ['dihedral'], _calc,
            energy_cutoff_kcal=energy_cutoff_kcal, min_weight=min_weight,
        )

    # ── Bond lengths ─────────────────────────────────────────────────────────

    def get_bond_length(self, atoms, energy_cutoff_kcal=None, min_weight=None):
        """
        Boltzmann-average a bond length across every conformer.

        Parameters
        ----------
        atoms : list[int]
            Exactly 2 atom indices, 1-based, e.g. [1, 5].

        Returns
        -------
        (per_conformer_df, boltzmann_row) -- `boltzmann_row` has a single
        value column 'bond_length' (Angstrom) plus n_conformers/n_conformers_used.
        """
        if len(atoms) != 2:
            raise ValueError(f"get_bond_length expects exactly 2 atom indices, got {atoms}")

        def _calc(conf):
            return {'bond_length': float(conf.get_bond_length_single(atoms).iloc[0, 0])}

        return self._boltzmann_average_descriptor(
            f"bond_length{atoms}", ['bond_length'], _calc,
            energy_cutoff_kcal=energy_cutoff_kcal, min_weight=min_weight,
        )

    # ── Buried volume ────────────────────────────────────────────────────────

    def get_buried_volume(self, metal_index, radii='bondi', radius=3.5, radii_scale=1.17,
                           include_hs=False, z_axis_atoms=None, xz_plane_atoms=None,
                           energy_cutoff_kcal=None, min_weight=None):
        """
        Boltzmann-average %V_bur (and related quantities) around `metal_index`
        (1-based) across every conformer, via morfeus `BuriedVolume`.

        Returns
        -------
        (per_conformer_df, boltzmann_row) -- `boltzmann_row` has columns
        percent_buried_volume, fraction_buried_volume, buried_volume,
        free_volume, plus n_conformers/n_conformers_used.
        """
        value_cols = ['percent_buried_volume', 'fraction_buried_volume', 'buried_volume', 'free_volume']

        def _calc(conf):
            df = conf.get_buried_volume(
                metal_index=metal_index, radii=radii, radius=radius, radii_scale=radii_scale,
                include_hs=include_hs, z_axis_atoms=z_axis_atoms, xz_plane_atoms=xz_plane_atoms,
            )
            return df.iloc[0].to_dict()

        return self._boltzmann_average_descriptor(
            f"buried_volume(metal={metal_index})", value_cols, _calc,
            energy_cutoff_kcal=energy_cutoff_kcal, min_weight=min_weight,
        )


class GoatEnsemble(ConformerEnsemble):
    """CREST-style Sterimol averaging on an ORCA GOAT multi-XYZ file.

    Discovers ``{name}.ens.xyz`` / ``{name}.finalensemble.xyz`` and reads the
    **last** float on each comment line as the energy. For metal-referenced
    CS3 descriptors (``fromM_*``, ``q_anc_sum``, …) use
    :class:`~MolFeatures.M2_data_extractor.metal_complex.MetalComplexEnsemble`
    instead — this class is the generic Sterimol path, not the CS3 metal frame.
    """

    CONFORMER_FILENAME_PRIORITY = (
        '{name}.ens.xyz',
        '{name}.finalensemble.xyz',
        os.path.join('results', '{name}.finalensemble.xyz'),
        '{name}_conformers.xyz',
        os.path.join('raw_crest_run', 'crest_conformers.xyz'),
        '{name}_best.xyz',
    )

    def __init__(self, ligand_dir=None, ensemble_file=None, molecule_name=None,
                 temperature=298.15, threshold=None, energy_convention='last'):
        super().__init__(
            ligand_dir=ligand_dir,
            ensemble_file=ensemble_file,
            molecule_name=molecule_name,
            temperature=temperature,
            threshold=threshold,
            energy_convention=energy_convention,
        )


class CrestLigandSet:
    """
    Batch driver over a directory of CREST ligand subfolders, e.g.::

        doyle_conformers/
          L1/  L1_conformers.xyz  L1_best.xyz  ...
          L2/  L2_conformers.xyz  L2_best.xyz  ...
          ...

    Every subfolder is treated as one ligand. Subfolders that don't contain a
    usable CREST ensemble (empty directories, leftover Gaussian .log-only
    folders, corrupted ensembles, etc.) are skipped rather than crashing the
    whole batch -- check `self.failed_ligands` / `self.warnings` afterward.
    """

    def __init__(self, root_dir, temperature=298.15, threshold=None):
        self.root_dir = os.path.abspath(root_dir)
        self.temperature = temperature
        self.threshold = threshold
        self.warnings = []
        self.failed_ligands = {}
        self.ensembles = {}

        def _looks_like_ligand_dir(d):
            # Skip subfolders with no .xyz anywhere (e.g. a `sterimol_results/`
            # output folder written by a previous run into this same root_dir) --
            # avoids noisy false "failures" on re-run rather than actually
            # signaling a broken ligand.
            full = os.path.join(self.root_dir, d)
            if not os.path.isdir(full):
                return False
            return bool(glob.glob(os.path.join(full, '**', '*.xyz'), recursive=True))

        subdirs = sorted(d for d in os.listdir(self.root_dir) if _looks_like_ligand_dir(d))
        for name in subdirs:
            ligand_dir = os.path.join(self.root_dir, name)
            try:
                ensemble = ConformerEnsemble(
                    ligand_dir=ligand_dir, molecule_name=name,
                    temperature=temperature, threshold=threshold,
                )
                self.ensembles[name] = ensemble
            except Exception as e:
                self.failed_ligands[name] = str(e)
                self.warnings.append(f"[{name}] skipped: {e}")

    @staticmethod
    def _natural_key(s):
        return [int(t) if t.isdigit() else t for t in re.split(r'(\d+)', str(s))]

    def get_sterimol_df(self, base_atoms=None, radii='CPK', sub_structure=False, drop_atoms=None,
                         energy_cutoff_kcal=None, min_weight=None, save_dir=None):
        """
        Compute Boltzmann-averaged Sterimol for every ligand in the set.

        Parameters mirror `ConformerEnsemble.get_sterimol`. `base_atoms` may
        also be a dict keyed by ligand name ({'L1': [1, 2], 'L2': [1, 3], ...})
        for datasets where the anchor-atom numbering isn't consistent across
        ligands; a plain list applies the same atom pair(s) to every ligand.

        Returns
        -------
        (summary_df, detail_df) : (pd.DataFrame, pd.DataFrame)
            `summary_df`: one row per ligand -- the Boltzmann-averaged Sterimol
            table you'd put in an SI. `detail_df`: one row per conformer per
            ligand, for auditing the averages before they go in the paper.
        If `save_dir` is given, both are written there as CSVs, along with
        `sterimol_run_warnings.txt` listing every skipped ligand/conformer.
        """
        summary_rows, detail_frames, run_warnings = [], [], list(self.warnings)

        for name in sorted(self.ensembles, key=self._natural_key):
            ensemble = self.ensembles[name]
            atoms_for_ligand = base_atoms.get(name) if isinstance(base_atoms, dict) else base_atoms
            per_conf, boltz = ensemble.get_sterimol(
                base_atoms=atoms_for_ligand, radii=radii, sub_structure=sub_structure,
                drop_atoms=drop_atoms, energy_cutoff_kcal=energy_cutoff_kcal, min_weight=min_weight,
            )
            run_warnings.extend(ensemble.warnings)
            per_conf.insert(0, 'ligand', name)
            detail_frames.append(per_conf)
            summary_rows.append(boltz)

        summary_df = pd.concat(summary_rows, axis=0) if summary_rows else pd.DataFrame()
        detail_df = pd.concat(detail_frames, axis=0) if detail_frames else pd.DataFrame()
        self.last_run_warnings = run_warnings

        if save_dir:
            os.makedirs(save_dir, exist_ok=True)
            summary_df.to_csv(os.path.join(save_dir, 'sterimol_boltzmann_avg.csv'))
            detail_df.to_csv(os.path.join(save_dir, 'sterimol_per_conformer.csv'))
            with open(os.path.join(save_dir, 'sterimol_run_warnings.txt'), 'w') as f:
                f.write('\n'.join(run_warnings) if run_warnings else 'No warnings.')

        return summary_df, detail_df


if __name__=='__main__':
    pass
