import pandas as pd
import numpy as np
import torch
import random
try:
    from .df_funs import *
    from torch.autograd import Variable
    from .rotations import *
    from .loss import *
    from .RDKIT_utils import filter_lone
except ImportError:
    from df_funs import *
    from torch.autograd import Variable
    from rotations import *
    from loss import *
    from RDKIT_utils import filter_lone
from sklearn.model_selection import ParameterGrid

def analyze_df (df, atoms_dictionary = {}):

    coords_tensor = torch.tensor(df.iloc[2:,1:].values.astype('float64'))

    atoms_list = []
    if atoms_dictionary == {}:
        target_flag = True
    else:
        target_flag = False
    
    for atom in df.iloc[2:,:1].values:
        atom = atom[0]
        if target_flag:
            if atom not in atoms_dictionary:            
                atoms_dictionary[atom] = len(atoms_dictionary)
            atoms_list.append(atoms_dictionary[atom])

        else:
            if atom not in atoms_dictionary:
                atoms_list.append(len(atoms_dictionary)) 
            else:
                atoms_list.append(atoms_dictionary[atom])
        
    atoms_tensor = torch.nn.functional.one_hot(torch.tensor(atoms_list), num_classes=len(atoms_dictionary)+1)         
    full_tensor = torch.concat([coords_tensor, atoms_tensor], dim = 1)

    return coords_tensor, atoms_tensor, full_tensor, atoms_dictionary

def init_random():

    initial_vars = torch.rand(6).to(torch.float64)
    initial_vars[:3] *= (2 * torch.pi)
    initial_vars[3:] *= 2
    
    return [Variable(initial_vars, requires_grad=True)]

def init_angles(Angles, init_vars):

    initial_vars = torch.zeros(6).to(torch.float64)
    initial_vars [3] = init_vars[2][0] + Angles['X_angle']
    initial_vars [4] = init_vars[2][1] + Angles['Y_angle']
    initial_vars [5] = init_vars[2][2] + Angles['Z_angle']

    return [Variable(initial_vars, requires_grad=True)]



def Coords_Mapping (target_df, sample_df, sample_mol, dist_threshold = 0.5):

    final_results = [[],{}]     
    _, _, target, atoms_dictionary = analyze_df(target_df)
    sample_coords, sample_atoms, sample, _ = analyze_df(sample_df, atoms_dictionary)
    atom_mapping = {}
    for idx, point in enumerate(sample):
        if torch.min(torch.pow(torch.sum(torch.pow (target - point, 2), 1),0.5)) < dist_threshold:
            atom_mapping [idx] = torch.argmin(torch.pow(torch.sum(torch.pow (target - point, 2), 1),0.5)).item()
    
    atom_mapping = filter_lone(atom_mapping, sample_mol)

    return atom_mapping

def optimize_numbered (target_df, sample_df, optimize_atoms, init_vars = None, lr=0.1):
 
    _, _, target, atoms_dictionary = analyze_df(target_df)
    sample_coords, sample_atoms, _, _ = analyze_df(sample_df, atoms_dictionary)
    
    initial_vars = torch.zeros(6).to(torch.float64)

    if init_vars != None:
        sample_coords -= init_vars [0]
        target[:,:3] -= init_vars [1]
    
        initial_vars[3] += init_vars[2][0] 
        initial_vars[4] += init_vars[2][1] 
        initial_vars[5] += init_vars[2][2] 

    var_tensor = [Variable(initial_vars, requires_grad=True)]

    optim = torch.optim.Adam(var_tensor, lr=lr)

    output = gen_rot(sample_coords, var_tensor[0])
    output = torch.concat([output, sample_atoms], dim = 1)

    it = 0
    loss = torch.sum(loss_function(output[optimize_atoms], target[optimize_atoms], 30))
    memory_list = [loss]

    while (memory_list[0] - loss) > 1e-6 or it < 500:

        it += 1
        optim.zero_grad()
        loss.backward()
        optim.step()
        
        output = gen_rot(sample_coords, var_tensor[0])
        output = torch.concat([output, sample_atoms], dim = 1)
        loss = torch.sum(loss_function(output[optimize_atoms], target[optimize_atoms], 30))


        memory_list.append(loss)
        if len(memory_list) == 501:
            memory_list = memory_list[1:]

    final_coords = gen_rot(sample_coords, var_tensor[0])

    if init_vars != None:
        final_coords += init_vars[1]

    res_df = df_w_symbols(sample_df, final_coords)
    
    return res_df


