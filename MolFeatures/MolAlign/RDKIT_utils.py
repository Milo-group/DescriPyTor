from rdkit import Chem
from rdkit.Chem import rdDetermineBonds, Draw, rdFMCS, MolFromSmiles, MolFromSmarts
from rdkit.Chem.rdmolops import AddHs, DeleteSubstructs, RenumberAtoms
from rdkit.Chem.rdmolops import GetMolFrags, FragmentOnBonds
import sys
import os

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
try:
    from .fmcs import fmcs
    from .df_funs import *
    from .rotations import *
except ImportError:
    from fmcs import fmcs
    from df_funs import *
    from rotations import *
from PIL import Image
import io
import torch
import torch.nn.functional as F

import itertools
import copy

from sklearn.decomposition import PCA

charges_list = [0, 1, -1, 2, -2, 3, -3, 4, -4]



def df_to_mol(df):

    temp_df = df.copy()
    original_atoms = temp_df.iloc[:, 0].copy()

    # Replace Pd atoms with C temporarily
    temp_df.iloc[:, 0] = temp_df.iloc[:, 0].replace('Pd', 'C')

    charges_list = [0, 1, -1, 2, -2, 3, -3, 4, -4]  # Example charge list

    mol = None
    for charge in charges_list:
        try:
            xyz = temp_df.to_string(header=False, index=False)
            mol = Chem.MolFromXYZBlock(xyz)
            rdDetermineBonds.DetermineBonds(mol, charge=charge)

            if mol is None or mol.GetNumBonds() == 0:
                raise ValueError("No bonds detected after DetermineBonds.")
            
            break
        except Exception as e:
            print(f"Charge {charge} failed: {e}")
            mol = None
            continue

    if mol is None:
        print("Molecule creation failed for all charges.")
        return None

    return mol

def fragmented_to_mol(mol, main_idx, dummies = True):

    frags = GetMolFrags(mol)
    mw = Chem.RWMol(mol)
    
    upt_idx = main_idx
    
    atoms_to_remove = [a for f in frags if main_idx not in f for a in f]
    atoms_to_remove.sort()
    atoms_to_keep = [a for f in frags if main_idx in f for a in f]

    upt_idx = {j:j - len ([i for i in atoms_to_remove if i < j]) for j in atoms_to_keep}
    
    for a_idx in atoms_to_remove[::-1]:
        for b in mw.GetAtoms()[a_idx].GetBonds():
            mw.RemoveBond(b.GetBeginAtomIdx(), b.GetEndAtomIdx())

    for a_idx in atoms_to_remove[::-1]:
        mw.RemoveAtom(a_idx) 
    
          
    for a in list(mw.GetAtoms())[::-1]:
        if a.GetSymbol() == "*":
            if dummies == False or mw.GetBondBetweenAtoms(upt_idx[main_idx], a.GetIdx()) is None:
                mw.RemoveAtom(a.GetIdx())
            else:
                mw.ReplaceAtom(a.GetIdx(), Chem.rdchem.Atom("*"))
    
    if dummies:
        mw.ReplaceAtom(upt_idx[main_idx], Chem.rdchem.Atom("*"))
    
    
    upt_idx = {upt_idx[i]:i for i in upt_idx}                    
                    
    res_mol = mw.GetMol()
    
    return res_mol, upt_idx

def produce_maps_list (ref, query, main_idx, ref_main_idx, bondCommpare = Chem.rdFMCS.BondCompare.CompareOrder,
                       atomCompare=Chem.rdFMCS.AtomCompare.CompareAny):

    
    maps_list = [{}]
    
    ref_len = len(ref.GetAtoms())
    query_len = len(query.GetAtoms())
    
    ref_org_map = {i:i for i in range(ref_len)}
    query_org_map = {i:i for i in range(query_len)}

    break_flag = False

    while True:

        new_maps_list = []

        mcs = fmcs([ref,query], maximize = "bonds_multiple").smarts

        if mcs == None:
            break
        
        for mcs_smarts in mcs.split("..."):
            mcs_mol = MolFromSmarts(mcs_smarts)  

            ref_mcs_maps = list(ref.GetSubstructMatches(mcs_mol, uniquify = False))
            query_mcs_maps = list(query.GetSubstructMatches(mcs_mol, uniquify = False))
        
            for q_map in query_mcs_maps:
                for r_map in ref_mcs_maps:
                    m = {query_org_map[q_map[i]]:ref_org_map[r_map[i]] for i in range (len(mcs_mol.GetAtoms()))}
                    if main_idx in m:
                        if m[main_idx] != ref_main_idx:
                            continue
                        else:
                            break_flag = True

                    elif ref_main_idx in m.values():
                        continue

                    temp_list = [copy.deepcopy(n) for n in maps_list]
                    for n in temp_list:
                        n.update(m)
                    new_maps_list += temp_list

        if new_maps_list != []:
            maps_list = new_maps_list

            ref = DeleteSubstructs(ref, mcs_mol) 
            ref_org_map = {len([1 for n in range (j) if n not in ref_mcs_maps[0]]):j for j in ref_org_map if j not in ref_mcs_maps[0]}

        if break_flag:
            break

        for mcs_smarts in mcs.split("..."):
            mcs_mol = MolFromSmarts(mcs_smarts)
            query = DeleteSubstructs(query, mcs_mol)

        deleted_values = []
        for q_map in query_mcs_maps:
            deleted_values += list (q_map)

        query_org_map = {len([1 for n in range (j) if n not in deleted_values] ):j for j in query_org_map if j not in query_mcs_maps[0]}

        
    return maps_list

def find_best_fitted_surface(points):
    pts = points.numpy()
    n_components = min(2, pts.shape[0], pts.shape[1])
    if n_components < 1:
        return torch.tensor([0.0, 0.0, 1.0], dtype=torch.float32)
    pca = PCA(n_components=n_components)
    pca.fit(pts)
    normal_vector = torch.tensor(pca.components_[0], dtype=torch.float32)
    return normal_vector

def rotation_angles_between_vectors(vec1, vec2):

    vec1 = vec1 / torch.norm(vec1)
    vec2 = vec2 / torch.norm(vec2)

    cross_product = torch.cross(vec1, vec2)

    yaw = torch.atan2(cross_product[0], cross_product[2])
    pitch = torch.asin(-cross_product[1])
    roll = torch.atan2(vec1[0]*vec2[1] - vec1[1]*vec2[0], vec1[0]*vec2[0] + vec1[1]*vec2[1])

    return torch.stack([pitch, yaw, roll])

def even_fragments(mol):
        
    mols_list = []
    gap_tensor = torch.tensor([])

    mol_atoms = mol.GetAtoms()

    for atom_idx in range(len(mol_atoms)):
    
        bonds = mol_atoms[atom_idx].GetBonds()
        bonds_to_break = [b.GetIdx() for b in bonds if b.GetIdx() > atom_idx] 

        if bonds_to_break == []:
            continue

        for bond in bonds_to_break:
            frag_mol = fragmented_to_mol(FragmentOnBonds (mol, [bond]), atom_idx, dummies = False)[0]
            if len(mol_atoms) - len(frag_mol.GetAtoms()) > 0:
                mols_list.append (bond)
                gap_tensor = torch.cat([gap_tensor, torch.abs(torch.tensor([len(mol_atoms) - 2*len(frag_mol.GetAtoms())]))])

    bond = mols_list[torch.argmin(gap_tensor)]
    return [fragmented_to_mol(FragmentOnBonds (mol, [bond]), mol.GetBonds()[bond].GetBeginAtomIdx(), dummies = False)[0],
            fragmented_to_mol(FragmentOnBonds (mol, [bond]), mol.GetBonds()[bond].GetEndAtomIdx(), dummies = False)[0]] 
        


def find_init_vars(ref, query, df_1, df_2, Frag_mcs = -1, anchors = None):

    query_len = len(query.GetAtoms())
    ref_len = len(ref.GetAtoms())

    if ref_len >= 10:

        query_to_match = Chem.rdmolops.RemoveHs(query)
        ref_to_match = Chem.rdmolops.RemoveHs(ref)
        
        ref_org_map = {}
        for i in range(ref_len):
            if ref.GetAtoms()[i].GetSymbol() != 'H':
                ref_org_map[len(ref_org_map)] = i

        
        query_org_map = {}
        for i in range(query_len):
            if query.GetAtoms()[i].GetSymbol() != 'H':
                query_org_map[len(query_org_map)] = i

    
    else:
        query_to_match = query
        ref_to_match = ref

        ref_org_map = {i:i for i in range(ref_len)}
        query_org_map = {i:i for i in range(query_len)}

    print ("search for MCS...")

    mcs = rdFMCS.FindMCS([ref_to_match, query_to_match], atomCompare=Chem.rdFMCS.AtomCompare.CompareAny, bondCompare=Chem.rdFMCS.BondCompare.CompareOrderExact).smartsString
    mcs_mol = MolFromSmarts(mcs)
    if Frag_mcs >= 0:
        mcs_mol = even_fragments(mcs_mol)[Frag_mcs]
    
    
    if anchors != None:
        mapping = []

        ref_mcs_maps = ref_to_match.GetSubstructMatches(mcs_mol, uniquify = False)
        query_mcs_maps = query_to_match.GetSubstructMatches(mcs_mol, uniquify = False)

        for r_map in ref_mcs_maps:
            for q_map in query_mcs_maps:
                    upt_query_map = [query_org_map[i] for i in q_map]
                    upt_ref_map = [ref_org_map[i] for i in r_map]
                    m = {upt_query_map[i]:upt_ref_map[i] for i in range (len(mcs_mol.GetAtoms()))}
                    if anchors[1] in m:
                        if m[anchors[1]] == anchors[0]:
                            
                            flag = True
                            for n in mapping:
                                if n == m:
                                    flag = False
                            if flag:
                                mapping.append(m)

        ring_info = [[upt_query_map[i] for i in ring] for ring in Chem.GetSSSR(mcs_mol)]
    
    else:
        ref_mcs_map = ref_to_match.GetSubstructMatch(mcs_mol)
        query_mcs_map = query_to_match.GetSubstructMatch(mcs_mol)

        
        upt_query_map = [query_org_map[i] for i in query_mcs_map]
        upt_ref_map = [ref_org_map[i] for i in ref_mcs_map]
                
                
        mapping = [{upt_query_map[i]:upt_ref_map[i] for i in range (len(mcs_mol.GetAtoms()))}]
        ring_info = [[upt_query_map[i] for i in ring] for ring in Chem.GetSSSR(mcs_mol)]

        if len(ref_to_match.GetSubstructMatches(mcs_mol)) == 1 and len (ref_to_match.GetSubstructMatches(mcs_mol, uniquify = False)) == len(ref_mcs_map) * 2 :
            mapping.append({upt_query_map[i]:upt_ref_map[i-1] for i in range (len(mcs_mol.GetAtoms()))})

    df_2_tensor = torch.tensor(df_2.iloc[2:,1:].values.astype('float64'))
    df_1_tensor = torch.tensor(df_1.iloc[2:,1:].values.astype('float64'))

    COM_2 = torch.sum(df_2_tensor [upt_query_map], 0) / len (upt_query_map)
    COM_1 = torch.sum(df_1_tensor [upt_ref_map], 0) / len (upt_ref_map)

    PCA2 = find_best_fitted_surface(df_2_tensor [upt_query_map])
    PCA1 = find_best_fitted_surface(df_1_tensor [upt_ref_map])

    PCA_angles = rotation_angles_between_vectors(PCA1, PCA2)

    return COM_2, COM_1, PCA_angles, mapping, ring_info 



def find_init_vars_ring(file_path1, file_path2, map, atom3):

    query_mol = df_to_mol(file_path1)
    ref_mol = df_to_mol(file_path2)
    
    ring_ref = [[i for i in ring] for ring in Chem.GetSSSR(ref_mol)]
    ring_query = [[i for i in ring] for ring in Chem.GetSSSR(query_mol)]

    
    for itm in ring_ref:
        if list(map.keys())[0] in itm and list(map.keys())[1] in itm and (atom3 is True or atom3[0] in itm):
            ring_ref = itm
            break
    
    for itm in ring_query:
        if list(map.values())[0] in itm and list(map.values())[1] in itm and (atom3 is True or atom3[1] in itm):
            ring_query = itm
            break

    df_2_tensor = torch.tensor(file_path2.iloc[2:,1:].values.astype('float64'))
    df_1_tensor = torch.tensor(file_path1.iloc[2:,1:].values.astype('float64'))

    COM_2 = torch.sum(df_2_tensor [ring_query], 0) / len (ring_query)
    COM_1 = torch.sum(df_1_tensor [ring_ref], 0) / len (ring_ref)


    PCA2 = find_best_fitted_surface(df_2_tensor [ring_query])
    PCA1 = find_best_fitted_surface(df_1_tensor [ring_ref])

    PCA_angles = rotation_angles_between_vectors(PCA1, PCA2)

    df_2_tensor_rotated = rot(df_2_tensor.type(torch.FloatTensor), PCA_angles[0], PCA_angles[1], PCA_angles[2])


    return df_1_tensor - COM_1, df_2_tensor_rotated - COM_2

def check_identities(mols, matched_bonded_atoms):
    
    unique_structures = len (mols)
    identity_list = []
    for i in range(len(mols)):
        if matched_bonded_atoms[i] in [d for n in identity_list for d in n]:
            continue
        identity_list.append({matched_bonded_atoms[i]:mols[i][1]})
        temp_ids = 1
        for j in range(i+1, len(mols)):
            if mols[i][0].HasSubstructMatch(mols[j][0], useChirality=True) and mols[j][0].HasSubstructMatch(mols[i][0], useChirality=True):
                temp_ids += 1
                identity_list[-1][matched_bonded_atoms[j]] = {idx:mols[j][1][atom] for (idx, atom) in enumerate(mols[j][0].GetSubstructMatch(mols[i][0]))}#mols[j][1]

        unique_structures -= temp_ids
        if unique_structures < 2:
            break
        
    return identity_list


def filter_lone(atom_mapping, sample_mol):
    
    atoms_to_keep = []
    for atom1 in atom_mapping:
        if atom1 in atoms_to_keep:
            continue

        for atom2 in atom_mapping:
            if atom1 == atom2: 
                continue

            if sample_mol.GetBondBetweenAtoms(atom1, atom2) != None:
                atoms_to_keep += [atom1, atom2]

    
    return {i:atom_mapping[i] for i in atom_mapping if i in atoms_to_keep}


def check_order(s, new_map, df2, df1):
    s = s.keys()
    ref_tensor = torch.tensor(df1.iloc[2:,1:].values.astype('float64')).to(torch.float64)
    ref_s = [new_map[i] for i in s]
    query_tensor = torch.tensor(df2.iloc[2:,1:].values.astype('float64')).to(torch.float64)
    
    permutations = [list(i) for i in itertools.permutations(s)]
    distances = torch.zeros(len(permutations))
    
    for idx, perm in enumerate(permutations):
        distances[idx] = torch.std(torch.pow(torch.sum(torch.pow (query_tensor[perm] - ref_tensor[ref_s], 2), 1),0.5))

    return (permutations[torch.argmin(distances)])

def choose_map(MCS_maps, frag_ref_mapping, frag_query_mapping, df2, df1):
    
    ref_tensor = torch.tensor(df1.iloc[2:,1:].values.astype('float64')).to(torch.float64)
    query_tensor = torch.tensor(df2.iloc[2:,1:].values.astype('float64')).to(torch.float64)
    
    distances = torch.zeros(len(MCS_maps))

    for idx, mapping in enumerate(MCS_maps):
        for i in mapping:
            if (frag_query_mapping[i] < query_tensor.shape[0] and frag_ref_mapping[mapping[i]] >= ref_tensor.shape[0]) \
            or (frag_query_mapping[i] >= query_tensor.shape[0] and frag_ref_mapping[mapping[i]] < ref_tensor.shape[0]):
                distances[idx] = 10000000
                break
        if distances[idx] == 10000000:
            continue
        query_set = [frag_query_mapping[i] for i in mapping if frag_query_mapping[i] < query_tensor.shape[0]]
        ref_set = [frag_ref_mapping[mapping[i]] for i in mapping if frag_ref_mapping[mapping[i]] < ref_tensor.shape[0]]
        distances[idx] = torch.sum(torch.pow(torch.sum(torch.pow (query_tensor[query_set] - ref_tensor[ref_set], 2), 1),0.5))

    return (MCS_maps[torch.argmin(distances)])

def compare_to_previous(MCS_map, MCSandCoords_map, coords_mapping, frag_query_mapping, frag_ref_mapping , df2, df1):

    common_new_frag = [n for n in MCS_map if frag_ref_mapping[MCS_map[n]] in MCSandCoords_map.values() and frag_ref_mapping[MCS_map[n]] not in coords_mapping.values()]
    common_new = [frag_query_mapping[n] for n in common_new_frag]
    common_previous = [n for n in MCSandCoords_map if MCSandCoords_map[n] in [frag_ref_mapping[MCS_map[i]] for i in MCS_map] and MCSandCoords_map[n] not in coords_mapping.values()]
    common_ref = [MCSandCoords_map[n] for n in MCSandCoords_map if MCSandCoords_map[n] in [frag_ref_mapping[MCS_map[i]] for i in MCS_map] and MCSandCoords_map[n] not in coords_mapping.values()]
    
    ref_tensor = torch.tensor(df1.iloc[2:,1:].values.astype('float64')).to(torch.float64)
    query_tensor = torch.tensor(df2.iloc[2:,1:].values.astype('float64')).to(torch.float64)

    dist_new = torch.sum(torch.pow(torch.sum(torch.pow (query_tensor[common_new] - ref_tensor[common_ref], 2), 1),0.5))
    dist_previous = torch.sum(torch.pow(torch.sum(torch.pow (query_tensor[common_previous] - ref_tensor[common_ref], 2), 1),0.5))
    
    if dist_new > dist_previous:
        for i in common_new_frag:
            del MCS_map[i]
    
    return MCS_map
            

def correct_order(new_order, existing_order, MCSandCoords_map, update_map, coords_mapping):

    for idx, n in enumerate(new_order):
        init_dic = existing_order[list(existing_order.keys())[idx]]
        target_dic = existing_order[n]


        if init_dic == target_dic or len (init_dic) != len (target_dic):
            continue
 
        changes = [init_dic[i] for i in init_dic if init_dic[i]!=target_dic[i]]
        changes += [target_dic[i] for i in target_dic if init_dic[i]!=target_dic[i]]

        if len (changes & coords_mapping.keys()) > 1:
            continue

        for i in init_dic:

            if init_dic[i] in update_map:
                continue
            if init_dic[i] == target_dic[i]:
                continue
            if init_dic[i] in MCSandCoords_map and target_dic[i] in MCSandCoords_map:
                update_map[init_dic[i]] = MCSandCoords_map[target_dic[i]]
            
    return update_map

def check_ring(map, ring1):
    leng = len(ring1) - 1
    indx_1_1 = 0
    indx_1_2 = 0
    if map[0] in ring1:
        indx_1_1 = ring1.index(map[0])
    if map[1] in ring1:
        indx_1_2 = ring1.index(map[1])
    return (abs(indx_1_1 - indx_1_2) == 1 or (indx_1_1 == 0 and indx_1_2 == leng) or (indx_1_1 == leng and indx_1_2 == 0))
    


def apply_MCS(coords_mapping, mol1, mol2, df_1, df_2):
    

    MCSandCoords_map = copy.deepcopy(coords_mapping)
    
    mol1_atoms = mol1.GetAtoms()
    mol2_atoms = mol2.GetAtoms()

    for atom_idx in coords_mapping:
        bonds = mol2_atoms[atom_idx].GetBonds()
        bonded_atoms = [bond.GetEndAtomIdx() if bond.GetEndAtomIdx()!=atom_idx else bond.GetBeginAtomIdx() for bond in bonds]
        matched_bonded_atoms = [a for a in bonded_atoms if a in MCSandCoords_map]

        if len(matched_bonded_atoms) == len(bonded_atoms):
            continue
        
        bonds_to_break = [b.GetIdx() for b in bonds if (b.GetEndAtomIdx() in matched_bonded_atoms or b.GetBeginAtomIdx() in matched_bonded_atoms)] 
        if bonds_to_break == []:
            continue
        
        ref_bonds_to_break = [mol1.GetBondBetweenAtoms(MCSandCoords_map[atom_idx], MCSandCoords_map[i]) for i in matched_bonded_atoms]
        ref_bonds_to_break = [b.GetIdx() for b in ref_bonds_to_break if b != None]

        if bonds_to_break == []:
            continue

        ref_f_mol, frag_ref_mapping = fragmented_to_mol(FragmentOnBonds (mol1, ref_bonds_to_break), coords_mapping[atom_idx], dummies = False)
        query_f_mol, frag_query_mapping = fragmented_to_mol(FragmentOnBonds (mol2, bonds_to_break), atom_idx, dummies = False)
        
        ref_frag_mapping = {frag_ref_mapping[i]:i for i in frag_ref_mapping}
        query_frag_mapping = {frag_query_mapping[i]:i for i in frag_query_mapping}
        
        MCS_maps = produce_maps_list(ref_f_mol, query_f_mol, query_frag_mapping[atom_idx], ref_frag_mapping[coords_mapping[atom_idx]])

        if MCS_maps == [{}]:
            continue

        if len (MCS_maps) == 1:
            MCS_map = MCS_maps[0]
        else:
            MCS_map = choose_map(MCS_maps, frag_ref_mapping, frag_query_mapping, df_2, df_1)

        MCS_map = compare_to_previous(MCS_map, MCSandCoords_map, coords_mapping, frag_query_mapping, frag_ref_mapping, df_2, df_1)

        for n in MCS_map:            
            if frag_query_mapping[n] not in coords_mapping and frag_ref_mapping[MCS_map[n]] not in coords_mapping.values() and frag_query_mapping[n] < len(mol2_atoms) and frag_ref_mapping[MCS_map[n]] < len(mol1_atoms):           
                MCSandCoords_map = {i:MCSandCoords_map[i] for i in MCSandCoords_map if (MCSandCoords_map[i]!= frag_ref_mapping[MCS_map[n]] or i in coords_mapping)}
                MCSandCoords_map[frag_query_mapping[n]] = frag_ref_mapping[MCS_map[n]]

    return MCSandCoords_map

def ring_mapping(coords_mapping, mol1, mol2, atom3):
    mol1_atoms = mol1.GetAtoms()
    mol2_atoms = mol2.GetAtoms()
    ring1 = [[i for i in ring] for ring in Chem.GetSSSR(mol1)]
    ring2 = [[i for i in ring] for ring in Chem.GetSSSR(mol2)]
    keys = list(coords_mapping.keys())
    vals = list(coords_mapping.values())
    for itm in ring1:
        if vals[0] in itm and vals[1] in itm and (atom3 is True or atom3[1] in itm):
            ring1 = itm
            break

    
    for itm in ring2:        
        if keys[0] in itm and keys[1] in itm and (atom3 is True or atom3[0] in itm):
            ring2 = itm
            break

    ring1_len = len(ring1)
    ring2_len = len(ring2)
    if ring1_len == ring2_len:

        while ring2[0] != list(coords_mapping.keys())[0]:
            ring2 = ring2[-1:] + ring2[:-1]
        
        while ring1[0] != list(coords_mapping.values())[0]:
            ring1 = ring1[-1:] + ring1[:-1]
        
        map = copy.deepcopy(coords_mapping)
        for i in range(len(ring1)//2):
            ring1i = ring1[i]
            ring2i = ring2[i]
            if ring2[i] not in list(map.keys()) and ring1[i] not in list(map.values()):
                map[ring2[i]] = ring1[i]
        for i in range(len(ring2) - 1, len(ring1)//2 - 1, -1):
            ring1i = ring1[i]
            ring2i = ring2[i]
            if ring2[i] not in list(map.keys()) and ring1[i] not in list(map.values()):
                map[ring2[i]] = ring1[i]
    else:
        if ring1_len < ring2_len:
            ring1, ring2 = ring2, ring1
            ring1_len, ring2_len = ring2_len, ring1_len
            mol1_atoms ,mol2_atoms = mol2_atoms, mol1_atoms

        while ring1[0] not in list(coords_mapping.keys()):
            ring1 = ring1[-1:] + ring1[:-1]
        
        if ring1[0] in ring2:
            while ring2[0] != ring1[0]:
               ring2 = ring2[-1:] + ring2[:-1]
        else:
            while ring2[0] not in list(coords_mapping.values()):
                ring2 = ring2[-1:] + ring2[:-1]

        map = copy.deepcopy(coords_mapping)
        ring1 = ring1[2:]
        ring2 = ring2[2:]
        for i in range(len(ring2)//2):
            if ring1[i] not in list(map.keys()) and ring2[i] not in list(map.values()):
                map[ring1[i]] = ring2[i]
        j = len(ring1) - 1
        for i in range(len(ring2) - 1, len(ring2)//2 - 1, -1):
            if ring1[i] not in list(map.keys()) and ring2[i] not in list(map.values()):
                map[ring1[j]] = ring2[i]
            j = j - 1
    temp_map = copy.deepcopy(map)
    print(ring1)
    print(ring2)
    bonded_atoms1 = []
    bonded_atoms2 = []
    for atom_idx in list(temp_map.keys()):
        bonds1 = mol1_atoms[atom_idx].GetBonds()
        bonded_atoms1.append([bond.GetEndAtomIdx() if bond.GetEndAtomIdx()!=atom_idx else bond.GetBeginAtomIdx() for bond in bonds1])
    for atom_idx in list(temp_map.values()):
        bonds2 = mol2_atoms[atom_idx].GetBonds()
        bonded_atoms2.append([bond.GetEndAtomIdx() if bond.GetEndAtomIdx()!=atom_idx else bond.GetBeginAtomIdx() for bond in bonds2])
    for k in range(len(bonded_atoms1)):
        for j in bonded_atoms2[k]:
            for i in bonded_atoms1[k]:
                if i not in list(map.keys()) and j not in list(map.values()) and i not in ring1:
                    map[i] = j
                    break
    return map


def spatial_corrections(coords_mapping, MCSandCoords_map, mol2, df_1, df_2):
    
    update_map = {}
    mol2_atoms = mol2.GetAtoms()

    for atom_idx in MCSandCoords_map: 
            
        bonds = mol2_atoms[atom_idx].GetBonds()
        if len(bonds) < 2:
            continue
        
        bonded_atoms = [bond.GetEndAtomIdx() if bond.GetEndAtomIdx()!=atom_idx else bond.GetBeginAtomIdx() for bond in bonds]
        matched_bonded_atoms = [a for a in bonded_atoms if a in MCSandCoords_map]

        if len(matched_bonded_atoms) < 2:
            continue


        mols_to_compare = []
        for r in matched_bonded_atoms:
            bonds_to_break = [b.GetIdx() for b in bonds if (b.GetEndAtomIdx()!= r and b.GetBeginAtomIdx()!=r)]
            mols_to_compare.append(fragmented_to_mol(FragmentOnBonds (mol2, bonds_to_break), atom_idx))

        identical_subtituents = check_identities(mols_to_compare, matched_bonded_atoms)
        ord_identical_subtituents = [check_order(s, MCSandCoords_map, df_2, df_1) if len(s)>1 else [] for s in identical_subtituents]
        for idx in range(len(ord_identical_subtituents)):
            
            if ord_identical_subtituents[idx] != list(identical_subtituents[idx].keys()) and ord_identical_subtituents[idx] != []:
                update_map = correct_order(ord_identical_subtituents[idx], identical_subtituents[idx], MCSandCoords_map, update_map, coords_mapping)

    return update_map