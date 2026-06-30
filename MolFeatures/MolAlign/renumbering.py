
import os
os.chdir(r'C:\Users\edens\Desktop\DescriPytor\DescriPyTor-main\MolFeatures\MolAlign')
import pandas as pd
from tqdm import tqdm
import argparse



from gaussian_handler import *



from df_funs import xyz_to_df
from RDKIT_utils import *
from MCS_optimization import *

from itertools import combinations
import math


from datetime import datetime


print(os.getcwd())
def add_unumbered(MCSandCoords_map, mol2):

    mol2_atoms = mol2.GetAtoms()

    unnumbered = {}
    for idx2 in range (len(mol2_atoms)):
        if idx2 not in MCSandCoords_map:
            for idx in range (len(mol2_atoms)):
                if idx not in MCSandCoords_map.values() and idx not in unnumbered.values():
                    unnumbered[idx2] = idx
                    break

    return unnumbered

def optimize_based_on_current_mapping(df_1, df_2, temp_map, mol2, init_vars = None):

    optimize_atoms = list(temp_map.values())
    temp_map.update(add_unumbered(temp_map, mol2))

    final_temp_map = {temp_map[i]:i for i in range(len(temp_map))}
    final_temp_order = [final_temp_map[i] for i in range(len(final_temp_map))] 
    final_reversed_map = {i:temp_map[i] for i in range(len(temp_map))}
    final_reversed_order = [final_reversed_map[i] for i in range(len(final_reversed_map))] 

    df_2.iloc[2:, :] = df_2.iloc[2:, :].iloc[final_temp_order, :]
    df_2 = optimize_numbered (df_1, df_2, optimize_atoms, init_vars = None)  
    df_2.iloc[2:, :] = df_2.iloc[2:, :].iloc[final_reversed_order, :]

    return df_2, [torch.zeros(3)]*3 + list(init_vars)[3:]


def spatial_matching(MCSandCoords_map, mol2, mol1, df2, df1):
    
    unnumbered = {}
    mol1_atoms = len(mol1.GetAtoms())
    mol2_atoms = len(mol2.GetAtoms())

    ref_tensor = torch.tensor(df1.iloc[2:,1:].values.astype('float64')).to(torch.float64)
    query_tensor = torch.tensor(df2.iloc[2:,1:].values.astype('float64')).to(torch.float64)
    
    for threshold in torch.arange(0.3, 5.0, 0.005):

        target_unmactched = [i for i in range (mol1_atoms) if i not in MCSandCoords_map.values() and i not in unnumbered.values()]
        query_unmactched = [i for i in range (mol2_atoms) if i not in MCSandCoords_map and i not in unnumbered]
        
        if target_unmactched == []:
            break

        for i in query_unmactched:
            distances = torch.zeros(len(target_unmactched))
            for idx, j in enumerate(target_unmactched):
                distances[idx] = torch.pow(torch.sum(torch.pow (query_tensor[i] - ref_tensor[j], 2)),0.5)
            if torch.min(distances) < threshold and target_unmactched[torch.argmin(distances)] not in unnumbered.values():
                unnumbered[i] = target_unmactched[torch.argmin(distances)]
        
   
    return (unnumbered)

# def gen_output_folder(files_list):
#
#     run_name = f"molecules/results/{files_list[0].split('/')[-2]}-{datetime.now().strftime('%D-%H-%M').replace('/', '-')}"
#     os.mkdir(run_name)
#
#     return run_name

def gen_output_folder(files_list):
    if not files_list:
        return f"molecules/results/run-{datetime.now().strftime('%D-%H-%M').replace('/', '-')}"

    # Get the full path of the first file and normalize it
    file_path = files_list[0].replace('\\', '/')

    # Extract the containing folder name
    # If path is like "molecules/feather_files/file.ext", we want "feather_files"
    path_parts = file_path.split('/')

    # Find the folder name - it's the parent directory of the file
    if len(path_parts) >= 2:
        folder_name = path_parts[-2]  # Second-to-last part is the folder name
    else:
        folder_name = "run"  # Default if we can't determine folder

    # Create timestamped output folder
    timestamp = datetime.now().strftime('%D-%H-%M').replace('/', '-')
    run_name = f"molecules/results/{folder_name}-{timestamp}"

    # Create the directory
    os.makedirs(run_name, exist_ok=True)

    return run_name


def gen_out_path(path, comment="", output_directory="molecules/results"):
    # Normalize path separators (convert backslashes to forward slashes)
    path = path.replace('\\', '/')

    # Get just the filename without path
    filename = path.split('/')[-1]

    # Handle different file extensions
    if filename.endswith('.feather'):
        basename = filename[:-8]  # Remove '.feather'
    elif filename.endswith('.xyz'):
        basename = filename[:-4]  # Remove '.xyz'
    else:
        # For any other extension, remove what's after the last dot
        basename = filename.rsplit('.', 1)[0]

    return f"{output_directory}/{basename}_{comment}.xyz"

def gen_out_path1 (path, comment = "", output_directory = "molecules/results"):

    return output_directory + "/" + path.split("/")[-1][10:-4] + "_" + comment + ".xyz"


def batch_renumbering(folder_path = None, log_folder_path = '', feather_path = '', basic_path = "",log_basic_path = " ",feather_basic_path = "", anchors_file = None):
    files_list = []
    min_size = 0
    min_file = ""
    error_files = []  # Added missing declaration
    invalid_files = []
    # os.chdir(folder_path)
    # if folder_path[-1] != '/':
    #     folder_path += '/'



    if anchors_file != None:
        # anchors_df = pd.read_excel(folder_path + "/" + anchors_file)
        anchors_df = pd.read_excel(anchors_file)
        anchors_dic = {folder_path + anchors_df['Mol'][i] + ".xyz": anchors_df['idx'][i] - 1 for i in
                       range(anchors_df.shape[0])}

    # Initialize mol as None to check later if it was assigned
    mol = None
    basic_size = 0  # Initialize basic_size with a default value

    if basic_path != "":  # if the user defines a basic structure file
        if not os.path.isfile(basic_path):
            raise Exception("basic file does not exist")

        try:
            df = xyz_to_df(basic_path)
        except:
            raise Exception("basic file is not valid")

        if log_basic_path != "":  # if the user defines a basic structure and wants to use log file
            mol = log_xyz_to_mol(basic_path, log_basic_path)
        else:
            mol = df_to_mol(df)
            if mol == False:
                raise Exception("basic file charge is not valid")

        basic_size = len(mol.GetAtoms())

    elif feather_basic_path != "":
        df = xyz_from_feather(feather_basic_path)
        mol = df_to_mol(df)
        basic_size = len(mol.GetAtoms())

    # Only calculate basic_size if mol has been assigned
    if mol is not None:
        basic_size = len(mol.GetAtoms())


    # --log files--
    invalid_files = []
    files_list = []
    min_size = 0
    min_path = None

    if log_folder_path is not None and log_folder_path != "":  # Check if log_folder_path is provided
        xyz_files = [f for f in os.listdir(folder_path) if os.path.isfile(os.path.join(folder_path, f))]
        log_files = [f for f in os.listdir(log_folder_path) if os.path.isfile(os.path.join(log_folder_path, f))]

        if len(xyz_files) != len(log_files):
            raise ValueError(
                f"Error: {len(xyz_files)} xyz files, while {len(log_files)} log files. They must have the same number of files.")

        print("Folders have the same number of files. Proceeding with processing...")

        for file_name in log_files:
            file_path = os.path.join(log_folder_path, file_name)
            # Fixed: need to provide both xyz and log path
            xyz_path = os.path.join(folder_path, os.path.splitext(file_name)[0] + ".xyz")
            mol = log_xyz_to_mol(xyz_path, file_path)

            if basic_path == "":
                if min_size == 0 or len(mol.GetAtoms()) < min_size:
                    min_path_log = file_path
                    min_size = len(mol.GetAtoms())
            else:
                if len(mol.GetAtoms()) < basic_size:
                    raise Exception(file_name + " size is smaller than the basic file")

            files_list.append(file_path)

        if basic_path == "":
            print(f"The basic structure path is {min_path_log}")
            # Get corresponding xyz file
            base_name = os.path.splitext(os.path.basename(min_path_log))[0]
            xyz_file = os.path.join(folder_path, base_name + ".xyz")

            if not os.path.exists(xyz_file):
                raise ValueError(f"Matching XYZ file not found for {min_path_log}")
            basic_path = xyz_file
            # Now set mol and basic_size
            df = xyz_to_df(basic_path)
            mol = log_xyz_to_mol(basic_path, min_path_log)
            basic_size = len(mol.GetAtoms())

    # --feather file--
    elif feather_path is not None and feather_path != "":
        min_size = 0
        for file_name in os.listdir(feather_path):  # Changed to iterate through directory
            file_path = os.path.join(feather_path, file_name)
            df = xyz_from_feather(file_path)  # Fixed to use file_path
            mol = df_to_mol(df)

            if feather_basic_path == "":
                if min_size == 0 or len(mol.GetAtoms()) < min_size:
                    min_path = file_path
                    min_size = len(mol.GetAtoms())
            else:
                if len(mol.GetAtoms()) < basic_size:
                    raise Exception(file_name + " size is smaller than the basic file")

            files_list.append(file_path)

        if feather_basic_path == "":
            feather_basic_path = min_path
            basic_path = feather_basic_path
            # Now set basic_size
            df = xyz_from_feather(basic_path)
            mol = df_to_mol(df)
            basic_size = len(mol.GetAtoms())

    # --xyz files only--
    else:
        print("No basic structure file provided, using the smallest file in the folder as the basic structure.")
        for file_name in os.listdir(folder_path):
            file_path = os.path.join(folder_path, file_name)
            try:
                df = xyz_to_df(file_path)
            except:
                print(f"{file_name} is not valid")
                invalid_files.append(file_name)
                continue

            mol = df_to_mol(df)
            if not mol:
                print(f"{file_name} charge is not valid")
                invalid_files.append(file_name)
                continue

            if basic_path == "":
                if min_size == 0 or len(mol.GetAtoms()) < min_size:
                    min_path = file_path
                    min_size = len(mol.GetAtoms())
            else:
                if len(mol.GetAtoms()) < basic_size:
                    raise Exception(f"{file_name} size is smaller than the basic file")

            files_list.append(file_path)

        if basic_path == "":
            basic_path = min_path
            # Now set basic_size
            df = xyz_to_df(basic_path)
            mol = df_to_mol(df)
            basic_size = len(mol.GetAtoms())

    if not files_list:
        raise Exception("No valid files found.")

    # Create output directory
    output_folder = gen_output_folder(files_list)


    if feather_path is not None and feather_path != "":
        df = xyz_from_feather(basic_path)
        df.to_csv(
            gen_out_path(basic_path, "target", output_folder),
            header=None, index=None, sep=' ', mode='w+'
        )
    else:
        # For regular xyz files
        xyz_to_df(basic_path).to_csv(
            gen_out_path(basic_path, "target", output_folder),
            header=None, index=None, sep=' ', mode='w+'
        )
    print(f"The basic structure path is {basic_path}")


    parts = basic_path.split("_")
    target_idx = 0
    for part in parts:
        try:
            target_idx = int(part)
            break
        except ValueError:
            continue
    target_idx = target_idx - 1
    print(f"Target index for renumbering is {target_idx}")
    renumbering_dict_list=[]
    # organize files_list by the number appearing in the name of each xyz file
    def extract_first_int_from_filename(x):
        for part in x.split("_"):
            try:
                return int(part)
            except ValueError:
                continue
        return float('inf')  # Files without a number will be sorted last

    files_list.sort(key=extract_first_int_from_filename)

    print(f"Files to optimize: {len(files_list)}")
    for idx, file_path in enumerate(tqdm(files_list, desc="Optimizing Molecules")):
        
        if file_path == basic_path:
            print("Skipping basic file")
            continue

        print(f"Optimizing {file_path}... ({idx + 1}/{len(files_list)})")

        # try:
            
        if anchors_file is not None:
            anchors = [anchors_dic[basic_path], anchors_dic[file_path]]
        else:
            anchors = None
        try:
        # Determine input mode and call renumbering accordingly
            if feather_path != "":
                
                numbering_dict = renumbering(
                    feather_file1=feather_basic_path,
                    feather_file2=file_path,
                    output_directory=output_folder,
                    anchors=anchors,
                    batch=True,
                    mode="feather"
                    )
            elif log_folder_path != "":
                numbering_dict = renumbering(
                    os.path.basename(basic_path), os.path.basename(file_path), output_folder,
                    anchors=anchors, batch=True, mode="log"
                )
            else:
                print(f'idx is {idx}, target_idx is {target_idx}')
                if idx == target_idx:
                    numbering_dict = renumbering(
                        file1_path=basic_path, file2_path=file_path, output_directory=output_folder,
                        anchors=anchors, batch=True, mode="xyz", one_run_before_target=True
                    )
                    
                else:
                    numbering_dict = renumbering(
                        basic_path, file_path, output_folder,
                        anchors=anchors, batch=True, mode="xyz"
                    )
                
            if numbering_dict is None:
                print(f"Failed to optimize {file_path} due to invalid mapping.")
                error_files.append(file_path)
                # renumbering_dict_list.append({})
                continue
            print(f"Successfully optimized {file_path}. Mapping: {numbering_dict}")
            renumbering_dict_list.append(numbering_dict)
        except FileNotFoundError as e:
            print(f"File not found: {file_path}. Error: {e}")
            error_files.append(file_path)
            continue

   
    print ("Optimization Ended!")
    print (f"{len(invalid_files)} files were invalid:")
    print (invalid_files)
    print (f"{len(error_files)} failed to optimize:")
    print (error_files)
    print(f"Results saved to: {output_folder}")
    print(f'Current working directory: {os.getcwd()}')
    return renumbering_dict_list , target_idx

def choose_init_map(init_vars, df_1, df_2, mol2):

    if len (init_vars[3]) == 1:
        return copy.deepcopy(init_vars[3][0])
    
    else:
        matches = torch.zeros(len(init_vars[3]))

        for idx, m in enumerate(init_vars[3]):
            temp_map = copy.deepcopy(m)
            df_2_temp = copy.deepcopy(df_2)
            df_2_temp, _ = optimize_based_on_current_mapping(df_1, df_2_temp, temp_map, mol2, init_vars)
            coords_mapping = Coords_Mapping(df_1, df_2_temp, mol2)
            matches[idx] = len (coords_mapping)

        return copy.deepcopy(init_vars[3][torch.argmax(matches)])

def complete_rings(corrds_mapping, temp_map, init_vars, mol1):
    
    update_map = {}
    for ring in init_vars[4]:
        alignment = 0
        for i in ring:
            if i >= len(mol1.GetAtoms()):
                alignment = 0
                break
            if i in corrds_mapping:
                if temp_map[i] == corrds_mapping[i]:
                    alignment += 1
        if alignment >= 3:
            update_map.update({i:temp_map[i] for i in ring})

            reversed_corrds_mapping = {corrds_mapping[i]:i for i in corrds_mapping}
            for i in update_map:
                if update_map[i] in reversed_corrds_mapping:
                    if i != reversed_corrds_mapping[update_map[i]]:
                        return {}

    return update_map

def center_coords(df):
    df = df.iloc[2:].reset_index(drop=True)
    for coord in ["X", "Y", "Z"]:
        df[coord] = pd.to_numeric(df[coord], errors="coerce")
        df[coord] -= df[coord].mean()
    return df


from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdmolops
from rdkit.Chem import AllChem
from PIL import Image

from rdkit import Chem


def get_smiles_with_indices(mol):
    # Generate the standard SMILES representation of the molecule
    mol = Chem.rdmolops.RemoveHs(mol)
    smiles = Chem.MolToSmiles(mol)

    # Create a dictionary to hold atom indices
    atom_labels = {}
    for atom in mol.GetAtoms():
        atom_idx = atom.GetIdx()
        atom_symbol = atom.GetSymbol()
        atom_charge = atom.GetFormalCharge()

        # If the atom has a charge, append the charge to the atom's symbol
        if atom_charge != 0:
            atom_symbol += f"{atom_charge:+d}"  # Adding charge with '+' or '-' sign

        # Store the atom symbol + index
        atom_labels[atom_idx] = atom_symbol + str(atom_idx)

    # Modify the SMILES to include the atom indices
    smiles_with_indices = ''.join([atom_labels[atom.GetIdx()] for atom in mol.GetAtoms()])

    # Print the SMILES string with atom indices
    print("SMILES with indices:", smiles_with_indices)

    return smiles_with_indices


def renumbering(
    file1_path=None,
    file2_path=None,
    log_file1=None,
    log_file2=None,
    feather_file1=None,
    feather_file2=None,
    output_directory="molecules/results",
    anchors=None,
    shuffle=False,
    batch=False,
    ring_pref=True,
    mode="xyz",
    one_run_before_target=False):

    if mode == "feather":
        df_1 = xyz_from_feather(feather_file1)
        df_2 = xyz_from_feather(feather_file2)
        mol1 = df_to_mol(df_1)
        mol2 = df_to_mol(df_2)



    elif mode == "log":

        df_1 = xyz_to_df(file1_path)
        df_2 = xyz_to_df(file2_path)
        mol1 = log_xyz_to_mol(file1_path, log_file1)
        mol2 = log_xyz_to_mol(file2_path, log_file2)
    else:  # default: xyz
        df_1 = xyz_to_df(file1_path)
        df_2 = xyz_to_df(file2_path)
        mol1 = df_to_mol(df_1)
        mol2 = df_to_mol(df_2)

    #df_1 = center_coords(df_1)
    #df_2 = center_coords(df_2)
    #print(df_1)

    if df_1.shape[0] > df_2.shape[0]:
        df_1, df_2 = df_2, df_1
        if mode == "feather":
            feather_file1, feather_file2 = feather_file2, feather_file1
        elif mode == "log":
            file1_path, file2_path = file2_path, file1_path
            log_file1, log_file2 = log_file2, log_file1
        else:  # xyz mode
            file1_path, file2_path = file2_path, file1_path
    if shuffle:
        df_1 = shuffle_df(df_1)


    init_vars = find_init_vars(mol1, mol2, df_1, df_2, anchors = anchors)
    temp_map = choose_init_map(init_vars, df_1, df_2, mol2)
    #try:
    cycle_num = 1
    cycle = 0
    mcs_cycles = 0
    while cycle != cycle_num:
        df_2, init_vars = optimize_based_on_current_mapping(df_1, df_2, temp_map, mol2, init_vars)
        coords_mapping = Coords_Mapping(df_1, df_2, mol2)

        if not ring_pref:
            temp_map = apply_MCS(coords_mapping, mol1, mol2, df_1, df_2)


            while type(temp_map) is not dict:
                temp_map = apply_MCS(coords_mapping, mol1, mol2, df_1, df_2)

            coords_mapping.update(complete_rings(coords_mapping, temp_map, init_vars, mol1))

        coords_mapping.update(complete_rings(coords_mapping, temp_map, init_vars, mol1))


        if coords_mapping == {} and mcs_cycles < 2:
             init_vars = find_init_vars(mol1, mol2, df_1, df_2, Frag_mcs = mcs_cycles, anchors = anchors)
             temp_map = choose_init_map(init_vars, df_1, df_2, mol2)
             mcs_cycles += 1
             continue


        mcs_cycles = 0
        if ring_pref:
            temp_map = apply_MCS(coords_mapping, mol1, mol2, df_1, df_2)

            while type(temp_map) is not dict:
                    temp_map = apply_MCS(coords_mapping, mol1, mol2, df_1, df_2)

        cycle += 1
    print(f"coords mapping is {coords_mapping}")
    MCSandCoords_map = copy.deepcopy(temp_map)
    MCSandCoords_map.update(spatial_corrections(coords_mapping, MCSandCoords_map, mol2, df_1, df_2))
    MCSandCoords_map.update(spatial_matching(MCSandCoords_map, mol2, mol1, df_2, df_1))
    MCSandCoords_map.update(add_unumbered(MCSandCoords_map, mol2))

    final_map = {MCSandCoords_map[i]:i for i in range(len(MCSandCoords_map))}
    final_order = [final_map[i] for i in range(len(final_map))]
    df_2.iloc[2:, :] = df_2.iloc[2:, :].iloc[final_order, :]

    #rearange the rest of the data by the new numberig
    if mode == "feather":
        nbo_charge_df = nbo_charge(feather_file2)
        nbo_charge_df.iloc[:] = nbo_charge_df.iloc[final_order, :]  # rearrange the nbo df

    if batch == False:
        if shuffle:
            if mode == "feather":
                df_1.to_csv(gen_out_path(feather_file1, "shuffled_target", output_directory), header=None, index=None,
                            sep=' ', mode='w+')
            elif mode == "log":
                df_1.to_csv(gen_out_path(file1_path, "shuffled_target", output_directory), header=None, index=None,
                            sep=' ', mode='w+')
            else:  # xyz mode
                df_1.to_csv(gen_out_path(file1_path, "shuffled_target", output_directory), header=None, index=None,
                            sep=' ', mode='w+')
        else:
            if mode == "feather":
                df_1.to_csv(gen_out_path(feather_file1, "target", output_directory), header=None, index=None, sep=' ',
                            mode='w+')
            elif mode == "log":
                df_1.to_csv(gen_out_path(file1_path, "target", output_directory), header=None, index=None, sep=' ',
                            mode='w+')
            else:  # xyz mode
                print(f"Saving {file1_path} to {gen_out_path(file1_path, 'target', output_directory)}")
                print(f"Saving {file2_path} to {gen_out_path(file2_path, 'optimized', output_directory)}")
    df_1.to_csv(gen_out_path(file1_path, "target", output_directory), header=None, index=None, sep=' ',
                mode='w+')
    df_2.to_csv(gen_out_path (file2_path, "optimized", output_directory), header=None, index=None, sep=' ', mode='w+')
    if one_run_before_target:
        df1_identity_map = {idx: idx for idx in df_1.index}
        xyz_name1 = os.path.splitext(os.path.basename(file1_path))[0]
        xyz_name2 = os.path.splitext(os.path.basename(file2_path))[0]
        return {xyz_name2:coords_mapping, xyz_name1:df1_identity_map}
    xyz_name = os.path.splitext(os.path.basename(file1_path))[0]
    return {xyz_name:coords_mapping}

    # except:
    #     print (f"{file2_path} wasn't optimized")
    #     return file2_path


def mapping_by_bond():
             
    print("Enter the bonds (use - to separate the atoms): ")
    bond_1 = input('Bond from molecule #1: ')
    bond_2 = input('Bond from molecule #2: ')
    atoms1 = bond_1.split("-")
    atoms2 = bond_2.split("-")
    mapping = {int(atoms1[0]) -1 : int(atoms2[0]) - 1, int(atoms1[1]) - 1 : int(atoms2[1]) - 1}
    return mapping

def check_mul_ring(map, ring1):
    for ring in ring1:
        if check_ring(map, ring):
            return True
    return False

def check_bond(map, ring1, ring2):
    if check_mul_ring(list(map.keys()), ring1) and check_mul_ring(list(map.values()), ring2):
        flag = 0
        keys = list(map.keys())
        vals = list(map.values())
        for itm in ring1:
            if (keys[0] in itm and keys[1] in itm) or (vals[0] in itm and vals[1] in itm) :
                flag +=1
                for itm in ring2:
                    if (keys[0] in itm and keys[1] in itm) or (vals[0] in itm and vals[1] in itm) :
                        flag +=1
        if flag == 2:
            return True
        elif flag >= 2:
            # molecule #1 will be the one with more atoms
            atom1 = input('Please enter a different atom on the ring from molecule #1: ')
            atom2 = input('Please enter a different atom on the ring from molecule #2: ')
            if atom1.isnumeric() and atom2.isnumeric():
                return [int(atom1) - 1, int(atom2) - 1]
            else:
                return False
        else:
            return False


def renumbering_ring(file1_path, file2_path, output_directory = "molecules/results", shuffle = False, batch = False, temp_map = {0:0, 1:1}, by_bond = False):

    df_1 = xyz_to_df(file1_path)
    df_2 = xyz_to_df(file2_path)
    
    flag = False

    if df_1.shape[0] > df_2.shape[0]:
        flag = True
        df_1, df_2 = df_2, df_1
        file1_path, file2_path = file2_path, file1_path
        keys = list(temp_map.keys())
        vals = list(temp_map.values())
        temp_map = {vals[0]:keys[0], vals[1]:keys[1]}
        

    if shuffle:
        df_1 = shuffle_df(df_1)

    mol1 = df_to_mol(df_1)
    mol2 = df_to_mol(df_2)

    ring1 = [[i for i in ring] for ring in Chem.GetSSSR(mol1)]
    ring2 = [[i for i in ring] for ring in Chem.GetSSSR(mol2)]
    
    # if not by_bond:
    #     temp_map = {}
    # else:
    #     temp_map = mapping_by_bond()

    checkbond = check_bond(temp_map, ring1, ring2)

    if not checkbond:
        print('Error, please choose a different bond')
        return 0

    if flag:
        keys = list(temp_map.keys())
        vals = list(temp_map.values())
        temp_map = {vals[0]:keys[0], vals[1]:keys[1]}
   
    init_vars = find_init_vars_ring(df_1, df_2, temp_map,checkbond)
    # try:
    tmap = copy.deepcopy(temp_map)
    df_2, init_vars = optimize_based_on_current_mapping(df_1, df_2, temp_map, mol2, init_vars)
    temp_map = ring_mapping(tmap, mol1, mol2, checkbond)
    map = copy.deepcopy(temp_map)
    map.update(spatial_corrections(temp_map, map, mol2, df_1, df_2))
    map.update(spatial_matching(map, mol2, mol1, df_2, df_1))
    map.update(add_unumbered(map, mol2))

    final_map = {map[i]:i for i in range(len(map))}
    final_order = [final_map[i] for i in range(len(final_map))]
    for i in range(len(final_order)):
        if i not in final_order:
            final_order.append(i)
    df_t = df_2.iloc[:2, :]
    df_2 = df_2.iloc[2:, :].iloc[final_order, :]
    dfs = [df_t, df_2]
    df_2 = pd.concat(dfs)

    if batch == False:

        if shuffle:
            df_1.to_csv(gen_out_path1 (file1_path, "shuffled_target", output_directory), header=None, index=None, sep=' ', mode='w+')

        else:
            df_1.to_csv(gen_out_path1 (file1_path, "target", output_directory), header=None, index=None, sep=' ', mode='w+')
    
    df_2.to_csv(gen_out_path1 (file2_path, "optimized", output_directory), header=None, index=None, sep=' ', mode='w+')  

    return



def main():
    parser = argparse.ArgumentParser(description="Renumber atoms in molecules.")

    # Arguments
    parser.add_argument("--file1", type=str,
                        help="Path to the first molecule file (basic structure for batch (optional) or regular renumbering).",
                        default="")
    parser.add_argument("--file2", type=str, help="Path to the second molecule file (used for regular renumbering).",
                        default=None)
    parser.add_argument("--output", default="molecules/results", help="Output directory.")
    parser.add_argument("--anchors", nargs='+', type=int, help="List of anchor indices.", default=None)
    parser.add_argument("--shuffle", action="store_false", help="Disable shuffling.", default=False)
    parser.add_argument("--batch", action="store_true", help="Enable batch mode.", default=True)
    parser.add_argument("--ring_pref", action="store_true", help="Enable ring preference.", default=True)
    parser.add_argument("--folder_path", type=str, help="Path to the folder for batch renumbering.", default=None)
    parser.add_argument("--anchors_file", type=str, help="Path to the Excel file containing anchor indices.",default=None)


    args = parser.parse_args()
    
    if args.batch:  # If batch mode is enabled
        # Ensure folder_path is provided for batch renumbering
        if not args.folder_path:
            raise Exception("Folder path must be provided for batch renumbering.")
        # organize args.folder_path folder by the number appearing in the name of each xyz file
        
        # Ensure anchors_file is provided if anchors are not specified
        if args.anchors_file:
            renumbering_dict_list=batch_renumbering(
                folder_path=args.folder_path,
                basic_path=args.file1,  #  file1 as the basic structure file
                anchors_file=args.anchors_file,
    
    
            )
        else:
            # No anchors_file, so set anchors to None
            renumbering_dict_list=batch_renumbering(
                folder_path=args.folder_path,
                basic_path=args.file1,
                anchors_file=None,  # No anchors file
    
            )
    else:  # Otherwise, run regular renumbering
        renumbering(
            file1_path=args.file1,
            file2_path=args.file2,
            output_directory=args.output,
            anchors=args.anchors,
            shuffle=args.shuffle,
            ring_pref=args.ring_pref
        )
        # if renumbering_dict_list write it to json file
    if renumbering_dict_list is not None:
        print(renumbering_dict_list)
        import json
        # files= os.listdir(args.output)
        # xyz_files_names = [file for file in files if file.endswith('.xyz')]
        
        
        # renumbering_dict_list = [{"file_name": file_name, "renumbering_dict": renumbering_dict} for file_name, renumbering_dict in zip(xyz_files_names, renumbering_dict_list)]
        json_path = os.path.join(args.output, 'renumbering_dict_list.json').replace('\\', '/')

        with open(json_path, 'w') as f:
            json.dump(renumbering_dict_list, f, indent=4)


if __name__ == "__main__":
    main()
