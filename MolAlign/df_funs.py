import pandas as pd
import torch
import random


def xyz_to_df(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    num_atoms = int(lines[0].strip())

    # Extract atomic symbols and coordinates
    symbols = []
    coordinates = []
    for line in lines[2:]:
        parts = line.strip().split()
        if len(parts) <= 1:
            break
        symbols.append(parts[0])
        coordinates.append([float(coord) for coord in parts[1:4]])

    # Create DataFrame
    df = pd.DataFrame({
        'Symbol': symbols,
        'X': [coord[0] for coord in coordinates],
        'Y': [coord[1] for coord in coordinates],
        'Z': [coord[2] for coord in coordinates]
    })
    
    new_row = pd.DataFrame({'Symbol': len(symbols),'X':'','Y': "",'Z': ''}, index=[0])
    new_gap = pd.DataFrame({'Symbol': '','X':'','Y': "",'Z': ''}, index=[0])
    df = pd.concat([new_row, new_gap, df])

    return df

def shuffle_df (df):
    
    coords = torch.tensor(df.iloc[2:,1:].values.astype('float64'))
    atoms = df.iloc[2:,:1].values
    
    val = list(range(coords.shape[0]))
    random.shuffle (val)
    coords = coords[val]
    atoms = atoms[val]
    
    df.iloc[2:,1:] = coords
    df.iloc[2:,:1] = atoms
    
    return df


def df_w_symbols(df, tensor, save = None):

    symbols = df.iloc[2:, [0]]
    numpy_array = tensor.detach().numpy()
    df = pd.DataFrame(numpy_array)
    df =  pd.concat([symbols, df], axis=1)
    df.columns = ['Symbol','X','Y','Z']
    new_row = pd.DataFrame({'Symbol': symbols.values.shape[0],'X':'','Y': "",'Z': ''}, index=[0])
    new_gap = pd.DataFrame({'Symbol': '','X':'','Y': "",'Z': ''}, index=[0])
    df = pd.concat([new_row, new_gap, df])

    if save != None:
        print (save)
        df.to_csv(save + ".xyz", header=None, index=None, sep=' ', mode='w+')

    return (df)

def read_big_xyz_file(file_path):

    lst = []
    with open(file_path, 'r') as file:
        lines = file.readlines()
    num_atoms = int(lines[0].strip())
    length = len(lines)
    i = 2
    end_index = num_atoms
    while i < length:
        symbols = []
        coordinates = []
        for line in lines[i:end_index+2]:
            parts = line.strip().split()
            symbols.append(parts[0])
            coordinates.append([float(coord) for coord in parts[1:4]])
        df = pd.DataFrame({
            'Symbol': symbols,
            'X': [coord[0] for coord in coordinates],
            'Y': [coord[1] for coord in coordinates],
            'Z': [coord[2] for coord in coordinates]
            })
        i = i+num_atoms+2
        end_index = end_index + int(num_atoms) + 2
        new_row = pd.DataFrame({'Symbol': len(symbols),'X':'','Y': "",'Z': ''}, index=[0])
        new_gap = pd.DataFrame({'Symbol': '','X':'','Y': "",'Z': ''}, index=[0])
        df = pd.concat([new_row, new_gap, df])
        lst.append(df)
    return lst
        
   
