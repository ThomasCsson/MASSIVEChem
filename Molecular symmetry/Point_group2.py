import pandas as pd


PG_dict = {
    'Point Group':['sigma','inversion','C2','C3','C4','C5','C6','C7','C8','S3','S4','S5','S6','S7','S8','S9','S10','S11','S12','S13','S14','S15','S16','Cinf','Sinf'],
    'C1': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    'Cs': [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    'Ci': [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    'C2': [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    'C3': [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    'C4': [0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    'C5': [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    'C6': [0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    'C7': [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    'C8': [0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    'D2': [0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    'D3': [0, 0, 3, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    'D4': [0, 0, 5, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    'D5': [0, 0, 5, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    'D6': [0, 0, 7, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    'D7': [0, 0, 7, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    'D8': [0, 0, 8, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    'C2v': [2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    'C3v': [3, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    'C4v': [4, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    'C5v': [5, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    'C6v': [6, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    'C2h': [1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    'C3h': [1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    'C4h': [1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    'C5h': [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    'C6h': [1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    'D2h': [3, 1, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    'D3h': [4, 0, 3, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    'D4h': [5, 1, 5, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    'D5h': [6, 0, 5, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    'D6h': [7, 1, 7, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    'D2d': [2, 0, 3, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    'D3d': [3, 1, 3, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    'D4d': [4, 0, 5, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    'D5d': [5, 1, 5, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
    'D6d': [6, 0, 7, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
    'S4': [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    'S6': [0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    'S8': [0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    'Td': [6, 0, 3, 8, 0, 0, 0, 0, 0, 0, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    'Oh': [9, 1, 9, 8, 6, 0, 0, 0, 0, 0, 6, 0, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    'Ih': [15, 1, 15, 20, 0, 12, 0, 0, 0, 0, 0, 0, 20, 0, 0, 0, 12, 0, 0, 0, 0, 0, 0, 0, 0],
    'Cinfv': [1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
    'Dinfh': [1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1]
}

pg_df = pd.DataFrame(PG_dict)
print(pg_df)

from rdkit import Chem
import numpy as np
from rdkit.Chem import AllChem
from rdkit.Chem import Draw

# Define a function to generate 3D coordinates for a molecule using RDKit
def generate_3d_coordinates(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    return mol

# Define a function to generate 3D coordinates for a molecule using RDKit
def generate_3d_coordinates(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    return mol

# Define a function to show the atomic coordinates of a molecule
def show_coordinates(mol):
    conformer = mol.GetConformer()
    num_atoms = mol.GetNumAtoms()

    # Create an empty list to store the coordinates
    coordinates_list = []

    for i in range(num_atoms):
        pos = conformer.GetAtomPosition(i)
        atom = mol.GetAtomWithIdx(i)
        symbol = atom.GetSymbol()
        x, y, z = pos.x, pos.y, pos.z

        # Append the coordinates to the list
        coordinates_list.append({'Atom': symbol, 'X': x, 'Y': y, 'Z': z})

    # Convert the list of dictionaries to a DataFrame
    df = pd.DataFrame(coordinates_list)
    print(df)
    return df


# Define a function to determine the principal axis of a molecule
"""def determine_principal_axis(mol):
    # Get the atomic coordinates
    coords = mol.GetConformer().GetPositions()

    # Calculate the center of mass
    center_of_mass = np.mean(coords, axis=0)

    # Translate the molecule to center it at the origin
    coords -= center_of_mass

    # Calculate the covariance matrix
    covariance_matrix = np.cov(coords.T)

    # Use singular value decomposition (SVD) to find principal axes
    _, _, v = np.linalg.svd(covariance_matrix)

    # The first principal axis corresponds to the direction of highest variance
    principal_axis = v[0]

    return principal_axis"""


# Define a function to visualize a molecule from its SMILES representation
def visualize_molecule(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    if mol is not None:
        return Draw.MolToImage(mol)
    else:
        return None

# Define the SMILES string for the molecule
smiles = 'c1ccccc1P(c1ccccc1)c1ccccc1'

# Generate 3D coordinates for the molecule
mol = generate_3d_coordinates(smiles)

""""# Determine the principal axis of the molecule
principal_axis = determine_principal_axis(mol)
print("Principal axis:", principal_axis)"""

# Visualize the molecule
image = visualize_molecule(smiles)
if image is not None:
    image.show()
else:
    print("Unable to visualize molecule from the given SMILES string.")

# Generate 3D coordinates for the molecule
mol = generate_3d_coordinates(smiles)

# Show the atomic coordinates of the molecule
show_coordinates(mol)

import pandas as pd
import numpy as np

df = show_coordinates(mol)

# Define atomic masses
atomic_masses = {'H': 1.0079, 'C': 12.0107, 'P': 30.9738}

# Calculate the mass of each atom and add a 'Mass' column to the DataFrame
df['Mass'] = df['Atom'].map(atomic_masses)

# Calculate the total mass of the molecule
total_mass = df['Mass'].sum()

# Calculate the center of mass along each axis
center_of_mass_x = np.dot(df['Mass'], df['X']) / total_mass
center_of_mass_y = np.dot(df['Mass'], df['Y']) / total_mass
center_of_mass_z = np.dot(df['Mass'], df['Z']) / total_mass

# Output the center of mass
print("Center of Mass (X):", center_of_mass_x)
print("Center of Mass (Y):", center_of_mass_y)
print("Center of Mass (Z):", center_of_mass_z)
