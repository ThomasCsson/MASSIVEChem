import pandas as pd
from rdkit import Chem

Symmetry_specis_PG = {
    'Point Group': ['C∞v','C1','C2','C2h','C2v','C3v','C4v','Ci','Cs','D∞h','D2','D2d','D2h','D3','D3d','D3h','D4d','D4h','D5d','D5h','D6h','D7h','Ih','Kh','Oh','Td'],
    'Symmetry Plane': [1000,0,0,1,2,3,4,0,1,1000,0,2,3,0,3,4,4,5,5,6,8,8,15,1000,9,6],
    'Proper Rotation Axis':[1000,0,1,1,1,2,3,0,0,1000,3,3,3,5,5,5,5,7,9,9,13,13,59,1000,23,11],
    'Improper Rotation Axis':[0,0,0,0,0,0,0,0,0,2,0,2,0,0,2,2,4,2,4,2,4,6,44,1000,14,6],
    'Inversion Center':[0,0,0,1,0,0,0,1,0,1,0,0,1,0,1,0,0,1,1,0,1,0,1,1,1,0]
}

symmetry_species_DF= pd.DataFrame(Symmetry_specis_PG)
print(symmetry_species_DF)

Chiral_PG = ['C1','C2','D2','D3']
Dipole_PG = ['C1']


with open('Molecular_point_group', 'r', encoding='ISO-8859-1') as file:
    # Sauter la première ligne
    next(file)
    lines = file.readlines()

colonne1 = []
colonne2 = []


# Parcours des lignes et extraction des valeurs
for line in lines:
    values = line.split()  # Séparation des valeurs par espace

    # Vérifier si la ligne contient le nombre attendu de valeurs
    if len(values) >= 2:
        # Remplacement des virgules par des points pour les décimales et conversion en float
        colonne1.append(values[0])
        colonne2.append(values[1])



Point_group_database = {
    'Point Group': colonne1,
    'Molecular Formula': colonne2,

}

Point_group_database_DF = pd.DataFrame(Point_group_database)
print(Point_group_database_DF)

from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np

# Define a molecule from a SMILES string
mol = Chem.MolFromSmiles('CCO')
mol = Chem.AddHs(mol)
# Generate 3D coordinates using the MMFF94 force field method
AllChem.EmbedMolecule(mol, AllChem.ETKDG())

# Get the atomic coordinates
coords = mol.GetConformer().GetPositions()

# Calculate the center of mass
center_of_mass = np.mean(coords, axis=0)

# Translate the molecule to center it at the origin
coords -= center_of_mass

# Output the 3D coordinates
for i, coord in enumerate(coords):
    print(f"{mol.GetAtomWithIdx(i).GetSymbol()} {coord[0]:.4f} {coord[1]:.4f} {coord[2]:.4f}")


from rdkit import Chem
from rdkit.Chem import AllChem

# Define a molecule from a SMILES string
mol = Chem.MolFromSmiles('CCO')
mol = Chem.AddHs(mol)
# Generate 3D coordinates using the MMFF94 force field method
AllChem.EmbedMolecule(mol, AllChem.ETKDG())

# Output the 3D coordinates
print(Chem.MolToMolBlock(mol))

