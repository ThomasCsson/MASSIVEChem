#the goal of this file is to simulate the shifting of a hydrogen on a benzene cycle

import rdkit
from rdkit import Chem
import pandas as pd

with open('Benzene_shifting.txt', 'r', encoding='ISO-8859-1') as file:
    # Sauter la première ligne
    next(file)
    lines = file.readlines()

colonne11 = []
colonne22 = []
colonne33 = []
colonne44 = []

# Parcours des lignes et extraction des valeurs
for line in lines:
    values = line.split()  # Séparation des valeurs par espace

    # Vérifier si la ligne contient le nombre attendu de valeurs
    if len(values) >= 2:
        # Remplacement des virgules par des points pour les décimales et conversion en float
        colonne11.append(values[0])
        colonne22.append(float(values[1]))
        colonne33.append(float(values[2]))
        colonne44.append(float(values[3]))

Benzene = {
    'Substituent': colonne11,
    'Z1': colonne22,
    'Z2': colonne33,
    'Z3': colonne44
}

Benzene_DF = pd.DataFrame(Benzene)

Benzene_DF = pd.DataFrame(Benzene)



def benezene_shifting(mol_smi,DF_benz):
    mol = Chem.MolFromSmiles(mol_smi)
    Ben_smi = 'c1ccccc1'
    Ben_mol = Chem.MolFromSmiles(Ben_smi)
    if mol.HasSubstructMatch(Ben_mol):
        l = mol.GetSubstructMatch(Ben_mol)
        l_0_neighbor = (mol.GetAtomWithIdx(l[0])).GetNeighbors()
        l_1_neighbor = (mol.GetAtomWithIdx(l[1])).GetNeighbors()



        return l_0_neighbor
        z2 = 1
        z3 = 1
        z4 = 1
        Hydrogen_shift = 7.26 + z2 + z3 + z4
    else:
        return 'No benzene in the Smile'

print(benezene_shifting('c1ccccc1',1))

from rdkit import Chem

# Load your molecule (replace 'your_molecule.smiles' with your file name or SMILES string)
mol = Chem.MolFromSmiles('c1ccccc1')

# Index of the carbon atom you're interested in
carbon_index = 5

# Get the atom with index i
carbon_atom = mol.GetAtomWithIdx(carbon_index)

# Get the neighboring atoms
neighboring_atoms = [neighbor.GetIdx() for neighbor in carbon_atom.GetNeighbors()]

print(f"Neighboring atoms of carbon {carbon_index}: {neighboring_atoms}")