from rdkit import Chem
from rdkit.Chem import Draw

def calculate_unsaturation(mol_smile) -> int:
    #---------------------------------------------------------------------------------------------#
    '''
    calculate_unsaturation(mol_smile)

    Input: molecule under SMILEs representation

    Output: unsaturation of the input molecule (integer value)
    '''
    #---------------------------------------------------------------------------------------------# 

    C, N, HX, halogens_hydrogen = 0, 0, 0, ['F', 'Cl', 'Br', 'I', 'At', 'H']

    mol = Chem.AddHs(Chem.MolFromSmiles(mol_smile))

    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':
            C += 1
        elif atom.GetSymbol() == 'N':
            N += 1
        elif atom.GetSymbol() in halogens_hydrogen:
            HX += 1

    unsaturation = C + 1 + (N - HX) / 2

    return unsaturation
