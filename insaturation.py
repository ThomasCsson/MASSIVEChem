from rdkit import Chem
from rdkit.Chem import Draw


mol_smi = input('Enter the SMILES of a molecule: ')
mol_1 = Chem.MolFromSmiles(mol_smi)
mol = Chem.AddHs(mol_1)
def insaturation_level (mol):
    C,N,HX, = 0,0,0
    halogens_hydrogen = ['F','Cl','Br','I','At','H']
    for atom in mol.GetAtoms():
            atom_sym = atom.GetSymbol()
            if atom_sym == 'C':
                  C = C+1
            elif atom_sym == 'N':
                  N = N + 1
            elif atom_sym in halogens_hydrogen:
                  HX = HX + 1
    insaturation = C + 1 + (N-HX)/2
    return insaturation
print(f'The insaturation of the given molecule is {insaturation_level(mol)}')
img = Draw.MolToImage(mol)
img.show()
