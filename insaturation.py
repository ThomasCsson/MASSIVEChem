from rdkit import Chem
from rdkit.Chem import Draw
import time

#This calculates the insaturation of the molecule.
#Input: SMILES
#Output: Insaturation degree & Image of the molecule
mol_smi = input('Enter the SMILES of a molecule: ')
start = time.time()
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
end = time.time()
print(f'The insaturation of the given molecule is {insaturation_level(mol)}')
print(f'This program took{end - start} seconds to complete')
img = Draw.MolToImage(mol_1)
img.show()

def subgroup_finder(mol):
      print('')
