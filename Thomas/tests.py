from rdkit import Chem
from rdkit.Chem import Draw
import numpy as np

mol_smi = input('SMILEs: ')
mol = Chem.MolFromSmiles(mol_smi)
def list_atoms_present(mol):
    atoms_present = []
    for atom in mol.GetAtoms():
        AtomSymbol = atom.GetSymbol()
        atoms_present.append(AtomSymbol)
    return atoms_present
def mass_main_peak(mol):
    mass_main = 0
    for atom in mol.GetAtoms():
        mass_main = mass_main + atom.GetAtomicNum()
    return mass_main
def list_number_atoms_present(mol):
    pass

print(mass_main_peak(mol))



