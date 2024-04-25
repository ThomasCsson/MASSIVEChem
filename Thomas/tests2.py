from rdkit import Chem
from rdkit.Chem import Draw
import numpy as np

mol_smi = input('SMILEs: ')
mol = Chem.MolFromSmiles(mol_smi)

atoms_present = []
for atom in mol.GetAtoms():
    AtomSymbol = atom.GetSymbol()
    atoms_present.append(AtomSymbol)
print(atoms_present)
