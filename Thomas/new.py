from rdkit import Chem


mol_smi = input('Enter SMILEs: ')
mol_without_Hs = Chem.MolFromSmiles(mol_smi)
mol = Chem.AddHs(mol_without_Hs)

list_atoms = []
for atom in mol.GetAtoms():
    list_atoms.append(atom.GetSymbol())
list_atoms.remove('H')
print(list_atoms)