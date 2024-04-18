from rdkit import Chem

l = []
mol_smi = 'CC=O'
match_smi = 'C=O'
mol = Chem.MolFromSmiles(mol_smi)
match = Chem.MolFromSmiles(match_smi)

idx = (mol.GetSubstructMatch(match))
lowidx = min(idx)
print(lowidx)

