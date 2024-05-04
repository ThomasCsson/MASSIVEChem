from rdkit import Chem

mol = Chem.MolFromSmarts('*[Sh0]*')
mol1 = Chem.MolFromSmiles('CSC')
print(mol1.HasSubstructMatch(mol))