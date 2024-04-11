from rdkit import Chem
from rdkit.Chem import Draw

mol_smi = input('SMILES: ')
mol = Chem.MolFromSmiles(mol_smi)
molH = Chem.AddHs(mol)
img = Draw.MolToImage(molH)
img.show()
