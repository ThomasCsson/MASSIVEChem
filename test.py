from rdkit import Chem
from rdkit.Chem import Draw

mol_smi = input('SMILES: ')
mol = Chem.MolFromSmiles(mol_smi)
img = Draw.MolToImage(mol)
img.show()

