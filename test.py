from rdkit import Chem
from rdkit.Chem import Draw

mol_smi = input('Add SMILES here: ')
mol = Chem.MolFromSmiles(mol_smi)
img = Draw.MolToImage(mol)
img.show()