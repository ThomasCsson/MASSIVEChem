from rdkit import Chem
from rdkit.Chem import Draw

mol_smi = 'CCC(CC)(CC)CC'
mol = Chem.MolFromSmiles(mol_smi)
img = Draw.MolToImage(mol)
img.show()