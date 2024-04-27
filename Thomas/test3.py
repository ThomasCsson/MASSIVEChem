from rdkit import Chem
from rdkit.Chem import Draw

mol_smi = input('SMILEs: ')
mol = Chem.MolFromSmiles(mol_smi)
molHs = Chem.AddHs(mol)

img = Draw.MolToImage(molHs)
img.show()