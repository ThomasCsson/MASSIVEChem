from rdkit import Chem
from rdkit.Chem import Draw

image = Draw.MolToImage(Chem.MolFromSmiles('C=O'))

image.show()