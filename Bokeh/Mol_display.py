from rdkit import Chem
from bokeh.plotting import figure, show
import numpy as np
from rdkit.Chem import Draw



def validate_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    return True

input_mol = input('SMILES: ')


# Validate the input SMILES string
if not validate_smiles(input_mol):
    print("Invalid SMILES input.")
    exit()
mol = Chem.MolFromSmiles(input_mol)

image = Draw.MolToImage(mol)
image.show()

from rdkit import Chem
from PIL import Image
from io import BytesIO

# Function to generate RDKit molecule image and save to file
def save_molecule_image_to_file(mol, file_path):
    # Generate the image from the molecule
    image = Chem.Draw.MolToImage(mol)

    # Save the image to a file
    image.save(file_path)

if mol is None:
    print("Invalid SMILES input.")
    exit()

output_file_path = "molecule_image.png"
save_molecule_image_to_file(mol, output_file_path)

