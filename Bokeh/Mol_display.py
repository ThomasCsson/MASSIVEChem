
from bokeh.plotting import figure, show
from rdkit.Chem import Draw, AllChem
from rdkit import Chem
from bokeh.plotting import row
import numpy as np
import matplotlib as plt


input_mol = input('SMILES: ')

# Function to generate RDKit molecule image and save to file
def save_molecule_image_to_file(smi, file_path, show_Hs=False, show_3D = False):
    # Generate the image from the molecule
    mol = Chem.MolFromSmiles(smi)

    # Adds the hydrogens to the molecule if specified
    if show_Hs:
        mol = Chem.AddHs(mol)

    # Show the molecule in 3D if specified
    if show_3D:
        mol_3D = AllChem.EmbedMolecule(mol)

    image = Draw.MolToImage(mol)

    # Save the image to a file
    image.save(file_path)

output_file_path = "Bokeh/molecule_image.png"
save_molecule_image_to_file(input_mol, output_file_path, False, False)

def mol_web_show(image_url):

    # Creating a Bokeh figure to display the molecule
    p = figure(width=400, height=400,toolbar_location=None, x_range=(0, 1), y_range=(0, 1))
    p.image_url(url=[image_url], x=0, y=1, w=1, h=1)

    # Hide grid lines and axes
    p.xgrid.grid_line_color = None
    p.ygrid.grid_line_color = None
    p.xaxis.visible = False
    p.yaxis.visible = False

    return p



image_url = "Bokeh/molecule_image.png"

show(mol_web_show(image_url))

def final_layout(p1, p2):
    layout = row(p1,p2)
    show(layout)

