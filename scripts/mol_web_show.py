from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from bokeh.plotting import figure, show
import os

def mol_web_show(mol_smi,show_Hs= False, show_3D= False):
        # name of the file name to create
    filename = 'molecule_image.png'

    #finds the current directory
    current_directory = os.getcwd()

    #creates a file path to the current directory
    filepath = os.path.join(current_directory, filename)

    #checks if the file already exists
    if not os.path.exists(filepath):

        #if no, creates the path
        with open(filepath, 'a'):
            pass
    else:

        #else pass
        pass

    # Generate the image from the molecule
    mol = Chem.MolFromSmiles(mol_smi)
    # Adds the hydrogens to the molecule if specified
    if show_Hs:
        mol = Chem.AddHs(mol)

    # Show the molecule in 3D if specified
    if show_3D:
        mol = AllChem.EmbedMolecule(mol)

    image = Draw.MolToImage(mol)

    # Save the image to a file
    image.save(filepath)
    print(filepath)

    # Creating a Bokeh figure to display the molecule
    p = figure(width=350, height=350,toolbar_location=None, x_range=(0, 1), y_range=(0, 1))
    p.image_url(url=[filepath], x=0, y=1, w=1, h=1)

    # Hide grid lines and axes
    p.xgrid.grid_line_color = None
    p.ygrid.grid_line_color = None
    p.xaxis.visible = False
    p.yaxis.visible = False
    
    return p


mol_smi = input('MOL:  ')
show(mol_web_show(mol_smi))