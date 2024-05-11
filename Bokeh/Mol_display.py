
from bokeh.plotting import figure, show
from rdkit.Chem import Draw, AllChem
from rdkit import Chem
from bokeh.plotting import row



input_mol = input('SMILES: ')


# Function to generate RDKit molecule image and save to file
def save_molecule_image_to_file(mol_smi, file_path, show_Hs=False, show_3D = False):

    #---------------------------------------------------------------------------------------------#
    '''
    save_molecule_image_to_file(mol_smi)
    
    Input: molecule under SMILEs representation
    
    Output: image of the molecule in the molecule_image.png file

    Functionnality: - if show_Hs= True:
                        shows all the hydrogens of the molecule
                    - if show_3D= True:
                        shows the molecule in 3D and the chirality
    '''
    #---------------------------------------------------------------------------------------------#


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
    image.save(file_path)

    

output_file_path = "Bokeh/molecule_image.png"
save_molecule_image_to_file(input_mol, output_file_path, False, False)

def mol_web_show(image_url):

    #---------------------------------------------------------------------------------------------#
    '''
    mol_web_show(image_url)
    
    Input: path of the molecule image
    
    Output: image of the molecule in bokeh

    '''
    #---------------------------------------------------------------------------------------------#


    # Creating a Bokeh figure to display the molecule
    p = figure(width=400, height=400,toolbar_location=None, x_range=(0, 1), y_range=(0, 1))
    p.image_url(url=[image_url], x=0, y=1, w=1, h=1)

    # Hide grid lines and axes
    p.xgrid.grid_line_color = None
    p.ygrid.grid_line_color = None
    p.xaxis.visible = False
    p.yaxis.visible = False

    show(p)

    return p

image_url = "molecule_image.png"

mol_web_show(image_url)



