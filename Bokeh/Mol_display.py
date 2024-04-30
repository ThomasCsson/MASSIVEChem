
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

output_file_path = "molecule_image.png"
save_molecule_image_to_file(input_mol, output_file_path, False, False)

def mol_web_show(image_url):

    # Creating a Bokeh figure
    p = figure(width=400, height=400,toolbar_location=None, x_range=(0, 1), y_range=(0, 1))
    p.image_url(url=[image_url], x=0, y=1, w=1, h=1)

    # Hide grid lines and axes
    p.xgrid.grid_line_color = None
    p.ygrid.grid_line_color = None
    p.xaxis.visible = False
    p.yaxis.visible = False

    return p



image_url = "/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Bokeh/molecule_image.png"

show(mol_web_show(image_url))

def final_layout(p1, p2):
    layout = row(p1,p2)
    show(layout)


def subgroup_nmr_value (mol_smi, dict_functional_groups):

    list_contained_subgroups, list_contained_subgroups_values = [],[]
    mol = Chem.MolFromSmiles(mol_smi)

    for SMILES, value in dict_functional_groups.items():
        substruct = Chem.MolFromSmiles(SMILES)

        if mol.HasSubstructMatch(substruct):

            list_contained_subgroups.append(SMILES)
            list_contained_subgroups_values.append(value)
            indices = mol.GetSubstructMatch(substruct)
            print(indices)

        else:

            list_contained_subgroups.append(f'NOT {SMILES}')

    return list_contained_subgroups,list_contained_subgroups_values


a = 1.0
x1,x2=-0.1+a,0.1+a

def delta(x,eps):
    return  1.0/(2.0*eps*np.cosh(x/eps)**2)
delta=np.vectorize(delta)
# plotting delta functions for different eps values
eps=0.0002
x=np.linspace(-100,100,1000)
y=delta(x,eps)


# Create a new plot with a title and axis labels
p = figure(title="Sharp delta Distribution", x_axis_label='x', y_axis_label='Probability Density')



# Add a line renderer
p.line(x, y, line_width=2)

show(p)