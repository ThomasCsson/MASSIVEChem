from xyz2graph import MolGraph, to_networkx_graph, to_plotly_figure
from plotly.offline import offline

from rdkit import Chem
from rdkit.Chem import AllChem

import plotly.graph_objects as go
from bokeh.io import output_file, show
from bokeh.layouts import layout
from bokeh.models import Div
from plotly.io import to_html
import pandas as pd

def smiles_to_xyz(smiles, output_file):
    # Generate 3D coordinates from SMILES string
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)  # Add hydrogens for better geometry optimization
    AllChem.EmbedMolecule(mol, randomSeed=42)  # Embed the molecule in 3D space
    AllChem.MMFFOptimizeMolecule(mol)  # Optimize the geometry using MMFF94 force field

    # Write coordinates to XYZ file
    with open(output_file, 'w') as f:
        f.write(str(mol.GetNumAtoms()) + '\n\n')  # Write number of atoms
        for atom in mol.GetAtoms():
            pos = mol.GetConformer().GetAtomPosition(atom.GetIdx())
            f.write('{} {} {} {}\n'.format(atom.GetSymbol(), pos.x, pos.y, pos.z))

# Example usage:
smiles = "CCC(=O)OOCCC"  # Ethanol molecule
output_file = "ethanol.xyz"
smiles_to_xyz(smiles, output_file)


# Create the MolGraph object
mg = MolGraph()

# Read the data from the .xyz file
mg.read_xyz('ethanol.xyz')

# Create the Plotly figure object
fig = to_plotly_figure(mg)

# Plot the figure
#offline.plot(fig)
import panel as pn

plotly_pane = pn.pane.Plotly(fig)


from rdkit import Chem
from rdkit.Chem import Draw, AllChem

from io import BytesIO
import base64

from bokeh.plotting import show
from bokeh.io import show
from bokeh.models import Div


def mol_web_show(mol_smi, show_Hs=False, show_3D = False):

    #---------------------------------------------------------------------------------------------#
    '''
    mol_web_show(mol_smi, show_Hs=False, show_3D = False)
    
    Input: SMILEs of a molecule. Also specify if want the function to show the hydrogens explicitely or the 3D
    
    Output: image of the molecule as a bokeh plot
    '''
    #---------------------------------------------------------------------------------------------#

    # Generate the image from the molecule
    mol = Chem.MolFromSmiles(mol_smi)

    # Show the molecule in 3D if specified
    if show_3D:
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
    
    # Adds the hydrogens to the molecule if specified
    if show_Hs:
        mol = Chem.AddHs(mol)

    #Draws the image
    image = Draw.MolToImage(mol)

    #stocks the image in a base64 format
    buffered = BytesIO()
    image.save(buffered, format="PNG")
    image_base64 = base64.b64encode(buffered.getvalue()).decode("utf-8")
    image_url = f"data:image/png;base64,{image_base64}"

    # Create a Div element to display the image
    img_div = Div(text=f'<img src="{image_url}" style="width:350px;height:350px;">')
   
    return img_div
    
input_mol = input('MOL:  ')

bokeh =mol_web_show(input_mol,False,True)

bokeh_pane = pn.pane.Bokeh(bokeh)

layout = pn.Row(bokeh_pane,plotly_pane)



