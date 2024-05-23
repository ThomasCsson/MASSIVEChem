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

import unittest


class TestMolWebShow(unittest.TestCase):

    def test_mol_web_show_no_Hs_no_3D(self):
        # Test with no explicit hydrogens and no 3D
        mol_smi = "CCO"  # Example SMILES string
        result = mol_web_show(mol_smi, show_Hs=False, show_3D=False)
        self.assertIsInstance(result, Div)  # Check if the output is a Div element

    def test_mol_web_show_with_Hs_no_3D(self):
        # Test with explicit hydrogens but no 3D
        mol_smi = "CCO"  # Example SMILES string
        result = mol_web_show(mol_smi, show_Hs=True, show_3D=False)
        self.assertIsInstance(result, Div)  # Check if the output is a Div element

    def test_mol_web_show_no_Hs_with_3D(self):
        # Test with no explicit hydrogens but with 3D
        mol_smi = "CCO"  # Example SMILES string
        result = mol_web_show(mol_smi, show_Hs=False, show_3D=True)
        self.assertIsInstance(result, Div)  # Check if the output is a Div element

    def test_mol_web_show_with_Hs_with_3D(self):
        # Test with explicit hydrogens and with 3D
        mol_smi = "CCO"  # Example SMILES string
        result = mol_web_show(mol_smi, show_Hs=True, show_3D=True)
        self.assertIsInstance(result, Div)  # Check if the output is a Div element

if __name__ == '__main__':
    unittest.main()
