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
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from io import BytesIO
import base64
from bokeh.models import Div

class TestMolWebShow(unittest.TestCase):

    def test_default_parameters(self):
        mol_smi = "CCO"  # Ethanol
        result = mol_web_show(mol_smi)
        self.assertIsInstance(result, Div)

    def test_show_Hs(self):
        mol_smi = "CCO"  # Ethanol
        result = mol_web_show(mol_smi, show_Hs=True)
        self.assertIsInstance(result, Div)
        mol = Chem.MolFromSmiles(mol_smi)
        mol_with_hs = Chem.AddHs(mol)
        image_with_hs = Draw.MolToImage(mol_with_hs)
        buffered_with_hs = BytesIO()
        image_with_hs.save(buffered_with_hs, format="PNG")
        image_base64_with_hs = base64.b64encode(buffered_with_hs.getvalue()).decode("utf-8")
        self.assertIn(image_base64_with_hs[:100], result.text)  # Partial match for simplicity

    def test_invalid_smiles(self):
        mol_smi = "InvalidSMILES"
        with self.assertRaises(ValueError):
            mol_web_show(mol_smi)

    def test_complex_smiles(self):
        mol_smi = "C1=CC=C(C=C1)C2=CC=CC=C2"  # Biphenyl
        result = mol_web_show(mol_smi)
        self.assertIsInstance(result, Div)
    
    def test_benzene_show_Hs(self):
        mol_smi = "C1=CC=CC=C1"  # Benzene
        result = mol_web_show(mol_smi, show_Hs=True)
        self.assertIsInstance(result, Div)
        mol = Chem.MolFromSmiles(mol_smi)
        mol_with_hs = Chem.AddHs(mol)
        image_with_hs = Draw.MolToImage(mol_with_hs)
        buffered_with_hs = BytesIO()
        image_with_hs.save(buffered_with_hs, format="PNG")
        image_base64_with_hs = base64.b64encode(buffered_with_hs.getvalue()).decode("utf-8")
        self.assertIn(image_base64_with_hs[:100], result.text)  # Partial match for simplicity


if __name__ == '__main__':
    unittest.main()
