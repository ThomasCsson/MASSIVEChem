from rdkit import Chem
from rdkit.Chem import Draw, AllChem

from io import BytesIO
import base64

from bokeh.plotting import show
from bokeh.io import show
from bokeh.models import Div


def mol_web_show(mol_smi):

    #---------------------------------------------------------------------------------------------#
    '''
    mol_web_show(mol_smi)
    
    Input: SMILEs of a molecule
    
    Output: image of the molecule as a bokeh plot
    '''
    #---------------------------------------------------------------------------------------------#

    # Generate the image from the molecule
    mol = Chem.MolFromSmiles(mol_smi)

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

        mol_smi = "CCO"
         
        result = mol_web_show(mol_smi)

        self.assertIsInstance(result, Div)

    def test_small_smiles(self):
        mol_smi = "CCO"  

        result = mol_web_show(mol_smi)

        self.assertIsInstance(result, Div)

        mol = Chem.MolFromSmiles(mol_smi)
        
        image = Draw.MolToImage(mol)
        buffered = BytesIO()
        image.save(buffered, format="PNG")
        image_base64_with_hs = base64.b64encode(buffered.getvalue()).decode("utf-8")

        self.assertIn(image_base64_with_hs[:100], result.text)  

    def test_invalid_smiles(self):

        mol_smi = "InvalidSMILES"

        with self.assertRaises(ValueError):
            mol_web_show(mol_smi)

    def test_complex_smiles(self):

        mol_smi = "C1=CC=C(C=C1)C2=CC=CC=C2" 

        result = mol_web_show(mol_smi)

        self.assertIsInstance(result, Div)
    
    def test_benzene(self):

        mol_smi = "C1=CC=CC=C1" 

        result = mol_web_show(mol_smi)

        self.assertIsInstance(result, Div)

        mol = Chem.MolFromSmiles(mol_smi)

        image = Draw.MolToImage(mol)
        buffered = BytesIO()
        image.save(buffered, format="PNG")
        image_base64_with_hs = base64.b64encode(buffered.getvalue()).decode("utf-8")

        self.assertIn(image_base64_with_hs[:100], result.text) 


if __name__ == '__main__':
    unittest.main()
