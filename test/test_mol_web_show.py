from rdkit import Chem
from rdkit.Chem import Draw
from io import BytesIO
import base64
from bokeh.models import Div
import MASSiveChem.MASSiveChem as MC
import unittest


class TestMolWebShow(unittest.TestCase):

    def test_default_parameters(self):

        mol_smi = "CCO"
         
        result = MC.mol_web_show(mol_smi)

        self.assertIsInstance(result, Div)

    def test_small_smiles(self):
        mol_smi = "CCO"  

        result = MC.mol_web_show(mol_smi)

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
            MC.mol_web_show(mol_smi)

    def test_empty_smiles(self):

        mol_smi = ""

        with self.assertRaises(ValueError):
            MC.mol_web_show(mol_smi)

    def test_complex_smiles(self):

        mol_smi = "C1=CC=C(C=C1)C2=CC=CC=C2" 

        result = MC.mol_web_show(mol_smi)

        self.assertIsInstance(result, Div)
    
    def test_benzene(self):

        mol_smi = "C1=CC=CC=C1" 

        result = MC.mol_web_show(mol_smi)

        self.assertIsInstance(result, Div)

        mol = Chem.MolFromSmiles(mol_smi)

        image = Draw.MolToImage(mol)
        buffered = BytesIO()
        image.save(buffered, format="PNG")
        image_base64_with_hs = base64.b64encode(buffered.getvalue()).decode("utf-8")

        self.assertIn(image_base64_with_hs[:100], result.text) 


if __name__ == '__main__':
    unittest.main()
