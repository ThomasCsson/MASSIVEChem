import MASSiveChem.MASSiveChem as MC
import unittest

class TestSmilesTo3DPlot(unittest.TestCase):
    def test_valid_smiles(self):
        
        mol_smi = "CCO" 

        result = MC.smiles_to_3D_plot(mol_smi)

        self.assertIsNotNone(result)

    def test_empty_smiles(self):
        
        with self.assertRaises(ValueError):
            MC.smiles_to_3D_plot("")  

    def test_invalid_smiles(self):
        
        with self.assertRaises(ValueError):
            MC.smiles_to_3D_plot("invalid_smiles")    

    def test_plot_data(self):
        
        mol_smi = "CCO"  

        result = MC.smiles_to_3D_plot(mol_smi)

        self.assertIn('data', result)  

    def test_plot_layout(self):
        
        mol_smi = "CCO"  

        result = MC.smiles_to_3D_plot(mol_smi)

        self.assertIn('layout', result)  

    def test_plot_layout_attributes(self):
        
        mol_smi = "CCO" 

        result = MC.smiles_to_3D_plot(mol_smi)
        layout = result['layout']
        
        self.assertIn('title', layout)  
        self.assertIn('scene', layout)
if __name__ == '__main__':
    unittest.main()


