import MASSiveChem.MASSiveChem as MC
import unittest

class TestSpectrum3D(unittest.TestCase):
    
    def setUp(self):
        self.apparatus_resolution = 0.01

    def test_butane_False(self):
        mol_smi = "CCCC"
        imprecision_True_False = False
        result = MC.spectrum_3D(mol_smi, imprecision_True_False, self.apparatus_resolution)
        self.assertIsNotNone(result)
    
    def test_butane_True(self):
        mol_smi = "CCCC"
        imprecision_True_False = True
        result = MC.spectrum_3D(mol_smi, imprecision_True_False, self.apparatus_resolution)
        self.assertIsNotNone(result)

    def test_empty_input(self):
        
        with self.assertRaises(ValueError):
            imprecision_True_False= False
            MC.spectrum_3D("",imprecision_True_False, self.apparatus_resolution)

    def test_invalid_input(self):
        
        with self.assertRaises(ValueError):
            imprecision_True_False = False
            MC.spectrum_3D("invalid_smiles",imprecision_True_False, self.apparatus_resolution)  
        

if __name__ == '__main__':
    unittest.main()


