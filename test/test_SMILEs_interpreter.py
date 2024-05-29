import MASSiveChem.MASSiveChem as MC
from rdkit import Chem
import unittest

class TestSMILEsInterpreter(unittest.TestCase):

    def test_valid_smiles(self):

        smiles = "CCO"
        mol = MC.SMILEs_interpreter(smiles)

        self.assertIsNotNone(mol)
        self.assertEqual(Chem.MolToSmiles(Chem.RemoveHs(mol)), smiles)
    
    def test_invalid_smiles(self):

        smiles = "C1CCxVCC"

        with self.assertRaises(SystemExit):
            MC.SMILEs_interpreter(smiles)
    
    def test_simple_molecule(self):

        smiles = "C"
        mol = MC.SMILEs_interpreter(smiles)

        self.assertIsNotNone(mol)
        self.assertEqual(Chem.MolToSmiles(Chem.RemoveHs(mol)), smiles)
    
    def test_molecule_with_explicit_hydrogens(self):

        smiles = "C"
        mol = MC.SMILEs_interpreter(smiles)

        self.assertIsNotNone(mol)
        self.assertEqual(Chem.MolToSmiles(mol), "[H]C([H])([H])[H]")
    

if __name__ == "__main__":
    unittest.main()