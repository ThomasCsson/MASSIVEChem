import MASSiveChem.MASSiveChem as MC
import unittest
from rdkit import Chem

class TestMoleculeListGenerator(unittest.TestCase):

    def test_benzene(self):

        mol_smi = "C1=CC=CC=C1"  
        mol = Chem.MolFromSmiles(mol_smi)

        result = MC.molecule_list_generator(mol)

        expected = ['C', 'C', 'C', 'C', 'C', 'C']

        self.assertEqual(result, expected)

    def test_water(self):

        mol_smi = "O"  
        mol = Chem.MolFromSmiles(mol_smi)

        result = MC.molecule_list_generator(mol)

        expected = ['O']
        
        self.assertEqual(result, expected)

    def test_no_molecule(self):

        with self.assertRaises(ValueError):
            MC.molecule_list_generator('') 
    

if __name__ == '__main__':
    unittest.main()
