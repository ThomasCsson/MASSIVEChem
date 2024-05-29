import MASSiveChem.MASSiveChem as MC
import unittest

class TestFunctionalGroupFinder(unittest.TestCase):
    def test_functional_group_finder_empty(self):

        smiles = ""

        with self.assertRaises(ValueError):
            MC.functional_group_finder(smiles)

    def test_functional_group_finder_incorrect_smiles(self):

        smiles = "CCXrtCC"

        with self.assertRaises(ValueError):
            MC.functional_group_finder(smiles)

    def test_functional_group_finder_no_functional_groups(self):
        
        mol_smi = "CC"

        functional_groups = MC.functional_group_finder(mol_smi)

        self.assertEqual(functional_groups, [])

    def test_functional_group_finder_single_functional_group(self):
        
        mol_smi = "CCO"

        functional_groups = MC.functional_group_finder(mol_smi)

        self.assertEqual(functional_groups, ['Alcohol'])

    def test_functional_group_finder_multiple_functional_groups(self):
        
        mol_smi = "CCOCC(=O)OC"

        functional_groups = MC.functional_group_finder(mol_smi)

        self.assertEqual(functional_groups, ['Ester', 'Ether'])

    def test_functional_group_finder_duplicate_functional_groups(self):
        
        mol_smi = "NCCCOCC1CC(C=O)CC(CCC=NC)C1"

        functional_groups = MC.functional_group_finder(mol_smi)

        self.assertEqual(functional_groups, ['Aldehyde', 'Ether', 'Amine', 'Amine', 'Imine'])

if __name__ == '__main__':
    unittest.main()
