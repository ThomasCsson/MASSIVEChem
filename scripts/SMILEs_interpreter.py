from rdkit import Chem

def SMILEs_interpreter(mol_smi):
    #---------------------------------------------------------------------------------------------#
    '''
    SMILEs_interpreter(mol_smi)

    Input: molecule under SMILEs representation
    
    Output: molecule under MOL representation
    '''
    #---------------------------------------------------------------------------------------------#

    #Checks that the SMILEs input is correct/allowed

    mol_without_Hs = Chem.MolFromSmiles(mol_smi)

    if mol_without_Hs is None:
        print('')
        print("Invalid SMILEs input.")
        print('Please try again with a different SMILEs.')
        exit()

    mol = Chem.AddHs(mol_without_Hs)

    return mol

import unittest

class TestSMILEsInterpreter(unittest.TestCase):
    def test_valid_smiles(self):
        smiles = "CCO"
        mol = SMILEs_interpreter(smiles)
        self.assertIsNotNone(mol)
        self.assertEqual(Chem.MolToSmiles(Chem.RemoveHs(mol)), smiles)
    
    """def test_invalid_smiles(self):
        smiles = "C1CCC1"
        with self.assertRaises(SystemExit):
            SMILEs_interpreter(smiles)"""
    
    def test_simple_molecule(self):
        smiles = "C"
        mol = SMILEs_interpreter(smiles)
        self.assertIsNotNone(mol)
        self.assertEqual(Chem.MolToSmiles(Chem.RemoveHs(mol)), smiles)
    
    def test_molecule_with_explicit_hydrogens(self):
        smiles = "C"
        mol = SMILEs_interpreter(smiles)
        self.assertIsNotNone(mol)
        self.assertEqual(Chem.MolToSmiles(mol), "[H]C([H])([H])[H]")
    
    def test_complex_molecule(self):
        smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"
        mol = SMILEs_interpreter(smiles)
        self.assertIsNotNone(mol)
        self.assertEqual(Chem.MolToSmiles(mol), smiles)

if __name__ == "__main__":
    unittest.main()