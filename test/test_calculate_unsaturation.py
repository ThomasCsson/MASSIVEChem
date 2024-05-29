from rdkit import Chem
import MASSiveChem.MASSiveChem as MC
"""
def calculate_unsaturation(mol_smile) -> int:
    #---------------------------------------------------------------------------------------------#
    '''
    calculate_unsaturation(mol_smile)

    Input: molecule under SMILEs representation

    Output: unsaturation of the input molecule (integer value)
    '''
    #---------------------------------------------------------------------------------------------# 
    if mol_smile == None:
        raise ValueError('Enter a non-empty input')
    
    C, N, HX, halogens_hydrogen = 0, 0, 0, ['F', 'Cl', 'Br', 'I', 'At', 'H']

    mol = Chem.AddHs(Chem.MolFromSmiles(mol_smile))
    
    if mol == None:
        raise ValueError('Invalid SMILEs representation')
    
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':
            C += 1
        elif atom.GetSymbol() == 'N':
            N += 1
        elif atom.GetSymbol() in halogens_hydrogen:
            HX += 1

    unsaturation = C + 1 + (N - HX) / 2

    return unsaturation
"""
import unittest

class TestCalculateUnsaturation(unittest.TestCase):

    def test_one_unsaturation(self):

        result = MC.calculate_unsaturation("C=C")  
        self.assertEqual(result, 1)

    def test_nitrogen_and_no_halogens(self):

        result = MC.calculate_unsaturation("CC=O")  
        self.assertEqual(result, 1)

    def test_nitrogen_and_halogens(self):

        result = MC.calculate_unsaturation("CNC(Cl)Br")  
        self.assertEqual(result, 0)

    def test_multiple_nitrogen_and_halogens(self):

        result = MC.calculate_unsaturation("BrC=CC(Br)(Cl)CC(N)CCN") 
        self.assertEqual(result, 1)

    def test_multiple_unsaturations(self):

        result = MC.calculate_unsaturation("C1=CC=CC=C1") 
        self.assertEqual(result, 4)

    def test_no_atoms(self):

        result = MC.calculate_unsaturation("") 
        self.assertEqual(result, 1)

    def test_no_nitrogen_and_halogens(self):

        result = MC.calculate_unsaturation("CCCl")
        self.assertEqual(result, 0)

if __name__ == '__main__':
    unittest.main()
