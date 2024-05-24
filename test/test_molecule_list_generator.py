import unittest
from rdkit import Chem
def molecule_list_generator(mol) -> list[str]:
    #---------------------------------------------------------------------------------------------#
    '''
    molecule_list_generator(mol)
    
    Input: molecule under MOL representation
    
    Output: list containing the atomic symbol of each atom in the input molecule
    '''
    #---------------------------------------------------------------------------------------------#
    if not mol:
        raise ValueError('Enter a non-empty input')
    if not Chem.MolToSmiles(mol):
        raise ValueError('Wrong format or mol type')
    
    list_atoms = []
    for atom in mol.GetAtoms():
        list_atoms.append(atom.GetSymbol())
    return list_atoms



class TestMoleculeListGenerator(unittest.TestCase):

    def test_benzene(self):

        mol_smi = "C1=CC=CC=C1"  
        mol = Chem.MolFromSmiles(mol_smi)

        result = molecule_list_generator(mol)

        expected = ['C', 'C', 'C', 'C', 'C', 'C']

        self.assertEqual(result, expected)

    def test_water(self):

        mol_smi = "O"  
        mol = Chem.MolFromSmiles(mol_smi)

        result = molecule_list_generator(mol)

        expected = ['O']
        
        self.assertEqual(result, expected)

    def test_no_molecule(self):

        with self.assertRaises(ValueError):
            molecule_list_generator('') 
    

if __name__ == '__main__':
    unittest.main()
