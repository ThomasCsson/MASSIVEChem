from rdkit import Chem

def atom_present_list(mol_smi):
    #---------------------------------------------------------------------------------------------#
    '''
    atom_present_list(mol_smi)
    
    Input: molecule under MOL representation
    
    Output: list of atoms present in the molecule with their respective count
    '''
    #---------------------------------------------------------------------------------------------#
    mol, list_atoms,list_atoms = Chem.MolFromSmiles(mol_smi),[],[]
    for atom in mol.GetAtoms():
        list_atoms.append(atom.GetSymbol())
    for atom in mol.GetAtoms():
        element = (f'{atom.GetSymbol()} : {list_atoms.count(atom.GetSymbol())}')   
        if element not in list_atoms:
            list_atoms.append(element)
    return list_atoms

import unittest
from rdkit import Chem

class TestAtomPresentList(unittest.TestCase):

    def test_single_atom(self):
        self.assertEqual(atom_present_list('C'), ['C', 'C : 1'])

    def test_multiple_same_atoms(self):
        self.assertEqual(atom_present_list('CC'),['C', 'C', 'C : 2'])

    def test_different_atoms(self):
        self.assertEqual(atom_present_list('CO'), ['C', 'O', 'C : 1', 'O : 1'])

    def test_more_complex_molecule(self):
        self.assertEqual(atom_present_list('CCO'), ['C', 'C', 'O', 'C : 2', 'O : 1'])

    def test_benzene(self):
        self.assertEqual(atom_present_list('c1ccccc1'), ['C', 'C', 'C', 'C', 'C', 'C', 'C : 6'])

    def test_with_nitrogen(self):
        self.assertEqual(atom_present_list('CCN'), ['C', 'C', 'N', 'C : 2', 'N : 1'])

    def test_with_halogen(self):
        self.assertEqual(atom_present_list('CCCl'), ['C', 'C', 'Cl', 'C : 2', 'Cl : 1'])

if __name__ == '__main__':
    unittest.main()
