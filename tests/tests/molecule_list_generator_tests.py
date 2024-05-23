import unittest
from rdkit import Chem

def molecule_list_generator(mol) -> list[str]:
    list_atoms = []
    for atom in mol.GetAtoms():
        list_atoms.append(atom.GetSymbol())
    return list_atoms

class TestMoleculeListGenerator(unittest.TestCase):

    def test_benzene(self):
        mol_smi = "C1=CC=CC=C1"  # Benzene
        mol = Chem.MolFromSmiles(mol_smi)
        result = molecule_list_generator(mol)
        expected = ['C', 'C', 'C', 'C', 'C', 'C']
        self.assertEqual(result, expected)

    def test_water(self):
        mol_smi = "O"  # Water
        mol = Chem.MolFromSmiles(mol_smi)
        result = molecule_list_generator(mol)
        expected = ['O']
        self.assertEqual(result, expected)

if __name__ == '__main__':
    unittest.main()
