from rdkit import Chem
from rdkit.Chem import Draw
def functional_group_finder(mol_smi):

    #---------------------------------------------------------------------------------------------#
    '''
    functional_group_finder(mol_smi)
    
    Input: molecule under SMILEs representation
    
    Output: list containing every functionl group contained (if a functional group is contained twice in the molecule, it will appear twice in this list)
    '''
    #---------------------------------------------------------------------------------------------#

    # initiate variables
    functional_groups_contained, mol_in = [], Chem.MolFromSmiles(mol_smi)

    # dictionnary of all the considered functional groups to check (some might be missing)

    functional_groups_smarts = {
        'Alcohol': 'C[Oh1+0]',
        'Aldehyde': 'C[Ch1]=O',
        'Ketone': 'CC(=O)C',
        'Carboxylic Acid': 'CC(=O)[Oh1]',
        'Ester': 'CC(=O)[Oh0]',
        'Ether': '*[Oh0]*',
        'Amide': 'C(=O)N',
        'Amine': '[C][Nh2]',
        'Nitrile': 'C#N',
        'Chloride': 'Cl',
        'Bromide': 'Br',
        'Fluoride': 'F',
        'Iodide': 'I',
        'Alkene': 'C=C',
        'Alkyne': 'C#C',
        'Imine': 'C=N*',
        'Amino acid': '[Nh2][Ch1*]C(=O)O',
        'Proline': '[Nh1][Ch1*]C(=O)O',
        'Thiol': '[Sh1]',
        'Sulfide': '*[Sh0]*',
        'Acyl Chloride': 'CC(=O)Cl',
        'Anhydride': '*[Ch0](=O)O[Ch0](=O)*',
        'Nitro': 'C[N+](=O)[O-]',
        'Enamine': 'C=C[Nh0]',
        'Enamine2': 'C=C[Nh1]',
        'Enamine3': 'C=C[Nh2]',
        'Imide': 'C(=O)NC(=O)*',
        'Azide': 'CNNN',
        'Enol': 'C=C([Oh1])C',
        'Hemiacetal': 'CC(O)(O)C',
        'Carbonate': '[Oh0]C(=O)[Oh0]',
        'Carbonate2': '[Oh1]C(=O)[Oh1]',
        'Disulfide': 'CSSC',
        'Sulfoxide': 'CS(=O)C',
        'Sulfone': '*[So2](=O)(=O)*',
        'Sulfonic acid': '*S(=O)(=O)[Oh1]',
        'Thioester': 'C(=O)S*',
        'Phosphine': '*[Po0](*)*',
        'Phosphate': '*OP(=O)(O)O',
        'Benzene': 'c1ccccc1',
        'Peroxide':'C[Oh0][Oh0]C'
    }

    # check that the substructure from functional_groups_smarts are contained in mol_smi

    for name, smarts in functional_groups_smarts.items():
        if mol_in.HasSubstructMatch(Chem.MolFromSmarts(smarts)):
            for _ in range(len(mol_in.GetSubstructMatches(Chem.MolFromSmarts(smarts)))):
                functional_groups_contained.append(name)

    # exceptions for conflicts during the iteration of functional groups

    if 'Ester' and 'Anhydride' in functional_groups_contained:
        functional_groups_contained.remove('Ether')
    if 'Carboxylic Acid' in functional_groups_contained:
        functional_groups_contained.remove('Alcohol')
    if 'Ester' in functional_groups_contained:
        functional_groups_contained.remove('Ether')
    if 'Phosphate' in functional_groups_contained:
        functional_groups_contained.remove('Ether')
    if 'Thioester' in functional_groups_contained:
        functional_groups_contained.remove('Sulfide')
    if 'Sulfonic acid' in functional_groups_contained:
        functional_groups_contained.remove('Sulfide')
    if 'Sulfoxide' in functional_groups_contained:
        functional_groups_contained.remove('Sulfide')
    if 'Acyl Chloride' in functional_groups_contained:
        functional_groups_contained.remove('Chloride')
    if 'Anhydride' in functional_groups_contained:
        functional_groups_contained.remove('Ester')
        functional_groups_contained.remove('Ester')
    if 'Enamine2' in functional_groups_contained:
        functional_groups_contained.remove('Enamine2')
        functional_groups_contained.append('Enamine')
    if 'Enamine3' in functional_groups_contained:
        functional_groups_contained.remove('Enamine3')
        functional_groups_contained.remove('Amine')
        functional_groups_contained.append('Enamine')
    if 'Imide' in functional_groups_contained:
        functional_groups_contained.remove('Amide')
        functional_groups_contained.remove('Amide')
    if 'Enol' in functional_groups_contained:
        functional_groups_contained.remove('Alkene')
        functional_groups_contained.remove('Alcohol')
    if 'Hemiacetal' in functional_groups_contained:
        functional_groups_contained.remove('Alcohol')
        functional_groups_contained.remove('Alcohol')
    if 'Carbonate2' in functional_groups_contained:
        functional_groups_contained.remove('Alcohol')
        functional_groups_contained.remove('Alcohol')
        functional_groups_contained.remove('Carbonate2')
        functional_groups_contained.append('Carbonate')
    if 'Disulfide' in functional_groups_contained:
        functional_groups_contained.remove('Sulfide')
        functional_groups_contained.remove('Sulfide')
    if 'Amine2' in functional_groups_contained:
        functional_groups_contained.remove('Amine2')
        functional_groups_contained.append('Amine')
    if 'Peroxide' in functional_groups_contained:
        functional_groups_contained.remove('Ether')
        functional_groups_contained.remove('Ether')
    
    return functional_groups_contained


<<<<<<< HEAD

import unittest

class TestFunctionalGroupFinder(unittest.TestCase):
    def test_functional_group_finder_empty(self):
        # Test the function with an empty SMILES string
        mol_smi = ""
        functional_groups = functional_group_finder(mol_smi)
        self.assertEqual(functional_groups, [])

    def test_functional_group_finder_no_functional_groups(self):
        # Test the function with a molecule without any functional groups
        mol_smi = "CC"
        functional_groups = functional_group_finder(mol_smi)
        self.assertEqual(functional_groups, [])

    def test_functional_group_finder_single_functional_group(self):
        # Test the function with a molecule containing a single functional group
        mol_smi = "CCO"
        functional_groups = functional_group_finder(mol_smi)
        self.assertEqual(functional_groups, ['Alcohol'])

    def test_functional_group_finder_multiple_functional_groups(self):
        # Test the function with a molecule containing multiple functional groups
        mol_smi = "CCOCC(=O)OC"
        functional_groups = functional_group_finder(mol_smi)
        self.assertEqual(functional_groups, ['Ester', 'Ether'])

    def test_functional_group_finder_duplicate_functional_groups(self):
        # Test the function with a molecule containing duplicate functional groups
        mol_smi = "CCOCC(=O)OCCC(=O)OC(=O)C"
        functional_groups = functional_group_finder(mol_smi)
        self.assertEqual(functional_groups, ['Ether', 'Ester', 'Anhydryde'])

if __name__ == '__main__':
    unittest.main()
=======
print(functional_group_finder('CCCC(=O)OCCC(=O)O'))
>>>>>>> b878f9398c1dbb1cf1a7b59c79fd9441e9f76da2
