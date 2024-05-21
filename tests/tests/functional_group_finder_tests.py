from rdkit import Chem

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
        'Amine': '[C][N]',
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
    for functional_group in functional_groups_contained:
        if 'Ester' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                functional_groups_contained.remove('Ether')
        elif functional_group == 'Carboxylic Acid':
            for _ in range (functional_groups_contained.count(functional_group)):
                functional_groups_contained.remove('Alcohol')
        elif 'Ester' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                functional_groups_contained.remove('Ether')
        elif 'Phosphate' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                functional_groups_contained.remove('Ether')
        elif 'Thioester' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                functional_groups_contained.remove('Sulfide')
        elif 'Sulfonic acid' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                functional_groups_contained.remove('Sulfide')
        elif 'Sulfoxide' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                functional_groups_contained.remove('Sulfide')
        elif 'Acyl Chloride' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                functional_groups_contained.remove('Chloride')
        elif 'Anhydride' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                functional_groups_contained.remove('Ester')
                functional_groups_contained.remove('Ester')
                functional_groups_contained.append('Ether')
        elif 'Enamine2' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                functional_groups_contained.remove('Enamine2')
                functional_groups_contained.append('Enamine')
        elif 'Enamine3' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                functional_groups_contained.remove('Enamine3')
                functional_groups_contained.remove('Amine')
                functional_groups_contained.append('Enamine')
        elif 'Imide' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                functional_groups_contained.remove('Amide')
                functional_groups_contained.remove('Amide')
        elif 'Enol' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                functional_groups_contained.remove('Alkene')
                functional_groups_contained.remove('Alcohol')
        elif 'Hemiacetal' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                functional_groups_contained.remove('Alcohol')
                functional_groups_contained.remove('Alcohol')
        elif 'Carbonate2' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                functional_groups_contained.remove('Alcohol')
                functional_groups_contained.remove('Alcohol')
                functional_groups_contained.remove('Carbonate2')
                functional_groups_contained.append('Carbonate')
        elif 'Disulfide' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                functional_groups_contained.remove('Sulfide')
                functional_groups_contained.remove('Sulfide')
        elif 'Peroxide' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                functional_groups_contained.remove('Ether')
                functional_groups_contained.remove('Ether')
    
    
    return functional_groups_contained

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
        self.assertEqual(functional_groups, ['Ester', 'Ether', 'Anhydride'])

if __name__ == '__main__':
    unittest.main()





print(functional_group_finder('CCCC(=O)OCCC(=O)O'))