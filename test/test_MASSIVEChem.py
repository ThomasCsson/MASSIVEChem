
#This file contains the tests for the entire package 

import MASSiveChem.MASSiveChem as MC
import unittest
from rdkit import Chem
from rdkit.Chem import Draw
from io import BytesIO
import base64
from bokeh.models import Div

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

class TestDeltaFunctionPlotter(unittest.TestCase):

    def test_single_peak(self):
        
        x_in = [1.0]
        y_in = [0.5]

        expected_x = [0.5, 1.0 - 10**(-100), 1.0, 1.0 + 10**(-100), 2.0]
        expected_y = [0, 0, 0.5, 0, 0]

        x_axis, y_axis = MC.delta_function_plotter(x_in, y_in)

        self.assertEqual(x_axis, expected_x)
        self.assertEqual(y_axis, expected_y)

    def test_two_peaks(self):

        x_in = [1.0, 2.0]
        y_in = [0.5, 0.8]

        expected_x = [
            0.5, 1.0 - 10**(-100), 1.0, 1.0 + 10**(-100), 
            2.0 - 10**(-100), 2.0, 2.0 + 10**(-100), 3.0
        ]
        expected_y = [0, 0, 0.5, 0, 0, 0.8, 0, 0]

        x_axis, y_axis = MC.delta_function_plotter(x_in, y_in)

        self.assertEqual(x_axis, expected_x)
        self.assertEqual(y_axis, expected_y)

    def test_three_peaks(self):

        x_in = [1.0, 2.0, 3.0]
        y_in = [0.5, 0.8, 0.2]

        expected_x = [
            0.5, 1.0 - 10**(-100), 1.0, 1.0 + 10**(-100),
            2.0 - 10**(-100), 2.0, 2.0 + 10**(-100),
            3.0 - 10**(-100), 3.0, 3.0 + 10**(-100), 4.0
        ]
        expected_y = [0, 0, 0.5, 0, 0, 0.8, 0, 0, 0.2, 0, 0]

        x_axis, y_axis = MC.delta_function_plotter(x_in, y_in)

        self.assertEqual(x_axis, expected_x)
        self.assertEqual(y_axis, expected_y)

    def test_empty_input_raises_value_error(self):

        with self.assertRaises(ValueError):

            MC.delta_function_plotter([], [])
        
        with self.assertRaises(ValueError):

            MC.delta_function_plotter([1.0], [])
        
        with self.assertRaises(ValueError):

            MC.delta_function_plotter([], [0.5])

    def test_invalid_input_length_mismatch(self):

        with self.assertRaises(ValueError):
            MC.delta_function_plotter([1.0, 2.0], [0.5])

class TestDoublePlot(unittest.TestCase):
    def test_double_plot(self):

        x_in = [100, 200, 300, 400, 500]
        y_in = [0.1, 0.2, 0.3, 0.4, 0.5]

        layout = MC.double_plot(x_in, y_in)

        self.assertIsNotNone(layout)
        self.assertEqual(len(layout.children), 2)  # Check if there are two plots in the layout

    def test_double_plot_empty_data(self):

        x_in = []
        y_in = []

        with self.assertRaises(ValueError):
            MC.double_plot(x_in, y_in)
    
    def test_double_plot_1_empty_data(self):

        x_in = [1]
        y_in = []

        with self.assertRaises(ValueError):
            MC.double_plot(x_in, y_in)
    
    def test_double_plot_1_empty_data_2(self):

        x_in = []
        y_in = [1]

        with self.assertRaises(ValueError):
            MC.double_plot(x_in, y_in)

    def test_double_plot_single_point(self):
        
        x_in = [100]
        y_in = [0.1]

        layout = MC.double_plot(x_in, y_in)

        self.assertIsNotNone(layout)
        self.assertEqual(len(layout.children), 2)  # Check if there are two plots in the layout

    def test_double_plot_identical_data(self):
        
        x_in = [100, 200, 300, 400, 500]
        y_in = [0.1, 0.1, 0.1, 0.1, 0.1]

        layout = MC.double_plot(x_in, y_in)

        self.assertIsNotNone(layout)
        self.assertEqual(len(layout.children), 2)  # Check if there are two plots in the layout

class TestFunctionalGroupDisplay(unittest.TestCase):
    def test_functional_group_display(self):
        groups_list = ['Alcohol', 'Aldehyde', 'Ketone']
        table = MC.functional_group_display(groups_list)

        self.assertIsNotNone(table)
        self.assertEqual(len(table.columns), 2)  # Check if there are two columns
        self.assertEqual(table.columns[0].field, "groups")  # Check if the first column is for functional groups
        self.assertEqual(table.columns[1].field, "images")  # Check if the second column is for images

    def test_functional_group_display_empty_list(self):
        table = MC.functional_group_display([])

        self.assertIsNotNone(table)
        self.assertEqual(len(table.columns), 2)  # Check if there are two columns
        self.assertEqual(len(table.source.data['groups']), 0)  # Check if the number of groups is 0

    def test_functional_group_display_invalid_group(self):
        groups_list = ['Alcohol', 'InvalidGroup']
        table = MC.functional_group_display(groups_list)

        self.assertIsNotNone(table)
        self.assertEqual(len(table.columns), 2)  # Check if there are two columns
        self.assertEqual(len(table.source.data['groups']), 2)  # Check if the number of groups is 1
        self.assertEqual(table.source.data['groups'][0], 'Alcohol')  # Check if the valid group is displayed

    # Additional tests
    def test_functional_group_display_multiple_valid_groups(self):
        groups_list = ['Alcohol', 'Ketone', 'Amine', 'Benzene']
        table = MC.functional_group_display(groups_list)

        self.assertIsNotNone(table)
        self.assertEqual(len(table.columns), 2)
        self.assertEqual(len(table.source.data['groups']), 4)
        self.assertIn('Alcohol', table.source.data['groups'])
        self.assertIn('Ketone', table.source.data['groups'])
        self.assertIn('Amine', table.source.data['groups'])
        self.assertIn('Benzene', table.source.data['groups'])

    def test_functional_group_display_all_invalid_groups(self):
        groups_list = ['InvalidGroup1', 'InvalidGroup2']
        table = MC.functional_group_display(groups_list)

        self.assertIsNotNone(table)
        self.assertEqual(len(table.columns), 2)
        self.assertEqual(len(table.source.data['groups']), 2)

    def test_functional_group_display_single_group(self):
        groups_list = ['Alcohol']
        table = MC.functional_group_display(groups_list)

        self.assertIsNotNone(table)
        self.assertEqual(len(table.columns), 2)
        self.assertEqual(len(table.source.data['groups']), 1)
        self.assertEqual(table.source.data['groups'][0], 'Alcohol')

    def test_functional_group_display_mixed_valid_and_invalid_groups(self):
        groups_list = ['Alcohol', 'InvalidGroup', 'Benzene', 'AnotherInvalidGroup']
        table = MC.functional_group_display(groups_list)

        self.assertIsNotNone(table)
        self.assertEqual(len(table.columns), 2)
        self.assertEqual(len(table.source.data['groups']), 4)
        self.assertIn('Alcohol', table.source.data['groups'])
        self.assertIn('Benzene', table.source.data['groups'])

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

class TestIonisationMethod(unittest.TestCase):

    def test_remove_single_h(self):

        self.assertEqual(MC.ionisation_method(['H', 'C', 'O']), ['C', 'O'])

    def test_remove_multiple_h(self):

        self.assertEqual(MC.ionisation_method(['H', 'H', 'C', 'O']), ['H', 'C', 'O'])

    def test_no_h_in_list(self):

        self.assertEqual(MC.ionisation_method(['C', 'O', 'N']), ['C', 'O', 'N'])

    def test_only_h_in_list(self):

        self.assertEqual(MC.ionisation_method(['H']), [])

    def test_empty_list(self):

        self.assertEqual(MC.ionisation_method([]), [])

    def test_no_removal_for_other_atoms(self):

        self.assertEqual(MC.ionisation_method(['C', 'O', 'N']), ['C', 'O', 'N'])

    def test_remove_h_from_complex_list(self):
        
        self.assertEqual(MC.ionisation_method(['H', 'H', 'C', 'O', 'H']), ['H', 'C', 'O', 'H'])

class TestMainFunction(unittest.TestCase):

    def test_single_atom(self):
        
        list_atoms = ['H']
        imprecision_True_False = False
        
        result1, result2 = MC.main_function(list_atoms, imprecision_True_False)
        
        expected_masses = [2.014, 1.008]
        expected_probabilities = [0.00015, 0.99985]
        
        self.assertEqual((result1, result2), (expected_masses, expected_probabilities))

    def test_single_atom_multiple(self):
        
        list_atoms = ['C', 'C']
        imprecision_True_False = False
        
        result1, result2 = MC.main_function(list_atoms, imprecision_True_False)
        
        expected_masses = [26.007, 25.003, 24.0]
        expected_probabilities = [0.00012100000000000003, 0.021758000000000003, 0.9781210000000002]
        
        self.assertEqual((result1, result2), (expected_masses, expected_probabilities))

    def test_multiple_atoms(self):
        
        list_atoms = ['C', 'H', 'O']
        imprecision_True_False = False
        
        result1, result2 = MC.main_function(list_atoms, imprecision_True_False)
        
        expected_masses = [33.017, 32.016, 31.012, 32.01, 31.01, 30.006, 32.013, 31.013, 30.009, 31.007, 30.007, 29.003]
        expected_probabilities =  [3.3e-09, 6.600000000000001e-10, 1.6460400000000003e-06, 2.1996700000000002e-05, 4.399340000000001e-06, 0.010971953960000001, 2.967e-07, 5.934e-08, 0.00014799396, 0.0019777033, 0.00039554066000000007, 0.9864784060400001]
        
        self.assertEqual((result1, result2), (expected_masses, expected_probabilities))

    def test_odd_count_nitrogen(self):

        list_atoms = ['N', 'N', 'N']
        imprecision_True_False = False
        
        result1, result2 = MC.main_function(list_atoms, imprecision_True_False)
        
        expected_masses = [45.0, 44.003, 43.006, 42.009, 43.0]
        expected_probabilities = [5.0653000000000004e-08, 4.0918041e-05, 0.011018011958999999, 0.9889410193469999, 0.0103838807031435]
        
        self.assertEqual((result1, result2), (expected_masses, expected_probabilities))

    def test_even_count_nitrogen(self):

        list_atoms = ['N', 'N']
        imprecision_True_False = False
        
        result1, result2 = MC.main_function(list_atoms, imprecision_True_False)
        
        expected_masses = [30.0, 29.003, 28.006]
        expected_probabilities = [1.3690000000000001e-05, 0.00737262, 0.9926136899999999]
        
        self.assertEqual((result1, result2), (expected_masses, expected_probabilities))

    def test_odd_count_sulfur(self):
        list_atoms = ['S', 'S', 'S']
        imprecision_True_False = False
        
        result1, result2 = MC.main_function(list_atoms, imprecision_True_False)
        
        expected_masses = [107.901, 105.902, 104.906, 103.906, 103.903, 102.906, 101.907, 101.91, 100.911, 99.911, 101.904, 100.907, 99.908, 98.911, 97.912, 98.914, 97.915, 96.916, 95.916, 97.908]
        expected_probabilities = [8e-12, 5.052e-09, 9e-10, 1.1402399999999999e-07, 1.063446e-06, 3.789e-07, 4.800410400000001e-05, 3.375e-08, 8.5518e-06, 0.0005488323989999999, 7.4618461e-05, 3.9879225e-05, 0.005052431946, 0.0018001539, 0.11403374905199998, 4.21875e-07, 0.00016034625, 0.020314800899999996, 0.8579166140079998, 0.020589998736191994]
        
        self.assertEqual((result1, result2), (expected_masses, expected_probabilities))

    def test_imprecision_filter(self):

        list_atoms = ['C', 'C', 'H']
        imprecision_True_False = True

        result1, result2 = MC.main_function(list_atoms, imprecision_True_False)

        expected_masses = [27.015, 26.011, 26.014, 25.008]
        expected_probabilities = [0.00012098185000000003, 0.021754736300000004, 0.00014671815000000002, 0.9779742818500002]

        self.assertEqual((result1, result2), (expected_masses, expected_probabilities))

class TestMolWebShow(unittest.TestCase):

    def test_default_parameters(self):

        mol_smi = "CCO"
         
        result = MC.mol_web_show(mol_smi)

        self.assertIsInstance(result, Div)

    def test_small_smiles(self):
        mol_smi = "CCO"  

        result = MC.mol_web_show(mol_smi)

        self.assertIsInstance(result, Div)

        mol = Chem.MolFromSmiles(mol_smi)
        
        image = Draw.MolToImage(mol)
        buffered = BytesIO()
        image.save(buffered, format="PNG")
        image_base64_with_hs = base64.b64encode(buffered.getvalue()).decode("utf-8")

        self.assertIn(image_base64_with_hs[:100], result.text)  

    def test_invalid_smiles(self):

        mol_smi = "InvalidSMILES"

        with self.assertRaises(ValueError):
            MC.mol_web_show(mol_smi)

    def test_empty_smiles(self):

        mol_smi = ""

        with self.assertRaises(ValueError):
            MC.mol_web_show(mol_smi)

    def test_complex_smiles(self):

        mol_smi = "C1=CC=C(C=C1)C2=CC=CC=C2" 

        result = MC.mol_web_show(mol_smi)

        self.assertIsInstance(result, Div)
    
    def test_benzene(self):

        mol_smi = "C1=CC=CC=C1" 

        result = MC.mol_web_show(mol_smi)

        self.assertIsInstance(result, Div)

        mol = Chem.MolFromSmiles(mol_smi)

        image = Draw.MolToImage(mol)
        buffered = BytesIO()
        image.save(buffered, format="PNG")
        image_base64_with_hs = base64.b64encode(buffered.getvalue()).decode("utf-8")

        self.assertIn(image_base64_with_hs[:100], result.text) 

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

class TestPeakMerger(unittest.TestCase):

    def test_basic_merge(self):

        x_in = [1.0, 1.1, 2.0]
        y_in = [0.5, 0.5, 1.0]

        result = MC.peak_merger(x_in, y_in, 0.15)

        expected_x = [1.05, 2.0]
        expected_y = [1.0, 1.0]

        self.assertEqual(result, (expected_x, expected_y))

    def test_no_merge(self):

        x_in = [1.0, 2.0, 3.0]
        y_in = [0.5, 0.5, 1.0]

        result = MC.peak_merger(x_in, y_in, 0.5)

        expected_x = [1.0, 2.0, 3.0]
        expected_y = [0.5, 0.5, 1.0]

        self.assertEqual(result, (expected_x, expected_y))

    def test_all_merge(self):

        x_in = [1.0, 1.1, 1.2]
        y_in = [0.5, 0.3, 0.2]

        result = MC.peak_merger(x_in, y_in, 0.15)

        expected_x = [1.05, 1.2]
        expected_y = [0.8, 0.2]

        self.assertEqual(result, (expected_x, expected_y))

    def test_no_peaks(self):

        with self.assertRaises(ValueError):
            MC.peak_merger([], [],1)
        
        with self.assertRaises(ValueError):
            MC.peak_merger([1.0], [],1)
        
        with self.assertRaises(ValueError):
            MC.peak_merger([], [0.5],1)

    def test_close_but_unmerged(self):

        x_in = [1.0, 1.16, 2.0]
        y_in = [0.5, 0.5, 1.0]

        result = MC.peak_merger(x_in, y_in, 0.15)

        expected_x = [1.0, 1.16, 2.0]
        expected_y = [0.5, 0.5, 1.0]

        self.assertEqual(result, (expected_x, expected_y))

    def test_mixed_merging(self):

        x_in = [1.0, 1.1, 2.0, 2.05]
        y_in = [0.5, 0.5, 0.7, 0.3]

        result = MC.peak_merger(x_in, y_in, 0.1)

        expected_x = [1.0, 1.1, 2.025]
        expected_y = [0.5, 0.5, 1.0]

        self.assertEqual(result, (expected_x, expected_y))

    def test_edge_case_resolution(self):

        x_in = [1.0, 1.15, 2.0]
        y_in = [0.5, 0.5, 1.0]

        result = MC.peak_merger(x_in, y_in, 0.15)

        expected_x = [1.075, 2.0]
        expected_y = [1.0, 1.0]
        
        self.assertEqual(result, (expected_x, expected_y))

class TestPeakSorter(unittest.TestCase):

    def test_basic_unsorted(self):

        x_in = [3.0, 1.0, 2.0]
        y_in = [0.3, 0.1, 0.2]

        result = MC.peak_sorter(x_in, y_in)

        expected_x = [1.0, 2.0, 3.0]
        expected_y = [0.1, 0.2, 0.3]

        self.assertEqual(result, (expected_x, expected_y))

    def test_already_sorted(self):

        x_in = [1.0, 2.0, 3.0]
        y_in = [0.1, 0.2, 0.3]

        result = MC.peak_sorter(x_in, y_in)

        expected_x = [1.0, 2.0, 3.0]
        expected_y = [0.1, 0.2, 0.3]

        self.assertEqual(result, (expected_x, expected_y))

    def test_empty_lists(self):

        with self.assertRaises(ValueError):
            MC.peak_sorter("","") 
    
    def test_one_empty_list_1(self):

        with self.assertRaises(ValueError):
            MC.peak_sorter("",['y']) 
    
    def test_one_empty_list_2(self):

        with self.assertRaises(ValueError):
            MC.peak_sorter(['x'],"") 

    def test_single_element_lists(self):
        x_in = [2.0]
        y_in = [0.2]

        result = MC.peak_sorter(x_in, y_in)

        expected_x = [2.0]
        expected_y = [0.2]

        self.assertEqual(result, (expected_x, expected_y))

    def test_negative_masses(self):

        x_in = [-1.0, -2.0, 3.0]
        y_in = [0.1, 0.2, 0.3]

        result = MC.peak_sorter(x_in, y_in)

        expected_x = [-2.0, -1.0, 3.0]
        expected_y = [0.2, 0.1, 0.3]

        self.assertEqual(result, (expected_x, expected_y))

    def test_different_size_lists(self):

        x_in = [1.0, 2.0]
        y_in = [0.1]
        
        with self.assertRaises(ValueError) as context:
            MC.peak_sorter(x_in, y_in)
        self.assertEqual(str(context.exception), 'Lists should be of the same size')

class TestSMILEsInterpreter(unittest.TestCase):

    def test_valid_smiles(self):

        smiles = "CCO"
        mol = MC.SMILEs_interpreter(smiles)

        self.assertIsNotNone(mol)
        self.assertEqual(Chem.MolToSmiles(Chem.RemoveHs(mol)), smiles)
    
    def test_invalid_smiles(self):

        smiles = "C1CCxVCC"

        with self.assertRaises(ValueError):
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
    
class TestSmilesTo3DPlot(unittest.TestCase):
    def test_valid_smiles(self):
        
        mol_smi = "CCO" 

        result = MC.smiles_to_3D_plot(mol_smi)

        self.assertIsNotNone(result)

    def test_empty_smiles(self):
        
        with self.assertRaises(ValueError):
            MC.smiles_to_3D_plot("")  

    def test_invalid_smiles(self):
        
        with self.assertRaises(ValueError):
            MC.smiles_to_3D_plot("invalid_smiles")    

    def test_plot_data(self):
        
        mol_smi = "CCO"  

        result = MC.smiles_to_3D_plot(mol_smi)

        self.assertIn('data', result)  

    def test_plot_layout(self):
        
        mol_smi = "CCO"  

        result = MC.smiles_to_3D_plot(mol_smi)

        self.assertIn('layout', result)  

    def test_plot_layout_attributes(self):
        
        mol_smi = "CCO" 

        result = MC.smiles_to_3D_plot(mol_smi)
        layout = result['layout']
        
        self.assertIn('title', layout)  
        self.assertIn('scene', layout)

class TestSpectrum3D(unittest.TestCase):
    
    def setUp(self):
        self.apparatus_resolution = 0.01

    def test_butane_False(self):
        mol_smi = "CCCC"
        imprecision_True_False = False
        result = MC.spectrum_3D(mol_smi, imprecision_True_False, self.apparatus_resolution)
        self.assertIsNotNone(result)
    
    def test_butane_True(self):
        mol_smi = "CCCC"
        imprecision_True_False = True
        result = MC.spectrum_3D(mol_smi, imprecision_True_False, self.apparatus_resolution)
        self.assertIsNotNone(result)

    def test_invalid_input(self):
        
        with self.assertRaises(ValueError):
            imprecision_True_False = False
            MC.spectrum_3D("invalid_smiles",imprecision_True_False, self.apparatus_resolution)

class TestSpectrum(unittest.TestCase):
    
    def setUp(self):
        self.apparatus_resolution = 0.01

    def test_butane_False(self):
        mol_smi = "CCCC"
        imprecision_True_False = False
        result = MC.spectrum(mol_smi, imprecision_True_False, self.apparatus_resolution)
        self.assertIsNotNone(result)
    
    def test_butane_True(self):
        mol_smi = "CCCC"
        imprecision_True_False = True
        result = MC.spectrum(mol_smi, imprecision_True_False, self.apparatus_resolution)
        self.assertIsNotNone(result)

    def test_empty_input(self):
        
        with self.assertRaises(ValueError):
            imprecision_True_False= False
            MC.spectrum_3D("",imprecision_True_False, self.apparatus_resolution)

    def test_invalid_input(self):
        
        with self.assertRaises(ValueError):
            imprecision_True_False = False
            MC.spectrum("invalid_smiles",imprecision_True_False, self.apparatus_resolution)  

class TestSulphurNitrogenAdder(unittest.TestCase):

    def test_no_nitrogen_or_sulphur(self):

        x_in = [1.0, 2.0]
        y_in = [0.1, 0.2]
        
        result = MC.sulphur_nitrogen_adder(x_in, y_in, False, False, 0, 0)

        expected_x = [1.0, 2.0]
        expected_y = [0.1, 0.2]

        self.assertEqual(result, (expected_x, expected_y))

    def test_only_nitrogen(self):

        x_in = [1.0, 2.0]
        y_in = [0.1, 0.2]

        result = MC.sulphur_nitrogen_adder(x_in, y_in, True, False, 1, 0)

        expected_x = [1.0, 2.0, 1.994]
        expected_y = [0.1, 0.2, 0.0007000000000000001]

        self.assertEqual(result, (expected_x, expected_y))

    def test_only_sulphur(self):

        x_in = [1.0, 2.0]
        y_in = [0.1, 0.2]

        result = MC.sulphur_nitrogen_adder(x_in, y_in, False, True, 0, 1)

        expected_x = [1.0, 2.0, 1.996]
        expected_y = [0.1, 0.2, 0.0016]

        self.assertEqual(result, (expected_x, expected_y))

    def test_both_nitrogen_and_sulphur(self):

        x_in = [1.0, 2.0]
        y_in = [0.1, 0.2]
        result = MC.sulphur_nitrogen_adder(x_in, y_in, True, True, 1, 1)

        expected_x = [1.0, 2.0, 1.994, 1.996]
        expected_y = [0.1, 0.2, 0.0007000000000000001, 0.0016]

        self.assertEqual(result, (expected_x, expected_y))

    def test_multiple_nitrogen_and_sulphur(self):

        x_in = [1.0, 2.0]
        y_in = [0.1, 0.2]

        result = MC.sulphur_nitrogen_adder(x_in, y_in, True, True, 2, 2)

        expected_x = [1.0, 2.0, 1.994, 1.996]
        expected_y = [0.1, 0.2, 0.0014000000000000002, 0.0032]

        self.assertEqual(result, (expected_x, expected_y))

    def test_no_atoms(self):
        with self.assertRaises(ValueError):
            MC.sulphur_nitrogen_adder([],[], True, True, 2, 2)
    
    def test_not_a_list(self):
        with self.assertRaises(ValueError):
            MC.sulphur_nitrogen_adder('NNS', 'NNS', True, True, 2, 2)

class TestSulphurNitrogenFinder(unittest.TestCase):
    def test_no_atoms(self):
        with self.assertRaises(ValueError):
            MC.sulphur_nitrogen_finder([])
    
    def test_not_a_list(self):
        with self.assertRaises(ValueError):
            MC.sulphur_nitrogen_finder("NSS")
    
    def test_odd_number_of_nitrogen(self):
        self.assertEqual(MC.sulphur_nitrogen_finder(['N', 'N', 'N']), (True, False, 3, 0))
    
    def test_even_number_of_nitrogen(self):
        self.assertEqual(MC.sulphur_nitrogen_finder(['N', 'N']), (False, False, 0, 0))
    
    def test_odd_number_of_sulphur(self):
        self.assertEqual(MC.sulphur_nitrogen_finder(['S', 'S', 'S']), (False, True, 0, 3))
    
    def test_even_number_of_sulphur(self):
        self.assertEqual(MC.sulphur_nitrogen_finder(['S', 'S']), (False, False, 0, 0))
    
    def test_mixed_atoms(self):
        self.assertEqual(MC.sulphur_nitrogen_finder(['N', 'O', 'S', 'N', 'S', 'S', 'C']), (False, True, 0, 3))

if __name__ == '__main__':
    unittest.main()



