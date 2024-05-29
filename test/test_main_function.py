import MASSiveChem.MASSiveChem as MC
import unittest

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

if __name__ == '__main__':
    unittest.main()
