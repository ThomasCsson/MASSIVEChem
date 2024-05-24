def sulphur_nitrogen_adder(x_in, y_in, has_N, has_S, count_N, count_S) -> list[float]:
    #---------------------------------------------------------------------------------------------#
    '''
    sulphur_nitrogen_adder(x_in, y_in, has_N, has_S, count_N, count_S)
    
    Input: two lists + 4 booleans:
    1. ordered list of the masses with values on y merged together if peaks within precision of apparatus
    2. ordered list of the probabilities of apparation of each of the molecules
    3. boolean that tells if the molecule contains Nitrogen
    4. boolean that tells if the molecule contains Sulphur
    5. number of Nitrogen atoms in the molecule
    6. number of Sulphur atoms in the molecule
    
    Output: two lists:
    1. ordered list of the masses with values on y merged together if peaks within precision of apparatus
    2. ordered list of the probabilities of apparation of each of the molecules
    
    (the mass in list 1 at index i is associated to the probability at index i in list 2)
    '''
    #---------------------------------------------------------------------------------------------#

    if not x_in:
        raise ValueError("Enter a non-empty list")
    if not y_in:
        raise ValueError("Enter a non-empty list")
    if len(x_in) != len(y_in):
        raise ValueError("The two lists must have the same length")
    if has_N not in [True, False]:
        raise ValueError("The third argument must be a boolean")
    if has_S not in [True, False]:
        raise ValueError("The fourth argument must be a boolean")
    if count_N < 0:
        raise ValueError("The fifth argument must be a positive integer")
    if count_S < 0:
        raise ValueError("The sixth argument must be a positive integer")
    if type(count_N) != int:
        raise ValueError("The fifth argument must be an integer")
    if type(count_S) != int:
        raise ValueError("The sixth argument must be an integer")
    if x_in != sorted(x_in):
        raise ValueError("The first list must be ordered")

    maximum = max(y_in)

    if has_N:
        x_in.append(x_in[1] - 0.006)  
        y_in.append(0.0035*count_N*maximum)  
    
    if has_S:
        x_in.append(x_in[1]-0.004)  
        y_in.append(0.008*count_S*maximum)


    return x_in, y_in


import unittest

class TestSulphurNitrogenAdder(unittest.TestCase):

    def test_no_nitrogen_or_sulphur(self):

        x_in = [1.0, 2.0]
        y_in = [0.1, 0.2]
        
        result = sulphur_nitrogen_adder(x_in, y_in, False, False, 0, 0)

        expected_x = [1.0, 2.0]
        expected_y = [0.1, 0.2]

        self.assertEqual(result, (expected_x, expected_y))

    def test_only_nitrogen(self):

        x_in = [1.0, 2.0]
        y_in = [0.1, 0.2]

        result = sulphur_nitrogen_adder(x_in, y_in, True, False, 1, 0)

        expected_x = [1.0, 2.0, 1.994]
        expected_y = [0.1, 0.2, 0.0007000000000000001]

        self.assertEqual(result, (expected_x, expected_y))

    def test_only_sulphur(self):

        x_in = [1.0, 2.0]
        y_in = [0.1, 0.2]

        result = sulphur_nitrogen_adder(x_in, y_in, False, True, 0, 1)

        expected_x = [1.0, 2.0, 1.996]
        expected_y = [0.1, 0.2, 0.0016]

        self.assertEqual(result, (expected_x, expected_y))

    def test_both_nitrogen_and_sulphur(self):

        x_in = [1.0, 2.0]
        y_in = [0.1, 0.2]
        result = sulphur_nitrogen_adder(x_in, y_in, True, True, 1, 1)

        expected_x = [1.0, 2.0, 1.994, 1.996]
        expected_y = [0.1, 0.2, 0.0007000000000000001, 0.0016]

        self.assertEqual(result, (expected_x, expected_y))

    def test_multiple_nitrogen_and_sulphur(self):

        x_in = [1.0, 2.0]
        y_in = [0.1, 0.2]

        result = sulphur_nitrogen_adder(x_in, y_in, True, True, 2, 2)

        expected_x = [1.0, 2.0, 1.994, 1.996]
        expected_y = [0.1, 0.2, 0.0014000000000000002, 0.0032]

        self.assertEqual(result, (expected_x, expected_y))

    def test_no_atoms(self):
        with self.assertRaises(ValueError):
            sulphur_nitrogen_adder([],[], True, True, 2, 2)
    
    def test_not_a_list(self):
        with self.assertRaises(ValueError):
            sulphur_nitrogen_adder("NSS", "NNS", True, True, 2, 2)

if __name__ == '__main__':
    unittest.main()


#python -m unittest discover test