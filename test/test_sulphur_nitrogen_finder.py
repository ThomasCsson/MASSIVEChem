def sulphur_nitrogen_finder(list_atoms):
    #---------------------------------------------------------------------------------------------#
    '''
    sulpher_nitrogen_finder(list_atoms)

    Input: list of atoms in the molecule

    Output: 2 booleans and 2 integers:
    1. boolean that tells if the molecule contains an odd number of Nitrogen atoms
    2. boolean that tells if the molecule contains an odd number of Sulphur atoms
    3. number of Nitrogen atoms in the molecule (in case where there is an odd number)
    4. number of Sulphur atoms in the molecule (in case where there is an odd number)
    '''
    #---------------------------------------------------------------------------------------------#
    if not list_atoms:
        raise ValueError("Enter a non-empty list")
    if type(list_atoms) != list:
        raise ValueError("Enter a list as argument")
    

    count_N = 0
    has_N = False
    if 'N' in list_atoms and list_atoms.count('N')%2 == 1:
        has_N = True
        count_N = list_atoms.count('N')
    count_S = 0
    has_S = False
    if 'S' in list_atoms and list_atoms.count('S')%2 == 1:
        has_S = True
        count_S = list_atoms.count('S')
    return has_N, has_S, count_N, count_S

import unittest

class TestSulphurNitrogenFinder(unittest.TestCase):
    def test_no_atoms(self):
        with self.assertRaises(ValueError):
            sulphur_nitrogen_finder([])
    
    def test_not_a_list(self):
        with self.assertRaises(ValueError):
            sulphur_nitrogen_finder("NSS")
    
    def test_odd_number_of_nitrogen(self):
        self.assertEqual(sulphur_nitrogen_finder(['N', 'N', 'N']), (True, False, 3, 0))
    
    def test_even_number_of_nitrogen(self):
        self.assertEqual(sulphur_nitrogen_finder(['N', 'N']), (False, False, 0, 0))
    
    def test_odd_number_of_sulphur(self):
        self.assertEqual(sulphur_nitrogen_finder(['S', 'S', 'S']), (False, True, 0, 3))
    
    def test_even_number_of_sulphur(self):
        self.assertEqual(sulphur_nitrogen_finder(['S', 'S']), (False, False, 0, 0))
    
    def test_mixed_atoms(self):
        self.assertEqual(sulphur_nitrogen_finder(['N', 'O', 'S', 'N', 'S', 'S', 'C']), (False, True, 0, 3))

if __name__ == '__main__':
    unittest.main()
