def ionisation_method (list_atoms) -> list[str]:
    #---------------------------------------------------------------------------------------------#
    '''
    ionisation_method (list_atoms)
    
    Input: list of atomic symbols of atoms in a given molecule
    
    Output: corrected list of atoms (i.e. list of atoms that enter spectrometry apparatus)
    '''
    #---------------------------------------------------------------------------------------------#

    '''In the case of ionisation by proton, we need to add a H+ ion, which is done in the following'''
    if 'H' in list_atoms:

        #Check that there is in fact a proton to remove
        list_atoms.remove('H')
    return list_atoms


import unittest

class TestIonisationMethod(unittest.TestCase):

    def test_remove_single_h(self):
        self.assertEqual(ionisation_method(['H', 'C', 'O']), ['C', 'O'])

    def test_remove_multiple_h(self):
        self.assertEqual(ionisation_method(['H', 'H', 'C', 'O']), ['H', 'C', 'O'])

    def test_no_h_in_list(self):
        self.assertEqual(ionisation_method(['C', 'O', 'N']), ['C', 'O', 'N'])

    def test_only_h_in_list(self):
        self.assertEqual(ionisation_method(['H']), [])

    def test_empty_list(self):
        self.assertEqual(ionisation_method([]), [])

    def test_no_removal_for_other_atoms(self):
        self.assertEqual(ionisation_method(['C', 'O', 'N']), ['C', 'O', 'N'])

    def test_remove_h_from_complex_list(self):
        self.assertEqual(ionisation_method(['H', 'H', 'C', 'O', 'H']), ['H', 'C', 'O', 'H'])

if __name__ == '__main__':
    unittest.main()
