import MASSiveChem.MASSiveChem as MC
import unittest

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

if __name__ == '__main__':
    unittest.main()
