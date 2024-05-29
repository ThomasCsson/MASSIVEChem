import MASSiveChem.MASSiveChem as MC
import unittest

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
