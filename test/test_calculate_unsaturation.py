import MASSiveChem.MASSiveChem as MC
import unittest

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

if __name__ == '__main__':
    unittest.main()
