import MASSiveChem.MASSiveChem as MC
import unittest

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
            MC.sulphur_nitrogen_adder("NSS", "NNS", True, True, 2, 2)

if __name__ == '__main__':
    unittest.main()
