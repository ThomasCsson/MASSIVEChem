import MASSiveChem.MASSiveChem as MC
import unittest

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

if __name__ == '__main__':
    unittest.main()
