def delta_function_plotter(x_in, y_in) -> list[float]:
    #---------------------------------------------------------------------------------------------#
    '''
    delta_function_plotter(x_in, y_in)
    
    Input: two lists:
    1. list of the masses (of individual molecules) of each possible combination of isotopes (ordered)
    2. list of the probabilities of apparation of each of the molecules (ordered)
    
    Output: two lists:
    1. list of the masses (of individual molecules) of each possible combination of isotopes (ordered) with values of 0 (y_axis) added on eiter side of the "peak"
    2. list of the probabilities of apparation of each of the molecules (ordered)
    
    (the mass in list 1 at index i is associated to the probability at index i in list 2)
    '''
    #---------------------------------------------------------------------------------------------#

    if not x_in:
        raise ValueError('Enter a valid input')
    if not y_in:
        raise ValueError('Enter a valid input')
    if len(x_in) != len(y_in):
        raise ValueError('Both entry should be of the same length')

    min_x , max_x = min(x_in), max(x_in)

    x_axis, y_axis = [min_x-0.5],[0]
    for i in range (len(x_in)):
        x_axis.append(x_in[i]-10**(-100))
        x_axis.append(x_in[i])
        x_axis.append(x_in[i]+10**(-100))
        y_axis.append(0)
        y_axis.append(y_in[i])
        y_axis.append(0)

    x_axis.append(max_x+1)
    y_axis.append(0)

    return x_axis, y_axis

import unittest

class TestDeltaFunctionPlotter(unittest.TestCase):

    def test_single_peak(self):
        
        x_in = [1.0]
        y_in = [0.5]

        expected_x = [0.5, 1.0 - 10**(-100), 1.0, 1.0 + 10**(-100), 2.0]
        expected_y = [0, 0, 0.5, 0, 0]

        x_axis, y_axis = delta_function_plotter(x_in, y_in)

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

        x_axis, y_axis = delta_function_plotter(x_in, y_in)

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

        x_axis, y_axis = delta_function_plotter(x_in, y_in)

        self.assertEqual(x_axis, expected_x)
        self.assertEqual(y_axis, expected_y)

    def test_empty_input_raises_value_error(self):

        with self.assertRaises(ValueError):

            delta_function_plotter([], [])
        
        with self.assertRaises(ValueError):

            delta_function_plotter([1.0], [])
        
        with self.assertRaises(ValueError):

            delta_function_plotter([], [0.5])

    def test_invalid_input_length_mismatch(self):

        with self.assertRaises(ValueError):
            delta_function_plotter([1.0, 2.0], [0.5])

if __name__ == '__main__':
    unittest.main()
