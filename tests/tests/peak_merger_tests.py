import unittest

def peak_merger(x_in, y_in, apparatus_resolution) -> list[float]:
    #---------------------------------------------------------------------------------------------#
    '''
    delta_function_plotter(x_in, y_in)
    
    Input: two lists + float:
    1. ordered list of the masses (of individual molecules) of each possible combination of isotopes
    2. ordered list of the probabilities of apparation of each of the molecules
    3. apparatus precision (float representation of number which is the limit between two peaks above which they appear merged together)
    
    Output: two lists:
    1. ordered list of the masses with values on y merged together if peaks within precision of apparatus
    2. ordered list of the probabilities of apparation of each of the molecules
    
    (the mass in list 1 at index i is associated to the probability at index i in list 2)
    '''
    #---------------------------------------------------------------------------------------------#
    
    if not x_in:
        raise ValueError('Empty list')
    if not y_in:
        raise ValueError('Empty list')
    if len(x_in) != len(y_in):
        raise ValueError('Lists should be of the same size')
    
    x_out, y_out = [],[]
    while len(x_in)>1:
        if x_in[0]>x_in[1]-apparatus_resolution:
            y_in[1] = y_in[0] + y_in[1]
            x_in[1] = (x_in[0] + x_in[1])/2
            x_in.pop(0)
            y_in.pop(0)
        else:
            x_out.append(x_in[0])
            y_out.append(y_in[0])
            x_in.pop(0)
            y_in.pop(0)
    x_out.append(x_in[0])
    y_out.append(y_in[0])

    return x_out, y_out



class TestPeakMerger(unittest.TestCase):

    def test_basic_merge(self):

        x_in = [1.0, 1.1, 2.0]
        y_in = [0.5, 0.5, 1.0]

        result = peak_merger(x_in, y_in, 0.15)

        expected_x = [1.05, 2.0]
        expected_y = [1.0, 1.0]

        self.assertEqual(result, (expected_x, expected_y))

    def test_no_merge(self):

        x_in = [1.0, 2.0, 3.0]
        y_in = [0.5, 0.5, 1.0]

        result = peak_merger(x_in, y_in, 0.5)

        expected_x = [1.0, 2.0, 3.0]
        expected_y = [0.5, 0.5, 1.0]

        self.assertEqual(result, (expected_x, expected_y))

    def test_all_merge(self):

        x_in = [1.0, 1.1, 1.2]
        y_in = [0.5, 0.3, 0.2]

        result = peak_merger(x_in, y_in, 0.15)

        expected_x = [1.05, 1.2]
        expected_y = [0.8, 0.2]

        self.assertEqual(result, (expected_x, expected_y))

    def test_no_peaks(self):

        with self.assertRaises(ValueError):
            peak_merger([], [],1)
        
        with self.assertRaises(ValueError):
            peak_merger([1.0], [],1)
        
        with self.assertRaises(ValueError):
            peak_merger([], [0.5],1)

    def test_close_but_unmerged(self):

        x_in = [1.0, 1.16, 2.0]
        y_in = [0.5, 0.5, 1.0]

        result = peak_merger(x_in, y_in, 0.15)

        expected_x = [1.0, 1.16, 2.0]
        expected_y = [0.5, 0.5, 1.0]

        self.assertEqual(result, (expected_x, expected_y))

    def test_mixed_merging(self):

        x_in = [1.0, 1.1, 2.0, 2.05]
        y_in = [0.5, 0.5, 0.7, 0.3]

        result = peak_merger(x_in, y_in, 0.1)

        expected_x = [1.0, 1.1, 2.025]
        expected_y = [0.5, 0.5, 1.0]

        self.assertEqual(result, (expected_x, expected_y))

    def test_edge_case_resolution(self):

        x_in = [1.0, 1.15, 2.0]
        y_in = [0.5, 0.5, 1.0]

        result = peak_merger(x_in, y_in, 0.15)

        expected_x = [1.075, 2.0]
        expected_y = [1.0, 1.0]
        
        self.assertEqual(result, (expected_x, expected_y))

if __name__ == '__main__':
    unittest.main()
