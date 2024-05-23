def peak_sorter(x_in, y_in) -> list[float]:
    #---------------------------------------------------------------------------------------------#
    '''
    peak_sorter(x_in, y_in)
    
    Input: two lists:
    1. ordered list of the masses with values on y merged together if peaks within precision of apparatus
    2. ordered list of the probabilities of apparation of each of the molecules
    
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
    while len(x_in)>0:
        min_x = min(x_in)
        index_min = x_in.index(min_x)
        x_out.append(min_x)
        y_out.append(y_in[index_min])
        x_in.pop(index_min)
        y_in.pop(index_min)

    return x_out, y_out

import unittest

class TestPeakSorter(unittest.TestCase):

    def test_basic_unsorted(self):
        x_in = [3.0, 1.0, 2.0]
        y_in = [0.3, 0.1, 0.2]
        result = peak_sorter(x_in, y_in)
        expected_x = [1.0, 2.0, 3.0]
        expected_y = [0.1, 0.2, 0.3]
        self.assertEqual(result, (expected_x, expected_y))

    def test_already_sorted(self):
        x_in = [1.0, 2.0, 3.0]
        y_in = [0.1, 0.2, 0.3]
        result = peak_sorter(x_in, y_in)
        expected_x = [1.0, 2.0, 3.0]
        expected_y = [0.1, 0.2, 0.3]
        self.assertEqual(result, (expected_x, expected_y))

    def test_empty_lists(self):
        x_in = []
        y_in = []
        with self.assertRaises(ValueError) as context:
            peak_sorter(x_in, y_in)
        self.assertEqual(str(context.exception), 'Empty list')

    def test_single_element_lists(self):
        x_in = [2.0]
        y_in = [0.2]
        result = peak_sorter(x_in, y_in)
        expected_x = [2.0]
        expected_y = [0.2]
        self.assertEqual(result, (expected_x, expected_y))

    def test_negative_masses(self):
        x_in = [-1.0, -2.0, 3.0]
        y_in = [0.1, 0.2, 0.3]
        result = peak_sorter(x_in, y_in)
        expected_x = [-2.0, -1.0, 3.0]
        expected_y = [0.2, 0.1, 0.3]
        self.assertEqual(result, (expected_x, expected_y))

    def test_different_size_lists(self):
        x_in = [1.0, 2.0]
        y_in = [0.1]
        with self.assertRaises(ValueError) as context:
            peak_sorter(x_in, y_in)
        self.assertEqual(str(context.exception), 'Lists should be of the same size')

if __name__ == '__main__':
    unittest.main()