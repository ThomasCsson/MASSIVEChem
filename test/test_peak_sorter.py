import MASSiveChem.MASSiveChem as MC
import unittest

class TestPeakSorter(unittest.TestCase):

    def test_basic_unsorted(self):

        x_in = [3.0, 1.0, 2.0]
        y_in = [0.3, 0.1, 0.2]

        result = MC.peak_sorter(x_in, y_in)

        expected_x = [1.0, 2.0, 3.0]
        expected_y = [0.1, 0.2, 0.3]

        self.assertEqual(result, (expected_x, expected_y))

    def test_already_sorted(self):

        x_in = [1.0, 2.0, 3.0]
        y_in = [0.1, 0.2, 0.3]

        result = MC.peak_sorter(x_in, y_in)

        expected_x = [1.0, 2.0, 3.0]
        expected_y = [0.1, 0.2, 0.3]

        self.assertEqual(result, (expected_x, expected_y))

    def test_empty_lists(self):

        with self.assertRaises(ValueError):
            MC.peak_sorter("","") 
    
    def test_one_empty_list_1(self):

        with self.assertRaises(ValueError):
            MC.peak_sorter("",['y']) 
    
    def test_one_empty_list_2(self):

        with self.assertRaises(ValueError):
            MC.peak_sorter(['x'],"") 

    def test_single_element_lists(self):
        x_in = [2.0]
        y_in = [0.2]

        result = MC.peak_sorter(x_in, y_in)

        expected_x = [2.0]
        expected_y = [0.2]

        self.assertEqual(result, (expected_x, expected_y))

    def test_negative_masses(self):

        x_in = [-1.0, -2.0, 3.0]
        y_in = [0.1, 0.2, 0.3]

        result = MC.peak_sorter(x_in, y_in)

        expected_x = [-2.0, -1.0, 3.0]
        expected_y = [0.2, 0.1, 0.3]

        self.assertEqual(result, (expected_x, expected_y))

    def test_different_size_lists(self):

        x_in = [1.0, 2.0]
        y_in = [0.1]
        
        with self.assertRaises(ValueError) as context:
            MC.peak_sorter(x_in, y_in)
        self.assertEqual(str(context.exception), 'Lists should be of the same size')

if __name__ == '__main__':
    unittest.main()