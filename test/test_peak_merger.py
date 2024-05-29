import MASSiveChem.MASSiveChem as MC
import unittest

class TestPeakMerger(unittest.TestCase):

    def test_basic_merge(self):

        x_in = [1.0, 1.1, 2.0]
        y_in = [0.5, 0.5, 1.0]

        result = MC.peak_merger(x_in, y_in, 0.15)

        expected_x = [1.05, 2.0]
        expected_y = [1.0, 1.0]

        self.assertEqual(result, (expected_x, expected_y))

    def test_no_merge(self):

        x_in = [1.0, 2.0, 3.0]
        y_in = [0.5, 0.5, 1.0]

        result = MC.peak_merger(x_in, y_in, 0.5)

        expected_x = [1.0, 2.0, 3.0]
        expected_y = [0.5, 0.5, 1.0]

        self.assertEqual(result, (expected_x, expected_y))

    def test_all_merge(self):

        x_in = [1.0, 1.1, 1.2]
        y_in = [0.5, 0.3, 0.2]

        result = MC.peak_merger(x_in, y_in, 0.15)

        expected_x = [1.05, 1.2]
        expected_y = [0.8, 0.2]

        self.assertEqual(result, (expected_x, expected_y))

    def test_no_peaks(self):

        with self.assertRaises(ValueError):
            MC.peak_merger([], [],1)
        
        with self.assertRaises(ValueError):
            MC.peak_merger([1.0], [],1)
        
        with self.assertRaises(ValueError):
            MC.peak_merger([], [0.5],1)

    def test_close_but_unmerged(self):

        x_in = [1.0, 1.16, 2.0]
        y_in = [0.5, 0.5, 1.0]

        result = MC.peak_merger(x_in, y_in, 0.15)

        expected_x = [1.0, 1.16, 2.0]
        expected_y = [0.5, 0.5, 1.0]

        self.assertEqual(result, (expected_x, expected_y))

    def test_mixed_merging(self):

        x_in = [1.0, 1.1, 2.0, 2.05]
        y_in = [0.5, 0.5, 0.7, 0.3]

        result = MC.peak_merger(x_in, y_in, 0.1)

        expected_x = [1.0, 1.1, 2.025]
        expected_y = [0.5, 0.5, 1.0]

        self.assertEqual(result, (expected_x, expected_y))

    def test_edge_case_resolution(self):

        x_in = [1.0, 1.15, 2.0]
        y_in = [0.5, 0.5, 1.0]

        result = MC.peak_merger(x_in, y_in, 0.15)

        expected_x = [1.075, 2.0]
        expected_y = [1.0, 1.0]
        
        self.assertEqual(result, (expected_x, expected_y))

if __name__ == '__main__':
    unittest.main()
