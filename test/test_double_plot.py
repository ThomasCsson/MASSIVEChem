import MASSiveChem.MASSiveChem as MC
import unittest

class TestDoublePlot(unittest.TestCase):
    def test_double_plot(self):

        x_in = [100, 200, 300, 400, 500]
        y_in = [0.1, 0.2, 0.3, 0.4, 0.5]

        layout = MC.double_plot(x_in, y_in)

        self.assertIsNotNone(layout)
        self.assertEqual(len(layout.children), 2)  # Check if there are two plots in the layout

    def test_double_plot_empty_data(self):

        x_in = []
        y_in = []

        with self.assertRaises(ValueError):
            MC.double_plot(x_in, y_in)
    
    def test_double_plot_1_empty_data(self):

        x_in = [1]
        y_in = []

        with self.assertRaises(ValueError):
            MC.double_plot(x_in, y_in)
    
    def test_double_plot_1_empty_data_2(self):

        x_in = []
        y_in = [1]

        with self.assertRaises(ValueError):
            MC.double_plot(x_in, y_in)

    def test_double_plot_single_point(self):
        
        x_in = [100]
        y_in = [0.1]

        layout = MC.double_plot(x_in, y_in)

        self.assertIsNotNone(layout)
        self.assertEqual(len(layout.children), 2)  # Check if there are two plots in the layout

    def test_double_plot_identical_data(self):
        
        x_in = [100, 200, 300, 400, 500]
        y_in = [0.1, 0.1, 0.1, 0.1, 0.1]

        layout = MC.double_plot(x_in, y_in)

        self.assertIsNotNone(layout)
        self.assertEqual(len(layout.children), 2)  # Check if there are two plots in the layout

if __name__ == '__main__':
    unittest.main()


