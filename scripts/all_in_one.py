from bokeh.layouts import row, column

def all_in_one(p1,p2,p3):

    #---------------------------------------------------------------------------------------------#
    '''
    all_in_one(p1,p2,p3)
    
    Input: 3 bokeh plots
            Usually used in this package:
                    - p1 : bokeh double plot of mass spectrometry
                    - p2 : image of the molecule
                    - p3 : table of functional groups
    
    Output: bokeh page with all 3 graphs well arranged

    '''
    #---------------------------------------------------------------------------------------------#

    #creates a layout in column with p2 and p3
    layout1 = column(p2, p3)

    #creates the final layout in row with layout1 and p1
    layout = row(p1, layout1)

    return layout

import unittest
from bokeh.plotting import figure
from bokeh.layouts import column

class TestAllInOne(unittest.TestCase):
    def test_all_in_one(self):
        # Create sample plots for testing with some content added to them
        p1 = figure()
        p1.circle([1, 2, 3], [4, 5, 6], size=10, color="navy", alpha=0.5)

        p2 = figure()
        p2.line([1, 2, 3], [4, 5, 6], line_width=2)

        p3 = figure()
        p3.square([1, 2, 3], [4, 5, 6], size=10, color="olive", alpha=0.5)

        # Call the function
        layout = all_in_one(p1, p2, p3)

        # Assertions
        self.assertIsInstance(layout, column)  # Assuming Bokeh's column layout is used
        self.assertEqual(len(layout.children), 2)  # Assuming there are two plots in the column layout
        self.assertEqual(len(layout.children[1].children), 2)  # Assuming there are two plots in the second column layout

if __name__ == '__main__':
    unittest.main()
