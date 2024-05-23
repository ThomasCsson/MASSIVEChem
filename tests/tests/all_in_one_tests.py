import unittest
from bokeh.plotting import figure
from bokeh.models import ColumnDataSource
from bokeh.layouts import column, row
from bokeh.models.widgets import DataTable, TableColumn

def all_in_one(p1, p2, p3, p4):

    #---------------------------------------------------------------------------------------------#
    '''
    all_in_one(p1,p2,p3,p4)
    
    Input: 4 bokeh plots
            Usually used in this package:
                    - p1 : bokeh double plot of mass spectrometry
                    - p2 : image of the molecule
                    - p3 : table of functional groups
                    - p4 : molecule in 3D interactive
    
    Output: bokeh page with all 4 graphs well arranged

    '''
    #---------------------------------------------------------------------------------------------#
    
    #creates a layout in row with p3 and p4
    layout1 = row(p3, p4)

    #creates a layout in column with p2 and layout1
    layout2 = column(p2, layout1)

    #creates the final layout in row with layout1 and p1
    layout = row(p1, layout2)

    return layout

import unittest
from bokeh.plotting import figure
from bokeh.layouts import row, column, layout
from bokeh.models import LayoutDOM
from bokeh.layouts import column 

from bokeh.models import LayoutDOM

class TestAllInOneFunction(unittest.TestCase):
    
    def setUp(self):
        # Create dummy Bokeh plots
        self.p1 = figure(title="Plot 1")
        self.p2 = figure(title="Plot 2")
        self.p3 = figure(title="Plot 3")
        self.p4 = figure(title="Plot 4")
    
    def test_all_in_one_output_type(self):
        result = all_in_one(self.p1, self.p2, self.p3, self.p4)
        self.assertIsInstance(result, LayoutDOM)
    
    def test_all_in_one_structure(self):
        result = all_in_one(self.p1, self.p2, self.p3, self.p4)
        
        # Check the structure of the layout
        self.assertEqual(len(result.children), 2)
        self.assertIs(result.children[0], self.p1)
        
        col_layout = result.children[1]
        self.assertIsInstance(col_layout, column)  # Use column as a type
        self.assertEqual(len(col_layout.children), 2)
        
        self.assertIs(col_layout.children[0], self.p2)
        
        row_layout = col_layout.children[1]
        self.assertIsInstance(row_layout, row)
        self.assertEqual(len(row_layout.children), 2)
        
        self.assertIs(row_layout.children[0], self.p3)
        self.assertIs(row_layout.children[1], self.p4)

if __name__ == '__main__':
    unittest.main()


