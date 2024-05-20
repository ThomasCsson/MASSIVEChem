import unittest
from bokeh.plotting import figure
from bokeh.models import ColumnDataSource
from bokeh.layouts import column, row
from bokeh.models.widgets import DataTable, TableColumn

def all_in_one(p1, p2, p3, p4):

    #---------------------------------------------------------------------------------------------#
    '''
    all_in_one(p1,p2,p3, p4)
    
    Input: 3 bokeh plots
            Usually used in this package:
                    - p1 : bokeh double plot of mass spectrometry
                    - p2 : image of the molecule
                    - p3 : table of functional groups
                    - p4 : buttons with info on the molecule
    
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


class TestAllInOneFunction(unittest.TestCase):
    def setUp(self):
        # Create sample bokeh plots
        self.p1 = figure(title="Mass Spectrometry")
        self.p2 = figure(title="Molecule Image")
        self.p3 = figure(title="Functional Groups Table")

        # Add some dummy data to the figures
        self.p1.line([1, 2, 3], [4, 5, 6])
        self.p2.image_url(url=["https://bokeh.org/static/img/logo.png"], x=0, y=1, w=1, h=1)
        
        # Create a simple table for p3
        data = dict(groups=["Alcohol", "Ketone"], images=["<img src='data:image/png;base64,some_base64_string'>", "<img src='data:image/png;base64,some_base64_string'>"])
        source = ColumnDataSource(data)
        columns = [
            TableColumn(field="groups", title="Functional Groups"),
            TableColumn(field="images", title="Images")
        ]
        self.p3 = DataTable(source=source, columns=columns, width=400, height=280)

    def test_output_type(self):
        layout = all_in_one(self.p1, self.p2, self.p3)
        self.assertIsInstance(layout, row)

    def test_layout_structure(self):
        layout = all_in_one(self.p1, self.p2, self.p3)
        self.assertEqual(len(layout.children), 2)
        self.assertIsInstance(layout.children[1], column)
        self.assertEqual(len(layout.children[1].children), 2)

    def test_correct_plots_in_layout(self):
        layout = all_in_one(self.p1, self.p2, self.p3)
        self.assertEqual(layout.children[0], self.p1)
        self.assertEqual(layout.children[1].children[0], self.p2)
        self.assertEqual(layout.children[1].children[1], self.p3)

if __name__ == "__main__":
    unittest.main()