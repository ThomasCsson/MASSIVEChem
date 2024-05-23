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
    
    if not p1:
        raise ValueError("Enter a non-empty plot")
    if not p2:
        raise ValueError("Enter a non-empty plot")
    if not p3:
        raise ValueError("Enter a non-empty plot")
    if not p4:
        raise ValueError("Enter a non-empty plot")
    
    #creates a layout in row with p3 and p4
    layout1 = row(p3, p4)

    #creates a layout in column with p2 and layout1
    layout2 = column(p2, layout1)

    #creates the final layout in row with layout1 and p1
    layout = row(p1, layout2)

    return layout
