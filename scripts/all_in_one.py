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

