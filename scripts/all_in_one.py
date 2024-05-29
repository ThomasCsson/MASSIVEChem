from bokeh.layouts import column, row

def all_in_one(p1, p2, p3):

    #---------------------------------------------------------------------------------------------#
    '''
    all_in_one(p1,p2,p3)
    
    Input: 3 bokeh plots
            Usually used in this package:
                    - p1 : bokeh double plot of mass spectrometry
                    - p2 : image of the molecule
                    - p3 : table of functional groups
    
    Output: bokeh page with all 4 graphs well arranged

    '''
    #---------------------------------------------------------------------------------------------#
    
    if not p1:
        raise ValueError("Enter a non-empty plot")
    if not p2:
        raise ValueError("Enter a non-empty plot")
    if not p3:
        raise ValueError("Enter a non-empty plot")

    #creates a layout in column with p2 and p3
    layout2 = column(p2, p3)

    #creates the final layout in row with layout1 and p1
    layout = row(p1, layout2)

    return layout
