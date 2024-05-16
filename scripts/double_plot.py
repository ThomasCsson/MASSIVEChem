from bokeh.plotting import figure
from bokeh.models.tickers import FixedTicker
from bokeh.layouts import column
from bokeh.models import WheelPanTool, WheelZoomTool, BoxAnnotation, CustomJS



def double_plot(x_in,y_in):

    #---------------------------------------------------------------------------------------------#
    '''
    double_plot(x_in,y_in)
    
    Input: list of masses (x_in) and intensities (y_in)
    
    Output: 2 Bokeh graphs:
            - One that shows the mass spectrum of the molecule to which the user can interact
            -Another that is the same graph but shows where the user is zooming on the first graph
    '''
    #---------------------------------------------------------------------------------------------#

    # tells where to put the graduation on the graph

    ticked_peaks = []
    for i in range(len(x_in)):
        if y_in[i] > 0.0001:
            ticked_peaks.append(x_in[i])

    #creates the principal graph, mass spectrum of the molecule (interactive)

    p1 = figure(width=700, title=f'Mass spectrum of molecule')
    p1 = figure(x_axis_label = '[m/z]')
    p1 = figure(y_axis_label = 'Abundance')
    p1.xaxis.axis_label = "[m/z]"
    p1.height = 500
    p1.xaxis.ticker = FixedTicker(ticks=ticked_peaks)
    p1.add_tools(WheelPanTool(dimension="height"))
    p1.add_tools(WheelZoomTool(dimensions="height"))
    p1.line(x_in, y_in, line_width=1)
    p1.xaxis.major_label_orientation = "horizontal"

     #creates the secondary graph, mass spectrum of the molecule (non-interactive)

    p2 = figure(title="Simulated Mass Spectrum", x_axis_label='Mass [Th]', y_axis_label='Intensity')
    p2 = figure(width=300, title=f'Mass spectrum of molecule')
    p2 = figure(toolbar_location=None)
    p2.height = 300
    p2.line(x_in, y_in, legend_label="Mass spectrum", line_width=1)

    #creates a tool in order that the second graph shows where the zoom is on the first one

    box = BoxAnnotation(left=0, right=0, bottom=0, top=0,
    fill_alpha=0.1, line_color='red', fill_color='cornflowerblue')

    jscode = """
        box[%r] = cb_obj.start
        box[%r] = cb_obj.end
    """

    xcb = CustomJS(args=dict(box=box), code=jscode % ('left', 'right'))
    ycb = CustomJS(args=dict(box=box), code=jscode % ('bottom', 'top'))

    p1.x_range.js_on_change('start', xcb)
    p1.x_range.js_on_change('end', xcb)
    p1.y_range.js_on_change('start', ycb)
    p1.y_range.js_on_change('end', ycb)

    # adds the functionnality to the second figure

    p2.add_layout(box)

    # creates a layout that displays the 2 graphs

    layout = column(p1, p2)
    print('here')
    return layout


import unittest

class TestDoublePlot(unittest.TestCase):
    def test_double_plot(self):
        # Test the function with sample data
        x_in = [100, 200, 300, 400, 500]
        y_in = [0.1, 0.2, 0.3, 0.4, 0.5]
        layout = double_plot(x_in, y_in)
        self.assertIsNotNone(layout)
        self.assertEqual(len(layout.children), 2)  # Check if there are two plots in the layout

    def test_double_plot_empty_data(self):
        # Test the function with empty data
        x_in = []
        y_in = []
        layout = double_plot(x_in, y_in)
        self.assertIsNotNone(layout)
        self.assertEqual(len(layout.children), 2)  # Check if there are two plots in the layout

    def test_double_plot_single_point(self):
        # Test the function with a single point
        x_in = [100]
        y_in = [0.1]
        layout = double_plot(x_in, y_in)
        self.assertIsNotNone(layout)
        self.assertEqual(len(layout.children), 2)  # Check if there are two plots in the layout

    def test_double_plot_identical_data(self):
        # Test the function with identical x and y data
        x_in = [100, 200, 300, 400, 500]
        y_in = [0.1, 0.1, 0.1, 0.1, 0.1]
        layout = double_plot(x_in, y_in)
        self.assertIsNotNone(layout)
        self.assertEqual(len(layout.children), 2)  # Check if there are two plots in the layout

if __name__ == '__main__':
    unittest.main()


