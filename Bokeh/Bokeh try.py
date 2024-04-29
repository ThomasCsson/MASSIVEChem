

from bokeh.models import CustomJS, ColumnDataSource, Slider
from bokeh.plotting import figure, output_file, show
from bokeh.layouts import column
from bokeh.models import WheelZoomTool, PanTool
import numpy as np

output_file("callback.html")

x = np.linspace(-5, 5, 500)
mu = 0
sigma = 1
y = 1 / (sigma * np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((x - mu) / sigma) ** 2)

source = ColumnDataSource(data=dict(x=x, y=y))

plot = figure(width=400, height=400)
line = plot.line('x', 'y', source=source, line_width=3, line_alpha=0.6)
plot.add_tools(PanTool(),WheelZoomTool(dimensions="width"))

callback = CustomJS(args=dict(source=source), code="""
        const data = source.data;
        const sigma = cb_obj.value;
        const x = data['x'];
        const y = data['y'];
        for (let i = 0; i < x.length; i++) {
            y[i] = 1 / (sigma * Math.sqrt(2 * Math.PI)) * Math.exp(-0.5 * Math.pow((x[i] - 0) / sigma, 2));
        }
        source.change.emit();
    """)

slider = Slider(start=0.0001, end=10, value=1, step=0.01, title="Variance")
slider.js_on_change('value', callback)


from bokeh.io import show
from bokeh.models import CustomJS, Dropdown

menu = [("Item 1", "item_1"), ("Item 2", "item_2"), None, ("Item 3", "item_3")]

dropdown = Dropdown(label="Dropdown button", button_type="warning", menu=menu)
dropdown.js_on_event("menu_item_click", CustomJS(code="console.log('dropdown: ' + this.item, this.toString())"))

from bokeh.io import show
from bokeh.models import CheckboxButtonGroup, CustomJS

LABELS = ["Option 1", "Option 2", "Option 3"]

checkbox_button_group = CheckboxButtonGroup(labels=LABELS, active=[0, 1])
checkbox_button_group.js_on_event("button_click", CustomJS(args=dict(btn=checkbox_button_group), code="""
    console.log('checkbox_button_group: active=' + btn.active, this.toString())
"""))

from bokeh.io import show
from bokeh.models import (ColumnDataSource, DataCube, GroupingInfo,
                          StringFormatter, SumAggregator, TableColumn)

source = ColumnDataSource(data=dict(
    d0=['A', 'E', 'E', 'E', 'J', 'L', 'M'],
    d1=['B', 'D', 'D', 'H', 'K', 'L', 'N'],
    d2=['C', 'F', 'G', 'H', 'K', 'L', 'O'],
    px=[10, 20, 30, 40, 50, 60, 70],
))

target = ColumnDataSource(data=dict(row_indices=[], labels=[]))

formatter = StringFormatter(font_style='bold')

columns = [
    TableColumn(field='d2', title='Name', width=80, sortable=False, formatter=formatter),
    TableColumn(field='px', title='Price', width=40, sortable=False),
]

grouping = [
    GroupingInfo(getter='d0', aggregators=[SumAggregator(field_='px')]),
    GroupingInfo(getter='d1', aggregators=[SumAggregator(field_='px')]),
]

cube = DataCube(source=source, columns=columns, grouping=grouping, target=target)

from bokeh.io import show
from bokeh.models import CustomJS, RangeSlider

range_slider = RangeSlider(start=0, end=10, value=(1,9), step=.1, title="Stuff")
range_slider.js_on_change("value", CustomJS(code="""
    console.log('range_slider: value=' + this.value, this.toString())
"""))

from bokeh.io import show
from bokeh.models import CustomJS, TextInput

text_input = TextInput(value="default", title="Label:")
text_input.js_on_change("value", CustomJS(code="""
    console.log('text_input: value=' + this.value, this.toString())
"""))


layout = column(slider, plot, checkbox_button_group, dropdown, cube, range_slider, text_input)


show(layout)



from bokeh.models import PanTool, WheelZoomTool
import numpy as np
from bokeh.plotting import figure, output_file, show, curdoc
from bokeh.models import Button, Div, CustomJS, TextInput
from bokeh.layouts import column

# Generate data for a Gaussian distribution with a very small standard deviation
initial_mu, initial_sigma = 0, 0.005 # Adjust sigma to make the Gaussian very sharp
x = np.linspace(-0.05, 0.05, 1000)  # Adjust x range accordingly
y = 1/(initial_sigma * np.sqrt(2 * np.pi)) * np.exp(-(x - initial_mu)**2 / (2 * initial_sigma**2))

# Output to a static HTML file
output_file("../../../Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Mass Spectrometry/interactive_gaussian_plot.html")

# Create a new plot with a title and axis labels
p = figure(title="Sharp Gaussian Distribution", x_axis_label='x', y_axis_label='Probability Density')

# Set the range of the y-axis
p.y_range.start = 0
p.y_range.end = 100

# Add a line renderer
p.line(x, y, line_width=2)

# Add tools for panning and zooming
p.add_tools(PanTool(),WheelZoomTool(dimensions="width"))

info = Div(text="<h2>Information:</h2><p>Hi my guy</p>")

# Create a button widget
button = Button(label="Information about the distribution")

button_callback = CustomJS(args=dict(info=info), code="""
    info.text = "<h2>Gaussian Distrib info</h2><p>mu = 0, sigma = 0.005</p>";
""")

# Add the callback to the button widget
button.js_on_click(button_callback)
renderer = p.line(x, y, line_width=2)

mu_input = TextInput(value=str(initial_mu), title="Enter mu:")
sigma_input = TextInput(value=str(initial_sigma), title="Enter sigma:")

def update_plot(attr, old, new):
    mu_val = float(mu_input.value)
    sigma_val = float(sigma_input.value)
    y = 1/(sigma_val * np.sqrt(2 * np.pi)) * np.exp(-(x - mu_val)**2 / (2 * sigma_val**2))
    renderer.data_source.data['y'] = y

# Attach the callback to the TextInput widgets
mu_input.on_change('value', update_plot)
sigma_input.on_change('value', update_plot)

# Create a layout with the plot, TextInput widgets, and information
layout = column(p, button, info, mu_input, sigma_input)


show(layout)

