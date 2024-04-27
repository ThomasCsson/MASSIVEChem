

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


