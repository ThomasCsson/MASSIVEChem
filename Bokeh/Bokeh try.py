

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

layout = column(slider, plot)

show(layout)


