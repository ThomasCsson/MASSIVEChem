from bokeh.plotting import figure, show
from bokeh.io.export import get_screenshot_as_png
import os

# Create your Bokeh plot
p = figure()

# Add your plot elements here...

# Show the plot
show(p)

# Get the path to the HTML file
html_file_path = show(p, notebook_handle=False)

# Delete the HTML file
os.remove(html_file_path)