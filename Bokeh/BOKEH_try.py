from bokeh.plotting import figure, show
from bokeh.embed import file_html
import os

# Create your Bokeh plot
p = figure()

# Add your plot elements here...

# Generate the HTML content
html_content = file_html(p, "my_plot")

# Specify the output file path
output_file_path = "my_plot.html"

# Write the HTML content to a file
with open(output_file_path, "w") as html_file:
    html_file.write(html_content)

# Show the plot (if necessary)
show(p)

# Remove the HTML file
if os.path.exists(output_file_path):
    os.remove(output_file_path)


