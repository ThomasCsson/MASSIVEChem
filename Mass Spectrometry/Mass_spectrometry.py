import pandas as pd
import numpy as np

with open('Natural abudance', 'r', encoding='ISO-8859-1') as file:
    # Sauter la première ligne
    next(file)
    lines = file.readlines()

colonne1 = []
colonne2 = []
colonne3 = []


# Parcours des lignes et extraction des valeurs
for line in lines:
    values = line.split()  # Séparation des valeurs par espace

    # Vérifier si la ligne contient le nombre attendu de valeurs
    if len(values) >= 2:
        # Remplacement des virgules par des points pour les décimales et conversion en float
        colonne1.append(values[0])
        colonne2.append(float(values[1]))
        colonne3.append(float(values[2]))


Mass_Abundance = {
    'Atom': colonne1,
    'Mass': colonne2,
    '%': colonne3,

}

Mass_abundance_DF = pd.DataFrame(Mass_Abundance)
print(Mass_abundance_DF)

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