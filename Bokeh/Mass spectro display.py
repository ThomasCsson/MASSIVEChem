#import of packages
from bokeh.plotting import figure, show
from bokeh.models import Range1d  # Import Range1d for specifying ranges
import numpy as np

def peak_posi(list_mass_abundance):

    masses = []
    intensities = []

    for i in range(0,len(list_mass_abundance)):

        masses.append(list_mass_abundance[i][0])
        intensities.append(list_mass_abundance[i][0])

    aggregated_intensities = {'Peak positions': [], 'Peak Intensities': []}

    for mass, intensity in zip(masses, intensities):

        if mass in aggregated_intensities['Peak positions']:
            idx = aggregated_intensities['Peak positions'].index(mass)
            aggregated_intensities['Peak Intensities'][idx] += intensity

        else:
            aggregated_intensities['Peak positions'].append(mass)
            aggregated_intensities['Peak Intensities'].append(intensity)

    return aggregated_intensities

    def mass_spectrum(aggregated_intensities)
mass_range = np.linspace(0, 500, 1000)


intensity = np.zeros_like(mass_range)


for peak_position in aggre:
    # Define the peak shape using Lorentzian distribution
    peak_shape = 1 / (1 + ((mass_range - peak_position) / 10) ** 2)  # Lorentzian example (width = 10)

    # Add the peak to the intensity array
    intensity += peak_shape

# Create a new plot with a title and axis labels
p = figure(title="Simulated Mass Spectrum", x_axis_label='Mass', y_axis_label='Intensity')

# Add a line renderer with legend and line thickness
p.line(mass_range, intensity, legend_label="Intensity", line_width=2)

# Set the range of visualization for the x-axis and y-axis using Range1d
p.x_range = Range1d(0, 500)  # Set the range of the x-axis from 0 to 500
p.y_range = Range1d(0, intensity.max())  # Set the range of the y-axis from 0 to the maximum intensity value

# Show the plot
show(p)
