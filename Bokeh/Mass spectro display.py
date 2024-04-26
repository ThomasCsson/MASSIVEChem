from bokeh.plotting import figure, show
import numpy as np


peak_positions= [78.918335, 80.91629]
peak_intensities= [0.5069, 0.49310000000000004]


mass_range = np.linspace(min(peak_positions)-1, max(peak_positions)+1, 1000)

intensity = np.zeros_like(mass_range)

for peak_position, peak_intensity in zip(peak_positions, peak_intensities):

    peak_shape = peak_intensity * np.exp(-((mass_range - peak_position) ** 2) / (2 * 0.02 ** 2))  # Gaussian example


    intensity += peak_shape

# Create a new plot with a title and axis labels
p = figure(title="Simulated Mass Spectrum", x_axis_label='Mass', y_axis_label='Intensity')

# Add a line renderer with legend and line thickness
p.line(mass_range, intensity, legend_label="Intensity", line_width=2)

# Show the plot
show(p)









