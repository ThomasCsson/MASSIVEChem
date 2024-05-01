#Imports

from rdkit import Chem
from rdkit.Chem import Draw
import time
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from bokeh.plotting import figure, show
from bokeh.models import WheelPanTool, WheelZoomTool
from bokeh.models.tickers import FixedTicker




def bokeh_plotter(x_axis_final, y_axis_final):

    #---------------------------------------------------------------------------------------------#
    '''
    bokeh_plotter(x_axis_final, y_axis_final)
    
    Input: two lists:
    1. list of the masses (of individual molecules) of each possible combination of isotopes
    2. list of the probabilities of apparation of each of the molecules 
    (the mass in list 1 at index i is associated to the probability at index i in list 2) 
    
    Output: none

    Functionality: plots graph with bokeh (html format)
    '''
    #---------------------------------------------------------------------------------------------#

    eps = 10**(-20)

    mass_range = np.linspace(min(x_axis_final)-1, max(x_axis_final)+1, 1000)

    intensity = np.zeros_like(mass_range)

    for peak_position, peak_intensity in zip(x_axis_final, y_axis_final):
        peak_shape = peak_intensity * lorentzian(mass_range, peak_position, eps)
        intensity += peak_shape



    ticked_peaks = []
    for i in range(len(x_axis_final)):
        if y_axis_final[i]>0.0001:
            ticked_peaks.append(round(x_axis_final[i],4))


    # Create a new plot with a title and axis labels
    p = figure(title="Simulated Mass Spectrum", x_axis_label='Mass [m/z]', y_axis_label='Intensity [AU]')
    p = figure(width=700 , title= f'Mass spectrum of molecule')
    p.height = 500
    p.xaxis.ticker = FixedTicker(ticks= ticked_peaks)
    p.toolbar.autohide = False
    p.add_tools(WheelPanTool(dimension="height"))
    p.add_tools(WheelZoomTool(dimensions="height"))

    # Add a line renderer with legend and line thickness
    '''p.line(mass_range, intensity, legend_label="Intensity", line_width=1)'''
    p.line(x_axis_final, y_axis_final, legend_label = "Intensity", line_width=1)

    # Show the plot
    show(p)
    print('')
    return  (f'Computation complete.')