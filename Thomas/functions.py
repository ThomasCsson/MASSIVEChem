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


def list_generator():
    #Turn data of (Symbol | Mass | Probability) into lists 

    df = pd.read_csv('Thomas/abundance.txt'
                    , sep='\t'
                    , header=None
                    , names=['Atom', 'Mass', 'Percentage'])


    #mass = [mass1, mass2, mass3,...]

    mass = df['Mass'].tolist()


    #change from elements (not sure) to floats

    mass = [float(m) for m in mass]


    #abundance = [ab1, ab2, ab3,...]

    abundance_percent = df['Percentage'].tolist()


    #change from elements (not sure) to floats

    abundance_percent = [float(ap) for ap in abundance_percent]


    #from percent to proba

    abundance = []
    for percent in abundance_percent:
        abundance.append(percent/100)


    #isotopes = [iso1, iso2, iso3,...]

    isotopes = df['Atom'].tolist()
    return mass, abundance, isotopes

def SMILEs_interpreter(mol_smi):
    #Turn SMILEs representation into list of atomic symbols

    mol_without_Hs = Chem.MolFromSmiles(mol_smi)
    mol = Chem.AddHs(mol_without_Hs)

    image = Draw.MolToImage(mol)
    image.show()
    return mol


def main_function (mol):
    list_atoms = []
    for atom in mol.GetAtoms():
        list_atoms.append(atom.GetSymbol())

    '''In the case of ionisation by proton, we need to add a H+ ion, which is done in the following'''
    
    if 'H' in list_atoms:
        #Check that there is in fact a proton to remove
        list_atoms.remove('H')



    #check for sulphur and nitrogen

    has_N = False
    count_N = 0
    has_S = False
    count_S = 0
    if 'N' in list_atoms:
        has_N = True
        count_N = list_atoms.count('N')
    elif 'S' in list_atoms:
        has_S = True
        count_S = list_atoms.count('S')



    

    print(list_atoms)

    list_output = []
    mass_copy = mass.copy()
    abundance_copy = abundance.copy()
    isotopes_copy = isotopes.copy()

    for i in range (isotopes.count(list_atoms[0])):
        
        index = isotopes.index(list_atoms[0])
        list_output.append([mass_copy[index],abundance_copy[index]])
        mass_copy.pop(index)
        abundance_copy.pop(index)
        isotopes_copy.pop(index)

    list_atoms = list_atoms[1:]

    while len(list_atoms)>0:


        #This runs over all atoms in molecule

        list_output_new = []

        for i in range (len(list_output)):


            #This makes us run over all lists in list obtained before

            mass_copy = mass.copy()
            abundance_copy = abundance.copy()
            isotopes_copy = isotopes.copy() 

            for _ in range (isotopes.count(list_atoms[0])):


                #This for-loop runs over all isotope types of the atom type in pos 0 in list_atoms (input list)

                index = isotopes.index(list_atoms[0])
                new_mass = list_output[i][0] + mass_copy[index]
                new_proba = list_output[i][1] * abundance_copy[index]
                list_output_new.append([new_mass,new_proba])
                mass_copy.pop(index)
                abundance_copy.pop(index)
                isotopes_copy.pop(index)
                
        list_output = list_output_new
        list_atoms.pop(0)


    #Conversion of list_output (which is a list of lists) to a combination of two lists (x_axis & y_axis)

    x_axis, y_axis = [],[]
    x_axis_final, y_axis_final = [],[]
    for j in range (len(list_output)):
        x_axis.append(list_output[j][0])
        y_axis.append(list_output[j][1])


    #Compression of lists x_axis & y_axis into x_axis_final & y_axis_final so that peaks corresponding to same mass will be represented together 
    
    for j in range (len(x_axis)):
            if x_axis.count(x_axis[j]) == 1:
                x_axis_final.append(x_axis[j])
                y_axis_final.append(y_axis[j])
            elif x_axis_final.count(x_axis[j]) == 0:
                x_axis_final.append(x_axis[j])
                y_axis_final.append(y_axis[j])
            else:
                index = x_axis_final.index(x_axis[j])
                y_axis_final[index] =y_axis_final[index] + y_axis[j]


    #if there is any, add peaks corresponding to Sulphur/Nitrogen presence
    maximum = max(y_axis_final)
    maximum_2 = 0
    for i in range (len(y_axis_final)):
        if y_axis_final[i]>maximum_2 and y_axis_final[i]< maximum:
            maximum_2 = y_axis_final[i]
    index = y_axis_final.index(maximum_2)

    
    

    if has_N:
        x_axis_final.append(x_axis_final[index] - 0.006)  
        y_axis_final.append(0.0035*count_N*maximum)  
    
    if has_S:
        x_axis_final.append(x_axis_final[index]-0.004)  
        y_axis_final.append(0.008*count_S*maximum)  
    return x_axis_final, y_axis_final

def matplotlib_plotter(x_axis_final, y_axis_final):
    max_x = max(x_axis_final)
    min_x = min(x_axis_final)
    diff = (max_x-min_x)
    x_axis_final_use = x_axis_final.copy()
    y_axis_final_use = y_axis_final.copy()

    for i in range(int(-(100*diff)/10),int(11*(100*diff)/10)):
        x_axis_final_use.append(i/100 + min_x)
        y_axis_final_use.append(0)
    x_final_final = []
    y_final_final = []
    while len(x_axis_final_use)>1 :
        minx = min(x_axis_final_use)
        index = x_axis_final_use.index(minx)
        x_final_final.append(minx)
        y_final_final.append(y_axis_final_use[index])
        x_axis_final_use.pop(index)
        y_axis_final_use.pop(index)


    #plotting with matpotlib
    plt.scatter(x_axis_final,y_axis_final)
    plt.show()
    return

def pyplot_plotter (x_axis_final, y_axis_final):
    x, y = x_axis_final, y_axis_final
    fig = go.Figure()
    fig.add_trace(go.Scatter(x = x,y = y, mode = 'markers'))

    fig.update_layout(xaxis=dict(rangeslider=dict(visible=True), type="linear"),yaxis=dict(range=[min(y)-1, max(y)], type="linear"),dragmode='zoom',)


    fig.show()
    return 

def bokeh_plotter(x_axis_final, y_axis_final):
    x = 0.004

    mass_range = np.linspace(min(x_axis_final)-1, max(x_axis_final)+1, 1000)

    intensity = np.zeros_like(mass_range)

    for peak_position, peak_intensity in zip(x_axis_final, y_axis_final):

        peak_shape = peak_intensity * np.exp(-((mass_range - peak_position) ** 2) / (2 * x ** 2))  # Gaussian example

        intensity += peak_shape


    ticked_peaks = []
    for i in range(len(x_axis_final)):
        if y_axis_final[i]>0.0001:
            ticked_peaks.append(x_axis_final[i])

    print(ticked_peaks)

    # Create a new plot with a title and axis labels
    p = figure(title="Simulated Mass Spectrum", x_axis_label='Mass [Th]', y_axis_label='Intensity')
    p = figure(width=700 , title= f'Mass spectrum of molecule')
    p.height = 500
    p.xaxis.ticker = FixedTicker(ticks= ticked_peaks)
    p.toolbar.autohide = True
    p.add_tools(WheelPanTool(dimension="height"))
    p.add_tools(WheelZoomTool(dimensions="height"))

    # Add a line renderer with legend and line thickness
    p.line(mass_range, intensity, legend_label="Intensity", line_width=1)

    # Show the plot
    show(p) 
    return 


mol_smi = input('Enter SMILEs: ')
mol = SMILEs_interpreter(mol_smi)
mass, abundance, isotopes = list_generator()
xvalues, yvalues = main_function(mol)
print(bokeh_plotter(xvalues,yvalues))
