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


def data_list_generator():

    #---------------------------------------------------------------------------------------------#
    '''
    data_list_generator ()

    Input: none

    Output: three lists:
    1. mass: [mass1, mass2, mass3,...]  
    2. abundance: [ab1, ab2, ab3,...]
    3. isotopes: [iso1, iso2, iso3,...]
    '''
    #---------------------------------------------------------------------------------------------#


    #Turn data of (Symbol | Mass | Probability) into lists 

    df = pd.read_csv('/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/THOMAS_IS_BEST/abundance.txt'
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
    #---------------------------------------------------------------------------------------------#
    '''
    SMILEs_interpreter(mol_smi)

    Input: molecule under SMILEs representation
    
    Output: molecule under MOL representation
    '''
    #---------------------------------------------------------------------------------------------#

    #Checks that the SMILEs input is correct/allowed

    mol_without_Hs = Chem.MolFromSmiles(mol_smi)

    if mol_without_Hs is None:
        print('')
        print("Invalid SMILEs input.")
        print('Please try again with a different SMILEs.')
        exit()

    mol = Chem.AddHs(mol_without_Hs)

    return mol

def molecule_list_generator(mol):
    #---------------------------------------------------------------------------------------------#
    '''
    molecule_list_generator(mol)
    
    Input: molecule under MOL representation
    
    Output: list containing the atomic symbol of each atom in the input molecule
    '''
    #---------------------------------------------------------------------------------------------#

    list_atoms = []
    for atom in mol.GetAtoms():
        list_atoms.append(atom.GetSymbol())
    return list_atoms


def ionisation_method (list_atoms):
    #---------------------------------------------------------------------------------------------#
    '''
    ionisation_method (list_atoms)
    
    Input: list of atomic symbols of atoms in a given molecule
    
    Output: corrected list of atoms (i.e. list of atoms that enter spectrometry apparatus)
    '''
    #---------------------------------------------------------------------------------------------#

    '''In the case of ionisation by proton, we need to add a H+ ion, which is done in the following'''
    if 'H' in list_atoms:

        #Check that there is in fact a proton to remove
        list_atoms.remove('H')
    return list_atoms






def main_function (list_atoms, imprecision_True_False):
    #---------------------------------------------------------------------------------------------#
    '''
    main_function(list_atoms)
    
    Input: list of atoms that enter the apparatus
    
    Output: two lists:
    1. list of the masses (of individual molecules) of each possible combination of isotopes
    2. list of the probabilities of apparation of each of the molecules 
    (the mass in list 1 at index i is associated to the probability at index i in list 2)
    '''
    #---------------------------------------------------------------------------------------------#
    render_imprecise_list = imprecision_True_False #Set arg to be True for long molecules, set arg to False for short molecules/if precision for minuscule peaks is important

    #check for sulphur and nitrogen

    has_N = False
    count_N = 0
    has_S = False
    count_S = 0
    if 'N' in list_atoms and list_atoms.count('N')%2 == 1:
        has_N = True
        count_N = list_atoms.count('N')
    elif 'S' in list_atoms and list_atoms.count('S')%2 == 1:
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


                #removes any molecule who's probability is below 0.0001

                if render_imprecise_list: #only removes low-probability arrangements if render_imprecise_list arg is True
                    if new_proba>0.0001:
                        list_output_new.append([new_mass,new_proba])

                else:
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
        '''x_axis.append(list_output[j][0])'''
        x_axis.append(round(list_output[j][0],3)) # Adds rounded value (should help with Python-limitations that render a diff of magnitude 10^(-7) to combinatorics
        
        y_axis.append(list_output[j][1]) #Adds the true value
        

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

def peak_merger (x_axis_final, y_axis_final):
    return

def matplotlib_plotter(x_axis_final, y_axis_final):

    #---------------------------------------------------------------------------------------------#
    '''
    matplotlib_plotter(x_axis_final, y_axis_final)
    
    Input: two lists:
    1. list of the masses (of individual molecules) of each possible combination of isotopes
    2. list of the probabilities of apparation of each of the molecules 
    (the mass in list 1 at index i is associated to the probability at index i in list 2)
    
    Output: none

    Functionality: plots a graph made of dots with matplotlib
    '''
    #---------------------------------------------------------------------------------------------#

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

    #---------------------------------------------------------------------------------------------#
    '''
    pyplot_plotter (x_axis_final, y_axis_final)
    
    Input: two lists:
    1. list of the masses (of individual molecules) of each possible combination of isotopes
    2. list of the probabilities of apparation of each of the molecules 
    (the mass in list 1 at index i is associated to the probability at index i in list 2) 
    
    Output: p (plot of mass spectrum on bokeh)

    Functionality: none
    '''
    #---------------------------------------------------------------------------------------------#

    x, y = x_axis_final, y_axis_final
    fig = go.Figure()
    fig.add_trace(go.Scatter(x = x,y = y, mode = 'markers'))

    fig.update_layout(xaxis=dict(rangeslider=dict(visible=True), 
                                type="linear"),
                                yaxis=dict(range=[min(y)-1, max(y)], 
                                type="linear"),
                                dragmode='zoom',
                                )


    fig.show()
    return 
def lorentzian(x, eps):
    return 1.0 / (np.pi * eps * (1 + (x / eps) ** 2))

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

    eps = 10**(-4)

    mass_range = np.linspace(min(x_axis_final)-1, max(x_axis_final)+1, 1000)

    intensity = np.zeros_like(mass_range)

    for peak_position, peak_intensity in zip(x_axis_final, y_axis_final):
        peak_shape = peak_intensity * lorentzian(mass_range - peak_position, eps)  # Lorentzian example

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
    p.line(mass_range, intensity, legend_label="Intensity", line_width=1)

    # Show the plot
    show(p)
    print('')
    return  (f'Computation complete.')

print('')

mol_smi = input('Enter SMILEs: ')

start_time = time.time()

mol = SMILEs_interpreter(mol_smi)
mass, abundance, isotopes = data_list_generator()
list_atoms_pre = molecule_list_generator(mol) 
list_atoms = ionisation_method(list_atoms_pre)
xvalues, yvalues = main_function(list_atoms, True)

end_time = time.time()

duration = end_time-start_time

print(bokeh_plotter(xvalues,yvalues))
print(f'Process took: {duration} s')
