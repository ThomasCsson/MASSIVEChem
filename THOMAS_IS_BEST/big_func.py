#Imports

from rdkit import Chem
from rdkit.Chem import Draw
import time
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from bokeh.plotting import figure, show
from bokeh.models import WheelPanTool, WheelZoomTool
from bokeh.models.tickers import FixedTicker

def big_function (mol_smi, imprecision_True_False, apparatus_resolution):


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

    df = pd.read_csv('data/abundance.txt'
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

    #---------------------------------------------------------------------------------------------#
    '''
    SMILEs_interpreter(mol_smi)

    Input: molecule under SMILEs representation
    
    Output: interface with mass spectrum of the given molecule as well as the contained functional groups
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
    print(list_output)


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

    

    x_in, y_in = x_axis_final, y_axis_final




    #---------------------------------------------------------------------------------------------#
    '''
    delta_function_plotter(x_in, y_in)
    
    Input: two lists:
    1. list of the masses (of individual molecules) of each possible combination of isotopes
    2. list of the probabilities of apparation of each of the molecules
    
    Output: two lists:
    1. ordered list of the masses (of individual molecules) of each possible combination of isotopes
    2. ordered list of the probabilities of apparation of each of the molecules

    (ordered to have increasing values along x)
    
    (the mass in list 1 at index i is associated to the probability at index i in list 2)
    '''
    #---------------------------------------------------------------------------------------------#
    
    x_out, y_out = [],[]

    while len(x_in)>0:
        min_x = min(x_in)
        index_min = x_in.index(min_x)
        x_out.append(min_x)
        y_out.append(y_in[index_min])
        x_in.pop(index_min)
        y_in.pop(index_min)



    x_in, y_in = x_out, y_out 



    #---------------------------------------------------------------------------------------------#
    '''
    delta_function_plotter(x_in, y_in)
    
    Input: two lists + float:
    1. ordered list of the masses (of individual molecules) of each possible combination of isotopes
    2. ordered list of the probabilities of apparation of each of the molecules
    3. apparatus precision (float representation of number which is the limit between two peaks above which they appear merged together)
    
    Output: two lists:
    1. ordered list of the masses with values on y merged together if peaks within precision of apparatus
    2. ordered list of the probabilities of apparation of each of the molecules
    
    (the mass in list 1 at index i is associated to the probability at index i in list 2)
    '''
    #---------------------------------------------------------------------------------------------#

    x_out, y_out = [],[]
    while len(x_in)>1:
        if x_in[0]>x_in[1]-apparatus_resolution:
            y_in[1] = y_in[0] + y_in[1]
            x_in[1] = (x_in[0] + x_in[1])/2
            x_in.pop(0)
            y_in.pop(0)
        else:
            x_out.append(x_in[0])
            y_out.append(y_in[0])
            x_in.pop(0)
            y_in.pop(0)
    x_out.append(x_in[0])
    y_out.append(y_in[0])

    x_in, y_in = x_out, y_out 
    #---------------------------------------------------------------------------------------------#
    '''
    delta_function_plotter(x_in, y_in)
    
    Input: two lists:
    1. list of the masses (of individual molecules) of each possible combination of isotopes (ordered)
    2. list of the probabilities of apparation of each of the molecules (ordered)
    
    Output: two lists:
    1. list of the masses (of individual molecules) of each possible combination of isotopes (ordered) with values of 0 (y_axis) added on eiter side of the "peak"
    2. list of the probabilities of apparation of each of the molecules (ordered)
    
    (the mass in list 1 at index i is associated to the probability at index i in list 2)
    '''
    #---------------------------------------------------------------------------------------------#

    min_x , max_x = min(x_in), max(x_in)

    x_axis, y_axis = [min_x-1],[0]
    for i in range (len(x_in)):
        x_axis.append(x_in[i]-10**(-10))
        x_axis.append(x_in[i])
        x_axis.append(x_in[i]+10**(-10))
        y_axis.append(0)
        y_axis.append(y_in[i])
        y_axis.append(0)

    x_axis.append(max_x+1)
    y_axis.append(0)


    #add graphing
    return x_axis, y_axis
    




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

    min_x = min(x_axis)
    max_x = max(x_axis)
    x_axis, y_axis = [min_x-1],[0]
    for i in range (len(x_axis_final)):
        x_axis.append(x_axis_final[i]-10**(-10))
        x_axis.append(x_axis_final[i])
        x_axis.append(x_axis_final[i]+10**(-10))
        y_axis.append(0)
        y_axis.append(y_axis_final[i])
        y_axis.append(0)
    x_axis.append(max_x+1)
    y_axis.append(0)
        

    #plotting with matpotlib
    plt.plot(x_axis,y_axis)
    
    return x_axis, y_axis

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
    fig.add_trace(go.Scatter(x = x,y = y, mode = 'lines'))

    for peak in x:
        if x.index(peak) > 0 and x.index(peak) < len(x)-1:
            if y_axis_final[x_axis_final.index(peak)] >0.001:
                fig.add_annotation(x=peak,y = 0,
                        text=round(peak,3), 
                        showarrow=True, arrowhead=1, ax=0, ay=30)
    for peak in x:
        if x.index(peak) > 0 and x.index(peak) < len(x)-1:
            if y_axis_final[x_axis_final.index(peak)] >0.001:
                fig.add_annotation(x=peak, y=y_axis_final[x_axis_final.index(peak)],
                        text=round(y_axis_final[x_axis_final.index(peak)],3), 
                        showarrow=True, arrowhead=1, ax=0, ay=-30)
                

    fig.update_layout(title = 'Mass spectrum of input molecule',
                    xaxis_title = '[m/z]',
                    yaxis_title = 'Abundance',   
                    xaxis=dict(range=[min(x), max(x)], type="linear"),
                    yaxis=dict(range=[min(y)-0.1, max(y)+0.1], type="linear"),
                    )
    

    fig.add_trace(go.Scatter(x=x, y=y, mode='lines', showlegend=False, xaxis='x2', yaxis='y2'))


    fig.update_layout(
        xaxis2=dict(domain=[0.7, 1], anchor='y2'),
        yaxis2=dict(domain=[0.7, 1], anchor='x2'),
        xaxis_title='[m/z]',
        yaxis_title='Abundance')
    
    fig.show()
    return 
def test(x,y):
    
    fig = go.Figure()
    fig.add_trace(go.Scatter(x = x,y = y, mode = 'lines'))
    initial_zoomed_x = x[5:15]  # Adjust the initial range as needed
    initial_zoomed_y = y[5:15]  # Adjust the initial range as needed

    # Add trace for the initial zoomed-in region
    fig.add_trace(go.Scatter(x=initial_zoomed_x, y=initial_zoomed_y, mode='lines', showlegend=False, xaxis='x2', yaxis='y2'))

    # Update layout for the main plot
    fig.update_layout(
        xaxis2=dict(domain=[0.7, 1], anchor='y2'),
        yaxis2=dict(domain=[0.7, 1], anchor='x2'),
        xaxis_title='[m/z]',
        yaxis_title='Abundance'
    )

    # Add rectangle to the smaller graph to indicate the initial zoomed-in region
    fig.add_shape(
        type="rect",
        xref="x2",
        yref="y2",
        x0=initial_zoomed_x[0],  # Initial coordinates based on the initial zoomed-in region
        y0=min(initial_zoomed_y),
        x1=initial_zoomed_x[-1],
        y1=max(initial_zoomed_y),
        line=dict(color="rgba(0,0,0,0.5)", width=2),
        fillcolor="rgba(0,0,0,0.1)",
        visible=True,  # Initially visible
        layer="below"
    )

    # Define a function to update the zoomed-in region in the smaller graph based on the main graph's zoom level
    def update_zoomed_region(layout, x_range, y_range):
        fig.layout.shapes[0].x0 = x_range[0]
        fig.layout.shapes[0].x1 = x_range[1]
        fig.layout.shapes[0].y0 = y_range[0]
        fig.layout.shapes[0].y1 = y_range[1]

    # Define the callback for relayout event
    fig.update_layout(
        xaxis=dict(domain=[0, 0.65]),
        yaxis=dict(domain=[0, 0.65]),
        xaxis2=dict(domain=[0.7, 1], anchor='y2'),
        yaxis2=dict(domain=[0.7, 1], anchor='x2'),
        xaxis_title='[m/z]',
        yaxis_title='Abundance',
        dragmode='zoom',
        autosize=True,
        margin=dict(l=0, r=0, t=0, b=0),
        showlegend=False,
        plot_bgcolor='white'
    )

    fig.update_xaxes(matches='x')
    fig.update_yaxes(matches='y')

    # Add event handler for relayout event
    fig.update_layout(
        scene=dict(onclick=lambda eventdata: update_zoomed_region(fig.layout, eventdata['xaxis.range'], eventdata['yaxis.range']))
    )

    fig.show()



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
    ticked_peaks = []
    for i in range(len(x_axis_final)):
        if y_axis_final[i]>0.0001:
            ticked_peaks.append(round(x_axis_final[i],4))
    p = figure(title="Simulated Mass Spectrum", x_axis_label='Mass [m/z]', y_axis_label='Intensity [AU]')
    p = figure(width=700 , title= f'Mass spectrum of molecule')
    p.height = 500
    p.xaxis.ticker = FixedTicker(ticks= ticked_peaks)
    p.toolbar.autohide = True
    p.add_tools(WheelPanTool(dimension="height"))
    p.add_tools(WheelZoomTool(dimensions="height"))
    p.line(x_axis_final, y_axis_final, legend_label = "Intensity", line_width=1)
    show(p)
    return  


def plotly_test(x_in, y_in):


    fig = make_subplots(
    rows=2, cols=2,
    specs=[[{}, {}],
           [{"colspan": 2}, None]],
    subplot_titles=("First Subplot","Second Subplot", "Third Subplot"))


    image = Draw.MolToImage(Chem.MolFromSmiles('C=O'))

    fig.add_trace(go.image(z = image),row = 1, col = 1)

    fig.add_trace(go.Scatter(x=x_in, y=y_in, mode = 'lines'),
                    row=1, col=2)
    fig.add_trace(go.Scatter(x=x_in, y=y_in, mode = 'lines'),
                    row=2, col=1)

    fig.update_layout(showlegend=False, title_text="Specs with Subplot Title")
    fig.show()

    return

print('')

def functional_group_finder(mol_smi):

    #---------------------------------------------------------------------------------------------#
    '''
    functional_group_finder(mol_smi)
    
    Input: molecule under SMILEs representation
    
    Output: list containing every functionl group contained (if a functional group is contained twice in the molecule, it will appear twice in this list)
    '''
    #---------------------------------------------------------------------------------------------#


    functional_groups_contained, mol_in = [], Chem.MolFromSmiles(mol_smi)
    functional_groups_smarts = {
        'Alcohol': 'C[Oh1+0]','Aldehyde': 'C[Ch1]=O','Ketone': 'CC(=O)C','Carboxylic Acid': 'CC(=O)[Oh1]',
        'Ester': 'CC(=O)[Oh0]','Ether': '*[Oh0]*','Amide': 'C(=O)N','Amine': '[C][Nh2]','Nitrile': 'C#N',
        'Chloride': 'Cl','Bromide': 'Br','Fluoride': 'F','Iodide': 'I','Alkene': 'C=C','Alkyne': 'C#C',
        'Imine': 'C=N*','Amino acid': '[Nh2][Ch1*]C(=O)O','Proline': '[Nh1][Ch1*]C(=O)O','Thiol': '[Sh1]',
        'Sulfide': '*[Sh0]*','Acyl Chloride': 'CC(=O)Cl','Anhydride': '*[Ch0](=O)O[Ch0](=O)*','Nitro': 'C[N+](=O)[O-]',
        'Enamine': 'C=C[Nh0]','Enamine2': 'C=C[Nh1]','Enamine3': 'C=C[Nh2]','Imide': 'C(=O)NC(=O)*','Azide': 'CNNN',
        'Enol': 'C=C([Oh1])C','Hemiacetal': 'CC(O)(O)C','Carbonate': '[Oh0]C(=O)[Oh0]','Carbonate2': '[Oh1]C(=O)[Oh1]',
        'Disulfide': 'CSSC','Sulfoxide': 'CS(=O)C','Sulfone': '*[So2](=O)(=O)*','Sulfonic acid': '*S(=O)(=O)[Oh1]',
        'Thioester': 'C(=O)S*','Phosphine': '*[Po0](*)*','Phosphate': '*OP(=O)(O)O','Benzene': 'c1ccccc1','Peroxide':'C[Oh0][Oh0]C'
    }

    
    for name, smarts in functional_groups_smarts.items():
        if mol_in.HasSubstructMatch(Chem.MolFromSmarts(smarts)):
            for _ in range(len(mol_in.GetSubstructMatches(Chem.MolFromSmarts(smarts)))):
                functional_groups_contained.append(name)
            
    if 'Carboxylic Acid' in functional_groups_contained:
        functional_groups_contained.remove('Alcohol')
    if 'Ester' in functional_groups_contained:
        functional_groups_contained.remove('Ether')
    if 'Phosphate' in functional_groups_contained:
        functional_groups_contained.remove('Ether')
    if 'Thioester' in functional_groups_contained:
        functional_groups_contained.remove('Sulfide')
    if 'Sulfonic acid' in functional_groups_contained:
        functional_groups_contained.remove('Sulfide')
    if 'Sulfoxide' in functional_groups_contained:
        functional_groups_contained.remove('Sulfide')
    if 'Acyl Chloride' in functional_groups_contained:
        functional_groups_contained.remove('Chloride')
    if 'Anhydride' in functional_groups_contained:
        functional_groups_contained.remove('Ester')
        functional_groups_contained.remove('Ester')
    if 'Enamine2' in functional_groups_contained:
        functional_groups_contained.remove('Enamine2')
        functional_groups_contained.append('Enamine')
    if 'Enamine3' in functional_groups_contained:
        functional_groups_contained.remove('Enamine3')
        functional_groups_contained.remove('Amine')
        functional_groups_contained.append('Enamine')
    if 'Imide' in functional_groups_contained:
        functional_groups_contained.remove('Amide')
        functional_groups_contained.remove('Amide')
    if 'Enol' in functional_groups_contained:
        functional_groups_contained.remove('Alkene')
        functional_groups_contained.remove('Alcohol')
    if 'Hemiacetal' in functional_groups_contained:
        functional_groups_contained.remove('Alcohol')
        functional_groups_contained.remove('Alcohol')
    if 'Carbonate2' in functional_groups_contained:
        functional_groups_contained.remove('Alcohol')
        functional_groups_contained.remove('Alcohol')
        functional_groups_contained.remove('Carbonate2')
        functional_groups_contained.append('Carbonate')
    if 'Disulfide' in functional_groups_contained:
        functional_groups_contained.remove('Sulfide')
        functional_groups_contained.remove('Sulfide')
    if 'Amine2' in functional_groups_contained:
        functional_groups_contained.remove('Amine2')
        functional_groups_contained.append('Amine')
    if 'Peroxide' in functional_groups_contained:
        functional_groups_contained.remove('Ether')
        functional_groups_contained.remove('Ether')
    
    return functional_groups_contained






#Actually makes code run


mol_smi = input('Enter SMILEs: ')






#Generate the data from data/abundance.txt




#Keeps two lists but adds zeros on y next to each point on x
x_axis, y_axis = big_function(mol_smi, True, 0.01)



print(x_axis, y_axis)

