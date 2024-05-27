from rdkit import Chem
from rdkit.Chem import Draw, AllChem

import base64

import os

from bokeh.plotting import figure, row
from bokeh.models import ColumnDataSource, HTMLTemplateFormatter, WheelPanTool, WheelZoomTool, BoxAnnotation, CustomJS, Div
from bokeh.models.tickers import FixedTicker
from bokeh.models.widgets import DataTable, TableColumn
from bokeh.layouts import row, column

from io import BytesIO

import tempfile

import panel as pn

#main functions:
#1
def calculate_unsaturation(mol_smile) -> int:
    #---------------------------------------------------------------------------------------------#
    '''
    calculate_unsaturation(mol_smile)

    Input: molecule under SMILEs representation

    Output: unsaturation of the input molecule (integer value)
    '''
    #---------------------------------------------------------------------------------------------# 
    if mol_smile == None:
        raise ValueError('Enter a non-empty input')
    
    C, N, HX, halogens_hydrogen = 0, 0, 0, ['F', 'Cl', 'Br', 'I', 'At', 'H']

    mol = Chem.AddHs(Chem.MolFromSmiles(mol_smile))
    
    if mol == None:
        raise ValueError('Invalid SMILEs representation')
    
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':
            C += 1
        elif atom.GetSymbol() == 'N':
            N += 1
        elif atom.GetSymbol() in halogens_hydrogen:
            HX += 1

    unsaturation = C + 1 + (N - HX) / 2

    return unsaturation

#2
def spectrum(mol_smi, imprecision_True_False, apparatus_resolution):

    #tests for wether the input is valid
    if mol_without_Hs is None:
        raise ValueError('\nInvalid SMILEs enterred.\nPlease enter a different SMILEs.')
    if imprecision_True_False not in [True, False]:
        raise ValueError('Enter a boolean value')
    if apparatus_resolution < 0:
        raise ValueError('Enter a positive value')
    if type(apparatus_resolution) != float:
        raise ValueError('Enter an integer value')


    #lists of the data to facilitise the pip-installability of the package

    mass = [108.904755, 106.90509, 26.981539, 39.962383, 37.96273, 35.967545, 74.92159, 196.96654, 11.009305, 10.012937, 137.90523, 
            136.9058, 135.90456, 134.90567, 133.90448, 131.90504, 129.90628, 9.012182, 208.98038, 78.918335, 80.91629, 13.003355, 
            12.0, 45.95369, 43.95548, 42.958767, 41.95862, 39.96259, 47.952534, 110.90418, 109.90301, 107.90418, 115.904755, 
            113.90336, 112.9044, 111.902756, 141.90924, 139.90543, 137.90599, 135.90714, 36.965904, 34.968853, 58.933197, 53.93888, 
            52.94065, 51.94051, 49.946045, 132.90543, 62.939598, 157.92441, 155.92528, 163.92917, 162.92873, 161.92679, 160.92693, 
            159.92519, 169.93546, 167.93237, 166.93205, 165.93028, 163.9292, 161.92877, 152.92122, 150.91985, 18.998404, 57.933277, 
            56.935394, 53.939613, 70.9247, 68.92558, 157.9241, 156.92395, 155.92212, 154.92262, 153.92087, 151.91978, 159.92705, 
            75.9214, 73.92118, 72.92346, 71.92208, 69.92425, 2.014, 1.007825, 4.0026, 3.01603, 173.94005, 179.94655, 178.94582, 
            177.9437, 176.94322, 203.97346, 201.97061, 200.97028, 199.9683, 198.96825, 197.96674, 195.9658, 164.93031, 126.90447, 
            114.90388, 112.90406, 190.96059, 192.96292, 40.961826, 39.964, 38.963707, 77.9204, 85.910614, 83.91151, 82.91414, 
            81.91348, 79.91638, 138.90634, 137.9071, 7.016003, 6.015121, 174.94077, 25.982594, 24.985838, 23.985043, 54.938046, 
            94.90584, 93.90508, 91.90681, 99.90748, 97.9054, 96.90602, 95.90468, 15.000108, 14.003074, 22.989767, 92.90638, 142.9098, 
            141.90771, 149.92088, 147.91689, 145.91312, 144.91257, 143.91008, 21.991383, 20.993843, 19.992435, 61.928345, 60.931057,
            59.930786, 57.935345, 63.927967, 237.0482, 17.99916, 16.99913, 15.994915, 189.95844, 188.95813, 187.95586, 186.95573, 
            185.95383, 183.95248, 191.96147, 30.973763, 206.97588, 205.97444, 203.97302, 207.97662, 109.90517, 107.90389, 105.90348, 
            104.90508, 103.90403, 101.90563, 140.90765, 189.95992, 197.96786, 195.96492, 194.96477, 193.96266, 191.96101, 86.90919, 
            84.9118, 186.95575, 184.95296, 102.9055, 103.905426, 101.90435, 100.90558, 99.90422, 98.90594, 97.90529, 95.9076, 
            35.96708, 33.967865, 32.971455, 31.97207, 122.90421, 120.903824, 44.95591, 77.917305, 76.919914, 75.91921, 73.92248, 
            81.916695, 79.91652, 29.97377, 28.976496, 27.976927, 153.92221, 151.91972, 149.91727, 148.91718, 147.91483, 146.9149, 
            143.912, 123.90527, 121.90344, 119.9022, 118.90331, 117.90161, 116.902954, 115.90175, 114.90335, 113.90279, 111.90482, 
            87.90562, 86.90888, 85.90926, 83.91343, 180.948, 179.94746, 158.92534, 125.90331, 124.904434, 123.902824, 122.904274, 
            121.90305, 119.904045, 129.90623, 127.904465, 232.03806, 46.951763, 45.95263, 49.944794, 48.947872, 47.94795, 204.9744, 
            202.97232, 168.93422, 50.943962, 49.947163, 185.95436, 183.95093, 182.95023, 181.9482, 179.9467, 173.93886, 172.9382, 
            171.93637, 170.93633, 169.93475, 167.9339, 175.94257, 69.92532, 67.92484, 66.92713, 65.92603, 63.929146, 93.90647, 
            91.90504, 90.90565, 89.9047, 95.90827]
    

    abundance = [0.48161000000000004, 0.51839, 1.0, 0.996, 0.00063, 0.00337, 1.0, 1.0, 0.8009999999999999, 0.19899999999999998, 
                0.7170000000000001, 0.11230000000000001, 0.0785, 0.06593, 0.0242, 0.00101, 0.00106, 1.0, 1.0, 0.5069, 
                0.49310000000000004, 0.011000000000000001, 0.9890000000000001, 4e-05, 0.02086, 0.00135, 0.00647, 0.96941, 0.00187, 
                0.128, 0.1249, 0.0089, 0.07490000000000001, 0.2873, 0.1222, 0.2413, 0.11130000000000001, 0.8843000000000001, 0.0025, 
                0.0019, 0.24230000000000002, 0.7576999999999999, 1.0, 0.02365, 0.095, 0.8379000000000001, 0.043449999999999996, 1.0, 
                0.6917, 0.001, 0.0006, 0.282, 0.249, 0.255, 0.18899999999999997, 0.023399999999999997, 0.149, 0.268, 
                0.22949999999999998, 0.336, 0.0161, 0.0014000000000000002, 0.522, 0.478, 1.0, 0.0028000000000000004, 0.021, 
                0.059000000000000004, 0.39892000000000005, 0.60108, 0.2484, 0.1565, 0.2047, 0.14800000000000002, 0.0218, 0.002, 
                0.2186, 0.07440000000000001, 0.3594, 0.07719999999999999, 0.2766, 0.21239999999999998, 0.00015, 0.99985, 1.0, 
                1.37e-08, 0.0016200000000000001, 0.35100000000000003, 0.13629, 0.27297, 0.18606, 0.0687, 0.2986, 0.1318, 0.231, 
                0.16870000000000002, 0.09970000000000001, 0.0015, 1.0, 1.0, 0.9570000000000001, 0.043, 0.373, 0.627, 0.067302, 
                0.000117, 0.932581, 0.0034999999999999996, 0.17300000000000001, 0.57, 0.115, 0.11599999999999999, 0.0225, 0.999088, 
                0.000902, 0.925, 0.075, 0.9741, 0.1101, 0.1, 0.7898999999999999, 1.0, 0.1592, 0.0925, 0.1484, 0.09630000000000001, 
                0.2413, 0.0955, 0.1668, 0.0037, 0.9963, 1.0, 1.0, 0.12179999999999999, 0.2713, 0.0564, 0.0576, 0.17190000000000003,
                0.083, 0.23800000000000002, 0.0925, 0.0027, 0.9048, 0.03634, 0.011399999999999999, 0.26222999999999996, 0.68077, 
                0.009260000000000001, 1.0, 0.002, 0.0004, 0.9976, 0.264, 0.161, 0.133, 0.016, 0.0158, 0.0002, 0.41, 1.0, 0.221, 
                0.24100000000000002, 0.013999999999999999, 0.524, 0.11720000000000001, 0.2646, 0.2733, 0.22329999999999997, 
                0.1114, 0.0102, 1.0, 0.0001, 0.07200000000000001, 0.253, 0.33799999999999997, 0.32899999999999996, 0.0079, 0.2783, 
                0.7217, 0.626, 0.374, 1.0, 0.18600000000000003, 0.316, 0.171, 0.126, 0.127, 0.018600000000000002, 0.0554, 0.0002, 
                0.0421, 0.0075, 0.9501999999999999, 0.4264, 0.5736, 1.0, 0.2377, 0.07629999999999999, 0.09359999999999999, 0.0089, 
                0.0874, 0.4961, 0.031, 0.0467, 0.9223, 0.22699999999999998, 0.267, 0.07400000000000001, 0.138, 0.113, 0.15, 0.031,
                0.0579, 0.0463, 0.3259, 0.0858, 0.2422, 0.0768, 0.14529999999999998, 0.0036, 0.006500000000000001, 0.0097,
                0.8258, 0.07, 0.0986, 0.005600000000000001, 0.9999800000000001, 0.00012, 1.0, 0.1893, 0.0712, 0.0479, 0.00905, 
                0.0259, 0.00095, 0.3387, 0.317, 1.0, 0.073, 0.08, 0.054000000000000006, 0.055, 0.738, 0.7047599999999999, 0.29524, 
                1.0, 0.9975, 0.0025, 0.28600000000000003, 0.307, 0.14279999999999998, 0.263, 0.0012, 0.318, 0.1612, 
                0.21899999999999997, 0.14300000000000002, 0.0305, 0.0013, 0.127, 0.006, 0.188, 0.040999999999999995, 
                0.27899999999999997, 0.486, 0.17379999999999998, 0.17149999999999999, 0.11220000000000001, 0.5145000000000001, 
                0.027999999999999997]
    

    isotopes = ['Ag', 'Ag', 'Al', 'Ar', 'Ar', 'Ar', 'As', 'Au', 'B', 'B', 'Ba', 'Ba', 'Ba', 'Ba', 'Ba', 'Ba', 'Ba', 'Be', 'Bi', 'Br', 
                'Br', 'C', 'C', 'Ca', 'Ca', 'Ca', 'Ca', 'Ca', 'Ca', 'Cd', 'Cd', 'Cd', 'Cd', 'Cd', 'Cd', 'Cd', 'Ce', 'Ce', 'Ce', 'Ce', 
                'Cl', 'Cl', 'Co', 'Cr', 'Cr', 'Cr', 'Cr', 'Cs', 'Cu', 'Dy', 'Dy', 'Dy', 'Dy', 'Dy', 'Dy', 'Dy', 'Er', 'Er', 'Er', 
                'Er', 'Er', 'Er', 'Eu', 'Eu', 'Fe ', 'Fe', 'Fe', 'Fe', 'Ga', 'Ga', 'Gd', 'Gd', 'Gd', 'Gd', 'Gd', 'Gd', 'Gd', 'Ge', 
                'Ge', 'Ge', 'Ge', 'Ge', 'H', 'H', 'He', 'He', 'Hf', 'Hf', 'Hf', 'Hf', 'Hf', 'Hg', 'Hg', 'Hg', 'Hg', 'Hg', 'Hg', 'Hg', 
                'Ho', 'I', 'In', 'In', 'Ir', 'Ir', 'K', 'K', 'K', 'Kr', 'Kr', 'Kr', 'Kr', 'Kr', 'Kr', 'La', 'La', 'Li', 'Li', 'Lu', 
                'Mg', 'Mg', 'Mg', 'Mn', 'Mo', 'Mo', 'Mo', 'Mo', 'Mo', 'Mo', 'Mo', 'N', 'N', 'Na', 'Nb', 'Nd', 'Nd', 'Nd', 'Nd', 'Nd', 
                'Nd', 'Nd', 'Ne', 'Ne', 'Ne', 'Ni', 'Ni', 'Ni', 'Ni', 'Ni', 'Np ', 'O', 'O', 'O', 'Os', 'Os', 'Os', 'Os', 'Os', 'Os', 
                'Os', 'P', 'Pb', 'Pb', 'Pb', 'Pb', 'Pd', 'Pd', 'Pd', 'Pd', 'Pd', 'Pd', 'Pr', 'Pt', 'Pt', 'Pt', 'Pt', 'Pt', 'Pt', 
                'Rb ', 'Rb', 'Re', 'Re', 'Rh', 'Ru', 'Ru', 'Ru', 'Ru', 'Ru', 'Ru', 'Ru', 'S', 'S', 'S', 'S', 'Sb', 'Sb', 'Sc', 'Se', 
                'Se', 'Se', 'Se', 'Se', 'Se', 'Si', 'Si', 'Si', 'Sm', 'Sm', 'Sm', 'Sm', 'Sm', 'Sm', 'Sm', 'Sn', 'Sn', 'Sn', 'Sn', 
                'Sn', 'Sn', 'Sn', 'Sn', 'Sn', 'Sn', 'Sr', 'Sr', 'Sr', 'Sr', 'Ta', 'Ta', 'Tb', 'Te', 'Te', 'Te', 'Te', 'Te', 'Te', 
                'Te', 'Te', 'Th', 'Ti', 'Ti', 'Ti', 'Ti', 'Ti', 'TI', 'TI', 'Tm', 'V', 'V', 'W', 'W', 'W', 'W', 'W', 'Yb', 'Yb', 'Yb', 
                'Yb', 'Yb', 'Yb', 'Yb', 'Zn', 'Zn', 'Zn', 'Zn', 'Zn', 'Zr', 'Zr', 'Zr', 'Zr', 'Zr ']
    

    mol_without_Hs = Chem.MolFromSmiles(mol_smi)


    
    mol = Chem.AddHs(mol_without_Hs)

    '''molecule list generator'''
    list_atoms = []
    for atom in mol.GetAtoms():
        list_atoms.append(atom.GetSymbol())

    #In the case of ionisation by proton, we need to add a H+ ion, which is done in the following

    '''Ionisation method'''
    if 'H' in list_atoms:

        #Check that there is in fact a proton to remove
        list_atoms.remove('H')

    render_imprecise_list = imprecision_True_False 
    #Set arg to be True for long molecules, set arg to False for short molecules/if precision for minuscule peaks is important

    
    has_N = False
    count_N = 0
    has_S = False
    count_S = 0

    '''sulphur_nitrogen_finder function'''
    if 'N' in list_atoms and list_atoms.count('N')%2 == 1:
        has_N = True
        count_N = list_atoms.count('N')
        
    if 'S' in list_atoms and list_atoms.count('S')%2 == 1:
        has_S = True
        count_S = list_atoms.count('S')



    list_output = []
    mass_copy = mass.copy()
    abundance_copy = abundance.copy()
    isotopes_copy = isotopes.copy()

    '''main function'''
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
                    if new_proba>0.00001:
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


    
    x_in, y_in = x_axis_final, y_axis_final
    
    x_out, y_out = [],[]

    '''list_sorter function'''
    while len(x_in)>0:
        min_x = min(x_in)
        index_min = x_in.index(min_x)
        x_out.append(min_x)
        y_out.append(y_in[index_min])
        x_in.pop(index_min)
        y_in.pop(index_min)



    x_in, y_in = x_out, y_out 
    
    '''peak_merger function'''
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

    x_in = x_out
    y_in = y_out
    maximum = max(y_in)

    '''sulphur_nitrogen_plotter function'''
    if has_N:
        x_in.append(x_in[1] - 0.006)  
        y_in.append(0.0035*count_N*maximum)  
        
    
    if has_S:
        x_in.append(x_in[1]-0.004)  
        y_in.append(0.008*count_S*maximum)
        

    x_out, y_out = [],[]

    while len(x_in)>0:
        min_x = min(x_in)
        index_min = x_in.index(min_x)
        x_out.append(min_x)
        y_out.append(y_in[index_min])
        x_in.pop(index_min)
        y_in.pop(index_min)

    x_in, y_in = x_out, y_out 

    min_x , max_x = min(x_in), max(x_in)

    '''delta_function_plotter function'''
    x_axis, y_axis = [min_x-0.2],[0]
    for i in range (len(x_in)):
        x_axis.append(x_in[i]-10**(-100))
        x_axis.append(x_in[i])
        x_axis.append(x_in[i]+10**(-100))
        y_axis.append(0)
        y_axis.append(y_in[i])
        y_axis.append(0)

    x_axis.append(max_x+0.2)
    y_axis.append(0)


    x_in = x_axis
    y_in = y_axis

    '''double_plot function'''
    ticked_peaks = []
    for i in range(len(x_in)):
        if imprecision_True_False:    
            if y_in[i] > 0.001:
                ticked_peaks.append(x_in[i])
        else:
            ticked_peaks.append(x_in[i])

    #creates the principal graph, mass spectrum of the molecule (interactive)

    princ_graph = figure(width=700, title=f'Mass spectrum of molecule')
    princ_graph = figure(x_axis_label = '[m/z]')
    princ_graph = figure(y_axis_label = 'Abundance')
    princ_graph.xaxis.axis_label = "[m/z]"
    princ_graph.yaxis.axis_label = "Abundance"
    princ_graph.title = 'Mass spectrum of molecule'
    princ_graph.height = 500
    princ_graph.xaxis.ticker = FixedTicker(ticks=ticked_peaks)
    princ_graph.add_tools(WheelPanTool(dimension="height"))
    princ_graph.add_tools(WheelZoomTool(dimensions="height"))
    princ_graph.line(x_in, y_in, line_width=1)
    princ_graph.xaxis.major_label_orientation = "horizontal"

     #creates the secondary graph, mass spectrum of the molecule (non-interactive)

    sec_graph = figure(title="Simulated Mass Spectrum", x_axis_label='Mass [Th]', y_axis_label='Intensity')
    sec_graph = figure(width=250, title=f'Mass spectrum of molecule')
    sec_graph = figure(toolbar_location=None)
    sec_graph.height = 250
    sec_graph.line(x_in, y_in, legend_label="Mass spectrum", line_width=1)

    #creates a tool in order that the second graph shows where the zoom is on the first one

    box = BoxAnnotation(left=0, right=0, bottom=0, top=0,
    fill_alpha=0.1, line_color='red', fill_color='cornflowerblue')

    jscode = """
        box[%r] = cb_obj.start
        box[%r] = cb_obj.end
    """

    xcb = CustomJS(args=dict(box=box), code=jscode % ('left', 'right'))
    ycb = CustomJS(args=dict(box=box), code=jscode % ('bottom', 'top'))

    princ_graph.x_range.js_on_change('start', xcb)
    princ_graph.x_range.js_on_change('end', xcb)
    princ_graph.y_range.js_on_change('start', ycb)
    princ_graph.y_range.js_on_change('end', ycb)

    # adds the functionnality to the second figure

    sec_graph.add_layout(box)

    # creates a layout that displays the 2 graphs

    double_graph = column(princ_graph, sec_graph)


    # Generate the image from the molecule
    mol = Chem.MolFromSmiles(mol_smi)

    #Draws the image
    image = Draw.MolToImage(mol)

    #stocks the image in a base64 format
    buffered = BytesIO()
    image.save(buffered, format="PNG")
    image_base64 = base64.b64encode(buffered.getvalue()).decode("utf-8")
    image_url = f"data:image/png;base64,{image_base64}"

    # Create a Div element to display the image
    img_div = Div(text=f'<img src="{image_url}" style="width:350px;height:350px;">')

    # initiate variables
    functional_groups_contained, mol_in = [], Chem.MolFromSmiles(mol_smi)

    # dictionnary of all the considered functional groups to check (some might be missing)

    '''functional_group_finder function'''
    functional_groups_smarts = {
        'Alcohol': 'C[Oh1+0]','Aldehyde': 'C[Ch1]=O','Ketone': 'CC(=O)C','Carboxylic Acid': 'CC(=O)[Oh1]','Ester': 'CC(=O)[Oh0]',
        'Ether': '*[Oh0]*','Amide': 'C(=O)N','Amine': '[C][N]','Nitrile': 'C#N','Chloride': 'Cl','Bromide': 'Br','Fluoride': 'F',
        'Iodide': 'I','Alkene': 'C=C','Alkyne': 'C#C','Imine': 'C=N*','Amino acid': '[Nh2][Ch1*]C(=O)O','Proline': '[Nh1][Ch1*]C(=O)O',
        'Thiol': '[Sh1]','Sulfide': '*[Sh0]*','Acyl Chloride': 'CC(=O)Cl','Anhydride': '*[Ch0](=O)O[Ch0](=O)*','Nitro': 'C[N+](=O)[O-]',
        'Enamine': 'C=C[Nh0]','Enamine2': 'C=C[Nh1]','Enamine3': 'C=C[Nh2]','Imide': 'C(=O)NC(=O)*','Azide': 'CNNN','Enol': 'C=C([Oh1])C',
        'Hemiacetal': 'CC(O)(O)C','Carbonate': '[Oh0]C(=O)[Oh0]','Carbonate2': '[Oh1]C(=O)[Oh1]','Disulfide': 'CSSC','Sulfoxide': 'CS(=O)C',
        'Sulfone': '*[So2](=O)(=O)*','Sulfonic acid': '*S(=O)(=O)[Oh1]','Thioester': 'C(=O)S*','Phosphine': '*[Po0](*)*',
        'Phosphate': '*OP(=O)(O)O','Benzene': 'c1ccccc1','Peroxide':'C[Oh0][Oh0]C'
    }

    # check that the substructure from functional_groups_smarts are contained in mol_smi

    for name, smarts in functional_groups_smarts.items():
        if mol_in.HasSubstructMatch(Chem.MolFromSmarts(smarts)):
            for _ in range(len(mol_in.GetSubstructMatches(Chem.MolFromSmarts(smarts)))):
                functional_groups_contained.append(name)

    # exceptions for conflicts during the iteration of functional groups
    for functional_group in functional_groups_contained:
        if 'Ester' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Ether' in functional_groups_contained:
                    functional_groups_contained.remove('Ether')
        elif functional_group == 'Carboxylic Acid':
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Alcohol' in functional_groups_contained:
                    functional_groups_contained.remove('Alcohol')
        elif 'Phosphate' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Ether' in functional_groups_contained:
                    functional_groups_contained.remove('Ether')
        elif 'Thioester' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Sulfide' in functional_groups_contained:
                    functional_groups_contained.remove('Sulfide')
        elif 'Sulfonic acid' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Sulfide' in functional_groups_contained:
                    functional_groups_contained.remove('Sulfide')
        elif 'Sulfoxide' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Sulfide' in functional_groups_contained:
                    functional_groups_contained.remove('Sulfide')
        elif 'Acyl Chloride' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Sulfide' in functional_groups_contained:
                    functional_groups_contained.remove('Sulfide')
        elif 'Anhydride' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Ester' in functional_groups_contained:
                    functional_groups_contained.remove('Ester')
                if 'Ester' in functional_groups_contained:
                    functional_groups_contained.remove('Ester')
                if 'Ether' in functional_groups_contained:
                    functional_groups_contained.append('Ether')
        elif 'Enamine2' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Enamine2' in functional_groups_contained:
                    functional_groups_contained.remove('Enamine2')
                if 'Enamine' in functional_groups_contained:
                    functional_groups_contained.append('Enamine')
        elif 'Enamine3' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Enamine3' in functional_groups_contained:
                    functional_groups_contained.remove('Enamine3')
                if 'Amine' in functional_groups_contained:
                    functional_groups_contained.remove('Amine')
                functional_groups_contained.append('Enamine')
        elif 'Imide' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Amide' in functional_groups_contained:
                    functional_groups_contained.remove('Amide')
                if 'Amide' in functional_groups_contained:
                    functional_groups_contained.remove('Amide')
        elif 'Enol' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Alkene' in functional_groups_contained:
                    functional_groups_contained.remove('Alkene')
                if 'Alcohol' in functional_groups_contained:
                    functional_groups_contained.remove('Alcohol')
        elif 'Hemiacetal' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Alcohol' in functional_groups_contained:
                    functional_groups_contained.remove('Alcohol')
                if 'Alcohol' in functional_groups_contained:
                    functional_groups_contained.remove('Alcohol')
        elif 'Carbonate2' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Alcohol' in functional_groups_contained:
                    functional_groups_contained.remove('Alcohol')
                if 'Alcohol' in functional_groups_contained:
                    functional_groups_contained.remove('Alcohol')
                if 'Carbonate2' in functional_groups_contained:
                    functional_groups_contained.remove('Carbonate2')
                if 'Carbonate' in functional_groups_contained:
                    functional_groups_contained.append('Carbonate')
        elif 'Disulfide' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Sulfide' in functional_groups_contained:
                    functional_groups_contained.remove('Sulfide')
                if 'Sulfide' in functional_groups_contained:
                    functional_groups_contained.remove('Sulfide')
        elif 'Peroxide' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Ether' in functional_groups_contained:
                    functional_groups_contained.remove('Ether')
                if 'Ether' in functional_groups_contained:
                    functional_groups_contained.remove('Ether')
        elif 'Amide' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Amine' in functional_groups_contained:
                    functional_groups_contained.remove('Amine')
                if 'Amine' in functional_groups_contained:
                    functional_groups_contained.remove('Amine')

    #initiate empty variables
    present_group_smarts = []    
    present_group_images_base64 = []
    
    #appends the smarts of the contained functional groups
    for i,j in functional_groups_smarts.items():
        for x in functional_groups_contained:
            if x == i:
                present_group_smarts.append(j)

    #converts the smarts to images in base64 format
    for x in present_group_smarts:

        #converts SMARTs to SMILEs for the images to be nicer
        mol_x = Chem.MolFromSmarts(x)
        mol_smi = Chem.MolToSmiles(mol_x)
        mol = Chem.MolFromSmiles(mol_smi)

        if mol:

            if x == 'CC(=O)[Oh0]':
                mol2 = Chem.MolFromSmiles('CC(=O)OC')

                #creates the image
                image = Draw.MolToImage(mol2)

                # Convert the image to base64 format
                buffered = BytesIO()
                image.save(buffered, format="PNG")
                image_base64 = base64.b64encode(buffered.getvalue()).decode("utf-8")

                #appends the image to an empty dictionnary
                present_group_images_base64.append(image_base64)
            elif x == '[N]':
                mol2 = Chem.MolFromSmiles('CN')

                #creates the image
                image = Draw.MolToImage(mol2)

                # Convert the image to base64 format
                buffered = BytesIO()
                image.save(buffered, format="PNG")
                image_base64 = base64.b64encode(buffered.getvalue()).decode("utf-8")

                #appends the image to an empty dictionnary
                present_group_images_base64.append(image_base64)
            else:
                #creates the image
                image = Draw.MolToImage(mol)

                # Convert the image to base64 format
                buffered = BytesIO()
                image.save(buffered, format="PNG")
                image_base64 = base64.b64encode(buffered.getvalue()).decode("utf-8")

                #appends the image to an empty dictionnary
                present_group_images_base64.append(image_base64)

    #creates a dictionnary that links the name of the functional group to its image
    data = dict(groups=functional_groups_contained,images=present_group_images_base64)
    source = ColumnDataSource(data)

    #template for the bokeh table that read the base64 format 
    template = """
    <div>
        <img src="data:image/png;base64, <%= value %>" style="width:50px;height:50px;">
    </div>
    """

    # initiallizing the bokeh figure using the previous template for each functional group
    columns = [TableColumn(field="groups", title="Functional Groups"),
    TableColumn(field="images", title="Images", width=200, formatter=HTMLTemplateFormatter(template=template))]
    num_groups = len(functional_groups_contained)

    table_height = min(200 + num_groups * 60, 800)

    func_group_table = DataTable(source=source, columns=columns, width=250, height=table_height, row_height=60)

    last = column(img_div, func_group_table)
    final = row(double_graph, last)
    return final

#3
def spectrum_3D(mol_smi, imprecision_True_False, apparatus_resolution):

    from xyz2graph import MolGraph, to_plotly_figure

    #lists of the data to facilitise the pip-installability of the package

    mol_smi_3D = mol_smi

    mass = [108.904755, 106.90509, 26.981539, 39.962383, 37.96273, 35.967545, 74.92159, 196.96654, 11.009305, 10.012937, 137.90523, 
            136.9058, 135.90456, 134.90567, 133.90448, 131.90504, 129.90628, 9.012182, 208.98038, 78.918335, 80.91629, 13.003355, 
            12.0, 45.95369, 43.95548, 42.958767, 41.95862, 39.96259, 47.952534, 110.90418, 109.90301, 107.90418, 115.904755, 
            113.90336, 112.9044, 111.902756, 141.90924, 139.90543, 137.90599, 135.90714, 36.965904, 34.968853, 58.933197, 53.93888, 
            52.94065, 51.94051, 49.946045, 132.90543, 62.939598, 157.92441, 155.92528, 163.92917, 162.92873, 161.92679, 160.92693, 
            159.92519, 169.93546, 167.93237, 166.93205, 165.93028, 163.9292, 161.92877, 152.92122, 150.91985, 18.998404, 57.933277, 
            56.935394, 53.939613, 70.9247, 68.92558, 157.9241, 156.92395, 155.92212, 154.92262, 153.92087, 151.91978, 159.92705, 
            75.9214, 73.92118, 72.92346, 71.92208, 69.92425, 2.014, 1.007825, 4.0026, 3.01603, 173.94005, 179.94655, 178.94582, 
            177.9437, 176.94322, 203.97346, 201.97061, 200.97028, 199.9683, 198.96825, 197.96674, 195.9658, 164.93031, 126.90447, 
            114.90388, 112.90406, 190.96059, 192.96292, 40.961826, 39.964, 38.963707, 77.9204, 85.910614, 83.91151, 82.91414, 
            81.91348, 79.91638, 138.90634, 137.9071, 7.016003, 6.015121, 174.94077, 25.982594, 24.985838, 23.985043, 54.938046, 
            94.90584, 93.90508, 91.90681, 99.90748, 97.9054, 96.90602, 95.90468, 15.000108, 14.003074, 22.989767, 92.90638, 142.9098, 
            141.90771, 149.92088, 147.91689, 145.91312, 144.91257, 143.91008, 21.991383, 20.993843, 19.992435, 61.928345, 60.931057,
            59.930786, 57.935345, 63.927967, 237.0482, 17.99916, 16.99913, 15.994915, 189.95844, 188.95813, 187.95586, 186.95573, 
            185.95383, 183.95248, 191.96147, 30.973763, 206.97588, 205.97444, 203.97302, 207.97662, 109.90517, 107.90389, 105.90348, 
            104.90508, 103.90403, 101.90563, 140.90765, 189.95992, 197.96786, 195.96492, 194.96477, 193.96266, 191.96101, 86.90919, 
            84.9118, 186.95575, 184.95296, 102.9055, 103.905426, 101.90435, 100.90558, 99.90422, 98.90594, 97.90529, 95.9076, 
            35.96708, 33.967865, 32.971455, 31.97207, 122.90421, 120.903824, 44.95591, 77.917305, 76.919914, 75.91921, 73.92248, 
            81.916695, 79.91652, 29.97377, 28.976496, 27.976927, 153.92221, 151.91972, 149.91727, 148.91718, 147.91483, 146.9149, 
            143.912, 123.90527, 121.90344, 119.9022, 118.90331, 117.90161, 116.902954, 115.90175, 114.90335, 113.90279, 111.90482, 
            87.90562, 86.90888, 85.90926, 83.91343, 180.948, 179.94746, 158.92534, 125.90331, 124.904434, 123.902824, 122.904274, 
            121.90305, 119.904045, 129.90623, 127.904465, 232.03806, 46.951763, 45.95263, 49.944794, 48.947872, 47.94795, 204.9744, 
            202.97232, 168.93422, 50.943962, 49.947163, 185.95436, 183.95093, 182.95023, 181.9482, 179.9467, 173.93886, 172.9382, 
            171.93637, 170.93633, 169.93475, 167.9339, 175.94257, 69.92532, 67.92484, 66.92713, 65.92603, 63.929146, 93.90647, 
            91.90504, 90.90565, 89.9047, 95.90827]
    

    abundance = [0.48161000000000004, 0.51839, 1.0, 0.996, 0.00063, 0.00337, 1.0, 1.0, 0.8009999999999999, 0.19899999999999998, 
                0.7170000000000001, 0.11230000000000001, 0.0785, 0.06593, 0.0242, 0.00101, 0.00106, 1.0, 1.0, 0.5069, 
                0.49310000000000004, 0.011000000000000001, 0.9890000000000001, 4e-05, 0.02086, 0.00135, 0.00647, 0.96941, 0.00187, 
                0.128, 0.1249, 0.0089, 0.07490000000000001, 0.2873, 0.1222, 0.2413, 0.11130000000000001, 0.8843000000000001, 0.0025, 
                0.0019, 0.24230000000000002, 0.7576999999999999, 1.0, 0.02365, 0.095, 0.8379000000000001, 0.043449999999999996, 1.0, 
                0.6917, 0.001, 0.0006, 0.282, 0.249, 0.255, 0.18899999999999997, 0.023399999999999997, 0.149, 0.268, 
                0.22949999999999998, 0.336, 0.0161, 0.0014000000000000002, 0.522, 0.478, 1.0, 0.0028000000000000004, 0.021, 
                0.059000000000000004, 0.39892000000000005, 0.60108, 0.2484, 0.1565, 0.2047, 0.14800000000000002, 0.0218, 0.002, 
                0.2186, 0.07440000000000001, 0.3594, 0.07719999999999999, 0.2766, 0.21239999999999998, 0.00015, 0.99985, 1.0, 
                1.37e-08, 0.0016200000000000001, 0.35100000000000003, 0.13629, 0.27297, 0.18606, 0.0687, 0.2986, 0.1318, 0.231, 
                0.16870000000000002, 0.09970000000000001, 0.0015, 1.0, 1.0, 0.9570000000000001, 0.043, 0.373, 0.627, 0.067302, 
                0.000117, 0.932581, 0.0034999999999999996, 0.17300000000000001, 0.57, 0.115, 0.11599999999999999, 0.0225, 0.999088, 
                0.000902, 0.925, 0.075, 0.9741, 0.1101, 0.1, 0.7898999999999999, 1.0, 0.1592, 0.0925, 0.1484, 0.09630000000000001, 
                0.2413, 0.0955, 0.1668, 0.0037, 0.9963, 1.0, 1.0, 0.12179999999999999, 0.2713, 0.0564, 0.0576, 0.17190000000000003,
                0.083, 0.23800000000000002, 0.0925, 0.0027, 0.9048, 0.03634, 0.011399999999999999, 0.26222999999999996, 0.68077, 
                0.009260000000000001, 1.0, 0.002, 0.0004, 0.9976, 0.264, 0.161, 0.133, 0.016, 0.0158, 0.0002, 0.41, 1.0, 0.221, 
                0.24100000000000002, 0.013999999999999999, 0.524, 0.11720000000000001, 0.2646, 0.2733, 0.22329999999999997, 
                0.1114, 0.0102, 1.0, 0.0001, 0.07200000000000001, 0.253, 0.33799999999999997, 0.32899999999999996, 0.0079, 0.2783, 
                0.7217, 0.626, 0.374, 1.0, 0.18600000000000003, 0.316, 0.171, 0.126, 0.127, 0.018600000000000002, 0.0554, 0.0002, 
                0.0421, 0.0075, 0.9501999999999999, 0.4264, 0.5736, 1.0, 0.2377, 0.07629999999999999, 0.09359999999999999, 0.0089, 
                0.0874, 0.4961, 0.031, 0.0467, 0.9223, 0.22699999999999998, 0.267, 0.07400000000000001, 0.138, 0.113, 0.15, 0.031,
                0.0579, 0.0463, 0.3259, 0.0858, 0.2422, 0.0768, 0.14529999999999998, 0.0036, 0.006500000000000001, 0.0097,
                0.8258, 0.07, 0.0986, 0.005600000000000001, 0.9999800000000001, 0.00012, 1.0, 0.1893, 0.0712, 0.0479, 0.00905, 
                0.0259, 0.00095, 0.3387, 0.317, 1.0, 0.073, 0.08, 0.054000000000000006, 0.055, 0.738, 0.7047599999999999, 0.29524, 
                1.0, 0.9975, 0.0025, 0.28600000000000003, 0.307, 0.14279999999999998, 0.263, 0.0012, 0.318, 0.1612, 
                0.21899999999999997, 0.14300000000000002, 0.0305, 0.0013, 0.127, 0.006, 0.188, 0.040999999999999995, 
                0.27899999999999997, 0.486, 0.17379999999999998, 0.17149999999999999, 0.11220000000000001, 0.5145000000000001, 
                0.027999999999999997]
    

    isotopes = ['Ag', 'Ag', 'Al', 'Ar', 'Ar', 'Ar', 'As', 'Au', 'B', 'B', 'Ba', 'Ba', 'Ba', 'Ba', 'Ba', 'Ba', 'Ba', 'Be', 'Bi', 'Br', 
                'Br', 'C', 'C', 'Ca', 'Ca', 'Ca', 'Ca', 'Ca', 'Ca', 'Cd', 'Cd', 'Cd', 'Cd', 'Cd', 'Cd', 'Cd', 'Ce', 'Ce', 'Ce', 'Ce', 
                'Cl', 'Cl', 'Co', 'Cr', 'Cr', 'Cr', 'Cr', 'Cs', 'Cu', 'Dy', 'Dy', 'Dy', 'Dy', 'Dy', 'Dy', 'Dy', 'Er', 'Er', 'Er', 
                'Er', 'Er', 'Er', 'Eu', 'Eu', 'Fe ', 'Fe', 'Fe', 'Fe', 'Ga', 'Ga', 'Gd', 'Gd', 'Gd', 'Gd', 'Gd', 'Gd', 'Gd', 'Ge', 
                'Ge', 'Ge', 'Ge', 'Ge', 'H', 'H', 'He', 'He', 'Hf', 'Hf', 'Hf', 'Hf', 'Hf', 'Hg', 'Hg', 'Hg', 'Hg', 'Hg', 'Hg', 'Hg', 
                'Ho', 'I', 'In', 'In', 'Ir', 'Ir', 'K', 'K', 'K', 'Kr', 'Kr', 'Kr', 'Kr', 'Kr', 'Kr', 'La', 'La', 'Li', 'Li', 'Lu', 
                'Mg', 'Mg', 'Mg', 'Mn', 'Mo', 'Mo', 'Mo', 'Mo', 'Mo', 'Mo', 'Mo', 'N', 'N', 'Na', 'Nb', 'Nd', 'Nd', 'Nd', 'Nd', 'Nd', 
                'Nd', 'Nd', 'Ne', 'Ne', 'Ne', 'Ni', 'Ni', 'Ni', 'Ni', 'Ni', 'Np ', 'O', 'O', 'O', 'Os', 'Os', 'Os', 'Os', 'Os', 'Os', 
                'Os', 'P', 'Pb', 'Pb', 'Pb', 'Pb', 'Pd', 'Pd', 'Pd', 'Pd', 'Pd', 'Pd', 'Pr', 'Pt', 'Pt', 'Pt', 'Pt', 'Pt', 'Pt', 
                'Rb ', 'Rb', 'Re', 'Re', 'Rh', 'Ru', 'Ru', 'Ru', 'Ru', 'Ru', 'Ru', 'Ru', 'S', 'S', 'S', 'S', 'Sb', 'Sb', 'Sc', 'Se', 
                'Se', 'Se', 'Se', 'Se', 'Se', 'Si', 'Si', 'Si', 'Sm', 'Sm', 'Sm', 'Sm', 'Sm', 'Sm', 'Sm', 'Sn', 'Sn', 'Sn', 'Sn', 
                'Sn', 'Sn', 'Sn', 'Sn', 'Sn', 'Sn', 'Sr', 'Sr', 'Sr', 'Sr', 'Ta', 'Ta', 'Tb', 'Te', 'Te', 'Te', 'Te', 'Te', 'Te', 
                'Te', 'Te', 'Th', 'Ti', 'Ti', 'Ti', 'Ti', 'Ti', 'TI', 'TI', 'Tm', 'V', 'V', 'W', 'W', 'W', 'W', 'W', 'Yb', 'Yb', 'Yb', 
                'Yb', 'Yb', 'Yb', 'Yb', 'Zn', 'Zn', 'Zn', 'Zn', 'Zn', 'Zr', 'Zr', 'Zr', 'Zr', 'Zr ']
    

    mol_without_Hs = Chem.MolFromSmiles(mol_smi)

    if mol_without_Hs is None:
        raise ValueError('\nInvalid SMILEs enterred.\nPlease enter a different SMILEs.')
    
    mol = Chem.AddHs(mol_without_Hs)

    list_atoms = []
    for atom in mol.GetAtoms():
        list_atoms.append(atom.GetSymbol())
    
    if not list_atoms:
        return None
    
    #In the case of ionisation by proton, we need to add a H+ ion, which is done in the following

    if 'H' in list_atoms:

        #Check that there is in fact a proton to remove
        list_atoms.remove('H')

    render_imprecise_list = imprecision_True_False 
    #Set arg to be True for long molecules, set arg to False for short molecules/if precision for minuscule peaks is important

    #check for sulphur and nitrogen

    has_N = False
    count_N = 0
    has_S = False
    count_S = 0
    if 'N' in list_atoms and list_atoms.count('N')%2 == 1:
        has_N = True
        count_N = list_atoms.count('N')
        
    if 'S' in list_atoms and list_atoms.count('S')%2 == 1:
        has_S = True
        count_S = list_atoms.count('S')



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


                #removes any molecule who's probability is below 0.00001

                if render_imprecise_list: #only removes low-probability arrangements if render_imprecise_list arg is True
                    if new_proba>0.00001:
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


    
    x_in, y_in = x_axis_final, y_axis_final
    
    x_out, y_out = [],[]

    while len(x_in)>0:
        min_x = min(x_in)
        index_min = x_in.index(min_x)
        x_out.append(min_x)
        y_out.append(y_in[index_min])
        x_in.pop(index_min)
        y_in.pop(index_min)



    x_in, y_in = x_out, y_out 
#resolution
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


    x_in = x_out
    y_in = y_out
    maximum = max(y_in)

    if has_N:
        x_in.append(x_in[1] - 0.006)  
        y_in.append(0.0035*count_N*maximum)  
        
    
    if has_S:
        x_in.append(x_in[1]-0.004)  
        y_in.append(0.008*count_S*maximum)
        

    x_out, y_out = [],[]

    while len(x_in)>0:
        min_x = min(x_in)
        index_min = x_in.index(min_x)
        x_out.append(min_x)
        y_out.append(y_in[index_min])
        x_in.pop(index_min)
        y_in.pop(index_min)

    x_in, y_in = x_out, y_out 

    min_x , max_x = min(x_in), max(x_in)

    x_axis, y_axis = [min_x-0.2],[0]
    for i in range (len(x_in)):
        x_axis.append(x_in[i]-10**(-100))
        x_axis.append(x_in[i])
        x_axis.append(x_in[i]+10**(-100))
        y_axis.append(0)
        y_axis.append(y_in[i])
        y_axis.append(0)

    x_axis.append(max_x+0.2)
    y_axis.append(0)


    x_in = x_axis
    y_in = y_axis

    ticked_peaks = []
    for i in range(len(x_in)):
        if imprecision_True_False:    
            if y_in[i] > 0.001:
                ticked_peaks.append(x_in[i])
        else:
            ticked_peaks.append(x_in[i])

    #creates the principal graph, mass spectrum of the molecule (interactive)

    princ_graph = figure(width=700, title=f'Mass spectrum of molecule')
    princ_graph = figure(x_axis_label = '[m/z]')
    princ_graph = figure(y_axis_label = 'Abundance')
    princ_graph.xaxis.axis_label = "[m/z]"
    princ_graph.title = 'Mass spectrum of molecule'
    princ_graph.height = 500
    princ_graph.xaxis.ticker = FixedTicker(ticks=ticked_peaks)
    princ_graph.add_tools(WheelPanTool(dimension="height"))
    princ_graph.add_tools(WheelZoomTool(dimensions="height"))
    princ_graph.line(x_in, y_in, line_width=1)
    princ_graph.xaxis.major_label_orientation = "horizontal"

     #creates the secondary graph, mass spectrum of the molecule (non-interactive)

    sec_graph = figure(title="Simulated Mass Spectrum", x_axis_label='Mass [Th]', y_axis_label='Intensity')
    sec_graph = figure(width=250, title=f'Mass spectrum of molecule')
    sec_graph = figure(toolbar_location=None)
    sec_graph.height = 250
    sec_graph.line(x_in, y_in, legend_label="Mass spectrum", line_width=1)

    #creates a tool in order that the second graph shows where the zoom is on the first one

    box = BoxAnnotation(left=0, right=0, bottom=0, top=0,
    fill_alpha=0.1, line_color='red', fill_color='cornflowerblue')

    jscode = """
        box[%r] = cb_obj.start
        box[%r] = cb_obj.end
    """

    xcb = CustomJS(args=dict(box=box), code=jscode % ('left', 'right'))
    ycb = CustomJS(args=dict(box=box), code=jscode % ('bottom', 'top'))

    princ_graph.x_range.js_on_change('start', xcb)
    princ_graph.x_range.js_on_change('end', xcb)
    princ_graph.y_range.js_on_change('start', ycb)
    princ_graph.y_range.js_on_change('end', ycb)

    # adds the functionnality to the second figure

    sec_graph.add_layout(box)

    # creates a layout that displays the 2 graphs

    double_graph = column(princ_graph, sec_graph)


    # Generate the image from the molecule
    mol = Chem.MolFromSmiles(mol_smi)

    #Draws the image
    image = Draw.MolToImage(mol)

    #stocks the image in a base64 format
    buffered = BytesIO()
    image.save(buffered, format="PNG")
    image_base64 = base64.b64encode(buffered.getvalue()).decode("utf-8")
    image_url = f"data:image/png;base64,{image_base64}"

    # Create a Div element to display the image
    img_div = Div(text=f'<img src="{image_url}" style="width:350px;height:350px;">')

    # initiate variables
    functional_groups_contained, mol_in = [], Chem.MolFromSmiles(mol_smi)

    # dictionnary of all the considered functional groups to check (some might be missing)

    functional_groups_smarts = {
        'Alcohol': 'C[Oh1+0]','Aldehyde': 'C[Ch1]=O','Ketone': 'CC(=O)C','Carboxylic Acid': '*C(=O)[Oh1]','Ester': 'CC(=O)[Oh0]',
        'Ether': '*[Oh0]*','Amide': 'C(=O)N','Amine': '[C][N]','Nitrile': 'C#N','Chloride': 'Cl','Bromide': 'Br','Fluoride': 'F',
        'Iodide': 'I','Alkene': 'C=C','Alkyne': 'C#C','Imine': 'C=N*','Amino acid': '[Nh2][Ch1*]C(=O)O','Proline': '[Nh1][Ch1*]C(=O)O',
        'Thiol': '[Sh1]','Sulfide': '*[Sh0]*','Acyl Chloride': 'CC(=O)Cl','Anhydride': '*[Ch0](=O)O[Ch0](=O)*','Nitro': 'C[N+](=O)[O-]',
        'Enamine': 'C=C[Nh0]','Enamine2': 'C=C[Nh1]','Enamine3': 'C=C[Nh2]','Imide': 'C(=O)NC(=O)*','Azide': 'CNNN','Enol': 'C=C([Oh1])C',
        'Hemiacetal': 'CC(O)(O)C','Carbonate': '[Oh0]C(=O)[Oh0]','Carbonate2': '[Oh1]C(=O)[Oh1]','Disulfide': 'CSSC','Sulfoxide': 'CS(=O)C',
        'Sulfone': '*[So2](=O)(=O)*','Sulfonic acid': '*S(=O)(=O)[Oh1]','Thioester': 'C(=O)S*','Phosphine': '*[Po0](*)*','Phosphate': '*OP(=O)(O)O',
        'Benzene': 'c1ccccc1','Peroxide':'C[Oh0][Oh0]C'
    }

    # check that the substructure from functional_groups_smarts are contained in mol_smi

    for name, smarts in functional_groups_smarts.items():
        if mol_in.HasSubstructMatch(Chem.MolFromSmarts(smarts)):
            for _ in range(len(mol_in.GetSubstructMatches(Chem.MolFromSmarts(smarts)))):
                functional_groups_contained.append(name)

    # exceptions for conflicts during the iteration of functional groups
    for functional_group in functional_groups_contained:
        if 'Ester' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Ether' in functional_groups_contained:
                    functional_groups_contained.remove('Ether')
        elif functional_group == 'Carboxylic Acid':
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Alcohol' in functional_groups_contained:
                    functional_groups_contained.remove('Alcohol')
        elif 'Phosphate' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Ether' in functional_groups_contained:
                    functional_groups_contained.remove('Ether')
        elif 'Thioester' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Sulfide' in functional_groups_contained:
                    functional_groups_contained.remove('Sulfide')
        elif 'Sulfonic acid' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Sulfide' in functional_groups_contained:
                    functional_groups_contained.remove('Sulfide')
        elif 'Sulfoxide' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Sulfide' in functional_groups_contained:
                    functional_groups_contained.remove('Sulfide')
        elif 'Acyl Chloride' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Sulfide' in functional_groups_contained:
                    functional_groups_contained.remove('Sulfide')
        elif 'Anhydride' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Ester' in functional_groups_contained:
                    functional_groups_contained.remove('Ester')
                if 'Ester' in functional_groups_contained:
                    functional_groups_contained.remove('Ester')
                if 'Ether' in functional_groups_contained:
                    functional_groups_contained.append('Ether')
        elif 'Enamine2' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Enamine2' in functional_groups_contained:
                    functional_groups_contained.remove('Enamine2')
                if 'Enamine' in functional_groups_contained:
                    functional_groups_contained.append('Enamine')
        elif 'Enamine3' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Enamine3' in functional_groups_contained:
                    functional_groups_contained.remove('Enamine3')
                if 'Amine' in functional_groups_contained:
                    functional_groups_contained.remove('Amine')
                functional_groups_contained.append('Enamine')
        elif 'Imide' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Amide' in functional_groups_contained:
                    functional_groups_contained.remove('Amide')
                if 'Amide' in functional_groups_contained:
                    functional_groups_contained.remove('Amide')
        elif 'Enol' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Alkene' in functional_groups_contained:
                    functional_groups_contained.remove('Alkene')
                if 'Alcohol' in functional_groups_contained:
                    functional_groups_contained.remove('Alcohol')
        elif 'Hemiacetal' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Alcohol' in functional_groups_contained:
                    functional_groups_contained.remove('Alcohol')
                if 'Alcohol' in functional_groups_contained:
                    functional_groups_contained.remove('Alcohol')
        elif 'Carbonate2' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Alcohol' in functional_groups_contained:
                    functional_groups_contained.remove('Alcohol')
                if 'Alcohol' in functional_groups_contained:
                    functional_groups_contained.remove('Alcohol')
                if 'Carbonate2' in functional_groups_contained:
                    functional_groups_contained.remove('Carbonate2')
                if 'Carbonate' in functional_groups_contained:
                    functional_groups_contained.append('Carbonate')
        elif 'Disulfide' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Sulfide' in functional_groups_contained:
                    functional_groups_contained.remove('Sulfide')
                if 'Sulfide' in functional_groups_contained:
                    functional_groups_contained.remove('Sulfide')
        elif 'Peroxide' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Ether' in functional_groups_contained:
                    functional_groups_contained.remove('Ether')
                if 'Ether' in functional_groups_contained:
                    functional_groups_contained.remove('Ether')
        elif 'Amide' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Amine' in functional_groups_contained:
                    functional_groups_contained.remove('Amine')
                if 'Amine' in functional_groups_contained:
                    functional_groups_contained.remove('Amine')

    #initiate empty variables
    present_group_smarts = []    
    present_group_images_base64 = []
    
    #appends the smarts of the contained functional groups
    for i,j in functional_groups_smarts.items():
        for x in functional_groups_contained:
            if x == i:
                present_group_smarts.append(j)

    #converts the smarts to images in base64 format
    for x in present_group_smarts:

        #converts SMARTs to SMILEs for the images to be nicer
        mol_x = Chem.MolFromSmarts(x)
        mol_smi = Chem.MolToSmiles(mol_x)
        mol = Chem.MolFromSmiles(mol_smi)

        if mol:

            if x == 'CC(=O)[Oh0]':
                mol2 = Chem.MolFromSmiles('CC(=O)OC')

                #creates the image
                image = Draw.MolToImage(mol2)

                # Convert the image to base64 format
                buffered = BytesIO()
                image.save(buffered, format="PNG")
                image_base64 = base64.b64encode(buffered.getvalue()).decode("utf-8")

                #appends the image to an empty dictionnary
                present_group_images_base64.append(image_base64)
            else:
                #creates the image
                image = Draw.MolToImage(mol)

                # Convert the image to base64 format
                buffered = BytesIO()
                image.save(buffered, format="PNG")
                image_base64 = base64.b64encode(buffered.getvalue()).decode("utf-8")

                #appends the image to an empty dictionnary
                present_group_images_base64.append(image_base64)

    #creates a dictionnary that links the name of the functional group to its image
    data = dict(
        groups=functional_groups_contained,
        images=present_group_images_base64
    )
    source = ColumnDataSource(data)

    #template for the bokeh table that read the base64 format 
    template = """
    <div>
        <img src="data:image/png;base64, <%= value %>" style="width:50px;height:50px;">
    </div>
    """

    # initiallizing the bokeh figure using the previous template for each functional group
    columns = [
        TableColumn(field="groups", title="Functional Groups"),
        TableColumn(field="images", title="Images", width=200, formatter=HTMLTemplateFormatter(template=template))
    ]
    num_groups = len(functional_groups_contained)

    table_height = min(200 + num_groups * 60, 800)

    func_group_table = DataTable(source=source, columns=columns, width=250, height=table_height, row_height=60)

    last = row(img_div, func_group_table)

    # Generate 3D coordinates from SMILES string
    mol_3D = Chem.MolFromSmiles(mol_smi_3D)
    mol_3D = Chem.AddHs(mol_3D)  # Add hydrogens for better geometry optimization
    AllChem.EmbedMolecule(mol_3D, randomSeed=42)  # Embed the molecule in 3D space
    AllChem.MMFFOptimizeMolecule(mol_3D)  # Optimize the geometry using MMFF94 force field

    # Create a temporary file to store the 
    with tempfile.NamedTemporaryFile(delete=False) as tmp:
        tmp.write(f"{mol_3D.GetNumAtoms()}\n\n".encode('utf-8'))  # Write number of atoms
        for atom in mol_3D.GetAtoms():
            pos = mol_3D.GetConformer().GetAtomPosition(atom.GetIdx())
            tmp.write(f"{atom.GetSymbol()} {pos.x} {pos.y} {pos.z}\n".encode('utf-8'))
        tmp_path = tmp.name

    # Create the MolGraph object
    mg = MolGraph()
    
    # Read the data from the temporary file
    mg.read_xyz(tmp_path)

    # Create the Plotly figure object
    fig = to_plotly_figure(mg)


    # Converts all the graphs to pane to combine Bokeh and Plotly
    fig_pane = pn.pane.Plotly(fig)

    last_pane = pn.pane.Bokeh(last)

    left_pane = pn.pane.Bokeh(double_graph)

    # Creates the final plot
    right_pane = pn.Column(fig_pane, last_pane)

    total_plot_pane = pn.Row(left_pane, right_pane)

    return total_plot_pane


#functions that compose the main functions:

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
        raise ValueError('Invalid SMILEs representation')

    mol = Chem.AddHs(mol_without_Hs)

    return mol

def molecule_list_generator(mol) -> list[str]:
    #---------------------------------------------------------------------------------------------#
    '''
    molecule_list_generator(mol)
    
    Input: molecule under MOL representation
    
    Output: list containing the atomic symbol of each atom in the input molecule
    '''
    #---------------------------------------------------------------------------------------------#
    if not mol:
        raise ValueError('Enter a non-empty input')
    list_atoms = []
    for atom in mol.GetAtoms():
        list_atoms.append(atom.GetSymbol())
    return list_atoms

def ionisation_method (list_atoms) -> list[str]:
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

def sulphur_nitrogen_finder(list_atoms):
    #---------------------------------------------------------------------------------------------#
    '''
    sulpher_nitrogen_finder(list_atoms)

    Input: list of atoms in the molecule

    Output: 2 booleans and 2 integers:
    1. boolean that tells if the molecule contains an odd number of Nitrogen atoms
    2. boolean that tells if the molecule contains an odd number of Sulphur atoms
    3. number of Nitrogen atoms in the molecule (in case where there is an odd number)
    4. number of Sulphur atoms in the molecule (in case where there is an odd number)
    '''
    #---------------------------------------------------------------------------------------------#
    if not list_atoms:
        raise ValueError("Enter a non-empty list")
    if type(list_atoms) != list:
        raise ValueError("Enter a list as argument")
    

    count_N = 0
    has_N = False
    if 'N' in list_atoms and list_atoms.count('N')%2 == 1:
        has_N = True
        count_N = list_atoms.count('N')
    count_S = 0
    has_S = False
    if 'S' in list_atoms and list_atoms.count('S')%2 == 1:
        has_S = True
        count_S = list_atoms.count('S')
    return has_N, has_S, count_N, count_S

def main_function(list_atoms, imprecision_True_False) -> list[float]:
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

    if not list_atoms:
        raise ValueError('Enter a valid input for the first argument')
    if imprecision_True_False not in [True, False]:
        raise ValueError('Enter a valid input for the second argument')
    
    render_imprecise_list = imprecision_True_False 
    #Set arg to be True for long molecules, set arg to False for short molecules/if precision for minuscule peaks is important
    
    mass = [108.904755, 106.90509, 26.981539, 39.962383, 37.96273, 35.967545, 74.92159, 196.96654, 11.009305, 10.012937, 137.90523, 
            136.9058, 135.90456, 134.90567, 133.90448, 131.90504, 129.90628, 9.012182, 208.98038, 78.918335, 80.91629, 13.003355, 
            12.0, 45.95369, 43.95548, 42.958767, 41.95862, 39.96259, 47.952534, 110.90418, 109.90301, 107.90418, 115.904755, 
            113.90336, 112.9044, 111.902756, 141.90924, 139.90543, 137.90599, 135.90714, 36.965904, 34.968853, 58.933197, 53.93888, 
            52.94065, 51.94051, 49.946045, 132.90543, 62.939598, 157.92441, 155.92528, 163.92917, 162.92873, 161.92679, 160.92693, 
            159.92519, 169.93546, 167.93237, 166.93205, 165.93028, 163.9292, 161.92877, 152.92122, 150.91985, 18.998404, 57.933277, 
            56.935394, 53.939613, 70.9247, 68.92558, 157.9241, 156.92395, 155.92212, 154.92262, 153.92087, 151.91978, 159.92705, 
            75.9214, 73.92118, 72.92346, 71.92208, 69.92425, 2.014, 1.007825, 4.0026, 3.01603, 173.94005, 179.94655, 178.94582, 
            177.9437, 176.94322, 203.97346, 201.97061, 200.97028, 199.9683, 198.96825, 197.96674, 195.9658, 164.93031, 126.90447, 
            114.90388, 112.90406, 190.96059, 192.96292, 40.961826, 39.964, 38.963707, 77.9204, 85.910614, 83.91151, 82.91414, 
            81.91348, 79.91638, 138.90634, 137.9071, 7.016003, 6.015121, 174.94077, 25.982594, 24.985838, 23.985043, 54.938046, 
            94.90584, 93.90508, 91.90681, 99.90748, 97.9054, 96.90602, 95.90468, 15.000108, 14.003074, 22.989767, 92.90638, 142.9098, 
            141.90771, 149.92088, 147.91689, 145.91312, 144.91257, 143.91008, 21.991383, 20.993843, 19.992435, 61.928345, 60.931057,
            59.930786, 57.935345, 63.927967, 237.0482, 17.99916, 16.99913, 15.994915, 189.95844, 188.95813, 187.95586, 186.95573, 
            185.95383, 183.95248, 191.96147, 30.973763, 206.97588, 205.97444, 203.97302, 207.97662, 109.90517, 107.90389, 105.90348, 
            104.90508, 103.90403, 101.90563, 140.90765, 189.95992, 197.96786, 195.96492, 194.96477, 193.96266, 191.96101, 86.90919, 
            84.9118, 186.95575, 184.95296, 102.9055, 103.905426, 101.90435, 100.90558, 99.90422, 98.90594, 97.90529, 95.9076, 
            35.96708, 33.967865, 32.971455, 31.97207, 122.90421, 120.903824, 44.95591, 77.917305, 76.919914, 75.91921, 73.92248, 
            81.916695, 79.91652, 29.97377, 28.976496, 27.976927, 153.92221, 151.91972, 149.91727, 148.91718, 147.91483, 146.9149, 
            143.912, 123.90527, 121.90344, 119.9022, 118.90331, 117.90161, 116.902954, 115.90175, 114.90335, 113.90279, 111.90482, 
            87.90562, 86.90888, 85.90926, 83.91343, 180.948, 179.94746, 158.92534, 125.90331, 124.904434, 123.902824, 122.904274, 
            121.90305, 119.904045, 129.90623, 127.904465, 232.03806, 46.951763, 45.95263, 49.944794, 48.947872, 47.94795, 204.9744, 
            202.97232, 168.93422, 50.943962, 49.947163, 185.95436, 183.95093, 182.95023, 181.9482, 179.9467, 173.93886, 172.9382, 
            171.93637, 170.93633, 169.93475, 167.9339, 175.94257, 69.92532, 67.92484, 66.92713, 65.92603, 63.929146, 93.90647, 
            91.90504, 90.90565, 89.9047, 95.90827]
    

    abundance = [0.48161000000000004, 0.51839, 1.0, 0.996, 0.00063, 0.00337, 1.0, 1.0, 0.8009999999999999, 0.19899999999999998, 
                0.7170000000000001, 0.11230000000000001, 0.0785, 0.06593, 0.0242, 0.00101, 0.00106, 1.0, 1.0, 0.5069, 
                0.49310000000000004, 0.011000000000000001, 0.9890000000000001, 4e-05, 0.02086, 0.00135, 0.00647, 0.96941, 0.00187, 
                0.128, 0.1249, 0.0089, 0.07490000000000001, 0.2873, 0.1222, 0.2413, 0.11130000000000001, 0.8843000000000001, 0.0025, 
                0.0019, 0.24230000000000002, 0.7576999999999999, 1.0, 0.02365, 0.095, 0.8379000000000001, 0.043449999999999996, 1.0, 
                0.6917, 0.001, 0.0006, 0.282, 0.249, 0.255, 0.18899999999999997, 0.023399999999999997, 0.149, 0.268, 
                0.22949999999999998, 0.336, 0.0161, 0.0014000000000000002, 0.522, 0.478, 1.0, 0.0028000000000000004, 0.021, 
                0.059000000000000004, 0.39892000000000005, 0.60108, 0.2484, 0.1565, 0.2047, 0.14800000000000002, 0.0218, 0.002, 
                0.2186, 0.07440000000000001, 0.3594, 0.07719999999999999, 0.2766, 0.21239999999999998, 0.00015, 0.99985, 1.0, 
                1.37e-08, 0.0016200000000000001, 0.35100000000000003, 0.13629, 0.27297, 0.18606, 0.0687, 0.2986, 0.1318, 0.231, 
                0.16870000000000002, 0.09970000000000001, 0.0015, 1.0, 1.0, 0.9570000000000001, 0.043, 0.373, 0.627, 0.067302, 
                0.000117, 0.932581, 0.0034999999999999996, 0.17300000000000001, 0.57, 0.115, 0.11599999999999999, 0.0225, 0.999088, 
                0.000902, 0.925, 0.075, 0.9741, 0.1101, 0.1, 0.7898999999999999, 1.0, 0.1592, 0.0925, 0.1484, 0.09630000000000001, 
                0.2413, 0.0955, 0.1668, 0.0037, 0.9963, 1.0, 1.0, 0.12179999999999999, 0.2713, 0.0564, 0.0576, 0.17190000000000003,
                0.083, 0.23800000000000002, 0.0925, 0.0027, 0.9048, 0.03634, 0.011399999999999999, 0.26222999999999996, 0.68077, 
                0.009260000000000001, 1.0, 0.002, 0.0004, 0.9976, 0.264, 0.161, 0.133, 0.016, 0.0158, 0.0002, 0.41, 1.0, 0.221, 
                0.24100000000000002, 0.013999999999999999, 0.524, 0.11720000000000001, 0.2646, 0.2733, 0.22329999999999997, 
                0.1114, 0.0102, 1.0, 0.0001, 0.07200000000000001, 0.253, 0.33799999999999997, 0.32899999999999996, 0.0079, 0.2783, 
                0.7217, 0.626, 0.374, 1.0, 0.18600000000000003, 0.316, 0.171, 0.126, 0.127, 0.018600000000000002, 0.0554, 0.0002, 
                0.0421, 0.0075, 0.9501999999999999, 0.4264, 0.5736, 1.0, 0.2377, 0.07629999999999999, 0.09359999999999999, 0.0089, 
                0.0874, 0.4961, 0.031, 0.0467, 0.9223, 0.22699999999999998, 0.267, 0.07400000000000001, 0.138, 0.113, 0.15, 0.031,
                0.0579, 0.0463, 0.3259, 0.0858, 0.2422, 0.0768, 0.14529999999999998, 0.0036, 0.006500000000000001, 0.0097,
                0.8258, 0.07, 0.0986, 0.005600000000000001, 0.9999800000000001, 0.00012, 1.0, 0.1893, 0.0712, 0.0479, 0.00905, 
                0.0259, 0.00095, 0.3387, 0.317, 1.0, 0.073, 0.08, 0.054000000000000006, 0.055, 0.738, 0.7047599999999999, 0.29524, 
                1.0, 0.9975, 0.0025, 0.28600000000000003, 0.307, 0.14279999999999998, 0.263, 0.0012, 0.318, 0.1612, 
                0.21899999999999997, 0.14300000000000002, 0.0305, 0.0013, 0.127, 0.006, 0.188, 0.040999999999999995, 
                0.27899999999999997, 0.486, 0.17379999999999998, 0.17149999999999999, 0.11220000000000001, 0.5145000000000001, 
                0.027999999999999997]
    

    isotopes = ['Ag', 'Ag', 'Al', 'Ar', 'Ar', 'Ar', 'As', 'Au', 'B', 'B', 'Ba', 'Ba', 'Ba', 'Ba', 'Ba', 'Ba', 'Ba', 'Be', 'Bi', 'Br', 
                'Br', 'C', 'C', 'Ca', 'Ca', 'Ca', 'Ca', 'Ca', 'Ca', 'Cd', 'Cd', 'Cd', 'Cd', 'Cd', 'Cd', 'Cd', 'Ce', 'Ce', 'Ce', 'Ce', 
                'Cl', 'Cl', 'Co', 'Cr', 'Cr', 'Cr', 'Cr', 'Cs', 'Cu', 'Dy', 'Dy', 'Dy', 'Dy', 'Dy', 'Dy', 'Dy', 'Er', 'Er', 'Er', 
                'Er', 'Er', 'Er', 'Eu', 'Eu', 'Fe ', 'Fe', 'Fe', 'Fe', 'Ga', 'Ga', 'Gd', 'Gd', 'Gd', 'Gd', 'Gd', 'Gd', 'Gd', 'Ge', 
                'Ge', 'Ge', 'Ge', 'Ge', 'H', 'H', 'He', 'He', 'Hf', 'Hf', 'Hf', 'Hf', 'Hf', 'Hg', 'Hg', 'Hg', 'Hg', 'Hg', 'Hg', 'Hg', 
                'Ho', 'I', 'In', 'In', 'Ir', 'Ir', 'K', 'K', 'K', 'Kr', 'Kr', 'Kr', 'Kr', 'Kr', 'Kr', 'La', 'La', 'Li', 'Li', 'Lu', 
                'Mg', 'Mg', 'Mg', 'Mn', 'Mo', 'Mo', 'Mo', 'Mo', 'Mo', 'Mo', 'Mo', 'N', 'N', 'Na', 'Nb', 'Nd', 'Nd', 'Nd', 'Nd', 'Nd', 
                'Nd', 'Nd', 'Ne', 'Ne', 'Ne', 'Ni', 'Ni', 'Ni', 'Ni', 'Ni', 'Np ', 'O', 'O', 'O', 'Os', 'Os', 'Os', 'Os', 'Os', 'Os', 
                'Os', 'P', 'Pb', 'Pb', 'Pb', 'Pb', 'Pd', 'Pd', 'Pd', 'Pd', 'Pd', 'Pd', 'Pr', 'Pt', 'Pt', 'Pt', 'Pt', 'Pt', 'Pt', 
                'Rb ', 'Rb', 'Re', 'Re', 'Rh', 'Ru', 'Ru', 'Ru', 'Ru', 'Ru', 'Ru', 'Ru', 'S', 'S', 'S', 'S', 'Sb', 'Sb', 'Sc', 'Se', 
                'Se', 'Se', 'Se', 'Se', 'Se', 'Si', 'Si', 'Si', 'Sm', 'Sm', 'Sm', 'Sm', 'Sm', 'Sm', 'Sm', 'Sn', 'Sn', 'Sn', 'Sn', 
                'Sn', 'Sn', 'Sn', 'Sn', 'Sn', 'Sn', 'Sr', 'Sr', 'Sr', 'Sr', 'Ta', 'Ta', 'Tb', 'Te', 'Te', 'Te', 'Te', 'Te', 'Te', 
                'Te', 'Te', 'Th', 'Ti', 'Ti', 'Ti', 'Ti', 'Ti', 'TI', 'TI', 'Tm', 'V', 'V', 'W', 'W', 'W', 'W', 'W', 'Yb', 'Yb', 'Yb', 
                'Yb', 'Yb', 'Yb', 'Yb', 'Zn', 'Zn', 'Zn', 'Zn', 'Zn', 'Zr', 'Zr', 'Zr', 'Zr', 'Zr ']
    
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


                #removes any molecule who's probability is below 0.00001

                if render_imprecise_list: #only removes low-probability arrangements if render_imprecise_list arg is True
                    if new_proba>0.00001:
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

def peak_merger(x_in, y_in, apparatus_resolution) -> list[float]:
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

    if not x_in:
        raise ValueError('Empty list entry')
    if not y_in:
        raise ValueError('Empty list entry')
    if len(x_in) != len(y_in):
        raise ValueError('Lists should be of the same size')

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

    return x_out, y_out

def sulphur_nitrogen_adder(x_in, y_in, has_N, has_S, count_N, count_S) -> list[float]:
    #---------------------------------------------------------------------------------------------#
    '''
    sulphur_nitrogen_adder(x_in, y_in, has_N, has_S, count_N, count_S)
    
    Input: two lists + 4 booleans:
    1. ordered list of the masses with values on y merged together if peaks within precision of apparatus
    2. ordered list of the probabilities of apparation of each of the molecules
    3. boolean that tells if the molecule contains Nitrogen
    4. boolean that tells if the molecule contains Sulphur
    5. number of Nitrogen atoms in the molecule
    6. number of Sulphur atoms in the molecule
    
    Output: two lists:
    1. ordered list of the masses with values on y merged together if peaks within precision of apparatus
    2. ordered list of the probabilities of apparation of each of the molecules
    
    (the mass in list 1 at index i is associated to the probability at index i in list 2)
    '''
    #---------------------------------------------------------------------------------------------#

    if not x_in:
        raise ValueError("Enter a non-empty list")
    if not y_in:
        raise ValueError("Enter a non-empty list")
    if len(x_in) != len(y_in):
        raise ValueError("The two lists must have the same length")
    if has_N not in [True, False]:
        raise ValueError("The third argument must be a boolean")
    if has_S not in [True, False]:
        raise ValueError("The fourth argument must be a boolean")
    if count_N < 0:
        raise ValueError("The fifth argument must be a positive integer")
    if count_S < 0:
        raise ValueError("The sixth argument must be a positive integer")
    if type(count_N) != int:
        raise ValueError("The fifth argument must be an integer")
    if type(count_S) != int:
        raise ValueError("The sixth argument must be an integer")
    if x_in != sorted(x_in):
        raise ValueError("The first list must be ordered")

    maximum = max(y_in)

    if has_N:
        x_in.append(x_in[1] - 0.006)  
        y_in.append(0.0035*count_N*maximum)  
    
    if has_S:
        x_in.append(x_in[1]-0.004)  
        y_in.append(0.008*count_S*maximum)


    return x_in, y_in

def peak_sorter(x_in, y_in) -> list[float]:
    #---------------------------------------------------------------------------------------------#
    '''
    peak_sorter(x_in, y_in)
    
    Input: two lists:
    1. ordered list of the masses with values on y merged together if peaks within precision of apparatus
    2. ordered list of the probabilities of apparation of each of the molecules
    
    Output: two lists:
    1. ordered list of the masses with values on y merged together if peaks within precision of apparatus
    2. ordered list of the probabilities of apparation of each of the molecules
    
    (the mass in list 1 at index i is associated to the probability at index i in list 2)
    '''
    #---------------------------------------------------------------------------------------------#
    
    if not x_in:
        raise ValueError('Empty list entry')
    if not y_in:
        raise ValueError('Empty list entry')
    if len(x_in) != len(y_in):
        raise ValueError('Lists should be of the same size')

    x_out, y_out = [],[]
    while len(x_in)>0:
        min_x = min(x_in)
        index_min = x_in.index(min_x)
        x_out.append(min_x)
        y_out.append(y_in[index_min])
        x_in.pop(index_min)
        y_in.pop(index_min)

    return x_out, y_out

def delta_function_plotter(x_in, y_in) -> list[float]:
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


    if not x_in:
        raise ValueError('Enter a valid input')
    if not y_in:
        raise ValueError('Enter a valid input')
    if len(x_in) != len(y_in):
        raise ValueError('Both entry should be of the same length')
    
    min_x , max_x = min(x_in), max(x_in)
    

    x_axis, y_axis = [min_x-0.5],[0]
    for i in range (len(x_in)):
        x_axis.append(x_in[i]-10**(-100))
        x_axis.append(x_in[i])
        x_axis.append(x_in[i]+10**(-100))
        y_axis.append(0)
        y_axis.append(y_in[i])
        y_axis.append(0)

    x_axis.append(max_x+1)
    y_axis.append(0)

    return x_axis, y_axis

def double_plot(x_in,y_in):

    #---------------------------------------------------------------------------------------------#
    '''
    double_plot(x_in,y_in)
    
    Input: list of masses (x_in) and intensities (y_in)
    
    Output: 2 Bokeh graphs:
            - One that shows the mass spectrum of the molecule to which the user can interact
            -Another that is the same graph but shows where the user is zooming on the first graph
    '''
    #---------------------------------------------------------------------------------------------#

    # tells where to put the graduation on the graph
    if not x_in:
        raise ValueError("Enter a non-empty list as first argument")
    if not y_in:
        raise ValueError("Enter a non-empty list as second argument")
    if len(x_in) != len(y_in):
        raise ValueError("The two lists must have the same length")
    
    ticked_peaks = []
    for i in range(len(x_in)):
        if y_in[i] > 0.0001:
            ticked_peaks.append(x_in[i])

    #creates the principal graph, mass spectrum of the molecule (interactive)

    p1 = figure(width=700, title=f'Mass spectrum of molecule')
    p1 = figure(x_axis_label = '[m/z]')
    p1 = figure(y_axis_label = 'Abundance')
    p1.xaxis.axis_label = "[m/z]"
    p1.height = 500
    p1.xaxis.ticker = FixedTicker(ticks=ticked_peaks)
    p1.add_tools(WheelPanTool(dimension="height"))
    p1.add_tools(WheelZoomTool(dimensions="height"))
    p1.line(x_in, y_in, line_width=1)
    p1.xaxis.major_label_orientation = "horizontal"

     #creates the secondary graph, mass spectrum of the molecule (non-interactive)

    p2 = figure(title="Simulated Mass Spectrum", x_axis_label='Mass [Th]', y_axis_label='Intensity')
    p2 = figure(width=300, title=f'Mass spectrum of molecule')
    p2 = figure(toolbar_location=None)
    p2.height = 300
    p2.line(x_in, y_in, legend_label="Mass spectrum", line_width=1)

    #creates a tool in order that the second graph shows where the zoom is on the first one

    box = BoxAnnotation(left=0, right=0, bottom=0, top=0,
    fill_alpha=0.1, line_color='red', fill_color='cornflowerblue')

    jscode = """
        box[%r] = cb_obj.start
        box[%r] = cb_obj.end
    """

    xcb = CustomJS(args=dict(box=box), code=jscode % ('left', 'right'))
    ycb = CustomJS(args=dict(box=box), code=jscode % ('bottom', 'top'))

    p1.x_range.js_on_change('start', xcb)
    p1.x_range.js_on_change('end', xcb)
    p1.y_range.js_on_change('start', ycb)
    p1.y_range.js_on_change('end', ycb)

    # adds the functionnality to the second figure

    p2.add_layout(box)

    # creates a layout that displays the 2 graphs

    layout = column(p1, p2)
    print('here')
    return layout

def functional_group_finder(mol_smi):

    #---------------------------------------------------------------------------------------------#
    '''
    functional_group_finder(mol_smi)
    
    Input: molecule under SMILEs representation
    
    Output: list containing every functionl group contained (if a functional group is contained twice in the molecule, it will appear twice in this list)
    '''
    #---------------------------------------------------------------------------------------------#

    if not mol_smi:
        raise ValueError("Enter a non-empty string as argument")
    # initiate variables
    functional_groups_contained, mol_in = [], Chem.MolFromSmiles(mol_smi)

    # dictionnary of all the considered functional groups to check (some might be missing)

    functional_groups_smarts = {'Alcohol': 'C[Oh1+0]','Aldehyde': 'C[Ch1]=O','Ketone': 'CC(=O)C','Carboxylic Acid': 'CC(=O)[Oh1]',
        'Ester': 'CC(=O)[Oh0]','Ether': '*[Oh0]*','Amide': 'C(=O)N','Amine': '[C][N]','Nitrile': 'C#N','Chloride': 'Cl',
        'Bromide': 'Br','Fluoride': 'F','Iodide': 'I','Alkene': 'C=C','Alkyne': 'C#C','Imine': 'C=N*','Amino acid': '[Nh2][Ch1*]C(=O)O',
        'Proline': '[Nh1][Ch1*]C(=O)O','Thiol': '[Sh1]','Sulfide': '*[Sh0]*','Acyl Chloride': 'CC(=O)Cl','Anhydride': '*[Ch0](=O)O[Ch0](=O)*',
        'Nitro': 'C[N+](=O)[O-]','Enamine': 'C=C[Nh0]','Enamine2': 'C=C[Nh1]','Enamine3': 'C=C[Nh2]','Imide': 'C(=O)NC(=O)*',
        'Azide': 'CNNN','Enol': 'C=C([Oh1])C','Hemiacetal': 'CC(O)(O)C','Carbonate': '[Oh0]C(=O)[Oh0]','Carbonate2': '[Oh1]C(=O)[Oh1]',
        'Disulfide': 'CSSC','Sulfoxide': 'CS(=O)C','Sulfone': '*[So2](=O)(=O)*','Sulfonic acid': '*S(=O)(=O)[Oh1]','Thioester': 'C(=O)S*',
        'Phosphine': '*[Po0](*)*','Phosphate': '*OP(=O)(O)O','Benzene': 'c1ccccc1','Peroxide':'C[Oh0][Oh0]C'
    }

    # check that the substructure from functional_groups_smarts are contained in mol_smi

    for name, smarts in functional_groups_smarts.items():
        if mol_in.HasSubstructMatch(Chem.MolFromSmarts(smarts)):
            for _ in range(len(mol_in.GetSubstructMatches(Chem.MolFromSmarts(smarts)))):
                functional_groups_contained.append(name)

    # exceptions for conflicts during the iteration of functional groups
    for functional_group in functional_groups_contained:
        if 'Ester' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Ether' in functional_groups_contained:
                    functional_groups_contained.remove('Ether')
        elif functional_group == 'Carboxylic Acid':
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Alcohol' in functional_groups_contained:
                    functional_groups_contained.remove('Alcohol')
        elif 'Phosphate' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Ether' in functional_groups_contained:
                    functional_groups_contained.remove('Ether')
        elif 'Thioester' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Sulfide' in functional_groups_contained:
                    functional_groups_contained.remove('Sulfide')
        elif 'Sulfonic acid' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Sulfide' in functional_groups_contained:
                    functional_groups_contained.remove('Sulfide')
        elif 'Sulfoxide' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Sulfide' in functional_groups_contained:
                    functional_groups_contained.remove('Sulfide')
        elif 'Acyl Chloride' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Sulfide' in functional_groups_contained:
                    functional_groups_contained.remove('Sulfide')
        elif 'Anhydride' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Ester' in functional_groups_contained:
                    functional_groups_contained.remove('Ester')
                if 'Ester' in functional_groups_contained:
                    functional_groups_contained.remove('Ester')
                if 'Ether' in functional_groups_contained:
                    functional_groups_contained.append('Ether')
        elif 'Enamine2' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Enamine2' in functional_groups_contained:
                    functional_groups_contained.remove('Enamine2')
                if 'Enamine' in functional_groups_contained:
                    functional_groups_contained.append('Enamine')
        elif 'Enamine3' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Enamine3' in functional_groups_contained:
                    functional_groups_contained.remove('Enamine3')
                if 'Amine' in functional_groups_contained:
                    functional_groups_contained.remove('Amine')
                functional_groups_contained.append('Enamine')
        elif 'Imide' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Amide' in functional_groups_contained:
                    functional_groups_contained.remove('Amide')
                if 'Amide' in functional_groups_contained:
                    functional_groups_contained.remove('Amide')
        elif 'Enol' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Alkene' in functional_groups_contained:
                    functional_groups_contained.remove('Alkene')
                if 'Alcohol' in functional_groups_contained:
                    functional_groups_contained.remove('Alcohol')
        elif 'Hemiacetal' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Alcohol' in functional_groups_contained:
                    functional_groups_contained.remove('Alcohol')
                if 'Alcohol' in functional_groups_contained:
                    functional_groups_contained.remove('Alcohol')
        elif 'Carbonate2' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Alcohol' in functional_groups_contained:
                    functional_groups_contained.remove('Alcohol')
                if 'Alcohol' in functional_groups_contained:
                    functional_groups_contained.remove('Alcohol')
                if 'Carbonate2' in functional_groups_contained:
                    functional_groups_contained.remove('Carbonate2')
                if 'Carbonate' in functional_groups_contained:
                    functional_groups_contained.append('Carbonate')
        elif 'Disulfide' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Sulfide' in functional_groups_contained:
                    functional_groups_contained.remove('Sulfide')
                if 'Sulfide' in functional_groups_contained:
                    functional_groups_contained.remove('Sulfide')
        elif 'Peroxide' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Ether' in functional_groups_contained:
                    functional_groups_contained.remove('Ether')
                if 'Ether' in functional_groups_contained:
                    functional_groups_contained.remove('Ether')
        elif 'Amide' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Amine' in functional_groups_contained:
                    functional_groups_contained.remove('Amine')
                if 'Amine' in functional_groups_contained:
                    functional_groups_contained.remove('Amine')
    
    
    return functional_groups_contained

def functional_group_display(contained_functional_groups):

    #---------------------------------------------------------------------------------------------#
    '''
    functional_group_display(contained_functional_groups)
    
    Input: list of all contained functional groups in the molecule

    Output: bokeh table with the contained functional groups and their image
    '''
    #---------------------------------------------------------------------------------------------#

    #dictionnary with functional groups and associated SMARTs
    functional_groups_smarts = {
        'Alcohol': 'C[Oh1+0]',
        'Aldehyde': 'C[Ch1]=O',
        'Ketone': 'CC(=O)C',
        'Carboxylic Acid': 'CC(=O)[Oh1]',
        'Ester': 'CC(=O)[Oh0]',
        'Ether': '*[Oh0]*',
        'Amide': 'C(=O)N',
        'Amine': '[C][N]',
        'Nitrile': 'C#N',
        'Chloride': 'Cl',
        'Bromide': 'Br',
        'Fluoride': 'F',
        'Iodide': 'I',
        'Alkene': 'C=C',
        'Alkyne': 'C#C',
        'Imine': 'C=N*',
        'Amino acid': '[Nh2][Ch1*]C(=O)O',
        'Proline': '[Nh1][Ch1*]C(=O)O',
        'Thiol': '[Sh1]',
        'Sulfide': '*[Sh0]*',
        'Acyl Chloride': 'CC(=O)Cl',
        'Anhydride': '*[Ch0](=O)O[Ch0](=O)*',
        'Nitro': 'C[N+](=O)[O-]',
        'Enamine': 'C=C[Nh0]',
        'Enamine2': 'C=C[Nh1]',
        'Enamine3': 'C=C[Nh2]',
        'Imide': 'C(=O)NC(=O)*',
        'Azide': 'CNNN',
        'Enol': 'C=C([Oh1])C',
        'Hemiacetal': 'CC(O)(O)C',
        'Carbonate': '[Oh0]C(=O)[Oh0]',
        'Carbonate2': '[Oh1]C(=O)[Oh1]',
        'Disulfide': 'CSSC',
        'Sulfoxide': 'CS(=O)C',
        'Sulfone': '*[So2](=O)(=O)*',
        'Sulfonic acid': '*S(=O)(=O)[Oh1]',
        'Thioester': 'C(=O)S*',
        'Phosphine': '*[Po0](*)*',
        'Phosphate': '*OP(=O)(O)O',
        'Benzene': 'c1ccccc1',
        'Peroxide':'C[Oh0][Oh0]C'
    }

    #initiate empty variables
    present_group_smarts = []    
    present_group_images_base64 = []
    
    #appends the smarts of the contained functional groups
    for i,j in functional_groups_smarts.items():
        for x in contained_functional_groups:
            if x == i:
                present_group_smarts.append(j)

    #converts the smarts to images in base64 format
    for x in present_group_smarts:

        #converts SMARTs to SMILEs for the images to be nicer
        mol_x = Chem.MolFromSmarts(x)
        mol_smi = Chem.MolToSmiles(mol_x)
        mol = Chem.MolFromSmiles(mol_smi)

        if mol:

            #creates the image
            image = Draw.MolToImage(mol)

            # Convert the image to base64 format
            buffered = BytesIO()
            image.save(buffered, format="PNG")
            image_base64 = base64.b64encode(buffered.getvalue()).decode("utf-8")

            #appends the image to an empty dictionnary
            present_group_images_base64.append(image_base64)

    #creates a dictionnary that links the name of the functional group to its image
    data = dict(
        groups=contained_functional_groups,
        images=present_group_images_base64
    )
    source = ColumnDataSource(data)

    #template for the bokeh table that read the base64 format 
    template = """
    <div>
        <img src="data:image/png;base64, <%= value %>" style="width:50px;height:50px;">
    </div>
    """

    # initiallizing the bokeh figure using the previous template for each functional group
    columns = [
        TableColumn(field="groups", title="Functional Groups"),
        TableColumn(field="images", title="Images", width=200, formatter=HTMLTemplateFormatter(template=template))
    ]
    num_groups = len(contained_functional_groups)

    table_height = min(200 + num_groups * 60, 800)

    data_table = DataTable(source=source, columns=columns, width=250, height=table_height, row_height=60)

    return data_table

def mol_web_show(mol_smi, show_Hs=False, show_3D = False):

    #---------------------------------------------------------------------------------------------#
    '''
    mol_web_show(mol_smi, show_Hs=False, show_3D = False)
    
    Input: SMILEs of a molecule. Also specify if want the function to show the hydrogens explicitely or the 3D
    
    Output: image of the molecule as a bokeh plot
    '''
    #---------------------------------------------------------------------------------------------#

    # name of the file name to create
    filename = 'molecule_image.png'

    #finds the current directory
    current_directory = os.getcwd()

    #creates a file path to the current directory
    filepath = os.path.join(current_directory, filename)

    #checks if the file already exists
    if not os.path.exists(filepath):

        #if no, creates the path
        with open(filepath, 'a'):
            pass
    else:

        #else pass
        pass

    # Generate the image from the molecule
    mol = Chem.MolFromSmiles(mol_smi)

    # Adds the hydrogens to the molecule if specified
    if show_Hs:
        mol = Chem.AddHs(mol)

    # Show the molecule in 3D if specified
    if show_3D:
        mol = AllChem.EmbedMolecule(mol)

    image = Draw.MolToImage(mol)

    # Save the image to a file
    image.save(filepath)

    # Creating a Bokeh figure to display the molecule
    p = figure(width=350, height=350,toolbar_location=None, x_range=(0, 1), y_range=(0, 1))
    p.image_url(url=[filepath], x=0, y=1, w=1, h=1)

    # Hide grid lines and axes
    p.xgrid.grid_line_color = None
    p.ygrid.grid_line_color = None
    p.xaxis.visible = False
    p.yaxis.visible = False

    return p
 
def all_in_one(p1, p2, p3, p4):

    #---------------------------------------------------------------------------------------------#
    '''
    all_in_one(p1,p2,p3, p4)
    
    Input: 3 bokeh plots
            Usually used in this package:
                    - p1 : bokeh double plot of mass spectrometry
                    - p2 : image of the molecule
                    - p3 : table of functional groups
                    - p4 : buttons with info on the molecule
    
    Output: bokeh page with all 4 graphs well arranged

    '''
    #---------------------------------------------------------------------------------------------#
    
    if not p1:
        raise ValueError("Enter a non-empty plot")
    if not p2:
        raise ValueError("Enter a non-empty plot")
    if not p3:
        raise ValueError("Enter a non-empty plot")
    if not p4:
        raise ValueError("Enter a non-empty plot")
    
    #creates a layout in row with p3 and p4
    layout1 = row(p3, p4)

    #creates a layout in column with p2 and layout1
    layout2 = column(p2, layout1)

    #creates the final layout in row with layout1 and p1
    layout = row(p1, layout2)

    return layout

def smiles_to_3D_plot(mol_smi):

     #---------------------------------------------------------------------------------------------#
    '''
    smiles_to_3D_plot(smiles)
    
    Input: molecule is a SMILEs format

    Output: panel graph of the molecule in 3D and interactive
    '''
    #---------------------------------------------------------------------------------------------#

    # Generate 3D coordinates from SMILES string
    mol2 = Chem.MolFromSmiles(mol_smi)
    mol2 = Chem.AddHs(mol2)  # Add hydrogens for better geometry optimization
    AllChem.EmbedMolecule(mol2, randomSeed=42)  # Embed the molecule in 3D space
    AllChem.MMFFOptimizeMolecule(mol2)  # Optimize the geometry using MMFF94 force field

    # Create a temporary file to store the 
    with tempfile.NamedTemporaryFile(delete=False) as tmp:
        tmp.write(f"{mol2.GetNumAtoms()}\n\n".encode('utf-8'))  # Write number of atoms
        for atom in mol2.GetAtoms():
            pos = mol2.GetConformer().GetAtomPosition(atom.GetIdx())
            tmp.write(f"{atom.GetSymbol()} {pos.x} {pos.y} {pos.z}\n".encode('utf-8'))
        tmp_path = tmp.name

    # Create the MolGraph object
    mg = MolGraph()
    
    # Read the data from the temporary file
    mg.read_xyz(tmp_path)

    # Create the Plotly figure object
    fig = to_plotly_figure(mg)

    return fig






