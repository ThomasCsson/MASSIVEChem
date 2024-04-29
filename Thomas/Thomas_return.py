#Imports

from rdkit import Chem
from rdkit.Chem import Draw
import time
import pandas as pd
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from bokeh.plotting import figure, show

#Turn data (of Symbol | Mass | Probability) into lists 

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


#Turn SMILEs representation into list of atomic symbols

mol_smi = input('Enter SMILEs: ')
mol_without_Hs = Chem.MolFromSmiles(mol_smi)
mol = Chem.AddHs(mol_without_Hs)
img = Draw.MolToImage(mol)
img.show()


#Timing element (not very useful but its pretty so it's staying)

start_time = time.time()


#Function that takes in list of atoms and then gives a list of list of shape [[mass,proba],[mass,proba],...]

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
    has_S = False
    if 'N' in list_atoms:
        has_N = True
    elif 'S' in list_atoms:
        has_S = True

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

    print(list_output)
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
    if has_N:
        x_axis_final.append()  #Add values
        y_axis_final.append()  #Add values
    
    if has_S:
        x_axis_final.append()  #Add values
        y_axis_final.append()  #Add values

    

    '''HERE, IF YOU WERE TO 'return x_axis_final, y_axis_final', THE OUTPUT IS TWO LISTS, THE FIRST OF THE VALUES OF THE X AXIS (COMBINED MASSES) AND THE SECOND OF THE VALUES ON Y '''
    
    #timing

    end_time = time.time()
    print(f'Runtime: {end_time-start_time}s')
    

    #plotting

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


    #plotting with pyplot

    


    #plotting with bokeh

    x, y = x_final_final, y_final_final

    p = figure(title="Simple Line Graph", x_axis_label='Mass of compound [g/mol]', y_axis_label='Abundance [%]')
    p.line(x,y,legend_label="Line", line_width=2)
    show(p)





    return x_axis_final,y_axis_final

print(main_function(mol))




end_time = time.time()
print(f'Runtime: {end_time-start_time}s')