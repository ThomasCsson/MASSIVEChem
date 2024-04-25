#Imports

from rdkit import Chem
import time
import pandas as pd




#Turn data into lists

df = pd.read_csv('/Users/thomaschristiansson/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Thomas/abundance.txt'
                 , sep='\t'
                 , header=None
                 , names=['Atom', 'Mass', 'Percentage'])

mass = df['Mass'].tolist()
mass = [float(m) for m in mass]

abundance_percent = df['Percentage'].tolist()
abundance_percent = [float(ap) for ap in abundance_percent]

isotopes = df['Atom'].tolist()

abundance = []
for percent in abundance_percent:
    abundance.append(percent/100)




#Turn SMILEs repre. into list of atomic symbols

mol_smi = input('Enter SMILEs: ')
mol = Chem.MolFromSmiles(mol_smi)
'mol = Chem.AddHs(mol_without_Hs)'
start_time = time.time()

#Function that takes in list of atoms and then gives a list of list of shape [[mass,proba],[mass,proba],...]

def main_function (mol):
    list_atoms = []
    for atom in mol.GetAtoms():
        list_atoms.append(atom.GetSymbol())
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

                #This for loop runs over all isotope types of the atom type in pos 0 in list_atoms (input list)
                index = isotopes.index(list_atoms[0])
                new_mass = list_output[i][0] + mass_copy[index]
                new_proba = list_output[i][1] * abundance_copy[index]
                list_output_new.append([new_mass,new_proba])
                mass_copy.pop(index)
                abundance_copy.pop(index)
                isotopes_copy.pop(index)
                
        list_output = list_output_new
        list_atoms.pop(0)
    return list_output


print(main_function(mol))




end_time = time.time()
print(f'Run: {end_time-start_time}s')