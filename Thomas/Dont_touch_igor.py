# Imports

from rdkit import Chem
import time
import pandas as pd
import matplotlib.pyplot as plt

# Turn data (of Symbol | Mass | Probability) into lists

df = pd.read_csv(
    '/Users/igorgonteri/Documents/GitHub/ppchem-project-Christiansson-Gonteri-Humery/Thomas/abundance.txt'
    , sep='\t'
    , header=None
    , names=['Atom', 'Mass', 'Percentage'])

# mass = [mass1, mass2, mass3,...]
mass = df['Mass'].tolist()
# change from elements (not sure) to floats
mass = [float(m) for m in mass]

# abundance = [ab1, ab2, ab3,...]
abundance_percent = df['Percentage'].tolist()
# change from elements (not sure) to floats
abundance_percent = [float(ap) for ap in abundance_percent]
# from percent to proba
abundance = []
for percent in abundance_percent:
    abundance.append(percent / 100)

# isotopes = [iso1, iso2, iso3,...]
isotopes = df['Atom'].tolist()

abundance = []
for percent in abundance_percent:
    abundance.append(percent / 100)

# Turn SMILEs representation into list of atomic symbols

mol_smi = input('Enter SMILEs: ')
mol_without_Hs = Chem.MolFromSmiles(mol_smi)
mol = Chem.AddHs(mol_without_Hs)

# Timing element (not very useful but its pretty so it's staying)
start_time = time.time()


# Function that takes in list of atoms and then gives a list of list of shape [[mass,proba],[mass,proba],...]

def main_function(mol):
    list_atoms = []
    for atom in mol.GetAtoms():
        list_atoms.append(atom.GetSymbol())

    '''In the case of ionisation by proton, we need to add a H+ ion, which is done in the following'''
    if 'H' in list_atoms:
        list_atoms.remove('H')
    print(list_atoms)

    list_output = []
    mass_copy = mass.copy()
    abundance_copy = abundance.copy()
    isotopes_copy = isotopes.copy()

    for i in range(isotopes.count(list_atoms[0])):
        index = isotopes.index(list_atoms[0])
        list_output.append([mass_copy[index], abundance_copy[index]])
        mass_copy.pop(index)
        abundance_copy.pop(index)
        isotopes_copy.pop(index)

    list_atoms = list_atoms[1:]

    while len(list_atoms) > 0:

        # This runs over all atoms in molecule
        list_output_new = []

        for i in range(len(list_output)):

            # This makes us run over all lists in list obtained before
            mass_copy = mass.copy()
            abundance_copy = abundance.copy()
            isotopes_copy = isotopes.copy()

            for _ in range(isotopes.count(list_atoms[0])):
                # This for-loop runs over all isotope types of the atom type in pos 0 in list_atoms (input list)
                index = isotopes.index(list_atoms[0])
                new_mass = list_output[i][0] + mass_copy[index]
                new_proba = list_output[i][1] * abundance_copy[index]
                list_output_new.append([new_mass, new_proba])
                mass_copy.pop(index)
                abundance_copy.pop(index)
                isotopes_copy.pop(index)

        list_output = list_output_new
        list_atoms.pop(0)
    # Conversion of list_output (which is a list of lists) to a combination of two lists (x_axis & y_axis)
    x_axis, y_axis = [], []
    x_axis_final, y_axis_final = [], []
    for j in range(len(list_output)):
        x_axis.append(list_output[j][0])
        y_axis.append(list_output[j][1])


    return x_axis, y_axis

print(main_function(mol))

masses, intensities = main_function(mol)

def peak_posi(masses,intensities):

    aggregated_intensities = {'Peak positions': [], 'Peak Intensities': []}

    for mass, intensity in zip(masses, intensities):
        if mass in aggregated_intensities['Peak positions']:
            idx = aggregated_intensities['Peak positions'].index(mass)
            aggregated_intensities['Peak Intensities'][idx] += intensity
        else:
            aggregated_intensities['Peak positions'].append(mass)
            aggregated_intensities['Peak Intensities'].append(intensity)

    return aggregated_intensities

aggregated_intensities = peak_posi(masses,intensities)
peak_positions = aggregated_intensities['Peak positions']
peak_intensities = aggregated_intensities['Peak Intensities']

print("Peak positions:", peak_positions)
print("Peak intensities:", peak_intensities)




end_time = time.time()
print(f'Runtime: {end_time-start_time}s')