from rdkit import Chem
import time








isotopes =  ['C','C','N','N']
mass = [12,13,14,15]
abundance_percent = [99,1,99.5,0.5]
abundance = []

for percent in abundance_percent:
    abundance.append(percent/100)





mol_smi = input('Enter SMILEs: ')
mol = Chem.MolFromSmiles(mol_smi)
start_time = time.time()

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
                list_output_new.append([(list_output[i][0] + mass_copy[index]),list_output[i][1] * abundance_copy[index]])
                mass_copy.pop(index)
                abundance_copy.pop(index)
                isotopes_copy.pop(index)
                
        list_output = list_output_new
        list_atoms.pop(0)
    return list_output


print(main_function(mol))




end_time = time.time()
print(f'Run: {end_time-start_time}s')