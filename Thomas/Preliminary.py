from rdkit import Chem

mol_smi = 'CC'
mol_pre = Chem.MolFromSmiles(mol_smi)
mol = Chem.AddHs(mol_pre)

atoms_present = []
for atom in mol.GetAtoms():
    AtomSymbol = atom.GetSymbol()
    atoms_present.append(AtomSymbol)

print(atoms_present)



isotopes =  ['C','C']
mass = [12,13]
abundance = [99,1]
dict = {}

list_atoms = ['C','C','N','C']
list_output = []
mass_copy = mass.copy()
abundance_copy = abundance.copy()
isotopes_copy = isotopes.copy()
for i in range (0,isotopes.count(list_atoms[0])):
    
    index = isotopes.index(list_atoms[0])
    list_output.append([mass_copy[index],abundance_copy[index]])
    mass_copy.pop(index)
    abundance_copy.pop(index)
    isotopes_copy.pop(index)

list_atoms = list_atoms[1:]
print(list_atoms)
print(isotopes_copy)

while len(list_atoms)>0:
    new_output = []
    output = isotopes.index(list_atoms[0])
    for i in range (0,isotopes.count(list_atoms[0])):
        for j in range(0,len(list_output)):
            mass_copy = mass
            abundance_copy = abundance
            isotopes_copy = isotopes
            new_output.append([mass[index]*output[j][0],mass[index]*output[j][0]])
            mass_copy.pop(index)
            abundance_copy.pop(index)
            isotopes_copy.pop(index)
            print(new_output)
        print(new_output)
    list_atoms[1:]

    
    

print(list_output)











        
    
    
