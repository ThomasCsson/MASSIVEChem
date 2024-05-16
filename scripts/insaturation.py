from rdkit import Chem
from rdkit.Chem import Draw

def calculate_insaturation(mol_smile) -> int:
    # This function takes as input a string representing the SMILES notation of a molecule. 
    # Using the RDKit library, it creates a representation of the molecule, adds hydrogens and then calculates the unsaturation level of the molecule. 
    # The level of unsaturation is determined by counting the number of carbon, nitrogen and halogen atoms, then applying a specific formula. 
    # Finally, the function displays the unsaturation level and the image of the molecule, and indicates the time required to perform the function.
    mol_1 = Chem.MolFromSmiles(mol_smile)
    mol = Chem.AddHs(mol_1)

    C, N, HX = 0, 0, 0
    halogens_hydrogen = ['F', 'Cl', 'Br', 'I', 'At', 'H']
    for atom in mol.GetAtoms():
        atom_sym = atom.GetSymbol()
        if atom_sym == 'C':
            C += 1
        elif atom_sym == 'N':
            N += 1
        elif atom_sym in halogens_hydrogen:
            HX += 1

    insaturation = C + 1 + (N - HX) / 2

    print(f'The insaturation of the given molecule is {insaturation}')
    img = Draw.MolToImage(mol_1)
    img.show()
