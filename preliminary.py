from rdkit import Chem
from rdkit.Chem import Draw
import time

#dictionary of all subgroups to check for alogn with their associated values for the NMR spectrum
list_subgroups_shift_dict = {
    'C1CC1': 0.5,
    'c1ccccc1': 2.5,
    'CNC': 2.8,
    'CSC': 2.5,
    'CS': 3.8,
    'CN': 4.2,
    'C(=O)O': 11,
    'C=O': 10,
}


#function that takes the SMILES of a molecule in entry and gives the subgroups contained along with the values for the NMR specturm

def subgroup_nmr_value (mol_smi):

    list_contained_subgroups, list_contained_subgroups_values = [],[]
    mol = Chem.MolFromSmiles(mol_smi)
    for SMILES, value in list_subgroups_shift_dict.items():
        substruct = Chem.MolFromSmiles(SMILES)
        if mol.HasSubstructMatch(substruct):
            list_contained_subgroups.append(SMILES)
            list_contained_subgroups_values.append(value)
            for atom in substruct.GetAtoms():
                 if atom.GetSymbol()=='C':
                    index = mol.GetIdx(atom)
                    print(index)
        else:
            list_contained_subgroups.append(f'NOT {SMILES}')

    return list_contained_subgroups,list_contained_subgroups_values

print(subgroup_nmr_value(input('Input a SMILES: ')))

#function that takes SMILES as input and gives insaturation as output
def insaturation_level (mol_smi):

    C,N,HX, = 0,0,0
    molwithoutHs = Chem.MolFromSmiles(mol_smi)
    mol = Chem.AddHs( molwithoutHs)
    halogens_hydrogen = ['F','Cl','Br','I','At','H']
    for atom in mol.GetAtoms():
            atom_sym = atom.GetSymbol()
            if atom_sym == 'C':
                  C = C+1
            elif atom_sym == 'N':
                  N = N + 1
            elif atom_sym in halogens_hydrogen:
                  HX = HX + 1
    insaturation = C + 1 + (N-HX)/2
    return insaturation

print(insaturation_level (input('SMILES: ')))
    




