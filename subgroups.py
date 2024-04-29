from rdkit import Chem
from rdkit.Chem import Draw
import time

list_subgroups_smi = ['C=O','C(=O)O','COC']
list_subgroups,contained_subgroup_list,contained_subgroup_list_smi = [],[],[]
#creates new list with functional groups in mol representation
for smile in list_subgroups_smi:
    list_subgroups.append(Chem.MolFromSmiles(smile))

def subgroup_finder(mol_smi):
    mol = Chem.MolFromSmiles(mol_smi)
    for subgroup in list_subgroups:
        substructmol = Chem.MolFromSmiles(subgroup)
        if mol.GetSubstructMatch(subgroup):
            contained_subgroup_list.append(subgroup)
            contained_subgroup_list_smi.append(substructmol)

            

            
            
        
    return contained_subgroup_list,contained_subgroup_list_smi #here choose either mol or SMILES format for output

print(subgroup_finder(input('Input a SMILES: ')))

def subgroug_hydrogen_counter(mol_smi):
    mol = Chem.MolFromSmiles(mol_smi)
