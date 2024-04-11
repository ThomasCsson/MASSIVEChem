from rdkit import Chem
from rdkit.Chem import Draw
import time

#dictionary of all subgroups to check for alogn with their associated values for the NMR spectrum
list_subgroups_dict = {'C=O':1.1,'C(=O)O':1.2,'COC':1.3}

#lists that will contain the subgroups & the values to add to final graph
list_contained_subgroups, list_contained_subgroups_values = [],[]

#function that takes the SMILES of a molecule in entry and gives the subgroupps contained along with the values for the NMR specturm
def subgroup_nmr_value (mol_smi):
    start_time = time.time()
    mol = Chem.MolFromSmiles(mol_smi)
    for SMILES, value in list_subgroups_dict.items():
        if mol.HasSubstructMatch(Chem.MolFromSmiles(SMILES)):
            list_contained_subgroups.append(SMILES)
            list_contained_subgroups_values.append(value)
        else:
            list_contained_subgroups.append(f'NOT {SMILES}')
    end_time = time.time()
    return list_contained_subgroups,list_contained_subgroups_values

print(subgroup_nmr_value(input('Input a SMILES: ')))
print(f'The program took {en}')
    


def timer(mol):
    start_time = time.time()
    def subgroup_nmr_value (mol_smi):
    
        mol = Chem.MolFromSmiles(mol_smi)
        for SMILES, value in list_subgroups_dict.items():
            if mol.HasSubstructMatch(Chem.MolFromSmiles(SMILES)):
                list_contained_subgroups.append(SMILES)
                list_contained_subgroups_values.append(value)
            else:
                list_contained_subgroups.append(f'NOT {SMILES}')
        end_time = time.time()
        return list_contained_subgroups,list_contained_subgroups_values
    end_time = time.time()
    return (end_time-start_time)


print(subgroup_nmr_value(input('Input a SMILES: ')))
print(timer(input('Input a SMILES: ')))


