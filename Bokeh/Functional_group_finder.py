from rdkit import Chem


#for now the code only works for 1 functionnal group by molecule
input_mol = input('SMILEs : ')

functional_groups_smiles = {
    'Alcohol': 'CO',
    'Aldehyde': 'CC=O',
    'Ketone': 'CC(=O)C',
    'Carboxylic Acid': 'CC(=O)O',
    'Ester': 'CC(=O)OC',
    'Ether': 'COC',
    'Amide': 'CC(=O)N',
    'Amine': 'CN',
    'Nitrile': 'C#N',
    'Chloride': 'CCl',
    'Bromide': 'CBr',
    'Fluoride': 'CF',
    'Iodide': 'CI',
    'Alkene': 'C=C',
    'Alkyne': 'C#C',
    'Imine': 'C=NC',
    'Amino acid': 'CC(N)C(=O)O',
    'Thiol': 'CS',
    'Sulfides': 'CSC',
    'Acyl Chloride': 'CC(=O)Cl',
    'Anhydride': 'CC(=O)OC(=O)C',
    'Nitro': 'C[N+](=O)[O-]',
    'Enamine': 'C=CN',
    'Imide': 'C(=O)NC(=O)C',
    'Azide': 'CNNN',
    'Enol': 'C=C(O)C',
    'Hemiacetal': 'CC(O)(O)C',
    'Carbonate': 'OC(=O)O',
    'Disulfide': 'CSSC',
    'Sulfoxide': 'CS(=O)C',
    'Sulfone': 'CS(=O)(=O)C',
    'Sulfonic  acid': 'CS(=O)(=O)O',
    'Thioester': 'C(=O)SC',
    'Phosphine': 'CP',
    'Phosphate': 'COP(=O)(O)O',
    'Benzene': 'C1=CC=CC=C1'
}

functional_groups_smarts ={
    'Alcohol': '[#6]-[#8]',
    'Aldehyde': '[#6]-[#6]=[#8]',
    'Ketone': '[#6]-[#6](=[#8])-[#6]',
    'Carboxylic Acid': '[#6]-[#6](=[#8])-[#8]',
    'Ester': '[#6]-[#6](=[#8])-[#8]-[#6]',
    'Ether': '[#6]-[#8]-[#6]',
    'Amide': '[#6]-[#6](=[#8])-[#7]',
    'Amine': '[#6]-[#7]',
    'Nitrile': '[#6]#[#7]',
    'Chloride': '[#6]-[#17]',
    'Bromide': '[#6]-[#35]',
    'Fluoride': '[#6]-[#9]',
    'Iodide': '[#6]-[#53]',
    'Alkene': '[#6]=[#6]',
    'Alkyne': '[#6]#[#6]',
    'Imine': '[#6]=[#7]-[#6]',
    'Amino acid': '[#6]-[#6](-[#7])-[#6](=[#8])-[#8]',
    'Thiol': '[#6]-[#16]',
    'Sulfides': '[#6]-[#16]-[#6]',
    'Acyl Chloride': '[#6]-[#6](=[#8])-[#17]',
    'Anhydride': '[#6]-[#6](=[#8])-[#8]-[#6](=[#8])-[#6]',
    'Nitro': '[#6]-[#7+](=[#8])-[#8-]',
    'Enamine': '[#6]=[#6]-[#7]',
    'Imide': '[#6](=[#8])-[#7]-[#6](=[#8])-[#6]',
    'Azide': '[#6]-[#7]-[#7]-[#7]',
    'Enol': '[#6]=[#6](-[#8])-[#6]',
    'Hemiacetal': '[#6]-[#6](-[#8])(-[#8])-[#6]',
    'Carbonate': '[#8]-[#6](=[#8])-[#8]',
    'Disulfide': '[#6]-[#16]-[#16]-[#6]',
    'Sulfoxide': '[#6]-[#16](=[#8])-[#6]',
    'Sulfone': '[#6]-[#16](=[#8])(=[#8])-[#6]',
    'Sulfonic  acid': '[#6]-[#16](=[#8])(=[#8])-[#8]',
    'Thioester': '[#6](=[#8])-[#16]-[#6]',
    'Phosphine': '[#6]-[#15]',
    'Phosphate': '[#6]-[#8]-[#15](=[#8])(-[#8])-[#8]',
    'Benzene': '[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1'
}



def subgroup_value(mol_smi, dict_functional_groups):
    list_contained_subgroups, list_contained_subgroups_values = [], []
    mol = Chem.MolFromSmiles(mol_smi)

    for SMARTS, value in dict_functional_groups.items():
        substruct = Chem.MolFromSmarts(dict_functional_groups[SMARTS])  # Get the SMILES string from the dictionary

        if mol.HasSubstructMatch(substruct):
            list_contained_subgroups.append(SMARTS)
            list_contained_subgroups_values.append(value)
            indices = mol.GetSubstructMatch(substruct)
            print(indices)
        else:
            list_contained_subgroups.append(f'NOT {SMARTS}')

    print(list_contained_subgroups)

    if 'Aldehyde' and 'Ketone' in list_contained_subgroups:
        mol_aldehyde = Chem.MolFromSmiles('CC=O')
        mol_ketone = Chem.MolFromSmiles('C(=O)C')

        list_aldehyde_index = mol.GetSubstructMatch(mol_aldehyde)
        list_ketone_index = mol.GetSubstructMatch(mol_ketone)

        for x in list_aldehyde_index:
            if x in list_ketone_index:
                list_contained_subgroups[1] = 'NOT aldehyde'

    if 'Carboxylic acid' and 'Aldehyde' and 'Alcohol' in list_contained_subgroups:

        mol_acid = Chem.MolFromSmiles('CC(=O)O')
        mol_aldehyde = Chem.MolFromSmiles('CC=O')
        mol_alcohol = Chem.MolFromSmiles('CO')

        list_aldehyde_index = mol.GetSubstructMatch(mol_aldehyde)
        list_acid_index = mol.GetSubstructMatch(mol_acid)
        list_alcohol_index = mol.GetSubstructMatch(mol_alcohol)

        for x in list_acid_index:
            if x in list_aldehyde_index:
                list_contained_subgroups[1] = 'NOT Aldehyde'
            if x in list_alcohol_index:
                list_contained_subgroups[0] = 'NOT Alcohol'

    if 'Carboxylic acid' and 'Ester' and 'Ether' in list_contained_subgroups:

        mol_acid = Chem.MolFromSmiles('CC(=O)O')
        mol_ester = Chem.MolFromSmiles('CC(=O)OC')
        mol_ether = Chem.MolFromSmiles('COC')

        list_acid_index = mol.GetSubstructMatch(mol_acid)
        list_ester_index = mol.GetSubstructMatch(mol_ester)
        list_ether_index = mol.GetSubstructMatch(mol_ether)

        for x in list_ester_index:
            if x in list_acid_index:
                list_contained_subgroups[3] = 'NOT Carboxylic Acid'
            if x in list_ether_index:
                list_contained_subgroups[5] = 'NOT Ether'

    if 'Ether' and 'Alcohol' in list_contained_subgroups:
        mol_alcohol = Chem.MolFromSmiles('CO')
        mol_ether = Chem.MolFromSmiles('COC')

        list_alcohol_index = mol.GetSubstructMatch(mol_alcohol)
        list_ether_index = mol.GetSubstructMatch(mol_ether)

        for x in list_ether_index:
            if x in list_alcohol_index:
                list_contained_subgroups[0] = 'NOT Alcohol'

    """if 'Amide' and 'Aldehyde' in list_contained_subgroups:
    if 'Acyl Chloride' and 'Aldehyde' in list_contained_subgroups:
        
    if 'Thioester' and 'Ester':
    
    if 'Imine' and 'Amine' in list_contained_subgroups:
    
    iif 'Thio' and 'Sulfide in list_contained_subgroups:
    
    if 'Azide' and amine in list_contained_subgroups: """

    return list_contained_subgroups, list_contained_subgroups_values


print(subgroup_value(input_mol,functional_groups_smiles))


def check_functional_groups(molecule, functional_groups):
    found_groups = []
    mol = Chem.MolFromSmiles(molecule)
    if mol is None:
        return found_groups
    groups_to_remove = []
    for group, smarts in functional_groups.items():
        pattern = Chem.MolFromSmarts(smarts)
        if mol.HasSubstructMatch(pattern):
            found_groups.append(group)
            # Add related functional groups to remove list
            """if group == 'Ketone':
                groups_to_remove.append('Aldehyde')
            if group == 'Carboxylic acid':
                groups_to_remove.append('Aldehyde')
                groups_to_remove.append('Alcohol')
            if group == 'Ester':
                groups_to_remove.append('Carboxylic acid')
                groups_to_remove.append('Ether')
                groups_to_remove.append('Aldehyde')
            if group == 'Ether':
                groups_to_remove.append('Alcohol')
            if group == 'Imine':
                groups_to_remove.append('Amine')
            if group == 'Imide':
                groups_to_remove.append('Amine')
                groups_to_remove.append('Amide')
            if group == 'Amide':
                groups_to_remove.append('Aldehyde')
                groups_to_remove.append('Amine')
            if group == 'Acyl Chloride':
                groups_to_remove.append('Aldehyde')
                groups_to_remove.append('Chloride')
            if group == 'Thioester':
                groups_to_remove.append('Ester')
                groups_to_remove.append('Thiol')
                groups_to_remove.append('Sulfides')
            if group == 'Amino acid':
                groups_to_remove.append('Alcohol')
                groups_to_remove.append('Aldehyde')
            if group == 'Sulfides':
                groups_to_remove.append('Thiol')
            if group == 'Azide':
                groups_to_remove.append('Amine')
            if group == 'Disulfide':
                groups_to_remove.append('Thiol')
            if group == 'Nitro':
                groups_to_remove.append('Amine')
            if group == 'Sulfone':
                groups_to_remove.append('Sulfides')
                groups_to_remove.append('Sulfoxide')
            if group == 'Sulfonic acid':
                groups_to_remove.append('Thiol')

    # Remove related functional groups after checking all groups
    print(found_groups)
    for group in groups_to_remove:
        if group not in found_groups:
            pass
        else:
            found_groups.remove(group)"""
    return found_groups

for name ,mol2 in functional_groups_smiles.items():
    print(name,mol2)
    print(check_functional_groups(mol2,functional_groups_smarts))


def smiles_to_smarts(dict):
    dict2 = {}
    for x,y in dict.items():
        mol = Chem.MolFromSmiles(y)
        smart = Chem.MolToSmarts(mol)
        dict2[x] = smart
    return dict2

print(smiles_to_smarts(functional_groups_smiles))
"""

Aldehyde CC=O
[]

Carboxylic Acid CC(=O)O
['Alcohol', 'Carboxylic Acid']
Ester CC(=O)OC
['Ester', 'Ether']

Amide CC(=O)N
['Amide', 'Amine']
Nitrile C#N
[]
Imine C=NC
[]
Amino acid CC(N)C(=O)O
['Alcohol', 'Carboxylic Acid', 'Amine']


Acyl Chloride CC(=O)Cl
['Chloride', 'Acyl Chloride']
Anhydride CC(=O)OC(=O)C
['Ester', 'Ether']
Nitro C[N+](=O)[O-]
[]
Enamine C=CN
['Amine', 'Alkene', 'Enamine']
Imide C(=O)NC(=O)C
[]
Azide CNNN
[]

Hemiacetal CC(O)(O)C
['Alcohol']
Carbonate OC(=O)O
['Alcohol']

Sulfoxide CS(=O)C
[]
Sulfone CS(=O)(=O)C
[]
Sulfonic  acid CS(=O)(=O)O
[]
Thioester C(=O)SC
['Sulfides']
Phosphine CP
[]
Phosphate COP(=O)(O)O
[]
"""

