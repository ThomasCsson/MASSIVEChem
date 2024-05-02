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
    'Sulfonic acid': 'CS(=O)(=O)O',
    'Sulfanote ester': 'CS(=O)(=O)OC',
    'Thioester': 'C(=O)SC',
    'Phosphine': 'CP',
    'Phosphate ester': 'COP(=O)(O)O',
    'Benzene ring': 'C1=CC=CC=C1'
}

def subgroup_value(mol_smi, dict_functional_groups):
    list_contained_subgroups, list_contained_subgroups_values = [], []
    mol = Chem.MolFromSmiles(mol_smi)

    for SMILES, value in dict_functional_groups.items():
        substruct = Chem.MolFromSmiles(dict_functional_groups[SMILES])  # Get the SMILES string from the dictionary

        if mol.HasSubstructMatch(substruct):
            list_contained_subgroups.append(SMILES)
            list_contained_subgroups_values.append(value)
            indices = mol.GetSubstructMatch(substruct)
            print(indices)
        else:
            list_contained_subgroups.append(f'NOT {SMILES}')

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
    for group, smiles in functional_groups.items():
        pattern = Chem.MolFromSmiles(smiles)
        if mol.HasSubstructMatch(pattern):
            found_groups.append(group)
            # Add related functional groups to remove list
            if group == 'Ketone':
                groups_to_remove.append('Aldehyde')
            if group == 'Carboxylic acid':
                groups_to_remove.append('Aldehyde')
                groups_to_remove.append('Alcohol')
            if group == 'Ester':
                groups_to_remove.append('Carboxylic acid')
                groups_to_remove.append('Ether')
            if group == 'Ether':
                groups_to_remove.append('Alcohol')
            if group == 'Imine':
                groups_to_remove.append('Amine')
            if group == 'Amide':
                groups_to_remove.append('Aldehyde')
            if group == 'Acyl Chloride':
                groups_to_remove.append('Aldehyde')
            if group == 'Thioester':
                groups_to_remove.append('Ester')
            if group == 'Sulfide':
                groups_to_remove.append('Thiol')
            if group == 'Azide':
                groups_to_remove.append('Amine')
    # Remove related functional groups after checking all groups
    print(found_groups)
    for group in groups_to_remove:
        found_groups.remove(group)
    return found_groups

print(check_functional_groups(input_mol,functional_groups_smiles))



