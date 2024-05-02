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
def check_functional_groups(molecule, functional_groups):
    found_groups = []
    mol = Chem.MolFromSmiles(molecule)
    if mol is None:
        return found_groups
    groups_to_remove = []
    for group, smarts in functional_groups.items():
        pattern = Chem.MolFromSmarts(smarts)
        if mol.HasSubstructMatch(pattern):
            if group == 'Alcohol':
                if mol.GetSubstructMatch(pattern) == mol.GetSubstructMatch(Chem.MolFromSmiles('CC(=O)O')):
                    pass
                if mol.GetSubstructMatch(pattern) == mol.GetSubstructMatch(Chem.MolFromSmiles('COC')):
                    pass
                if mol.GetSubstructMatch(pattern) == mol.GetSubstructMatch(Chem.MolFromSmiles('CC(N)C(=O)O')):
                    pass
                else:
                    found_groups.append(group)

            if group == 'Aldehyde':
                if mol.GetSubstructMatch(pattern) == mol.GetSubstructMatch(Chem.MolFromSmiles('CC(=O)C')):
                    pass
                if mol.GetSubstructMatch(pattern) == mol.GetSubstructMatch(Chem.MolFromSmiles('CC(=O)O')):
                    pass
                if  mol.GetSubstructMatch(pattern) == mol.GetSubstructMatch(Chem.MolFromSmiles('CC(=O)OC')):
                    pass
                else:
                    found_groups.append(group)
            if group == 'Carboxylic acid':
                if mol.GetSubstructMatch(pattern) == mol.GetSubstructMatch(Chem.MolFromSmiles('CC(=O)OC')):
                    pass
                else:
                    found_groups.append(group)
            if group == 'Ether':
                pass
            else:
                found_groups.append(group)
            print(found_groups)
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


    print(check_functional_groups(input_mol,functional_groups_smarts))


def smiles_to_smarts(dict):
    dict2 = {}
    for x,y in dict.items():
        mol = Chem.MolFromSmiles(y)
        smart = Chem.MolToSmarts(mol)
        dict2[x] = smart
    return dict2

"""


"""

