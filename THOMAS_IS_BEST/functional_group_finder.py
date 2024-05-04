from rdkit import Chem
def functional_group_finder(mol_smi):

    #---------------------------------------------------------------------------------------------#
    '''
    functional_group_finder(mol_smi)
    
    Input: molecule under SMILEs representation
    
    Output: list containing every functionl group contained (if a functional group is contained twice in the molecule, it will appear twice in this list)
    '''
    #---------------------------------------------------------------------------------------------#


    functional_groups_contained, mol_in = [], Chem.MolFromSmiles(mol_smi)
    functional_groups_smarts = {
        'Alcohol': 'C[Oh1+0]','Aldehyde': 'C[Ch1]=O','Ketone': 'CC(=O)C','Carboxylic Acid': 'CC(=O)[Oh1]',
        'Ester': 'CC(=O)[Oh0]','Ether': '*[Oh0]*','Amide': 'C(=O)N','Amine': '[C][Nh2]','Nitrile': 'C#N',
        'Chloride': 'Cl','Bromide': 'Br','Fluoride': 'F','Iodide': 'I','Alkene': 'C=C','Alkyne': 'C#C',
        'Imine': 'C=N*','Amino acid': '[Nh2][Ch1*]C(=O)O','Proline': '[Nh1][Ch1*]C(=O)O','Thiol': '[Sh1]',
        'Sulfide': '*[Sh0]*','Acyl Chloride': 'CC(=O)Cl','Anhydride': '*[Ch0](=O)O[Ch0](=O)*','Nitro': 'C[N+](=O)[O-]',
        'Enamine': 'C=C[Nh0]','Enamine2': 'C=C[Nh1]','Enamine3': 'C=C[Nh2]','Imide': 'C(=O)NC(=O)*','Azide': 'CNNN',
        'Enol': 'C=C([Oh1])C','Hemiacetal': 'CC(O)(O)C','Carbonate': '[Oh0]C(=O)[Oh0]','Carbonate2': '[Oh1]C(=O)[Oh1]',
        'Disulfide': 'CSSC','Sulfoxide': 'CS(=O)C','Sulfone': '*[So2](=O)(=O)*','Sulfonic acid': '*S(=O)(=O)[Oh1]',
        'Thioester': 'C(=O)S*','Phosphine': '*[Po0](*)*','Phosphate': '*OP(=O)(O)O','Benzene': 'c1ccccc1'
    }

    
    for name, smarts in functional_groups_smarts.items():
        functional_group = Chem.MolFromSmarts(smarts)
        if mol_in.HasSubstructMatch(functional_group):
            for _ in range(len(mol_in.GetSubstructMatches(functional_group))):
                functional_groups_contained.append(name)
            
    if 'Carboxylic Acid' in functional_groups_contained:
        functional_groups_contained.remove('Alcohol')
    if 'Ester' in functional_groups_contained:
        functional_groups_contained.remove('Ether')
    if 'Phosphate' in functional_groups_contained:
        functional_groups_contained.remove('Ether')
    if 'Thioester' in functional_groups_contained:
        functional_groups_contained.remove('Sulfide')
    if 'Sulfonic acid' in functional_groups_contained:
        functional_groups_contained.remove('Sulfide')
    if 'Sulfoxide' in functional_groups_contained:
        functional_groups_contained.remove('Sulfide')
    if 'Acyl Chloride' in functional_groups_contained:
        functional_groups_contained.remove('Chloride')
    if 'Anhydride' in functional_groups_contained:
        functional_groups_contained.remove('Ester')
        functional_groups_contained.remove('Ester')
    if 'Enamine2' in functional_groups_contained:
        functional_groups_contained.remove('Enamine2')
        functional_groups_contained.append('Enamine')
    if 'Enamine3' in functional_groups_contained:
        functional_groups_contained.remove('Enamine3')
        functional_groups_contained.remove('Amine')
        functional_groups_contained.append('Enamine')
    if 'Imide' in functional_groups_contained:
        functional_groups_contained.remove('Amide')
        functional_groups_contained.remove('Amide')
    if 'Enol' in functional_groups_contained:
        functional_groups_contained.remove('Alkene')
        functional_groups_contained.remove('Alcohol')
    if 'Hemiacetal' in functional_groups_contained:
        functional_groups_contained.remove('Alcohol')
        functional_groups_contained.remove('Alcohol')
    if 'Carbonate2' in functional_groups_contained:
        functional_groups_contained.remove('Alcohol')
        functional_groups_contained.remove('Alcohol')
        functional_groups_contained.remove('Carbonate2')
        functional_groups_contained.append('Carbonate')
    if 'Disulfide' in functional_groups_contained:
        functional_groups_contained.remove('Sulfide')
        functional_groups_contained.remove('Sulfide')
    if 'Amine2' in functional_groups_contained:
        functional_groups_contained.remove('Amine2')
        functional_groups_contained.append('Amine')
    
    return functional_groups_contained


mol_smi = input('SMILEs: ')
print(functional_group_finder(mol_smi))
