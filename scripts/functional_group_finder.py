from rdkit import Chem

def functional_group_finder(mol_smi):

    #---------------------------------------------------------------------------------------------#
    '''
    functional_group_finder(mol_smi)
    
    Input: molecule under SMILEs representation
    
    Output: list containing every functionl group contained (if a functional group is contained twice in the molecule, it will appear twice in this list)
    '''
    #---------------------------------------------------------------------------------------------#

    if not mol_smi:
        raise ValueError("Enter a non-empty string as argument")
    # initiate variables
    functional_groups_contained, mol_in = [], Chem.MolFromSmiles(mol_smi)

    # dictionnary of all the considered functional groups to check (some might be missing)

    functional_groups_smarts = {
        'Alcohol': 'C[Oh1+0]',
        'Aldehyde': 'C[Ch1]=O',
        'Ketone': 'CC(=O)C',
        'Carboxylic Acid': 'CC(=O)[Oh1]',
        'Ester': 'CC(=O)[Oh0]',
        'Ether': '*[Oh0]*',
        'Amide': 'C(=O)N',
        'Amine': '[C][N]',
        'Nitrile': 'C#N',
        'Chloride': 'Cl',
        'Bromide': 'Br',
        'Fluoride': 'F',
        'Iodide': 'I',
        'Alkene': 'C=C',
        'Alkyne': 'C#C',
        'Imine': 'C=N*',
        'Amino acid': '[Nh2][Ch1*]C(=O)O',
        'Proline': '[Nh1][Ch1*]C(=O)O',
        'Thiol': '[Sh1]',
        'Sulfide': '*[Sh0]*',
        'Acyl Chloride': 'CC(=O)Cl',
        'Anhydride': '*[Ch0](=O)O[Ch0](=O)*',
        'Nitro': 'C[N+](=O)[O-]',
        'Enamine': 'C=C[Nh0]',
        'Enamine2': 'C=C[Nh1]',
        'Enamine3': 'C=C[Nh2]',
        'Imide': 'C(=O)NC(=O)*',
        'Azide': 'CNNN',
        'Enol': 'C=C([Oh1])C',
        'Hemiacetal': 'CC(O)(O)C',
        'Carbonate': '[Oh0]C(=O)[Oh0]',
        'Carbonate2': '[Oh1]C(=O)[Oh1]',
        'Disulfide': 'CSSC',
        'Sulfoxide': 'CS(=O)C',
        'Sulfone': '*[So2](=O)(=O)*',
        'Sulfonic acid': '*S(=O)(=O)[Oh1]',
        'Thioester': 'C(=O)S*',
        'Phosphine': '*[Po0](*)*',
        'Phosphate': '*OP(=O)(O)O',
        'Benzene': 'c1ccccc1',
        'Peroxide':'C[Oh0][Oh0]C'
    }

    # check that the substructure from functional_groups_smarts are contained in mol_smi

    for name, smarts in functional_groups_smarts.items():
        if mol_in.HasSubstructMatch(Chem.MolFromSmarts(smarts)):
            for _ in range(len(mol_in.GetSubstructMatches(Chem.MolFromSmarts(smarts)))):
                functional_groups_contained.append(name)

    # exceptions for conflicts during the iteration of functional groups
    for functional_group in functional_groups_contained:
        if 'Ester' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Ether' in functional_groups_contained:
                    functional_groups_contained.remove('Ether')
        elif functional_group == 'Carboxylic Acid':
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Alcohol' in functional_groups_contained:
                    functional_groups_contained.remove('Alcohol')
        elif 'Phosphate' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Ether' in functional_groups_contained:
                    functional_groups_contained.remove('Ether')
        elif 'Thioester' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Sulfide' in functional_groups_contained:
                    functional_groups_contained.remove('Sulfide')
        elif 'Sulfonic acid' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Sulfide' in functional_groups_contained:
                    functional_groups_contained.remove('Sulfide')
        elif 'Sulfoxide' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Sulfide' in functional_groups_contained:
                    functional_groups_contained.remove('Sulfide')
        elif 'Acyl Chloride' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Sulfide' in functional_groups_contained:
                    functional_groups_contained.remove('Sulfide')
        elif 'Anhydride' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Ester' in functional_groups_contained:
                    functional_groups_contained.remove('Ester')
                if 'Ester' in functional_groups_contained:
                    functional_groups_contained.remove('Ester')
                if 'Ether' in functional_groups_contained:
                    functional_groups_contained.append('Ether')
        elif 'Enamine2' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Enamine2' in functional_groups_contained:
                    functional_groups_contained.remove('Enamine2')
                if 'Enamine' in functional_groups_contained:
                    functional_groups_contained.append('Enamine')
        elif 'Enamine3' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Enamine3' in functional_groups_contained:
                    functional_groups_contained.remove('Enamine3')
                if 'Amine' in functional_groups_contained:
                    functional_groups_contained.remove('Amine')
                functional_groups_contained.append('Enamine')
        elif 'Imide' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Amide' in functional_groups_contained:
                    functional_groups_contained.remove('Amide')
                if 'Amide' in functional_groups_contained:
                    functional_groups_contained.remove('Amide')
        elif 'Enol' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Alkene' in functional_groups_contained:
                    functional_groups_contained.remove('Alkene')
                if 'Alcohol' in functional_groups_contained:
                    functional_groups_contained.remove('Alcohol')
        elif 'Hemiacetal' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Alcohol' in functional_groups_contained:
                    functional_groups_contained.remove('Alcohol')
                if 'Alcohol' in functional_groups_contained:
                    functional_groups_contained.remove('Alcohol')
        elif 'Carbonate2' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Alcohol' in functional_groups_contained:
                    functional_groups_contained.remove('Alcohol')
                if 'Alcohol' in functional_groups_contained:
                    functional_groups_contained.remove('Alcohol')
                if 'Carbonate2' in functional_groups_contained:
                    functional_groups_contained.remove('Carbonate2')
                if 'Carbonate' in functional_groups_contained:
                    functional_groups_contained.append('Carbonate')
        elif 'Disulfide' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Sulfide' in functional_groups_contained:
                    functional_groups_contained.remove('Sulfide')
                if 'Sulfide' in functional_groups_contained:
                    functional_groups_contained.remove('Sulfide')
        elif 'Peroxide' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Ether' in functional_groups_contained:
                    functional_groups_contained.remove('Ether')
                if 'Ether' in functional_groups_contained:
                    functional_groups_contained.remove('Ether')
        elif 'Amide' == functional_group:
            for _ in range (functional_groups_contained.count(functional_group)):
                if 'Amine' in functional_groups_contained:
                    functional_groups_contained.remove('Amine')
                if 'Amine' in functional_groups_contained:
                    functional_groups_contained.remove('Amine')

    
    
    return functional_groups_contained
print(functional_group_finder('NCCC(=O)NCCCNCCCN'))
