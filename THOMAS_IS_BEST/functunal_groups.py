from rdkit import Chem


functional_groups_smarts = {
    'Alcohol': 'C[Oh1]',
    'Aldehyde': 'C[Ch1]=O',
    'Ketone': 'CC(=O)C',
    'Carboxylic Acid': 'CC(=O)[Oh1]',
    'Ester': 'CC(=O)OC',
    'Ether': 'COC',
    'Amide': 'CC(=O)N',
    'Amine': 'C[Nh2]',
    'Nitrile': 'C#N',
    'Chloride': 'CCl',
    'Bromide': 'CBr',
    'Fluoride': 'CF',
    'Iodide': 'CI',
    'Alkene': 'C=C',
    'Alkyne': 'C#C',
    'Imine': 'C=NC',
    'Amino acid': 'CC(N)C(=O)O',
    'Thiol': 'C[Sh1]',
    'Sulfides': 'CSC',
    'Acyl Chloride': 'CC(=O)Cl',
    'Anhydride': 'CC(=O)OC(=O)C',
    'Nitro': 'C[N+](=O)[O-]',
    'Enamine': 'C=CN',
    'Imide': 'C(=O)NC(=O)C',
    'Azide': 'CNNN',
    'Enol': 'C=C([Oh1])C',
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


for name, smarts in functional_groups_smarts.items():
    try:
        mol = Chem.MolFromSmarts(smarts)
    except:
        print(f"Invalid SMARTS pattern for {name}: {smarts}")
        continue
    for name2, smarts2 in functional_groups_smarts.items():
        try:
            mol2 = Chem.MolFromSmarts(smarts2)
        except:
            print(f"Invalid SMARTS pattern for {name2}: {smarts2}")
            continue
        if mol.HasSubstructMatch(mol2):
            print(name, name2)
    
