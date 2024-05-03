from rdkit import Chem

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


mol_smi = input('SMILEs: ')
mol = Chem.MolFromSmiles(mol_smi)
list_name = []
list_index = []




for name, smiles in functional_groups_smiles.items():
    print(smiles)
    pattern = Chem.MolFromSmiles(smiles)

    if mol.HasSubstructMatch(pattern):
        print('found')
        indexes = mol.GetSubstructMatch(pattern)
        list_index.append(indexes)
        list_name.append(name)

#list name
# list indexes
list_indexes_copy = list_index.copy()
list_name_copy = list_name.copy()

for i in range(len(list_index)):
    True

print(list_index)
            
    
