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
    for group, smarts in functional_groups.items():
        pattern = Chem.MolFromSmarts(smarts)
        if mol.HasSubstructMatch(pattern):
            if group == 'Alcohol':
                acid_pattern = Chem.MolFromSmarts(functional_groups['Carboxylic Acid'])
                if mol.HasSubstructMatch(acid_pattern):
                    continue
            else:
                found_groups.append(group)
    return found_groups



def smiles_to_smarts(dict):
    dict2 = {}
    for x,y in dict.items():
        mol = Chem.MolFromSmiles(y)
        smart = Chem.MolToSmarts(mol)
        dict2[x] = smart
    return dict2

functional_groups_smarts = {
    'Alcohol': 'CO',
    'Aldehyde': '[CX3H1](=O)',
    'Ketone': 'C(=O)C',
    'Carboxylic Acid': '[CX3](=O)[OX2H1]',
    'Ester': '[CX3](=O)[OX2H0][CX4]',
    'Ether': 'COC',
    'Amide': '[CX3](=O)[NX3H2]',
    'Amine': 'CN',
    'Nitrile': 'CC#N',
    'Chloride': '[Cl]',
    'Bromide': '[Br]',
    'Fluoride': '[F]',
    'Iodide': '[I]',
    'Alkene': '[CX3]=[CX3]',
    'Alkyne': '[CX2]#[CX2]',
    'Imine': 'C=NC',
    'Amino acid': 'CC(N)C(=O)O',
    'Thiol': '[SX2H1]',
    'Sulfides': 'CSC',
    'Acyl Chloride': '[CX3](=O)[Cl]',
    'Anhydride': '[CX3](=O)[OX2][CX3](=O)',
    'Nitro': 'C[N+](=O)[O-]',
    'Enamine': '[NX3][CX3]=[CX3]',
    'Imide': '[CX3](=[OX1])[NX3][CX3](=[OX1])',
    'Azide': 'CNNN',
    'Enol': 'C=C(O)C',
    'Hemiacetal': 'CC(O)(O)C',
    'Carbonate': 'OC(=O)O',
    'Disulfide': '[SX2][SX2]',
    'Sulfoxide': 'CS(=O)C',
    'Sulfone': 'CS(=O)(=O)C',
    'Sulfonic acid': 'CS(=O)(=O)O',
    'Thioester': 'C(=O)SC',
    'Phosphine': '[PX3]',
    'Phosphate': 'COP(=O)(O)O',
    'Benzene': '[c]',
}


for mol,smile in functional_groups_smiles.items():
    print(mol,smile)
    print(check_functional_groups(smile,functional_groups_smarts))


# Import necessary libraries
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np


# Function to convert SMILES to molecular fingerprints
def smiles_to_fingerprint(smiles):
    mol = Chem.MolFromSmiles(smiles)
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
    return np.array(fp)

# Function to predict functional groups
def predict_functional_groups(smiles, model):
    fingerprint = smiles_to_fingerprint(smiles)
    prediction = model.predict(np.expand_dims(fingerprint, axis=0))
    # Process prediction and return functional groups
    return prediction

# Load pre-trained deep learning model
# Replace 'model_path' with the path to your pre-trained model
model = tf.keras.models.load_model('model_path')

# Example usage
input_smiles = "CCO"
functional_groups = predict_functional_groups(input_smiles, model)
print("Predicted functional groups:", functional_groups)
